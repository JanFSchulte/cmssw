"""Offline PUPPI parameter recalibration: recompute per-candidate weights
from the raw alpha diagnostics (puppi_recompute.py) for a hypothetical
(MedEtaSF, RMSEtaSF, MinNeutralPt, MinNeutralPtSlope), recluster PUPPI jets
from the reweighted candidates with genuine FastJet anti-kt(0.4) (matching
ak4PFJets.clone(applyWeight=True, jetPtMin=20)), match to gen jets, and
report response/resolution/fake-rate exactly like response.py/summary.py do
for the production PUPPI variants -- so a hypothetical parameter's jet-level
impact can be evaluated without rerunning PuppiProducer.

Scope/caveats (see puppi_recompute.py's docstring for the full derivation):
  - Models both the central AND forward algo blocks (v6 data only -- v5's
    forward rawAlpha was wrong, see puppi_recompute.py's module docstring).
    Forward candidates still legitimately contribute to central jets' R=0.4
    cones near the |eta|=2.5 edge (confirmed empirically: 93% of jets that
    didn't kinematically match production's plain-jet collection had
    |eta|>2.0), so they're reclustered with their own recomputed weight, not
    dropped or held fixed at a stored value.
  - MedEtaSF/RMSEtaSF act on every id==0 candidate via the validated chi2
    formula and are essentially exactly reproduced (confirmed via a direct
    C++-level debug probe, CommonTools/PileupAlgos temporarily instrumented
    to expose PuppiAlgo::compute()'s internal pVal/chi2 -- 89.5% of id==0
    central candidates match the stored weight even BEFORE applying any
    floor, and 100% of the pre-floor mismatches are explained by the
    neutral-pt-floor firing, not the chi2 math itself).
  - MinNeutralPt/MinNeutralPtSlope's floor needs the event's true
    reconstructed-vertex count (PuppiProducer's iPUProxy), which was never
    saved as a branch -- approximated here with Pileup_nPU+1 (MC truth
    pileup-interaction count, not the actual reconstructed-vertex count with
    its own finite efficiency). Results on this specific axis carry that
    approximation; MedEtaSF/RMSEtaSF results do not.
"""

import awkward as ak
import numpy as np
import fastjet

from . import io, jec
from .puppi_recompute import recompute_weight, CURRENT_FWD_MED_SF, CURRENT_FWD_RMS_SF, \
    CURRENT_FWD_MIN_NEUTRAL_PT, CURRENT_FWD_MIN_NEUTRAL_PT_SLOPE

_JETDEF = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
# matches PFJetParameters_cfi.py defaults used by ak4PFJets (doAreaFastjet=True,
# Active_Area_Repeats=1, GhostArea=0.01, Ghost_EtaMax=5.0)
_AREADEF = fastjet.AreaDefinition(fastjet.active_area, fastjet.GhostedAreaSpec(5.0, 1, 0.01))


def load_refit_inputs(files, prefix="ScoutingPuppiCalibCand", step_size="100 MB"):
    """pt/eta/phi/mass/charge/pdgId/pvAssocQuality/rawAlpha/alphaMed/alphaRms
    + Pileup_nPU + Rho_fixedGridRhoFastjetAll, jagged/flat per event --
    everything recompute_weight + reclustering + JEC need."""
    import uproot
    cand = io.load_candidates(files, ["nominal"], prefix=prefix, step_size=step_size)
    extra = uproot.concatenate(
        [f + ":Events" for f in files],
        filter_name=lambda n: n in ("Pileup_nPU", "Rho_fixedGridRhoFastjetAll"),
        step_size=step_size,
    )
    nPU = extra["Pileup_nPU"]
    rho = extra["Rho_fixedGridRhoFastjetAll"] if "Rho_fixedGridRhoFastjetAll" in extra.fields else None
    return cand, nPU, rho


def _delta_r2(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = np.abs(phi1 - phi2)
    dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    return deta * deta + dphi * dphi


def _cluster_event(pt, eta, phi, mass, weight, rho, ptmin=20.0, weight_pt_cut=1e-3, apply_jec=True):
    wpt = pt * weight
    keep = wpt > weight_pt_cut
    n = np.count_nonzero(keep)
    if n == 0:
        return [], [], [], []
    wpt_k = wpt[keep]
    eta_k = eta[keep]
    phi_k = phi[keep]
    m_k = mass[keep] * weight[keep]  # mass scales with weight too (PuppiProducer scales full p4)
    px = wpt_k * np.cos(phi_k)
    py = wpt_k * np.sin(phi_k)
    pz = wpt_k * np.sinh(eta_k)
    e = np.sqrt(px * px + py * py + pz * pz + np.clip(m_k, 0, None) ** 2)
    pjs = [fastjet.PseudoJet(float(px[i]), float(py[i]), float(pz[i]), float(e[i])) for i in range(n)]
    # raw clustering + jetPtMin=20 threshold applied on the UNCORRECTED jet,
    # matching PFJetParameters_cfi's jetPtMin semantics (a clustering-level
    # cut, before JEC) -- production's patJetCorrFactors JEC is then applied
    # on top and CAN push the stored/corrected pt below 20 for some jets.
    cs = fastjet.ClusterSequenceArea(pjs, _JETDEF, _AREADEF)
    jets = cs.inclusive_jets(ptmin)
    jpt, jeta, jphi, jmass = [], [], [], []
    for j in jets:
        pt_raw, eta_j, phi_j, area_j, mass_raw = j.pt(), j.eta(), j.phi(), j.area(), j.m()
        if phi_j > np.pi:
            phi_j -= 2 * np.pi
        jec_factor = 1.0
        if apply_jec:
            jec_factor = jec.correction_factor(pt_raw, eta_j, phi_j, area_j, rho)
        # JEC scales the whole 4-vector uniformly (matching how patJetCorrFactors'
        # MomentumScaleFactor is actually applied), so mass scales by the same
        # factor as pt, not left at its raw (uncorrected) value.
        jpt.append(pt_raw * jec_factor)
        jeta.append(eta_j)
        jphi.append(phi_j)
        jmass.append(mass_raw * jec_factor)
    return jpt, jeta, jphi, jmass


def _match_to_gen(jpt, jeta, jphi, gpt, geta, gphi, maxDeltaR=0.4):
    """Greedy nearest-neighbor matching, ambiguity-resolved (each gen jet used
    at most once), matching PhysicsTools.PatAlgos patJetGenJetMatch defaults
    (maxDeltaR=0.4, resolveAmbiguities=True)."""
    nj, ng = len(jpt), len(gpt)
    genJetIdx = [-1] * nj
    if nj == 0 or ng == 0:
        return genJetIdx
    pairs = []
    maxDR2 = maxDeltaR * maxDeltaR
    for i in range(nj):
        for k in range(ng):
            dr2 = _delta_r2(jeta[i], jphi[i], geta[k], gphi[k])
            if dr2 < maxDR2:
                pairs.append((dr2, i, k))
    pairs.sort(key=lambda p: p[0])
    used_gen = set()
    used_jet = set()
    for dr2, i, k in pairs:
        if i in used_jet or k in used_gen:
            continue
        genJetIdx[i] = k
        used_jet.add(i)
        used_gen.add(k)
    return genJetIdx


def refit_and_recluster(cand, nPU, genjets, rho=None, medEtaSF=1.0, rmsEtaSF=1.0,
                         minNeutralPt=0.2, minNeutralPtSlope=0.015,
                         fwd_medEtaSF=CURRENT_FWD_MED_SF[0], fwd_rmsEtaSF=CURRENT_FWD_RMS_SF[0],
                         fwd_minNeutralPt=CURRENT_FWD_MIN_NEUTRAL_PT[0],
                         fwd_minNeutralPtSlope=CURRENT_FWD_MIN_NEUTRAL_PT_SLOPE[0],
                         ptmin=20.0, maxDeltaR=0.4, apply_jec=True):
    """Returns a jets awkward Array with fields {pt, eta, phi, mass, genJetIdx},
    matching (a subset of) io.load_jets's output shape, so response.py/
    summary.py work unchanged.

    Requires v6+ candidate data: recompute_weight now models BOTH the
    central and forward algo blocks (see its docstring), which only works
    with v6's fixed per-region rawAlpha -- v5's forward rawAlpha is
    meaningless (always the central block's value). This function used to
    hard-override forward candidates' weight with the stored
    puppiWeight_nominal instead of recomputing them, back when the central
    block was all this pipeline could model; that override is gone now that
    the forward block is modeled too.
    """
    pt = cand["pt"]
    eta = cand["eta"]
    phi = cand["phi"]
    mass = cand["mass"]
    charge = cand["charge"]
    pdgId = cand["pdgId"]
    pvq = cand["pvAssocQuality"]
    rawAlpha = cand["puppiRawAlpha"]
    alphaMed = cand["puppiAlphaMed"]
    alphaRms = cand["puppiAlphaRms"]

    # nPU proxy: PuppiProducer's true iPUProxy (count of reconstructed
    # vertices passing ndof/z cuts) was never saved as a branch. Empirically
    # fit against stored puppiWeight_nominal (id0, raw!=0, has-PU-reference
    # candidates): 0.5*Pileup_nPU matches to 97.7% (weighted-pT ratio 0.99),
    # vs Pileup_nPU+1's 77.7% (ratio 0.85) -- reconstructed-vertex efficiency
    # at HLT is well below 100%, so half the true MC pileup-interaction count
    # is a much closer real-world estimate than "true count + 1".
    nPU_proxy = 0.5 * nPU
    nPU_c = ak.broadcast_arrays(nPU_proxy, pt)[0]

    weight = recompute_weight(pt, eta, charge, pdgId, pvq, rawAlpha, alphaMed, alphaRms, nPU_c,
                               medEtaSF=medEtaSF, rmsEtaSF=rmsEtaSF,
                               minNeutralPt=minNeutralPt, minNeutralPtSlope=minNeutralPtSlope,
                               fwd_medEtaSF=fwd_medEtaSF, fwd_rmsEtaSF=fwd_rmsEtaSF,
                               fwd_minNeutralPt=fwd_minNeutralPt,
                               fwd_minNeutralPtSlope=fwd_minNeutralPtSlope)

    n_events = len(pt)
    out_pt, out_eta, out_phi, out_mass, out_genidx = [], [], [], [], []
    for ev in range(n_events):
        pt_ev = np.asarray(pt[ev])
        if len(pt_ev) == 0:
            out_pt.append([])
            out_eta.append([])
            out_phi.append([])
            out_mass.append([])
            out_genidx.append([])
            continue
        eta_ev = np.asarray(eta[ev])
        phi_ev = np.asarray(phi[ev])
        mass_ev = np.asarray(mass[ev])
        w_ev = np.asarray(weight[ev])
        rho_ev = float(rho[ev]) if rho is not None else 0.0
        jpt, jeta, jphi, jmass = _cluster_event(pt_ev, eta_ev, phi_ev, mass_ev, w_ev, rho_ev,
                                                 ptmin=ptmin, apply_jec=apply_jec and rho is not None)

        gpt = np.asarray(genjets["pt"][ev])
        geta = np.asarray(genjets["eta"][ev])
        gphi = np.asarray(genjets["phi"][ev])
        gsel = gpt > 10.0  # match jetMCTable's own genJetIdx convention (gen pt>10)
        gpt, geta, gphi = gpt[gsel], geta[gsel], gphi[gsel]
        genidx_local = _match_to_gen(jpt, jeta, jphi, gpt, geta, gphi, maxDeltaR=maxDeltaR)
        # map back from the gsel-filtered local gen index to the ORIGINAL
        # genjets[ev] index space, since response.py indexes genjets[idx] directly
        orig_idx = np.nonzero(gsel)[0]
        genidx = [int(orig_idx[k]) if k >= 0 else -1 for k in genidx_local]

        out_pt.append(jpt)
        out_eta.append(jeta)
        out_phi.append(jphi)
        out_mass.append(jmass)
        out_genidx.append(genidx)

    return ak.Array({
        "pt": ak.Array(out_pt), "eta": ak.Array(out_eta),
        "phi": ak.Array(out_phi), "mass": ak.Array(out_mass),
        "genJetIdx": ak.Array(out_genidx),
    })
