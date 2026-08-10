"""Offline reproduction of PuppiAlgo::compute() + PuppiProducer.cc's post-
processing, from the raw per-candidate diagnostics
(puppiRawAlpha/puppiAlphaMed/puppiAlphaRms) exposed by
PuppiRawAlphaToValueMapProducer, so MedEtaSF/RMSEtaSF/MinNeutralPt/
MinNeutralPtSlope can be re-scanned entirely offline without rerunning
PuppiProducer.

Traced directly from CommonTools/PileupAlgos/{src/PuppiAlgo.cc,
src/PuppiContainer.cc, plugins/PuppiProducer.cc} (CMSSW_16_1_0_pre4) and
Run3ScoutingParticleToPackedCandidateProducer.cc's baseline (non-improved)
vtxass/quality assignment.

REGION-AWARE (v6 data only): CRAB v5's PuppiRawAlphaToValueMapProducer had a
bug where every candidate's rawAlpha was read from the CENTRAL block
regardless of its own eta -- meaningless for forward (|eta|>=2.5) candidates
(wrong cone/rmsPtMin), which is why this module used to be central-only and
puppi_refit.py hard-overrode forward candidates' weight with the stored
puppiWeight_nominal instead of recomputing them. Fixed in
PuppiRawAlphaToValueMapProducer.cc (etaBoundaries parameter) for v6 --
rawAlpha now correctly reflects each candidate's own region. alphaMed/
alphaRms were ALWAYS correctly per-candidate-region (PuppiContainer::
calculatePuppiWeights fills them from fPuppiAlgo[pPupId].median()/.rms(),
already region-specific), so only rawAlpha needed the production fix; this
module's job is entirely on the analysis side, to correctly UN-scale that
already-SF-scaled alphaMed/alphaRms before applying a hypothetical new SF
(see recompute_weight's med/rms lines below).

Puppi_cff.py's production "algos" VPSet has 2 top-level entries: central
(|eta|<2.5, MedEtaSF=1.0/RMSEtaSF=1.0/MinNeutralPt=0.2/MinNeutralPtSlope=
0.015) and forward (a single PuppiAlgo object internally covering both
2.5<=|eta|<3.0, MedEtaSF=0.90/RMSEtaSF=1.20/MinNeutralPt=1.7/
MinNeutralPtSlope=0.08, and |eta|>=3.0, MedEtaSF=0.75/RMSEtaSF=0.95/
MinNeutralPt=2.0/MinNeutralPtSlope=0.08 -- see PuppiRawAlphaToValueMapProducer.cc's
own comment for why cone/rmsPtMin, and hence rawAlpha itself, don't
distinguish these two forward sub-bins even though MedEtaSF/RMSEtaSF/
MinNeutralPt do). In practice essentially all real scouting candidates that
reach the forward block sit in the first sub-bin (2.5-3.0) -- confirmed on
v6 data, max observed |eta| is exactly 3.0 -- so a forward parameter scan
only meaningfully probes that sub-bin; the second sub-bin's current values
are still applied correctly (for un-scaling only) to any candidate that
happens to have |eta|>=3.0, there just aren't real ones to calibrate against.

id-assignment reconstruction (PuppiProducer.cc:200-216, fUseVertexAssociation
branch, vertexAssociationQuality threshold default=0 so it never gates):
  charge==0                      -> id=0 (neutral, always goes through compute())
  pvAssocQuality==0  (NotReconstructedPrimary) -> id=0 (unassociated)
  pvAssocQuality==5  (CompatibilityDz)         -> id=2 (PU, baseline vtxass's
                                                   only "associated to a
                                                   non-leading vertex" tier)
  pvAssocQuality==7  (UsedInFitTight)          -> id=1 (PV, baseline vtxass's
                                                   only "associated to vtx 0" tier)
This 3-value-only correspondence was confirmed empirically against the CRAB
v4 candidate_weight_summary.txt (id=1/id=2 buckets have mean/median exactly
1.0/0.0 respectively, matching the fApplyCHS hard override) and is safe
because the baseline (pre-fix) vtxass producer only ever emits these three
quality codes for charged candidates.

weight formula (PuppiAlgo.cc:183-216, PuppiContainer.cc:298-330), central
block only (fNAlgos=1, fAlgoId=5, useExp=False so no external chi2 term,
puppiNoLep not used so id==3 never occurs):
  id==1                     -> weight = 1                         (CHS override)
  id==2                     -> weight = 0                         (CHS override)
  id==0:
    if event has NO id==2 candidate with |eta|<EtaMaxExtrap(2.0) that event:
      soft_w = 1            # PuppiAlgo::compute()'s "if (fNCount[i0]==0)
                             # return 1." -- the reference alpha distribution
                             # (built ONLY from id==2/PU-tagged candidates,
                             # PuppiAlgo.cc:100-103) is empty, so compute()
                             # bails out and returns 1 UNCONDITIONALLY for
                             # every id==0 candidate that event, regardless of
                             # its own alpha. Confirmed via direct C++ debug
                             # instrumentation: this triggers far more often
                             # for `nominal` specifically than one might
                             # expect, because nominal's OWN broken vtxass
                             # (see scouting_vertex_association_bug memory)
                             # starves the id==2 population that fNCount
                             # counts -- the same root bug motivating
                             # vtxAssocImproved also explains why nominal's
                             # weights are hard to reproduce from alpha alone.
    else:
      med = alphaMed_stored * (MedEtaSF_new / 1.0)   # stored already carries
      rms = alphaRms_stored * (RMSEtaSF_new / 1.0)   # the CURRENT central SF=1.0
      pVal = rawAlpha if rawAlpha != 0 else med       # algoId==5 zero-substitution
      chi2 = (pVal-med)*|pVal-med| / rms**2
      soft_w = chisquared_cdf(chi2, ndof=1)
    weight = soft_w
    if weight*pt < MinNeutralPt_new + MinNeutralPtSlope_new*nPUProxy:
      weight = 0                                    # "Basic Cuts" neutral-pt floor,
                                                       # applies to ANY id==0 candidate,
                                                       # not just true neutrals, and
                                                       # still applies even in the
                                                       # fNCount==0/soft_w=1 case above
    if pt>20 and |pdgId|==22 and |eta|<2.5: weight = 1   # PtMaxPhotons protection
    if weight < 0.01: weight = 0                    # MinPuppiWeight

KNOWN APPROXIMATION: nPUProxy should be the number of good reconstructed
vertices PuppiProducer counted that event (vtxNdofCut=4, vtxZCut=24cm) --
this was never saved as a branch (the alpha-diagnostics design in the prior
session missed this dependency). Empirically fit (puppi_refit.py) against
stored puppiWeight_nominal on the id0/raw!=0/has-PU-reference population:
0.5*Pileup_nPU matches to 97.7% (weighted-pT ratio 0.99) -- much better than
the originally-assumed Pileup_nPU+1 (77.7% match, ratio 0.85, systematically
too large -- HLT reconstructed-vertex efficiency is well below 100%, so
roughly half the true MC pileup-interaction count is a much closer estimate
of the actual reconstructed-vertex count than "true count + 1"). Still an
approximation, not exact -- treat fitted MinNeutralPt/MinNeutralPtSlope
values with appropriate caution.
"""

import awkward as ak
import numpy as np
from scipy.stats import chi2 as _chi2dist

PV_QUALITY_UNASSOCIATED = 0
PV_QUALITY_PU = 5   # CompatibilityDz
PV_QUALITY_PV = 7   # UsedInFitTight

MIN_PUPPI_WEIGHT = 0.01
PT_MAX_PHOTONS = 20.0
ETA_MAX_PHOTONS = 2.5

# Current production forward-region defaults (Puppi_cff.py), by sub-bin:
# [2.5<=|eta|<3.0, |eta|>=3.0]. Needed to UN-scale the already-SF-scaled
# alphaMed/alphaRms before applying a hypothetical new forward SF -- see
# module docstring.
FWD_ETA_BOUNDARY = 3.0
CURRENT_FWD_MED_SF = (0.90, 0.75)
CURRENT_FWD_RMS_SF = (1.20, 0.95)
CURRENT_FWD_MIN_NEUTRAL_PT = (1.7, 2.0)
CURRENT_FWD_MIN_NEUTRAL_PT_SLOPE = (0.08, 0.08)


def assign_id(charge, pvAssocQuality):
    """0=unassociated/neutral (soft alpha path), 1=PV (hard weight=1), 2=PU (hard weight=0)."""
    id_ = ak.zeros_like(charge, dtype=np.int32)
    charged = charge != 0
    id_ = ak.where(charged & (pvAssocQuality == PV_QUALITY_PV), 1, id_)
    id_ = ak.where(charged & (pvAssocQuality == PV_QUALITY_PU), 2, id_)
    # everything else (charge==0, or charged+quality==0/anything unrecognized) stays id=0
    return id_


def recompute_weight(pt, eta, charge, pdgId, pvAssocQuality, rawAlpha, alphaMed, alphaRms, nPUProxy,
                      medEtaSF=1.0, rmsEtaSF=1.0, minNeutralPt=0.2, minNeutralPtSlope=0.015,
                      fwd_medEtaSF=CURRENT_FWD_MED_SF[0], fwd_rmsEtaSF=CURRENT_FWD_RMS_SF[0],
                      fwd_minNeutralPt=CURRENT_FWD_MIN_NEUTRAL_PT[0],
                      fwd_minNeutralPtSlope=CURRENT_FWD_MIN_NEUTRAL_PT_SLOPE[0]):
    """Vectorized (awkward-array-safe) reproduction of the stored PUPPI
    weight for BOTH the central and forward eta regions, for arbitrary
    hypothetical (medEtaSF, rmsEtaSF, minNeutralPt, minNeutralPtSlope) in
    the central region and (fwd_medEtaSF, fwd_rmsEtaSF, fwd_minNeutralPt,
    fwd_minNeutralPtSlope) in the forward region (applied uniformly across
    both forward sub-bins -- see module docstring for why that's fine in
    practice). Pass the CURRENT production defaults (central: 1.0, 1.0, 0.2,
    0.015; forward: 0.90, 1.20, 1.7, 0.08) to reproduce puppiWeight_nominal
    exactly (modulo the nPUProxy approximation noted above) -- these are
    already every argument's default.

    Requires v6+ candidate data (correct per-region rawAlpha) -- on v5 data,
    forward candidates' rawAlpha is meaningless (always the central block's
    value), so any forward-region result from this function on v5 input is
    garbage.
    """
    id_ = assign_id(charge, pvAssocQuality)

    abs_eta = np.abs(eta)
    is_central = abs_eta < ETA_MAX_PHOTONS  # 2.5
    is_fwd2 = abs_eta >= FWD_ETA_BOUNDARY

    # UN-scale alphaMed/alphaRms by whatever SF is CURRENTLY baked into the
    # stored value for this candidate's region, then re-scale by the
    # hypothetical new SF. Central's current SF is 1.0 (no-op un-scale, kept
    # for symmetry/clarity).
    current_med_sf = ak.where(is_central, 1.0, ak.where(is_fwd2, CURRENT_FWD_MED_SF[1], CURRENT_FWD_MED_SF[0]))
    current_rms_sf = ak.where(is_central, 1.0, ak.where(is_fwd2, CURRENT_FWD_RMS_SF[1], CURRENT_FWD_RMS_SF[0]))
    new_med_sf = ak.where(is_central, medEtaSF, fwd_medEtaSF)
    new_rms_sf = ak.where(is_central, rmsEtaSF, fwd_rmsEtaSF)

    med = (alphaMed / current_med_sf) * new_med_sf
    rms = (alphaRms / current_rms_sf) * new_rms_sf
    pVal = ak.where(rawAlpha != 0, rawAlpha, med)
    diff = pVal - med
    # (diff)*|diff|/rms^2 = sign(diff)*diff^2/rms^2 -- can be negative when
    # pVal<med. ROOT::Math::chisquared_cdf(x<0,...)==0, and scipy's chi2.cdf
    # already returns 0 for negative x, so no special-casing needed.
    chi2val = diff * abs(diff) / (rms * rms)
    counts = ak.num(chi2val) if chi2val.ndim > 1 else None
    chi2val_np = ak.to_numpy(ak.flatten(chi2val, axis=None))
    soft_w = _chi2dist.cdf(chi2val_np, df=1)
    soft_w = ak.unflatten(soft_w, counts) if counts is not None else ak.Array(soft_w)

    # PuppiAlgo::compute()'s "if (fNCount[i0]==0) return 1." fallback: no
    # id==2/PU-tagged reference candidates within |eta|<2.0 this event ->
    # soft_w=1 unconditionally for every id==0 candidate that event.
    has_pu_ref = ak.any((id_ == 2) & (np.abs(eta) < 2.0), axis=1)
    has_pu_ref_c = ak.broadcast_arrays(has_pu_ref, pt)[0]
    soft_w = ak.where(has_pu_ref_c, soft_w, 1.0)

    weight = ak.where(id_ == 1, 1.0, ak.where(id_ == 2, 0.0, soft_w))

    region_minNeutralPt = ak.where(is_central, minNeutralPt, fwd_minNeutralPt)
    region_minNeutralPtSlope = ak.where(is_central, minNeutralPtSlope, fwd_minNeutralPtSlope)
    neutralPtCut = region_minNeutralPt + region_minNeutralPtSlope * nPUProxy
    below_floor = (weight * pt < neutralPtCut) & (id_ == 0)
    weight = ak.where(below_floor, 0.0, weight)

    is_photon_protected = (np.abs(pdgId) == 22) & (np.abs(eta) < ETA_MAX_PHOTONS) & (pt > PT_MAX_PHOTONS)
    weight = ak.where(is_photon_protected, 1.0, weight)

    weight = ak.where(weight < MIN_PUPPI_WEIGHT, 0.0, weight)
    return weight
