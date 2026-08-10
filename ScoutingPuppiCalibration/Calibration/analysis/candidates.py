"""Per-candidate PUPPI weight sanity checks.

IMPORTANT CAVEAT: pvAssocQuality comes from Run3ScoutingParticleToPackedCandidate
Producer's "vtxass" ValueMap, itself derived from particle.vertex() -- an
ONLINE proxy computed at HLT by matching each candidate's position to the
nearest hltPixelVertices vertex. It is NOT MC truth (no TrackingParticle/
TrackingVertex truth link exists for scouting candidates). Everything in this
module is a secondary sanity check on whether PUPPI weights behave sensibly
relative to that online proxy -- NOT a true-LV-vs-true-PU ROC/efficiency
study. Label plots accordingly; do not use these as the primary tuning metric
(use response.py's genjet-truth comparisons for that).

pat::PackedCandidate::PVAssociationQuality values seen here (see
DataFormats/PatCandidates/interface/PackedCandidate.h):
  0 = NotReconstructedPrimary (unassociated)
  5 = CompatibilityDz         (associated to a pileup vertex, per this producer)
  7 = UsedInFitTight          (associated to the leading/PV vertex, per this producer)
"""

import awkward as ak
import numpy as np

PV_ASSOC_LABELS = {
    0: "NotReconstructedPrimary (unassociated)",
    5: "CompatibilityDz (online proxy: PU vertex)",
    7: "UsedInFitTight (online proxy: PV)",
}


def weight_by_pvassoc(cand_table, variant, eta_bins=((0.0, 1.3), (1.3, 2.5))):
    """Per-eta-bin dict of {pvAssocQuality value: weight array} for one
    PUPPI variant's puppiWeight_<variant> column. Charged candidates only
    (weight is only meaningful/nontrivial for charged/neutral distinction
    when there is a pvAssocQuality label at all -- neutral candidates are
    always pvAssocQuality==0 by construction, see Run3ScoutingParticleTo
    PackedCandidateProducer.cc, so are reported separately).
    """
    key = "puppiWeight_%s" % variant
    if key not in cand_table.fields:
        raise KeyError("%s not found in candidate table -- check the variant name" % key)

    # cand_table's fields are jagged (variable candidates per event), same as
    # the jet tables -- flatten before converting to numpy (np.asarray on a
    # jagged awkward array raises "subarray lengths are not regular").
    w = np.asarray(ak.flatten(cand_table[key]))
    pvq = np.asarray(ak.flatten(cand_table["pvAssocQuality"]))
    eta = np.asarray(ak.flatten(cand_table["eta"]))
    pdgId = np.asarray(ak.flatten(cand_table["pdgId"]))
    is_neutral = np.isin(np.abs(pdgId), [22, 130, 1, 2, 0])

    out = {}
    for lo, hi in eta_bins:
        ebin = (np.abs(eta) >= lo) & (np.abs(eta) < hi)
        entry = {"neutral": w[ebin & is_neutral]}
        for q, label in PV_ASSOC_LABELS.items():
            sel = ebin & (~is_neutral) & (pvq == q)
            if np.count_nonzero(sel):
                entry[label] = w[sel]
        out[(lo, hi)] = entry
    return out


def summarize_weight_distributions(cand_table, variants):
    """Small text summary (mean/median weight per category, per variant) --
    a quick numeric sanity check before plotting: PV-associated charged
    candidates and neutrals near jet cores should trend toward weight~1,
    PU-vertex-associated charged candidates should be ~0 (PUPPI keeps them
    at their online-proxy-assigned value, since charged-PU is handled by
    hard removal, not reweighting).
    """
    lines = []
    for variant in variants:
        lines.append("variant=%s" % variant)
        by_eta = weight_by_pvassoc(cand_table, variant)
        for (lo, hi), cats in by_eta.items():
            lines.append("  |eta| in [%.1f, %.1f):" % (lo, hi))
            for label, arr in cats.items():
                if len(arr) == 0:
                    continue
                lines.append(
                    "    %-45s n=%8d  mean=%.3f  median=%.3f"
                    % (label, len(arr), np.mean(arr), np.median(arr))
                )
    return "\n".join(lines)
