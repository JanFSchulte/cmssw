"""Jet reconstruction efficiency: for each GEN jet, was it matched by ANY
reco jet? (the inverse direction of response.py's reco->gen matching, which
only asks "what did this reco jet match to" and says nothing about gen jets
that were missed entirely).

genJetIdx (jetMCTable's convention, see PhysicsTools/NanoAOD/python/
jetMC_cff.py) is a per-RECO-jet index into the per-event GenJet collection.
There's no stored reverse map, so it's rebuilt here via an event-local
cartesian match: for each (genjet, recojet) pair in the same event, flag a
hit if recojet.genJetIdx == genjet's own local index, then reduce over the
reco axis. No new information needed -- this is fully derivable from data
already being collected.
"""

import awkward as ak
import numpy as np


def genjet_matched_flags(jets, genjets):
    """Per-genjet boolean: True if some reco jet in the same event has
    genJetIdx pointing at it. Shape matches genjets (jagged, one entry per
    genjet per event).
    """
    local_idx = ak.local_index(genjets["pt"], axis=1)
    g, r = ak.unzip(ak.cartesian([local_idx, jets["genJetIdx"]], axis=1, nested=True))
    return ak.any(g == r, axis=2)


def _binned_efficiency(x, matched, bins):
    idx = np.digitize(x, bins) - 1
    centers, effs, counts = [], [], []
    for b in range(len(bins) - 1):
        sel = idx == b
        n = np.count_nonzero(sel)
        if n < 20:
            continue
        centers.append(0.5 * (bins[b] + bins[b + 1]))
        effs.append(np.mean(matched[sel]))
        counts.append(n)
    return np.array(centers), np.array(effs), np.array(counts)


def efficiency_vs_pt(jets, genjets, pt_bins=None, eta_max=2.5):
    """Fraction of gen jets (|eta|<eta_max) with a matching reco jet, binned
    in gen-jet pt.
    """
    if pt_bins is None:
        pt_bins = np.array([15, 20, 25, 30, 40, 50, 70, 100, 150, 200, 300, 500])
    matched = genjet_matched_flags(jets, genjets)
    sel = np.abs(ak.flatten(genjets["eta"])) < eta_max
    pt = np.asarray(ak.flatten(genjets["pt"]))[sel]
    m = np.asarray(ak.flatten(matched))[sel]
    return _binned_efficiency(pt, m, pt_bins)


def efficiency_vs_pu(jets, genjets, nTrueInt, pu_bins=None, genjet_pt_min=30.0, eta_max=2.5):
    """Fraction of gen jets (pt>genjet_pt_min, |eta|<eta_max) with a matching
    reco jet, binned in true pileup.
    """
    if pu_bins is None:
        pu_bins = np.arange(0, 80, 5)
    matched = genjet_matched_flags(jets, genjets)
    pu_per_gj = ak.broadcast_arrays(nTrueInt, genjets["pt"])[0]
    sel = (np.abs(genjets["eta"]) < eta_max) & (genjets["pt"] > genjet_pt_min)
    pu = np.asarray(ak.flatten(pu_per_gj[sel]))
    m = np.asarray(ak.flatten(matched[sel]))
    return _binned_efficiency(pu, m, pu_bins)
