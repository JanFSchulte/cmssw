"""Jet response/resolution vs gen-jet pt and vs true pileup.

response = reco_jet.pt / matched_gen_jet.pt (matching via genJetIdx, i.e. the
existing patJetGenJetMatch used throughout scoutingToMiniAODDerivedCollections_
cff.py -- no re-derivation of deltaR matching here). slimmedGenJets (the
GenJet table's source) is genuinely pileup-free truth since genParticles is
built solely from the signal HepMCProduct, so this is the primary, rigorous
calibration driver -- unlike the per-candidate pvAssocQuality proxy used in
candidates.py, which is only an online approximation.
"""

import numpy as np
import awkward as ak

from . import io


def response_table(jets, genjets, nTrueInt):
    """Flat per-matched-jet table: response, genjet_pt, jet_pt, jet_eta, nTrueInt."""
    response, gj_pt, j_pt, j_eta, has_match = io.matched_response(jets, genjets)
    pu = io.broadcast_pu(nTrueInt, has_match) if nTrueInt is not None else None
    return {
        "response": np.asarray(response),
        "genjet_pt": np.asarray(gj_pt),
        "jet_pt": np.asarray(j_pt),
        "jet_eta": np.asarray(j_eta),
        "nTrueInt": np.asarray(pu) if pu is not None else None,
    }


def _binned_stats(x, y, bins):
    """Median and IQR-based resolution of y in bins of x."""
    idx = np.digitize(x, bins) - 1
    centers, medians, resolutions, counts = [], [], [], []
    for b in range(len(bins) - 1):
        sel = idx == b
        n = np.count_nonzero(sel)
        if n < 20:
            continue
        yy = y[sel]
        q16, q50, q84 = np.percentile(yy, [16, 50, 84])
        centers.append(0.5 * (bins[b] + bins[b + 1]))
        medians.append(q50)
        resolutions.append(0.5 * (q84 - q16))  # IQR/2, robust to non-Gaussian tails
        counts.append(n)
    return np.array(centers), np.array(medians), np.array(resolutions), np.array(counts)


def response_vs_pt(table, pt_bins=None, eta_max=2.5):
    """Median response + resolution binned in gen-jet pt, |eta|<eta_max only
    (scouting jets/taggers elsewhere in this workflow are central-only)."""
    if pt_bins is None:
        pt_bins = np.array([15, 20, 25, 30, 40, 50, 70, 100, 150, 200, 300, 500])
    sel = np.abs(table["jet_eta"]) < eta_max
    return _binned_stats(table["genjet_pt"][sel], table["response"][sel], pt_bins)


def response_vs_pt_data(table, pt_bins=None, eta_max=2.5):
    """Same as response_vs_pt, but for the data/data comparison table from
    data_matching.build_comparison_table() (scout_pt/scout_eta/off_pt fields,
    dR-matched scouting-vs-offline jet pairs in real collisions -- no gen
    truth involved). response = scout_pt / off_pt, binned in off_pt (the
    fully-calibrated offline jet standing in for genjet_pt's role above).
    """
    if pt_bins is None:
        pt_bins = np.array([15, 20, 25, 30, 40, 50, 70, 100, 150, 200, 300, 500])
    sel = np.abs(table["scout_eta"]) < eta_max
    response = table["scout_pt"][sel] / table["off_pt"][sel]
    return _binned_stats(table["off_pt"][sel], response, pt_bins)


def response_vs_pu(table, pu_bins=None, genjet_pt_min=30.0, eta_max=2.5):
    """Median response + resolution binned in true pileup (Pileup_nTrueInt),
    for jets above a fixed gen-jet pt threshold so the pt-dependence of
    response doesn't leak into the PU-dependence being studied.
    """
    if table["nTrueInt"] is None:
        raise ValueError("Pileup_nTrueInt not available -- check puTable is in the input file")
    if pu_bins is None:
        pu_bins = np.arange(0, 80, 5)
    sel = (np.abs(table["jet_eta"]) < eta_max) & (table["genjet_pt"] > genjet_pt_min)
    return _binned_stats(table["nTrueInt"][sel], table["response"][sel], pu_bins)


def mass_table(jets, genjets, nTrueInt):
    """Flat per-matched-jet table for mass resolution: mass_diff (jet.mass -
    genjet.mass), genjet_pt, genjet_mass, jet_eta, nTrueInt. Kept separate
    from response_table since the matched population's mass_diff isn't
    derivable from response_table's pt-ratio fields alone.
    """
    mass_diff, gj_pt, gj_mass, j_eta, has_match = io.matched_mass(jets, genjets)
    pu = io.broadcast_pu(nTrueInt, has_match) if nTrueInt is not None else None
    return {
        "mass_diff": np.asarray(mass_diff),
        "genjet_pt": np.asarray(gj_pt),
        "genjet_mass": np.asarray(gj_mass),
        "jet_eta": np.asarray(j_eta),
        "nTrueInt": np.asarray(pu) if pu is not None else None,
    }


def mass_resolution_vs_pt(table, pt_bins=None, eta_max=2.5):
    """Median mass difference (bias) + resolution (IQR/2 of jet.mass -
    genjet.mass) binned in gen-jet pt, |eta|<eta_max only.
    """
    if pt_bins is None:
        pt_bins = np.array([15, 20, 25, 30, 40, 50, 70, 100, 150, 200, 300, 500])
    sel = np.abs(table["jet_eta"]) < eta_max
    return _binned_stats(table["genjet_pt"][sel], table["mass_diff"][sel], pt_bins)


def mass_resolution_vs_pu(table, pu_bins=None, genjet_pt_min=30.0, eta_max=2.5):
    """Median mass difference + resolution binned in true pileup, for jets
    above a fixed gen-jet pt threshold.
    """
    if table["nTrueInt"] is None:
        raise ValueError("Pileup_nTrueInt not available -- check puTable is in the input file")
    if pu_bins is None:
        pu_bins = np.arange(0, 80, 5)
    sel = (np.abs(table["jet_eta"]) < eta_max) & (table["genjet_pt"] > genjet_pt_min)
    return _binned_stats(table["nTrueInt"][sel], table["mass_diff"][sel], pu_bins)


def fake_rate_vs_pu(jets, nTrueInt, pt_cut=30.0, eta_max=2.5, pu_bins=None):
    """Fraction of reconstructed jets (pt>pt_cut, |eta|<eta_max) with no
    genJetIdx match, vs Pileup_nTrueInt -- the pileup-jet/fake-rate signal
    that MinNeutralPt/MinNeutralPtSlope tuning is meant to control.
    """
    if pu_bins is None:
        pu_bins = np.arange(0, 80, 5)
    sel = (jets["pt"] > pt_cut) & (np.abs(jets["eta"]) < eta_max)
    is_fake = jets["genJetIdx"] < 0
    pu_per_jet = ak.broadcast_arrays(nTrueInt, sel)[0]

    pu_flat = np.asarray(ak.flatten(pu_per_jet[sel]))
    fake_flat = np.asarray(ak.flatten(is_fake[sel]))

    idx = np.digitize(pu_flat, pu_bins) - 1
    centers, rates, counts = [], [], []
    for b in range(len(pu_bins) - 1):
        m = idx == b
        n = np.count_nonzero(m)
        if n < 20:
            continue
        centers.append(0.5 * (pu_bins[b] + pu_bins[b + 1]))
        rates.append(np.mean(fake_flat[m]))
        counts.append(n)
    return np.array(centers), np.array(rates), np.array(counts)
