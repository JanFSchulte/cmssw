#!/usr/bin/env python3
"""Validate puppi_recompute.recompute_weight against the stored
puppiWeight_nominal at full CRAB statistics, using default parameters
(medEtaSF=1, rmsEtaSF=1, minNeutralPt=0.2, minNeutralPtSlope=0.015 -- the
actual central-region production defaults, see Puppi_cff.py).

This closes the loop on the earlier ~78-94% match-rate caveat from the
50-event smoke test: that gap was mostly the missing MinNeutralPt/
MinNeutralPtSlope "Basic Cuts" floor (PuppiContainer.cc:319-320), which
applies to ANY id==0 candidate (not just true neutrals) and was not
modeled at all in the earlier quick check. useExp's external chi2 term was
also checked and ruled out (useExp defaults False and is never enabled here).
"""
import argparse
import glob
import os
import sys

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.puppi_recompute import recompute_weight


def _collect_files(patterns):
    files = []
    for p in patterns:
        matches = sorted(glob.glob(p))
        files.extend(matches if matches else [p])
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--prefix", default="ScoutingPuppiCalibCand")
    ap.add_argument("--variant", default="nominal")
    args = ap.parse_args()

    files = _collect_files(args.files)
    fields = ["pt", "eta", "charge", "pdgId", "pvAssocQuality", "puppiRawAlpha", "puppiAlphaMed", "puppiAlphaRms"]
    branches = ["%s_%s" % (args.prefix, f) for f in fields]
    branches.append("%s_puppiWeight_%s" % (args.prefix, args.variant))
    branches.append("Pileup_nPU")

    print("Loading %d file(s)..." % len(files))
    arrays = uproot.concatenate([f + ":Events" for f in files], filter_name=lambda n: n in branches, step_size="100 MB")

    pt = arrays["%s_pt" % args.prefix]
    eta = arrays["%s_eta" % args.prefix]
    charge = arrays["%s_charge" % args.prefix]
    pdgId = arrays["%s_pdgId" % args.prefix]
    pvq = arrays["%s_pvAssocQuality" % args.prefix]
    rawAlpha = arrays["%s_puppiRawAlpha" % args.prefix]
    alphaMed = arrays["%s_puppiAlphaMed" % args.prefix]
    alphaRms = arrays["%s_puppiAlphaRms" % args.prefix]
    stored = arrays["%s_puppiWeight_%s" % (args.prefix, args.variant)]
    nPU = arrays["Pileup_nPU"] + 1  # approximation, see puppi_recompute docstring

    nPU_per_cand = ak.broadcast_arrays(nPU, pt)[0]

    recomputed = recompute_weight(pt, eta, charge, pdgId, pvq, rawAlpha, alphaMed, alphaRms, nPU_per_cand,
                                   medEtaSF=1.0, rmsEtaSF=1.0, minNeutralPt=0.2, minNeutralPtSlope=0.015)

    stored_np = ak.to_numpy(ak.flatten(stored))
    recomputed_np = ak.to_numpy(ak.flatten(recomputed))
    id_np = ak.to_numpy(ak.flatten(ak.where(charge == 0, -1, pvq)))  # -1 marks neutral for reporting only

    diff = np.abs(stored_np - recomputed_np)
    n = len(stored_np)
    print("Total candidates: %d" % n)
    for tol in (1e-6, 1e-3, 1e-2, 0.05):
        match = np.sum(diff < tol)
        print("  match within %.g: %d / %d = %.4f%%" % (tol, match, n, 100.0 * match / n))

    # breakdown by id-like category (only id==0 candidates are non-trivial;
    # id==1/2 are hard overrides both sides and should match exactly)
    is_charged = ak.to_numpy(ak.flatten(charge)) != 0
    is_pv = is_charged & (id_np == 7)
    is_pu = is_charged & (id_np == 5)
    is_unassoc_or_neutral = ~is_pv & ~is_pu
    for name, mask in [("PV (id=1, hard=1)", is_pv), ("PU (id=2, hard=0)", is_pu),
                        ("unassoc+neutral (id=0, soft)", is_unassoc_or_neutral)]:
        d = diff[mask]
        if len(d) == 0:
            continue
        match = np.sum(d < 1e-3)
        print("  [%s] n=%d match(<1e-3)=%.4f%%  mean|diff|=%.5f  max|diff|=%.5f"
              % (name, len(d), 100.0 * match / len(d), d.mean(), d.max()))

    # show worst mismatches for id==0 to understand any residual gap
    mism = diff[is_unassoc_or_neutral]
    if len(mism) > 0:
        pct = np.percentile(mism, [50, 90, 99, 99.9])
        print("  id=0 |diff| percentiles [50,90,99,99.9]: %s" % pct)


if __name__ == "__main__":
    main()
