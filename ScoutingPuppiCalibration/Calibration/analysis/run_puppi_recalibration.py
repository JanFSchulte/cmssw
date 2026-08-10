#!/usr/bin/env python3
"""CLI entry point for the offline PUPPI parameter recalibration scan.

Recomputes candidate weights from the raw alpha diagnostics for one or more
hypothetical (MedEtaSF, RMSEtaSF, MinNeutralPt, MinNeutralPtSlope) tuples,
reclusters PUPPI jets with genuine FastJet anti-kt(0.4), matches to gen jets,
and reports response bias / resolution / fake rate -- the same metrics
summary.py already uses for the production PUPPI variants -- so a
hypothetical parameter choice's jet-level impact can be evaluated without a
new CRAB production pass.

Usage:
    python run_puppi_recalibration.py --files calib_smoke.root --sanity-check
    python run_puppi_recalibration.py --files "f1.root f2.root ..." \
        --medEtaSF 0.7 0.8 0.9 1.0 1.1 --rmsEtaSF 1.0
"""
import argparse
import glob
import itertools
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import io, puppi_refit, response as response_mod, summary


def _collect_files(patterns):
    files = []
    for p in patterns:
        matches = sorted(glob.glob(p))
        files.extend(matches if matches else [p])
    if not files:
        raise SystemExit("No input files found for patterns: %s" % patterns)
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--prefix", default="ScoutingPuppiCalibCand")
    ap.add_argument("--medEtaSF", nargs="*", type=float, default=[1.0])
    ap.add_argument("--rmsEtaSF", nargs="*", type=float, default=[1.0])
    ap.add_argument("--minNeutralPt", nargs="*", type=float, default=[0.2])
    ap.add_argument("--minNeutralPtSlope", nargs="*", type=float, default=[0.015])
    ap.add_argument("--sanity-check", action="store_true",
                     help="also recluster with the CURRENT production defaults and compare "
                          "against the actual stored nominal PUPPI jet collection, as an "
                          "end-to-end validation of the reclustering+matching pipeline itself.")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    files = _collect_files(args.files)
    print("Loading %d file(s)..." % len(files))
    t0 = time.time()
    cand, nPU, rho = puppi_refit.load_refit_inputs(files, prefix=args.prefix)
    genjets = io.load_genjets(files)
    print("  loaded in %.1fs, %d events" % (time.time() - t0, len(cand["pt"])))

    if args.sanity_check:
        print("\n=== Sanity check: reclustering pipeline vs actual production nominal jets ===")
        jets_nom, nTrueInt_nom = io.load_jets(files, io.puppi_collection_name("nominal"))
        table_nom = response_mod.response_table(jets_nom, genjets, nTrueInt_nom)
        rows_nom = summary.summarize_variant({"nominal (actual)": table_nom}, {"nominal (actual)": jets_nom},
                                              {"nominal (actual)": nTrueInt_nom})

        t0 = time.time()
        jets_refit = puppi_refit.refit_and_recluster(cand, nPU, genjets, rho=rho,
                                                       medEtaSF=1.0, rmsEtaSF=1.0,
                                                       minNeutralPt=0.2, minNeutralPtSlope=0.015)
        print("  reclustering took %.1fs" % (time.time() - t0))
        table_refit = response_mod.response_table(jets_refit, genjets, nTrueInt_nom)
        rows_refit = summary.summarize_variant({"nominal (refit, default params)": table_refit},
                                                 {"nominal (refit, default params)": jets_refit},
                                                 {"nominal (refit, default params)": nTrueInt_nom})
        print(summary.format_markdown_table(rows_nom + rows_refit))
        print()

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)

    print("\n=== Parameter scan ===")
    combos = list(itertools.product(args.medEtaSF, args.rmsEtaSF, args.minNeutralPt, args.minNeutralPtSlope))
    print("%d parameter combination(s)" % len(combos))
    _, nTrueInt_ref = io.load_jets(files, io.puppi_collection_name("nominal"))

    rows = []
    for (medSF, rmsSF, minPt, minPtSlope) in combos:
        name = "med%.2f_rms%.2f_minPt%.2f_slope%.3f" % (medSF, rmsSF, minPt, minPtSlope)
        t0 = time.time()
        jets = puppi_refit.refit_and_recluster(cand, nPU, genjets, rho=rho, medEtaSF=medSF, rmsEtaSF=rmsSF,
                                                minNeutralPt=minPt, minNeutralPtSlope=minPtSlope)
        table = response_mod.response_table(jets, genjets, nTrueInt_ref)
        row = summary.summarize_variant({name: table}, {name: jets}, {name: nTrueInt_ref})
        rows += row
        print("  %-40s bias=%+.4f res=%.4f fake=%.4f  (%.1fs)"
              % (name, row[0]["response_bias"], row[0]["resolution"], row[0]["fake_rate"], time.time() - t0))

    print()
    print(summary.format_markdown_table(rows))
    if args.outdir:
        summary.write_csv(rows, os.path.join(args.outdir, "recalibration_scan.csv"))
        with open(os.path.join(args.outdir, "recalibration_scan.md"), "w") as f:
            f.write(summary.format_markdown_table(rows) + "\n")
        print("\nWrote scan results to %s" % args.outdir)


if __name__ == "__main__":
    main()
