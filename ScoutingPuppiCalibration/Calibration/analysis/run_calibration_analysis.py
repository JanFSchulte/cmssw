#!/usr/bin/env python3
"""CLI entry point for the scouting-PUPPI calibration analysis.

Reads one or more NanoAOD-style ROOT files produced by
ScoutingPuppiCalibration/test/scoutingPuppiCalib_mc_cfg.py, computes jet
response/resolution/fake-rate for the plain/CHS baselines and every PUPPI
variant, a per-candidate PUPPI-weight sanity check (labeled as an online-
proxy check, not truth), and writes plots + a summary table to an output
directory.

Usage:
    python run_calibration_analysis.py --files calib_smoke.root --outdir plots/
    python run_calibration_analysis.py --files "crab_out/*.root" --outdir plots/
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis import io, response as response_mod, candidates, summary
# variants.py is deliberately free of any CMSSW/FWCore import (see that
# file), so the analysis side can read the variant names without needing a
# full CMSSW environment.
from python.variants import DEFAULT_VARIANTS


def _collect_files(patterns):
    files = []
    for p in patterns:
        matches = sorted(glob.glob(p))
        files.extend(matches if matches else [p])
    if not files:
        raise SystemExit("No input files found for patterns: %s" % patterns)
    return files


def _plot_response_vs_pt(results, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for name, table in results.items():
        centers, medians, res, counts = response_mod.response_vs_pt(table)
        if len(centers) == 0:
            continue
        axes[0].plot(centers, medians, marker="o", label=name)
        axes[1].plot(centers, res, marker="o", label=name)
    axes[0].axhline(1.0, color="gray", linestyle="--", linewidth=1)
    axes[0].set_xlabel("gen jet pt [GeV]")
    axes[0].set_ylabel("median response (reco pt / gen pt)")
    axes[0].set_title("Jet response vs pt  (|eta|<2.5)")
    axes[1].set_xlabel("gen jet pt [GeV]")
    axes[1].set_ylabel("resolution (IQR/2 of response)")
    axes[1].set_title("Jet resolution vs pt  (|eta|<2.5)")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "response_resolution_vs_pt.png"), dpi=150)
    plt.close(fig)


def _plot_response_vs_pu(results, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    any_pu = False
    for name, table in results.items():
        if table["nTrueInt"] is None:
            continue
        any_pu = True
        centers, medians, res, counts = response_mod.response_vs_pu(table)
        if len(centers) == 0:
            continue
        axes[0].plot(centers, medians, marker="o", label=name)
        axes[1].plot(centers, res, marker="o", label=name)
    if not any_pu:
        plt.close(fig)
        print("Pileup_nTrueInt not found in input files -- skipping response-vs-PU plot "
              "(check puTable is present in the ROOT file)")
        return
    axes[0].axhline(1.0, color="gray", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Pileup_nTrueInt")
    axes[0].set_ylabel("median response")
    axes[0].set_title("Jet response vs true pileup (gen pt>30 GeV, |eta|<2.5)")
    axes[1].set_xlabel("Pileup_nTrueInt")
    axes[1].set_ylabel("resolution (IQR/2)")
    axes[1].set_title("Jet resolution vs true pileup")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "response_resolution_vs_pu.png"), dpi=150)
    plt.close(fig)


def _plot_fake_rate(jets_by_name, nTrueInt_by_name, outdir):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    any_curve = False
    for name, jets in jets_by_name.items():
        nTrueInt = nTrueInt_by_name.get(name)
        if nTrueInt is None:
            continue
        centers, rates, counts = response_mod.fake_rate_vs_pu(jets, nTrueInt)
        if len(centers) == 0:
            continue
        any_curve = True
        ax.plot(centers, rates, marker="o", label=name)
    if not any_curve:
        plt.close(fig)
        return
    ax.set_xlabel("Pileup_nTrueInt")
    ax.set_ylabel("fraction of jets (pt>30 GeV) with no gen-jet match")
    ax.set_title("Fake-rate proxy vs true pileup")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fake_rate_vs_pu.png"), dpi=150)
    plt.close(fig)


def _plot_candidate_weights(cand_table, variants, outdir):
    for variant in variants:
        try:
            by_eta = candidates.weight_by_pvassoc(cand_table, variant)
        except KeyError:
            continue
        fig, axes = plt.subplots(1, len(by_eta), figsize=(5.5 * len(by_eta), 4.2), squeeze=False)
        for ax, ((lo, hi), cats) in zip(axes[0], by_eta.items()):
            for label, arr in cats.items():
                if len(arr) == 0:
                    continue
                ax.hist(arr, bins=40, range=(0, 1.2), histtype="step", label=label, density=True)
            ax.set_xlabel("PUPPI weight")
            ax.set_title("|eta| in [%.1f, %.1f)" % (lo, hi))
            ax.legend(fontsize=7)
        fig.suptitle(
            "variant=%s -- split by pvAssocQuality, an ONLINE HLT-vertex proxy, NOT MC truth" % variant,
            fontsize=9,
        )
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "puppi_weight_%s.png" % variant), dpi=150)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", nargs="+", required=True,
                         help="ROOT file path(s) or glob pattern(s)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--variants", nargs="*", default=list(DEFAULT_VARIANTS.keys()),
                         help="PUPPI variant names to include (default: all of DEFAULT_VARIANTS)")
    parser.add_argument("--cand-files", nargs="*", default=None,
                         help="Optional smaller file subset for the per-candidate weight sanity "
                              "check (candidates.py), which is O(300) rows/event vs O(1) for jets "
                              "and dominates memory/time at full dataset scale. Defaults to --files.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    files = _collect_files(args.files)
    cand_files = _collect_files(args.cand_files) if args.cand_files else files
    print("Loading %d file(s): %s" % (len(files), files))

    genjets = io.load_genjets(files)

    collections = dict(io.BASELINE_COLLECTIONS)
    for v in args.variants:
        collections[v] = io.puppi_collection_name(v)

    jets_by_name, nTrueInt_by_name, table_by_name = {}, {}, {}
    for name, coll in collections.items():
        jets, nTrueInt = io.load_jets(files, coll)
        jets_by_name[name] = jets
        nTrueInt_by_name[name] = nTrueInt
        table_by_name[name] = response_mod.response_table(jets, genjets, nTrueInt)
        print("  %-16s: %d jets loaded" % (name, len(jets)))

    _plot_response_vs_pt(table_by_name, args.outdir)
    _plot_response_vs_pu(table_by_name, args.outdir)
    _plot_fake_rate(jets_by_name, nTrueInt_by_name, args.outdir)

    # Variants with improvedVertexAssociation=True and/or
    # trackMatchedVertexAssociation=True (see variants.py) get their own
    # packedPFCandidates instance/table (different product, so their
    # puppiWeight column doesn't exist on the main cand_table below -- see
    # _addImprovedCandidateTable/_addTrackMatchedCandidateTable in
    # scoutingPuppiCalibration_cff.py). Must be excluded from the main
    # table's variant list, not just handled separately:
    # summarize_weight_distributions (unlike _plot_candidate_weights)
    # doesn't catch a missing puppiWeight_<variant> column and will raise
    # KeyError. trackMatchedVertexAssociation variants also set
    # improvedVertexAssociation=True (see variants.py), so that check must
    # come first, or they'd be routed to the wrong (ImprovedVtxAssoc) table.
    track_matched_variants = [v for v in args.variants
                               if DEFAULT_VARIANTS.get(v, {}).get("trackMatchedVertexAssociation")]
    improved_variants = [v for v in args.variants
                          if DEFAULT_VARIANTS.get(v, {}).get("improvedVertexAssociation")
                          and v not in track_matched_variants]
    main_variants = [v for v in args.variants if v not in improved_variants and v not in track_matched_variants]

    print("Loading candidate table from %d file(s): %s" % (len(cand_files), cand_files))
    cand_table = io.load_candidates(cand_files, main_variants)
    _plot_candidate_weights(cand_table, main_variants, args.outdir)
    summary_text = candidates.summarize_weight_distributions(cand_table, main_variants)

    if improved_variants:
        improved_cand_table = io.load_candidates(
            cand_files, improved_variants, prefix="ScoutingPuppiCalibImprovedVtxAssocCand")
        _plot_candidate_weights(improved_cand_table, improved_variants, args.outdir)
        summary_text += "\n" + candidates.summarize_weight_distributions(improved_cand_table, improved_variants)

    if track_matched_variants:
        track_matched_cand_table = io.load_candidates(
            cand_files, track_matched_variants, prefix="ScoutingPuppiCalibTrackMatchedVtxAssocCand")
        _plot_candidate_weights(track_matched_cand_table, track_matched_variants, args.outdir)
        summary_text += "\n" + candidates.summarize_weight_distributions(track_matched_cand_table, track_matched_variants)

    with open(os.path.join(args.outdir, "candidate_weight_summary.txt"), "w") as f:
        f.write(summary_text)

    rows = summary.summarize_variant(
        {v: table_by_name[v] for v in args.variants},
        {v: jets_by_name[v] for v in args.variants},
        {v: nTrueInt_by_name[v] for v in args.variants},
    )
    md = summary.format_markdown_table(rows)
    print("\n" + md)
    with open(os.path.join(args.outdir, "summary.md"), "w") as f:
        f.write(md + "\n")
    summary.write_csv(rows, os.path.join(args.outdir, "summary.csv"))

    print("\nWrote plots and summary to %s" % args.outdir)


if __name__ == "__main__":
    main()
