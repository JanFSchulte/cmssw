"""Match scouting NanoAOD events to existing *centrally-produced* offline
NanoAOD (JetMET0/1, Muon0/1, ...) by (run, luminosityBlock, event), then
dR-match the two jet collections within each matched event.

Scouting and offline reconstruction of data live in physically separate
primary datasets (unlike MC MiniAOD, where both are in one file), so there is
no in-file way to compare scouting jets to offline jets for data. This module
does NOT produce a companion offline sample -- central Run3 NanoAOD already
exists for the offline PDs and covers the same runs as the scouting stream
(verified directly: run 380945, the example run used throughout
ScoutingNanoProduction's standalone configs, is present in
/JetMET0/Run2024D-MINIv6NANOv15-v1/NANOAOD, and streaming that dataset's
run/luminosityBlock/event over xrootd found exact matches against an
already-produced local scouting NanoAOD file for the same run). So the whole
job here is: find the existing files, join by event id, dR-match jets.

Three stages, sized to avoid pulling full jet content for every event in a PD
that's far larger than what's needed:
  1. scouting_event_keys()   -- cheap run/lumi/event + provenance, scouting side
  2. find_offline_files() + match_offline_files() -- DAS-restricted offline
     key scan, joined against the scouting keys
  3. build_comparison_table() -- targeted jet extraction for matched events
     only, then dR matching within each event

Trigger note: DST_PFScouting_JetHT-triggered scouting events should be
matched against JetMET0+JetMET1 (both siblings of the trigger-hash-split PD,
since it isn't confirmed up front which half a given event lands in);
DST_PFScouting_SingleMuon-triggered events against Muon0+Muon1.
"""

import subprocess

import numpy as np
import pandas as pd
import uproot

from . import io

DEFAULT_REDIRECTOR = "root://cms-xrd-global.cern.ch/"


def das_query(query):
    """Same tiny dasgoclient wrapper as ScoutingNanoProduction/submit_scoutingNano.py's
    das_query -- duplicated here (not imported cross-package) to keep this
    module runnable standalone, outside the CRAB submission machinery.
    """
    cmd = ["dasgoclient", "-query", query]
    out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode().strip()
    return [line for line in out.splitlines() if line.strip()]


def scouting_event_keys(files):
    """Stage 1: run/luminosityBlock/event + (file, local entry) provenance
    for every scouting NanoAOD file, as one pandas DataFrame with columns
    run, luminosityBlock, event, scout_file, scout_entry.
    """
    frames = []
    for f in files:
        keys = io.load_event_keys(f)
        n = len(keys["run"])
        frames.append(pd.DataFrame({
            "run": keys["run"],
            "luminosityBlock": keys["luminosityBlock"],
            "event": keys["event"],
            "scout_file": f,
            "scout_entry": np.arange(n),
        }))
    if not frames:
        return pd.DataFrame(columns=["run", "luminosityBlock", "event", "scout_file", "scout_entry"])
    return pd.concat(frames, ignore_index=True)


def find_offline_files(dataset_names, runs):
    """Stage 2a: DAS file discovery for the offline PD(s), restricted to the
    given runs. dataset_names is a list of full DAS dataset names (pass both
    siblings of a trigger-hash-split PD -- see module docstring). Central
    NanoAOD files are NOT split per-run (a single file spans a range of
    runs), so this only narrows the *candidate* file list; the run/lumi/event
    join in match_offline_files() below still does the real filtering.
    Returns a deduplicated, sorted list of LFNs (no redirector prefix).
    """
    files = set()
    for dataset in dataset_names:
        for run in runs:
            files.update(das_query("file dataset=%s run=%d" % (dataset, int(run))))
    return sorted(files)


def match_offline_files(scout_df, offline_files, redirector=DEFAULT_REDIRECTOR):
    """Stage 2b: open each candidate offline file, mask to the runs present
    in scout_df, and inner-join on (run, luminosityBlock, event). Returns a
    concatenated DataFrame of matched provenance: run, luminosityBlock,
    event, scout_file, scout_entry, off_file, off_entry (off_file is the
    LFN, without redirector prefix; off_entry is the entry's position in the
    full file, matching how build_comparison_table() re-opens it).
    """
    if len(scout_df) == 0:
        return pd.DataFrame(columns=["run", "luminosityBlock", "event",
                                      "scout_file", "scout_entry", "off_file", "off_entry"])

    runs_of_interest = set(int(r) for r in scout_df["run"].unique())
    matched_frames = []
    for lfn in offline_files:
        url = lfn if lfn.startswith("root:") else redirector + lfn
        try:
            keys = io.load_event_keys(url)
        except Exception as exc:
            print("  WARNING: failed to open %s: %s" % (url, exc))
            continue
        mask = np.isin(keys["run"], list(runs_of_interest))
        if not mask.any():
            continue
        off_df = pd.DataFrame({
            "run": keys["run"][mask],
            "luminosityBlock": keys["luminosityBlock"][mask],
            "event": keys["event"][mask],
            "off_file": lfn,
            "off_entry": np.nonzero(mask)[0],
        })
        merged = scout_df.merge(off_df, on=["run", "luminosityBlock", "event"], how="inner")
        if len(merged):
            matched_frames.append(merged)

    if not matched_frames:
        return pd.DataFrame(columns=["run", "luminosityBlock", "event",
                                      "scout_file", "scout_entry", "off_file", "off_entry"])
    return pd.concat(matched_frames, ignore_index=True)


def _delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    return (dphi + np.pi) % (2 * np.pi) - np.pi


def dr_match_jets(eta_a, phi_a, eta_b, phi_b, max_dr=0.4):
    """Greedy nearest-neighbor dR matching between two jet collections from
    the SAME event (1D arrays, one entry per jet). Closest pairs are
    assigned first and each jet is used in at most one pair -- the standard
    approach for matching two independently-clustered jet collections (no
    pre-existing dR matcher exists anywhere in this repo: response.py's
    matching relies entirely on MC's pre-computed genJetIdx, which has no
    cross-dataset equivalent here). Returns a list of (idx_a, idx_b, dr)
    tuples, sorted by ascending dr.
    """
    eta_a = np.asarray(eta_a, dtype=float)
    phi_a = np.asarray(phi_a, dtype=float)
    eta_b = np.asarray(eta_b, dtype=float)
    phi_b = np.asarray(phi_b, dtype=float)
    na, nb = len(eta_a), len(eta_b)
    if na == 0 or nb == 0:
        return []

    deta = eta_a[:, None] - eta_b[None, :]
    dphi = _delta_phi(phi_a[:, None], phi_b[None, :])
    dr = np.sqrt(deta ** 2 + dphi ** 2)

    ia, ib = np.nonzero(dr < max_dr)
    candidates = sorted(zip(dr[ia, ib], ia.tolist(), ib.tolist()))

    used_a, used_b = set(), set()
    matches = []
    for d, i, j in candidates:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        matches.append((i, j, float(d)))
    return matches


_SCOUT_FIELDS = ["pt", "eta", "phi", "mass", "rawFactor"]
_OFF_FIELDS = ["pt", "eta", "phi", "mass", "rawFactor", "area"]


def build_comparison_table(matches_df, scout_collection="ScoutingPFJetRecluster2",
                            redirector=DEFAULT_REDIRECTOR, max_dr=0.4):
    """Stage 3: for the event pairs found by match_offline_files(), pull just
    the matched entries' jet branches (scouting: scout_collection, e.g.
    "ScoutingPFJetRecluster2" or "ScoutingPFJetReclusterCHS" -- see
    io.BASELINE_COLLECTIONS; offline: standard NanoAOD "Jet"), dR-match the
    two jet collections within each event, and return one flat table (dict
    of numpy arrays, one row per matched jet pair) shaped like response.py's
    response_table() output so it plugs into the same _binned_stats-based
    plotting.
    """
    out = {"run": [], "luminosityBlock": [], "event": [], "dr": []}
    for f in _SCOUT_FIELDS:
        out["scout_%s" % f] = []
    for f in _OFF_FIELDS:
        out["off_%s" % f] = []

    if len(matches_df) == 0:
        return {k: np.asarray(v) for k, v in out.items()}

    scout_cache = {}
    off_cache = {}

    for (scout_file, off_file), group in matches_df.groupby(["scout_file", "off_file"]):
        if scout_file not in scout_cache:
            branches = ["%s_%s" % (scout_collection, f) for f in _SCOUT_FIELDS]
            scout_cache[scout_file] = uproot.open(scout_file)["Events"].arrays(branches, library="ak")
        if off_file not in off_cache:
            url = off_file if off_file.startswith("root:") else redirector + off_file
            branches = ["Jet_%s" % f for f in _OFF_FIELDS]
            off_cache[off_file] = uproot.open(url)["Events"].arrays(branches, library="ak")

        s_arr = scout_cache[scout_file]
        o_arr = off_cache[off_file]

        for row in group.itertuples():
            s_jets = {f: np.asarray(s_arr["%s_%s" % (scout_collection, f)][row.scout_entry])
                      for f in _SCOUT_FIELDS}
            o_jets = {f: np.asarray(o_arr["Jet_%s" % f][row.off_entry])
                      for f in _OFF_FIELDS}

            pairs = dr_match_jets(s_jets["eta"], s_jets["phi"], o_jets["eta"], o_jets["phi"], max_dr=max_dr)
            for i, j, d in pairs:
                out["run"].append(row.run)
                out["luminosityBlock"].append(row.luminosityBlock)
                out["event"].append(row.event)
                out["dr"].append(d)
                for f in _SCOUT_FIELDS:
                    out["scout_%s" % f].append(s_jets[f][i])
                for f in _OFF_FIELDS:
                    out["off_%s" % f].append(o_jets[f][j])

    return {k: np.asarray(v) for k, v in out.items()}
