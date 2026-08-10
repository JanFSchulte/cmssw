"""Loading helpers for the scouting-PUPPI calibration NanoAOD files.

Reads the flat tables written by ScoutingPuppiCalibration/python/
scoutingPuppiCalibration_cff.py (and the pre-existing plain/CHS tables in
PhysicsTools.PatFromScouting.scoutingToMiniAODDerivedCollections_cff) with
uproot/awkward. No PyROOT/RDataFrame dependency.
"""

import awkward as ak
import numpy as np
import uproot

# branch-name prefix (NanoAOD table "name=") -> module label used for the
# jet collection, so callers can refer to either the plain/CHS baselines or
# a PUPPI variant by the same short key.
BASELINE_COLLECTIONS = {
    "plain": "ScoutingPFJetRecluster2",
    "chs": "ScoutingPFJetReclusterCHS",
}

# AK8: plain-only baseline (no CHS-AK8 collection exists), same PUPPI naming
# pattern as AK4 with "Fat" inserted -- see _addAK8Jets in
# scoutingPuppiCalibration_cff.py. Use with load_jets(..., genjet_idx_branch=
# "genJetAK8Idx") and load_genjets(..., prefix="GenJetAK8").
BASELINE_COLLECTIONS_AK8 = {
    "plain": "ScoutingFatPFJetRecluster2",
}


def puppi_collection_name(variant):
    return "ScoutingPFJetReclusterPUPPI_%s" % variant


def puppi_collection_name_ak8(variant):
    return "ScoutingFatPFJetReclusterPUPPI_%s" % variant


def _jet_branch_names(tree, prefix):
    """All branches for one jet table (kinematics table + its *_MCTable)."""
    names = [n for n in tree.keys() if n.startswith(prefix + "_")]
    return names


def load_jets(files, collection_name, step_size="100 MB", genjet_idx_branch=None):
    """Load one jet collection's kinematics + genJetIdx (from *_MCTable) into
    a single awkward record array, aligned by (event, jet) jaggedness.

    collection_name is the NanoAOD table "name=" (e.g. "ScoutingPFJetRecluster2"
    or puppi_collection_name("nominal")); the *_MCTable with the same name
    provides the gen-match index. AK4 collections use "<name>_genJetIdx"
    (jetMCTable's own field name, gated at pt>10 to match GenJetTable's own
    cut -- see PhysicsTools/NanoAOD/python/jetMC_cff.py); AK8 collections use
    "<name>_genJetAK8Idx" instead (fatJetMCTable, gated at pt>100 to match
    GenJetAK8Table's cut -- a DIFFERENT field name, not just a different
    collection, since the two use different pt gates baked into the Var
    expression itself). Pass genjet_idx_branch="genJetAK8Idx" for AK8
    collections; defaults to "genJetIdx" (AK4) for backward compatibility.
    """
    if genjet_idx_branch is None:
        genjet_idx_branch = "genJetIdx"
    fields = ["pt", "eta", "phi", "mass", "rawFactor", "charge", "nConstituents"]
    branches = ["%s_%s" % (collection_name, f) for f in fields]
    branches.append("%s_%s" % (collection_name, genjet_idx_branch))
    branches.append("Pileup_nTrueInt")

    arrays = uproot.concatenate(
        [f + ":Events" for f in files],
        filter_name=lambda n: n in branches,
        step_size=step_size,
    )
    out = {}
    for f in fields:
        key = "%s_%s" % (collection_name, f)
        if key in arrays.fields:
            out[f] = arrays[key]
    gj_key = "%s_%s" % (collection_name, genjet_idx_branch)
    out["genJetIdx"] = arrays[gj_key] if gj_key in arrays.fields else ak.Array([])
    pu_key = "Pileup_nTrueInt"
    nTrueInt = arrays[pu_key] if pu_key in arrays.fields else None
    return ak.Array(out), nTrueInt


def load_event_keys(file_or_url, branches=("run", "luminosityBlock", "event")):
    """Read event-id branches from a single file (local path or xrootd URL),
    in on-disk order, so the returned arrays' positions double as local entry
    indices. Used by data_matching.py to join scouting and offline files by
    (run, luminosityBlock, event) without touching any jet branches first.
    """
    tree = uproot.open(file_or_url)["Events"]
    arrays = tree.arrays(list(branches), library="np")
    return {b: arrays[b] for b in branches}


def load_genjets(files, step_size="100 MB", prefix="GenJet"):
    """Standard NanoAOD GenJet table (pt/eta/phi), the truth reference.
    Pass prefix="GenJetAK8" for the AK8 gen-truth collection.
    """
    fields = ["pt", "eta", "phi", "mass"]
    branches = ["%s_%s" % (prefix, f) for f in fields]
    arrays = uproot.concatenate(
        [f + ":Events" for f in files],
        filter_name=lambda n: n in branches,
        step_size=step_size,
    )
    return ak.Array({f: arrays["%s_%s" % (prefix, f)] for f in fields})


def load_candidates(files, variant_names, step_size="100 MB", prefix="ScoutingPuppiCalibCand"):
    """ScoutingPuppiCalibCand table: kinematics + pvAssocQuality + every
    variant's puppiWeight, as sibling columns keyed to the same candidate
    index. pvAssocQuality is the caveated ONLINE HLT-vertex proxy, not truth.

    prefix="ScoutingPuppiCalibImprovedVtxAssocCand" reads the separate table
    for variants with improvedVertexAssociation=True (see variants.py /
    _addImprovedCandidateTable in scoutingPuppiCalibration_cff.py) -- a
    different product instance, so its own pvAssocQuality/puppiWeight can't
    live as sibling columns on the main table.
    """
    fields = ["pt", "eta", "phi", "mass", "charge", "pdgId", "pvAssocQuality",
              "puppiRawAlpha", "puppiAlphaMed", "puppiAlphaRms"]
    branches = ["%s_%s" % (prefix, f) for f in fields]
    branches += ["%s_puppiWeight_%s" % (prefix, v) for v in variant_names]

    arrays = uproot.concatenate(
        [f + ":Events" for f in files],
        filter_name=lambda n: n in branches,
        step_size=step_size,
    )
    out = {}
    for f in fields:
        key = "%s_%s" % (prefix, f)
        if key in arrays.fields:
            out[f] = arrays[key]
    for v in variant_names:
        key = "%s_puppiWeight_%s" % (prefix, v)
        if key in arrays.fields:
            out["puppiWeight_%s" % v] = arrays[key]
    return ak.Array(out)


def matched_response(jets, genjets, genjet_pt_min=10.0):
    """response = jet.pt / matched_genjet.pt for jets with a valid genJetIdx
    (jetMCTable's own convention: genJetIdx>=0 already implies the matched
    gen jet passed pt>10, see PhysicsTools/NanoAOD/python/jetMC_cff.py).

    Returns (response, genjet_pt, jet_pt, jet_eta) flattened to 1D awkward
    arrays over all matched jets in all events, plus a companion nTrueInt
    array broadcast per matched jet if provided by the caller separately.
    """
    has_match = jets["genJetIdx"] >= 0
    idx = jets["genJetIdx"][has_match]
    gj_pt = genjets["pt"][idx]
    j_pt = jets["pt"][has_match]
    j_eta = jets["eta"][has_match]
    response = j_pt / gj_pt
    return ak.flatten(response), ak.flatten(gj_pt), ak.flatten(j_pt), ak.flatten(j_eta), has_match


def broadcast_pu(nTrueInt, has_match):
    """Broadcast the per-event Pileup_nTrueInt to the flattened per-matched-
    jet arrays produced by matched_response, using the same mask.
    """
    pu_per_jet = ak.broadcast_arrays(nTrueInt, has_match)[0]
    return ak.flatten(pu_per_jet[has_match])


def matched_mass(jets, genjets, genjet_pt_min=10.0):
    """mass_diff = jet.mass - matched_genjet.mass for the same matched-jet
    population matched_response uses (genJetIdx>=0). Returns (mass_diff,
    genjet_pt, genjet_mass, jet_eta) flattened to 1D, plus the has_match mask
    (for broadcast_pu reuse) -- a difference, not a ratio, since light-jet
    genjet mass runs close to zero and a ratio metric blows up in that tail.
    """
    has_match = jets["genJetIdx"] >= 0
    idx = jets["genJetIdx"][has_match]
    gj_pt = genjets["pt"][idx]
    gj_mass = genjets["mass"][idx]
    j_mass = jets["mass"][has_match]
    j_eta = jets["eta"][has_match]
    mass_diff = j_mass - gj_mass
    return (ak.flatten(mass_diff), ak.flatten(gj_pt), ak.flatten(gj_mass),
            ak.flatten(j_eta), has_match)
