"""Pure-python PUPPI parameter-variant grid, kept free of any CMSSW/FWCore
import so the analysis scripts (which don't need or want a full CMSSW
environment just to read this dict) can import it cheaply.

scoutingPuppiCalibration_cff.py imports DEFAULT_VARIANTS from here and turns
each entry into a puppi.clone() + _applyVariantParams(...) call; see that
module for how each key is interpreted.
"""

DEFAULT_VARIANTS = {
    "nominal":        {},
    #"coneSmall":      {"central_cone": 0.2, "forward_cone": 0.2},
    #"coneLarge":      {"central_cone": 0.6, "forward_cone": 0.6},
    #"rmsPtMinLoose":  {"central_rmsPtMin": 0.2, "forward_rmsPtMin": 1.0},
    #"neutralPtLoose": {"fwd_MinNeutralPt": [1.0, 1.2], "fwd_MinNeutralPtSlope": [0.05, 0.05]},
    #"neutralPtTight": {"fwd_MinNeutralPt": [2.5, 3.0], "fwd_MinNeutralPtSlope": [0.12, 0.12]},
    # All variants above keep useVertexAssociation=True, i.e. PuppiProducer
    # trusts the scouting candidate's pre-computed 3-way vtxass quality flag
    # (PVUsedInFit/PU/unassociated) and never reads DeltaZCut/
    # NumOfPUVtxsForCharged/UseFromPV2Recovery -- see PuppiProducer.cc
    # (fUseVertexAssociation branch, ~L200-216) vs the dz-threshold branch
    # those PUPPI-v15-tune parameters actually live in (~L297+, only reached
    # when useVertexAssociation=False). "nominal" et al. therefore run v15's
    # cone/rmsPtMin/MinNeutralPt parameters but NOT v15's actual charged-
    # candidate categorization logic. vtxAssocOff switches that off so
    # PuppiProducer derives its own categorization from the packed
    # candidate's real dz()/fromPV(), which Run3ScoutingParticleToPacked
    # CandidateProducer already populates.
    #
    # useFromPVLooseTight=True is required alongside it, not optional: a
    # 300-event local smoke test with useVertexAssociation=False alone showed
    # charged PUPPI weights collapsing to ~1 for nearly everything (mean 0.95
    # vs nominal's 0.23) and jet multiplicity roughly doubling -- i.e. charged
    # pileup rejection was effectively disabled. Root cause, traced in
    # PuppiProducer.cc: the dz/pT gates that would normally catch PU-vertex-
    # associated charged candidates only engage at |eta|>=EtaMinUseDeltaZ=2.4
    # or pT>PtMaxCharged=20 GeV (Puppi_cff.py defaults) -- both essentially
    # never true for scouting's |eta|<2.5, mostly-soft constituents -- so
    # candidates fall through to the same default (id=0) as neutrals instead
    # of being hard-vetoed. UseFromPVLooseTight=True closes that gap: it maps
    # fromPV()==PVLoose (candidates associated to a PU vertex, i.e. vtxIdx>0
    # in Run3ScoutingParticleToPackedCandidateProducer.cc) to a hard PU veto
    # (id=2), and fromPV()==PVTight (unassociated candidates -- pvKey defaults
    # to 0 there, which is what makes them read as PVTight rather than PVLoose)
    # to a hard LV keep (id=1), independent of eta/pT. That's the closest
    # CHS-like binary categorization obtainable from the 3-tier vtxass flag
    # via the useVertexAssociation=False code path -- still coarser than
    # v15's real multi-tier |dz|<0.2/0.3cm thresholds (no ranked-PU-vertex
    # info exists upstream to support that), but no longer structurally
    # broken.
    #"vtxAssocOff":    {"useVertexAssociation": False, "useFromPVLooseTight": True},
    # vtxAssocImproved probes a *different* fix at a *different* layer than
    # vtxAssocOff: instead of changing which PuppiProducer branch reads the
    # existing (broken) vtxass quality flag, it changes what goes into that
    # flag in the first place. Traced in Run3ScoutingParticleToPackedCandidate
    # Producer.cc (useImprovedVertexAssociation flag, default off/unused by
    # every other variant and by the production CHS/plain jets): the raw
    # Run3ScoutingParticle::vertex() index that vtxass is built from
    # (HLTScoutingPFProducer.cc) is a 100-micron 3D position match between the
    # PF candidate's *stored* vertex point and the scouting vertex positions --
    # for charged candidates that stored point is the track's helix reference
    # point near the beamline (PFAlgo.cc), which essentially never coincides
    # with a reconstructed vertex position. Result, measured on the full v2
    # CRAB sample: 56% of ALL charged candidates (17.2M/30.7M) get dumped into
    # "unassociated" (NotReconstructedPrimary) by this accident, not because
    # they're actually incompatible with every vertex. useImprovedVertexAssoc
    # -iation redoes the LV/PU/unassociated split from a genuine per-vertex dz
    # scan (particle.dz()/dzsig(), which -unlike vertex()- really is a
    # trk->dz() against the leading vertex) using the same linear
    # dz(otherVertex)=dz(vertex0)-Delta_z approximation pat::PackedCandidate::
    # dz(ipv) itself relies on, mirroring CommonTools/RecoAlgos/
    # PrimaryVertexAssignment's dz+dzSig window test. Left at
    # useVertexAssociation's True default deliberately: the point of this
    # variant is to test whether PuppiProducer's simple, standard
    # useVertexAssociation=True branch (PuppiProducer.cc ~L200-216) works fine
    # for scouting once the vtxass flag it trusts is actually trustworthy, so
    # no useFromPVLooseTight workaround should be needed here.
    #"vtxAssocImproved": {"improvedVertexAssociation": True},
    # vtxAssocTrackMatched probes a refinement of the *same* fix
    # vtxAssocImproved applies, replacing its one remaining approximation.
    # vtxAssocImproved's nearest-vertex-in-dz scan still has to approximate
    # dz to any vertex other than 0 with a z-only linear shift, because
    # Run3ScoutingParticle only ever persists dz()/dzsig() against vertex 0 --
    # no per-vertex info survives to the converter for that value alone.
    # But Run3ScoutingParticleToPackedCandidateProducer.cc *already* consumes
    # a second collection ("tracks", reco::Track, built by
    # Run3ScoutingTrackToRecoTrackProducer from Run3ScoutingTrackCollection)
    # to embed track details (hasTrackDetails()/covariance/hit counts) via a
    # kinematic (eta/phi/pT) nearest-neighbor match -- confirmed traced to the
    # same underlying tracks (hltPFMuonMerging feeds both hltLightPFTracks,
    # source of the PF candidates behind Run3ScoutingParticle, and
    # hltScoutingTrackPacker, source of Run3ScoutingTrack, in the production
    # HLT menu). Those matched reco::Track objects are full, real tracks
    # (pat::makeRecoTrack sets reference point + momentum + covariance from
    # the persisted helix), so useTrackMatchedVertexAssociation reuses that
    # existing match (no new join) to call the track's own
    # dz(Point)/dxy(Point)/dzError() -- CMS's standard closest-approach
    # formula (DataFormats/TrackReco/TrackBase.h), honoring the vertex's full
    # (x,y,z) position -- for every candidate vertex in the same nearest-
    # vertex-in-dz scan vtxAssocImproved already does, instead of the z-only
    # shift. Falls back to vtxAssocImproved's linear-shift behavior whenever
    # no confident track match exists for a given candidate.
    "vtxAssocTrackMatched": {"improvedVertexAssociation": True, "trackMatchedVertexAssociation": True},
    # "optimized": the recalibrated central-region operating point from the
    # v5/v6 offline recalibration scans (see the scouting_puppi_recalibration
    # memory / ScoutingPuppiCalibration/Calibration/analysis/results/
    # v5_detailed_scan_summary.md for the full derivation): MedEtaSF 1.0->0.5
    # recovers ~15% relative response bias and ~4% relative resolution for
    # only ~+6% relative fake-rate cost (confirmed at full v5 statistics,
    # 6.19M events); MinNeutralPtSlope 0.015->0.01 compounds additively on
    # top of that with no interaction. RMSEtaSF/MinNeutralPt are left at PF
    # defaults (1.0/0.2) -- neither was a useful lever in the scans. This is
    # now the production default (referenceVariant in
    # customiseScoutingPuppiCalibrationNano) -- "nominal" (true PF/stock
    # defaults) and "vtxAssocTrackMatched" are kept alongside it so future
    # validation can still compare against the un-recalibrated baseline.
    "optimized": {"central_MedEtaSF": 0.5, "central_MinNeutralPtSlope": 0.01},
    # "chargedV15": ports PUPPI v15's (DP-2021/001) charged-particle
    # protections into PuppiProducer.cc's fUseVertexAssociation branch
    # (patched this session -- see the "v15-style protection" comments at
    # PuppiProducer.cc's PU-associated/unassociated cases), instead of
    # collapsing to id=0 (neutral) whenever scouting's coarser HLT
    # track-to-vertex matching produces no association at all -- measured on
    # the v7 DY sample: 65% of ALL charged candidates, still 13-16% above
    # 20-50 GeV, have pvAssocQuality=NotReconstructedPrimary and were getting
    # PUPPI-suppressed (mean weight 0.43-0.56) despite often being genuine LV
    # tracks, not pileup. Two knobs, both reused directly from v15's own
    # config surface (no new C++ parameters):
    #  - UseFromPV2Recovery/PtMinForFromPV2Recovery (True/4. GeV) are
    #    *already* inherited unconditionally from _stockPuppi (Puppi_cff.py)
    #    by every variant above once PuppiProducer.cc is rebuilt with the
    #    patch -- this is the soft pT floor for unassociated candidates.
    #  - ptMaxCharged=20 (-> PtMaxCharged) is NOT inherited (stock default is
    #    -1/disabled) and applied only here, to isolate its incremental
    #    effect on top of the FromPV2Recovery floor everyone else already
    #    gets: unconditionally keep any charged candidate above 20 GeV as LV
    #    regardless of what the (possibly wrong) vertex fit concluded, for
    #    both the PU-associated and fully-unassociated cases.
    # Combined with "optimized"'s MedEtaSF/slope recalibration since that's
    # the production default and the two fixes target independent failure
    # modes (neutral suppression vs. charged mis-categorization).
    "chargedV15": {"central_MedEtaSF": 0.5, "central_MinNeutralPtSlope": 0.01, "ptMaxCharged": 20.},
}
