import FWCore.ParameterSet.Config as cms

from ..modules.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidatesLowPt_cfi import *
from ..modules.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracksLowPt_cfi import *
from ..modules.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracksLowPt_cfi import *
from ..modules.hltIter0Phase2L3FromL1TkMuonTrackCutClassifierLowPt_cfi import *
from ..modules.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurityLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidatesLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonClustersRefRemovalLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracksLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEventLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonMergedLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonPixelClusterCheckLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonPixelHitDoubletsLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonPixelHitTripletsLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonPixelLayerTripletsLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonPixelSeedsLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonTrackCutClassifierLowPt_cfi import *
from ..modules.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurityLowPt_cfi import *
from ..modules.hltPhase2L3GlbMuonLowPt_cfi import *
from ..modules.hltL1MuonsPt0_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelLayerQuadrupletsLowPt_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksLowPt_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksHitDoubletsLowPt_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksHitQuadrupletsLowPt_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegionsLowPt_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelVerticesLowPt_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonTrimmedPixelVerticesLowPt_cfi import *

from ..modules.hltPhase2L3MuonMergedLowPt_cfi import *
from ..modules.hltPhase2L3MuonPixelTracksFilterLowPt_cfi import *
from ..modules.hltPhase2L3MuonPixelTracksFitterLowPt_cfi import *
from ..modules.hltPhase2L3MuonsLowPt_cfi import *
from ..modules.hltPhase2L3MuonsNoIDLowPt_cfi import *
from ..modules.hltPhase2L3MuonTracks_cfi import *
from ..modules.hltL3MuonsPhase2L3LinksLowPt_cfi import *
from ..modules.hltPhase2L3OIL3MuonCandidatesLowPt_cfi import *
from ..modules.hltPhase2L3OIL3MuonsLowPt_cfi import *
from ..modules.hltPhase2L3OIL3MuonsLinksCombinationLowPt_cfi import *
from ..modules.hltPhase2L3OIMuCtfWithMaterialTracksLowPt_cfi import *
from ..modules.hltPhase2L3OIMuonTrackCutClassifierLowPt_cfi import *
from ..modules.hltPhase2L3OIMuonTrackSelectionHighPurityLowPt_cfi import *
from ..modules.hltPhase2L3OISeedsFromL2MuonsLowPt_cfi import *
from ..modules.hltPhase2L3OITrackCandidatesLowPt_cfi import *
from ..modules.hltL3MuonsPhase2L3OILowPt_cfi import *

HLTL3MuonLowPtTask = cms.Task(
    hltPhase2L3OISeedsFromL2MuonsLowPt,
    hltPhase2L3OITrackCandidatesLowPt,
    hltPhase2L3OIMuCtfWithMaterialTracksLowPt,
    hltPhase2L3OIMuonTrackCutClassifierLowPt,
    hltPhase2L3OIMuonTrackSelectionHighPurityLowPt,
    hltPhase2L3OIL3MuonsLowPt,
    hltL3MuonsPhase2L3OILowPt,
    hltPhase2L3OIL3MuonsLinksCombinationLowPt,
    hltPhase2L3OIL3MuonCandidatesLowPt,
    hltL1MuonsPt0,
    hltPhase2L3FromL1TkMuonPixelLayerQuadrupletsLowPt,
    hltPhase2L3FromL1TkMuonPixelTracksTrackingRegionsLowPt,
    hltPhase2L3FromL1TkMuonPixelTracksHitDoubletsLowPt,
    hltPhase2L3FromL1TkMuonPixelTracksHitQuadrupletsLowPt,
    hltPhase2L3FromL1TkMuonPixelTracksLowPt,
    hltPhase2L3FromL1TkMuonPixelVerticesLowPt,
    hltPhase2L3FromL1TkMuonTrimmedPixelVerticesLowPt,
    hltPhase2L3MuonPixelTracksFilterLowPt,
    hltPhase2L3MuonPixelTracksFitterLowPt,
    hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracksLowPt,
    hltIter0Phase2L3FromL1TkMuonCkfTrackCandidatesLowPt,
    hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracksLowPt,
    hltIter0Phase2L3FromL1TkMuonTrackCutClassifierLowPt,
    hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurityLowPt,
    hltIter2Phase2L3FromL1TkMuonClustersRefRemovalLowPt,
    hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEventLowPt,
    hltIter2Phase2L3FromL1TkMuonPixelLayerTripletsLowPt,
    hltIter2Phase2L3FromL1TkMuonPixelHitDoubletsLowPt,
    hltIter2Phase2L3FromL1TkMuonPixelHitTripletsLowPt,
    hltIter2Phase2L3FromL1TkMuonPixelClusterCheckLowPt,
    hltIter2Phase2L3FromL1TkMuonPixelSeedsLowPt,
    hltIter2Phase2L3FromL1TkMuonCkfTrackCandidatesLowPt,
    hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracksLowPt,
    hltIter2Phase2L3FromL1TkMuonTrackCutClassifierLowPt,
    hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurityLowPt,
    hltIter2Phase2L3FromL1TkMuonMergedLowPt,
    hltPhase2L3GlbMuonLowPt,
    hltPhase2L3MuonMergedLowPt,
    hltPhase2L3MuonsNoIDLowPt,   
    hltPhase2L3MuonsLowPt,
    hltL3MuonsPhase2L3LinksLowPt,
    hltPhase2L3MuonTracks,
)
