import FWCore.ParameterSet.Config as cms

hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurityLowPt = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter2Phase2L3FromL1TkMuonTrackCutClassifierLowPt","MVAValues"),
    originalQualVals = cms.InputTag("hltIter2Phase2L3FromL1TkMuonTrackCutClassifierLowPt","QualityMasks"),
    originalSource = cms.InputTag("hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracksLowPt")
)
