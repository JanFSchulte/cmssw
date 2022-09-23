import FWCore.ParameterSet.Config as cms

hltPhase2L3MuonTracksLowPt = cms.EDProducer("HLTMuonTrackSelector",
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    copyTrajectories = cms.untracked.bool(False),
    muon = cms.InputTag("hltPhase2L3MuonsLowPt"),
    originalMVAVals = cms.InputTag("none"),
    track = cms.InputTag("hltPhase2L3MuonMergedLowPt")
)
