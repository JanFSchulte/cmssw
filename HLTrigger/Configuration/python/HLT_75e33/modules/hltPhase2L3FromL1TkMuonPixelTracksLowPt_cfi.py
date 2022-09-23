import FWCore.ParameterSet.Config as cms

hltPhase2L3FromL1TkMuonPixelTracksLowPt = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltPhase2L3MuonPixelTracksFilterLowPt"),
    Fitter = cms.InputTag("hltPhase2L3MuonPixelTracksFitterLowPt"),
    SeedingHitSets = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksHitQuadrupletsLowPt"),
    passLabel = cms.string('')
)
