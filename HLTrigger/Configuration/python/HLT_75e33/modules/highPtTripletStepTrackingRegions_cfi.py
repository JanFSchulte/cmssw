import FWCore.ParameterSet.Config as cms

highPtTripletStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(4),
        originHalfLength = cms.double(0),
        originRadius = cms.double(2),
        precise = cms.bool(True),
        ptMin = cms.double(0.5),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)
