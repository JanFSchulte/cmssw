import FWCore.ParameterSet.Config as cms

hltPhase2L3FromL1TkMuonPixelTracksHitDoubletsLowPt = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag(""),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltPhase2L3FromL1TkMuonPixelLayerQuadrupletsLowPt"),
    trackingRegions = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksTrackingRegionsLowPt"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)
