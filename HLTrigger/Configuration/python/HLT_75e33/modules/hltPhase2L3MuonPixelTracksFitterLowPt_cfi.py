import FWCore.ParameterSet.Config as cms

hltPhase2L3MuonPixelTracksFitterLowPt = cms.EDProducer("PixelFitterByHelixProjectionsProducer",
    scaleErrorsForBPix1 = cms.bool(False),
    scaleFactor = cms.double(0.65)
)
