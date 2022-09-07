import FWCore.ParameterSet.Config as cms

hltL1MuonsPt0 = cms.EDProducer( "HLTL1TMuonSelector",
    L1MinPt = cms.double( -1.0 ),
    CentralBxOnly = cms.bool( True ),
    InputObjects = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    L1MinQuality = cms.uint32( 7 ),
    L1MaxEta = cms.double( 5.0 )
)

