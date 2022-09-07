import FWCore.ParameterSet.Config as cms

hltDimuonLinksTkMuMerge = cms.EDProducer( "MuonLinksProducerForHLT",
    InclusiveTrackerTrackCollection = cms.InputTag( "hltDimuonMergergingTkMu" ),
    LinkCollection = cms.InputTag( "hltL3MuonsPhase2L3Links" ),
    ptMin = cms.double( 2.5 ),
    pMin = cms.double( 2.5 ),
    shareHitFraction = cms.double( 0.19 )
)
