import FWCore.ParameterSet.Config as cms

hltdstau3muDisplaced3muVtxProducer = cms.EDProducer( "HLTDisplacedmumumuVtxProducer",
    Src = cms.InputTag( "hltPhase2L3MuonCandidatesLowPt" ),
    PreviousCandTag = cms.InputTag( "hltdstau3mumuontrkFltr" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 1.2 ),
#    MinPtTriplet = cms.double( 8.0 ),
#    MinInvMass = cms.double( 1.6 ),
#    MaxInvMass = cms.double( 2.1 ),
    MinPtTriplet = cms.double( 0.0 ),
    MinInvMass = cms.double( 0. ),
    MaxInvMass = cms.double( 100. ),
    ChargeOpt = cms.int32( -1 )
)

