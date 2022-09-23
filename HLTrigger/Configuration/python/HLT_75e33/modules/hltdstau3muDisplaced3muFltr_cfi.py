import FWCore.ParameterSet.Config as cms

hltdstau3muDisplaced3muFltr = cms.EDFilter( "HLTDisplacedmumumuFilter",
    saveTags = cms.bool( True ),
    FastAccept = cms.bool( False ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999.0 ),
    MinVtxProbability = cms.double( 0.0 ),
    MinCosinePointingAngle = cms.double( -999 ),
    DisplacedVertexTag = cms.InputTag( "hltdstau3muDisplaced3muVtxProducer" ),
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltPhase2L3MuonCandidatesLowPt" )
)

