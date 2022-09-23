import FWCore.ParameterSet.Config as cms

hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracksLowPt = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter2LowPt'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEventLowPt"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter2Phase2L3FromL1TkMuonCkfTrackCandidatesLowPt"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)
