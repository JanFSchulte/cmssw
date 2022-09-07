import FWCore.ParameterSet.Config as cms

hltPhase2L3OIL3MuonCandidatesLowPt = cms.EDProducer("L3MuonCandidateProducer",
    InputLinksObjects = cms.InputTag("hltPhase2L3OIL3MuonsLinksCombinationLowPt"),
    InputObjects = cms.InputTag("hltPhase2L3OIL3MuonsLowPt"),
    MuonPtOption = cms.string('Tracker')
)
