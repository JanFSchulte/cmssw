import FWCore.ParameterSet.Config as cms

hltPhase2L3MuonCandidatesLowPt = cms.EDProducer("L3MuonCandidateProducerFromMuons",
    InputObjects = cms.InputTag("hltGlbTrkMuonsLowPt")
)
