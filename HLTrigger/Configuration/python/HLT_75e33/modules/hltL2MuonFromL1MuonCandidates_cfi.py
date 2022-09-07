import FWCore.ParameterSet.Config as cms

hltL2MuonFromL1MuonCandidates = cms.EDProducer("L2MuonCandidateProducer",
    InputObjects = cms.InputTag("hltL2MuonsFromL1Muon","UpdatedAtVtx")
)
