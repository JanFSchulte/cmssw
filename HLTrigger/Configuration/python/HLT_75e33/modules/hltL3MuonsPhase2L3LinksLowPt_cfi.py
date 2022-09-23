import FWCore.ParameterSet.Config as cms

hltL3MuonsPhase2L3LinksLowPt = cms.EDProducer("MuonLinksProducer",
    inputCollection = cms.InputTag("hltPhase2L3MuonsLowPt")
)
