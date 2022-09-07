import FWCore.ParameterSet.Config as cms

hltL3fL1TkDoubleMu3TkMuDsTau3MuPreFiltered22 = cms.EDFilter("HLTMuonTrkL1TkMuFilter",
    inputCandCollection = cms.InputTag("hltPhase2L3MuonCandidatesLowPt"),
    inputMuonCollection = cms.InputTag("hltGlbTrkMuonsLowPt"),
    maxAbsEta = cms.double(2.5),
    maxNormalizedChi2 = cms.double(1e+99),
    minMuonHits = cms.int32(-1),
    minMuonStations = cms.int32(1),
    minN = cms.uint32(2),
    minPt = cms.double(2.0),
    minTrkHits = cms.int32(-1),
    previousCandTag = cms.InputTag("hltL1DoubleMuFiltered2"),
    saveTags = cms.bool(True)
)
