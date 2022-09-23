import FWCore.ParameterSet.Config as cms

hltL2MuonSeedsFromL1Muon = cms.EDProducer("L2MuonSeedGeneratorFromL1T",
        CentralBxOnly = cms.bool(True),
        EtaMatchingBins = cms.vdouble(0.0, 2.5),
        GMTReadoutCollection = cms.InputTag(""),
        InputObjects = cms.InputTag("hltGtStage2Digis","Muon"),
        L1MaxEta = cms.double(2.5),
        L1MinPt = cms.double(0.0),
        L1MinQuality = cms.uint32(7),
        MatchDR = cms.vdouble(0.3),
        MatchType = cms.uint32(0),
        OfflineSeedLabel = cms.untracked.InputTag("hltL2OfflineMuonSeeds"),
        Propagator = cms.string('SteppingHelixPropagatorAny'),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring('SteppingHelixPropagatorAny'),
            RPCLayers = cms.bool(True),
            UseMuonNavigation = cms.untracked.bool(True)
        ),
        SetMinPtBarrelTo = cms.double(3.5),
        SetMinPtEndcapTo = cms.double(1.0),
        SortType = cms.uint32(0),
        UseOfflineSeed = cms.untracked.bool(True),
        UseUnassociatedL1 = cms.bool(False)
)
