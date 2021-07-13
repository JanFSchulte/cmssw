import FWCore.ParameterSet.Config as cms
process= cms.Process("AN")
process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016C-23Sep2016-v1_cosmics/170519_144216/0000/output_1.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016C-23Sep2016-v1_cosmics/170519_144216/0000/output_3.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016C-23Sep2016-v1_cosmics/170519_144216/0000/output_4.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016D-23Sep2016-v1_cosmics/170519_144255/0000/output_1.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016D-23Sep2016-v1_cosmics/170519_144255/0000/output_2.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016D-23Sep2016-v1_cosmics/170519_144255/0000/output_3.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016D-23Sep2016-v1_cosmics/170519_144255/0000/output_4.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016D-23Sep2016-v1_cosmics/170519_144255/0000/output_5.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016D-23Sep2016-v1_cosmics/170519_144255/0000/output_6.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016E-23Sep2016-v1_cosmics/170519_144324/0000/output_1.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016E-23Sep2016-v1_cosmics/170519_144324/0000/output_2.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016E-23Sep2016-v1_cosmics/170519_144324/0000/output_3.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016E-23Sep2016-v1_cosmics/170519_144324/0000/output_4.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016E-23Sep2016-v1_cosmics/170519_144324/0000/output_5.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016F-23Sep2016-v1_cosmics/170519_144509/0000/output_1.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016F-23Sep2016-v1_cosmics/170519_144509/0000/output_2.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016F-23Sep2016-v1_cosmics/170519_144509/0000/output_3.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016F-23Sep2016-v1_cosmics/170519_144509/0000/output_4.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_1.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_2.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_3.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_4.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_5.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_6.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_7.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_8.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016G-23Sep2016-v1_cosmics/170519_144527/0000/output_9.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_1.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_10.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_2.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_3.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_4.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_5.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_6.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_7.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_8.root',
        '/store/user/dmytro/NoBPTX/crab_skim_aod_Run2016H-PromptReco-v2_cosmics/170519_114833/0000/output_9.root'
        ),
      inputCommands = cms.untracked.vstring("keep *", "drop *_lumiProducer_*_*")
)

process.analyzer = cms.EDAnalyzer("CosmicMuonAnalyzer",
    output_path = cms.string("/eos/user/d/dmytro/www/plots/cosmic_muon/"),
    muons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    min_pt = cms.double(100),
    max_d0 = cms.double(10),
    min_d0 = cms.double(1),
)

process.maxEvents = cms.untracked.PSet(
    # input= cms.untracked.int32(1000)
    input= cms.untracked.int32(-1)
)

process.path = cms.Path( process.analyzer )

process.schedule = cms.Schedule(
    process.path
)
