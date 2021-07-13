import FWCore.ParameterSet.Config as cms
process= cms.Process("AN")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_1.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_10.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_11.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_12.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_14.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_15.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_16.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_17.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_18.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_2.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_20.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_21.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_23.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_25.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_27.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_29.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_3.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_30.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_31.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_32.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_33.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_34.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_36.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_37.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_38.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_4.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_40.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_41.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_42.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_43.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_44.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_45.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_46.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_47.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_48.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_49.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_5.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_50.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_51.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_52.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_53.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_54.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_55.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_56.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_57.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_58.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_59.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_61.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_62.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_63.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_64.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_65.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_66.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_67.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_68.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_69.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_7.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_70.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_71.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_73.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_74.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_75.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_76.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_77.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_78.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_79.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_80.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_81.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_82.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_83.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_84.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_85.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_87.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_89.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_9.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_90.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_91.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_92.root',
        '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_smart_prescale/170617_084207/0000/output_93.root',
        ),
      inputCommands = cms.untracked.vstring("keep *", 
                                            "drop *_lumiProducer_*_*",
                                            "drop *_hltGtStage2ObjectMap__HLT",
                                            "drop *_cosmicDCTracks__RECO",
                                            "drop *_castorDigis__RECO")
)

process.analyzer = cms.EDFilter("CosmicMuonAnalyzer",
    output_path = cms.string("/eos/user/d/dmytro/www/plots/test/"),
    make_plots = cms.bool(False),
    muons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    min_pt = cms.double(100),
    max_d0 = cms.double(10),
    min_d0 = cms.double(-1),
    pp_mode = cms.bool(True),
    tk_hits = cms.bool(False),
    tight_id = cms.bool(False),
    debug = cms.bool(False),                            
    tag_location = cms.int32(0),
    min_pt_interesting = cms.double(600)                            
)

process.maxEvents = cms.untracked.PSet(
    # input= cms.untracked.int32(10)
    input= cms.untracked.int32(-1)
)

# output
process.output = cms.OutputModule("PoolOutputModule",
   SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filterPath')),
    fileName = cms.untracked.string('Run2016H-PromptReco_600GeV_failed_sta.root'),
)
process.output_step = cms.EndPath(process.output)


process.filterPath = cms.Path( process.analyzer )

process.schedule = cms.Schedule(
    process.filterPath, 
    process.output_step
)

