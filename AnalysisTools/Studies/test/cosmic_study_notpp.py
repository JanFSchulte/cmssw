import FWCore.ParameterSet.Config as cms
process= cms.Process("AN")
process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
        # 'file:test/Run2016-CosmicSP-PromptReco_750GeV_failed_sta.root'
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016A-CosmicSP-PromptReco-v1_cosmics/170530_235123/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016A-CosmicSP-PromptReco-v2_cosmics/170530_235144/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v1_cosmics/170530_235211/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v1_cosmics/170530_235211/0000/output_2.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v1_cosmics/170530_235211/0000/output_4.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v1_cosmics/170530_235211/0000/output_5.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v2_cosmics/170530_235237/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v2_cosmics/170530_235237/0000/output_10.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v2_cosmics/170530_235237/0000/output_11.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v2_cosmics/170530_235237/0000/output_2.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v2_cosmics/170530_235237/0000/output_4.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v2_cosmics/170530_235237/0000/output_5.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v2_cosmics/170530_235237/0000/output_7.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016B-CosmicSP-PromptReco-v2_cosmics/170530_235237/0000/output_9.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016C-CosmicSP-PromptReco-v2_cosmics/170530_235714/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016C-CosmicSP-PromptReco-v2_cosmics/170530_235714/0000/output_3.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016D-CosmicSP-PromptReco-v2_cosmics/170531_000134/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016E-CosmicSP-PromptReco-v2_cosmics/170531_000312/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016E-CosmicSP-PromptReco-v2_cosmics/170531_000312/0000/output_2.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016F-CosmicSP-PromptReco-v1_cosmics/170531_000618/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016F-CosmicSP-PromptReco-v1_cosmics/170531_000618/0000/output_2.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016F-CosmicSP-PromptReco-v1_cosmics/170531_000618/0000/output_3.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016F-CosmicSP-PromptReco-v1_cosmics/170531_000618/0000/output_4.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016F-CosmicSP-PromptReco-v1_cosmics/170531_000618/0000/output_5.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016G-CosmicSP-PromptReco-v1_cosmics/170531_000657/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016G-CosmicSP-PromptReco-v1_cosmics/170531_000657/0000/output_2.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016G-CosmicSP-PromptReco-v1_cosmics/170531_000657/0000/output_3.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016G-CosmicSP-PromptReco-v1_cosmics/170531_000657/0000/output_4.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016G-CosmicSP-PromptReco-v1_cosmics/170531_000657/0000/output_5.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016G-CosmicSP-PromptReco-v1_cosmics/170531_000657/0000/output_6.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016G-CosmicSP-PromptReco-v1_cosmics/170531_000657/0000/output_7.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v1_cosmics/170531_000745/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v1_cosmics/170531_000745/0000/output_2.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_1.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_10.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_11.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_12.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_13.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_14.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_15.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_16.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_17.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_18.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_19.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_2.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_20.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_21.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_22.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_23.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_24.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_25.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_26.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_27.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_28.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_29.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_3.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_30.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_31.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_32.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_33.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_34.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_35.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_36.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_37.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_38.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_39.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_4.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_40.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_41.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_42.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_43.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_44.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_5.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_6.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_7.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_8.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v2_cosmics/170530_185616/0000/output_9.root',
            '/store/user/dmytro/Cosmics/crab_skim_raw-reco_Run2016H-CosmicSP-PromptReco-v3_cosmics/170531_001004/0000/output_1.root',
        ),
      inputCommands = cms.untracked.vstring("keep *", 
                                            "drop *_lumiProducer_*_*",
                                            "drop *_hltGtStage2ObjectMap__HLT",
                                            "drop *_cosmicDCTracks__RECO",
                                            "drop *_castorDigis__RECO")
)

process.analyzer = cms.EDFilter("CosmicMuonAnalyzer",
    output_path = cms.string("/eos/user/d/dmytro/www/plots/cosmic_muon_notpp_Run2016_splitMuons_bottom/"),
    # output_path = cms.string("/eos/user/d/dmytro/www/plots/test/"),
    # muons = cms.InputTag("muons"),
    muons = cms.InputTag("splitMuons"),
    # tracks = cms.InputTag("generalTracks"),
    # tracks = cms.InputTag("ctfWithMaterialTracksP5LHCNavigation"),
    tracks = cms.InputTag("splittedTracksP5"),
    min_pt = cms.double(100),
    min_pt_interesting = cms.double(750),
    max_d0 = cms.double(10),
    min_d0 = cms.double(-1),
    pp_mode = cms.bool(False),
    tk_hits = cms.bool(False),
    tight_id = cms.bool(False),
    debug = cms.bool(False),
    tag_location = cms.int32(-1)                            
)

process.maxEvents = cms.untracked.PSet(
    # input= cms.untracked.int32(10)
    input= cms.untracked.int32(-1)
)

# output
process.output = cms.OutputModule("PoolOutputModule",
   SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filterPath')),
    fileName = cms.untracked.string('Run2016-CosmicSP-PromptReco_750GeV_failed_sta_bottom.root'),
)
process.output_step = cms.EndPath(process.output)


process.filterPath = cms.Path( process.analyzer )

process.schedule = cms.Schedule(
    process.filterPath, 
    # process.output_step
)

