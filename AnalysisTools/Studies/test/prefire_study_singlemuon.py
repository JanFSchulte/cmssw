import FWCore.ParameterSet.Config as cms
import commands
process= cms.Process("AN")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/d/dmytro/projects/CMSSW_9_2_6/src/AnalysisTools/TinyNtupleMaker/test/output.root'
            # '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_prefire/170608_173715/0000/output_1.root',
            # '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_prefire/170608_173715/0000/output_2.root',
            # '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_prefire/170608_173715/0000/output_4.root',
            # '/store/user/dmytro/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_prefire/170608_173715/0000/output_10.root',
        ),
      inputCommands = cms.untracked.vstring("keep *", 
                                            "drop *_lumiProducer_*_*",
                                            "drop *_hltGtStage2ObjectMap__HLT",
                                            "drop *_cosmicDCTracks__RECO",
                                            "drop *_castorDigis__RECO")
)

process.source.fileNames = cms.untracked.vstring()
## files = commands.getoutput("find /eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016G-23Sep2016-v1_singlemuon20_prefire/170805_140355/ | grep root | grep -v failed")
# files = commands.getoutput("find /eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_l1calo_prefire/170812_110837/ | grep root | grep -v failed")
# list_of_files = files.split("\n")
# for f in list_of_files:
#     process.source.fileNames.extend(["file:%s"%f])
# # print process.source.fileNames

process.analyzer = cms.EDFilter("GMTPrefireStudy",
    output_path = cms.string("/eos/user/d/dmytro/www/plots/gmt_prefire_SingleMuon_l1calo_Run2016/"),
    # output_path = cms.string("/eos/user/d/dmytro/www/plots/test/"),
    min_pt = cms.double(20),
    max_d0 = cms.double(10),
    max_dz = cms.double(20),
    min_deta = cms.double(0.2),
    nmuons = cms.uint32(1),
    hist_filename = cms.string("histograms.root"),
    make_plots = cms.bool(True)
)

## disable event-loop by reseting input list to empty string
# process.source.fileNames = cms.untracked.vstring()

process.maxEvents = cms.untracked.PSet(
    input= cms.untracked.int32(5000)
    # input= cms.untracked.int32(-1)
)

## output
process.output = cms.OutputModule("PoolOutputModule",
   SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filterPath')),
    fileName = cms.untracked.string('selected.root'),
)
process.output_step = cms.EndPath(process.output)

process.filterPath = cms.Path( process.analyzer )

process.schedule = cms.Schedule(
    process.filterPath, 
    process.output_step
)
