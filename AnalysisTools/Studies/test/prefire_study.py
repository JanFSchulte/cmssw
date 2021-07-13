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
      # inputCommands = cms.untracked.vstring("keep *", 
      #                                       "drop *_lumiProducer_*_*",
      #                                       "drop *_hltGtStage2ObjectMap__HLT",
      #                                       "drop *_cosmicDCTracks__RECO",
      #                                       "drop *_castorDigis__RECO")
)

process.source.fileNames = cms.untracked.vstring()
files = commands.getoutput("find /eos/cms/store/user/dmytro/DoubleMuon/crab_skim_aod_doublemuon_Run2017C-PromptReco-v1_prefire/170724_075536/ | grep root")
list_of_files = files.split("\n")
for f in list_of_files:
    process.source.fileNames.extend(["file:%s"%f])
files = commands.getoutput("find /eos/cms/store/user/dmytro/DoubleMuon/crab_skim_aod_doublemuon_Run2017C-PromptReco-v1_prefire_ext1/170726_070003/ | grep root")
list_of_files = files.split("\n")
for f in list_of_files:
    process.source.fileNames.extend(["file:%s"%f])
# print process.source.fileNames

process.analyzer = cms.EDFilter("GMTPrefireStudy",
    output_path = cms.string("/eos/user/d/dmytro/www/plots/gmt_prefire_DoubleMuon_Run2017C_1.5fb/"),
    min_pt = cms.double(20),
    max_d0 = cms.double(10),
    max_dz = cms.double(20),
    min_deta = cms.double(0.2),
)

process.maxEvents = cms.untracked.PSet(
    # input= cms.untracked.int32(5000)
    input= cms.untracked.int32(-1)
)

process.path = cms.Path( process.analyzer )

process.schedule = cms.Schedule(
    process.path
)
