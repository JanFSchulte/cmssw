import FWCore.ParameterSet.Config as cms
process= cms.Process("AN2")
process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
        '/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016H-PromptReco-v2_l1calo_prefire/170812_110837/0000/output_1.root '
        # '/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00933E2A-A0D5-E611-B2CD-00266CF89130.root'
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_0-selected.root',
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_1-selected.root',
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_2-selected.root',
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_3-selected.root',
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_4-selected.root',
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_5-selected.root',
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_6-selected.root',
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_7-selected.root',
        # 'file:/eos/user/d/dmytro/studies/170806-103510/prefire_study_singlemuon20_8-selected.root'
        ),
)

process.analyzer = cms.EDAnalyzer("L1TriggerInfoDumper",
)

process.maxEvents = cms.untracked.PSet(
    # input= cms.untracked.int32(1000)
    input= cms.untracked.int32(10)
)

process.path = cms.Path( process.analyzer )

process.schedule = cms.Schedule(
    process.path
)
