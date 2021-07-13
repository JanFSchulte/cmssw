import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/data/Run2018D/EGamma/AOD/12Nov2019_UL2018-v8/2710000/05BA7ACB-FA4A-9E4B-AAA5-2795B86BC68C.root',
'/store/data/Run2018D/EGamma/AOD/12Nov2019_UL2018-v8/2710000/1025EEDC-2E86-D04D-948D-C31EDC22B049.root',
'/store/data/Run2018D/EGamma/AOD/12Nov2019_UL2018-v8/2710000/3885B259-5939-A84C-AD5C-AE4E13993A28.root',
'/store/data/Run2018D/EGamma/AOD/12Nov2019_UL2018-v8/2710000/1B362507-4377-8F40-A19C-CE339AF01146.root',
'/store/data/Run2018D/EGamma/AOD/12Nov2019_UL2018-v8/2710000/7FC9BE5A-9DA5-F04F-A17A-DE28455E6A86.root',
'/store/data/Run2018D/EGamma/AOD/12Nov2019_UL2018-v8/2710000/BFD1D043-81EC-FA4A-9C14-58F87298E8CA.root',


    ),
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, options.globaltag,'')



# filter
process.ntuple = cms.EDAnalyzer("PrefiringMuAnaEGamma",
    triggerRule = cms.InputTag("prefireVetoFilter:ruleIndex"),
    muonSrc = cms.InputTag("muons"),
    vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    l1Src = cms.InputTag("gmtStage2Digis:Muon"),
    l1BMTFSrc = cms.InputTag("gmtStage2Digis","BMTF","RECO"),
    l1OMTFSrc = cms.InputTag("gmtStage2Digis","OMTF","RECO"),
    l1EMTFSrc = cms.InputTag("gmtStage2Digis","EMTF","RECO"),
    l1GtSrc = cms.InputTag("gtStage2Digis"),
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    triggerEvents = cms.InputTag('hltTriggerSummaryAOD','','HLT'),
)
process.singleMuL1Filter = cms.EDFilter(
    "SingleMuL1Filter",
    l1Src = cms.InputTag("gmtStage2Digis:Muon"),
    muonSrc = cms.InputTag("muons"),
)


process.skimSequence = cms.Sequence(
    process.singleMuL1Filter * 
    process.ntuple
    )

process.skimPath = cms.Path(process.skimSequence)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output.root"),
    closeFileFast = cms.untracked.bool(True)
)
# output

process.schedule = cms.Schedule(
    process.skimPath
)

# Spit out filter efficiency at the end.                                                                                                                                         
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
