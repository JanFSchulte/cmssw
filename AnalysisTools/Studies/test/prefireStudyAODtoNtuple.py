import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#Passing_Command_Line_Arguments_T
# https://github.com/cms-sw/cmssw/blob/master/FWCore/ParameterSet/python/VarParsing.py
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')
# defaults for testing
options.inputFiles = ['/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018D-22Jan2019-v2_unprefirable/190726_103909/0000/output_378.root']
options.register('globaltag',
                 '80X_dataRun2_2016SeptRepro_v7',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")
options.outputFile = 'test.root'
options.parseArguments()

# Usage example:
# cmsRun prefireStudyAODtoNtuple.py inputFiles=file1.root,file2.root globaltag='102X_dataRun2_v11' outputFile='result.root'


# Need a proper geometry to do L1 muon object matching at low energy
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')
process.GlobalTag = GlobalTag(process.GlobalTag, options.globaltag,'')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring([
            
    ]),
)

for file in options.inputFiles:
    process.source.fileNames.extend(cms.untracked.vstring('file:%s'%file))

process.ntuple = cms.EDAnalyzer("PrefiringMuAna",
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

process.skimPath = cms.Path(process.ntuple)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile),
    closeFileFast = cms.untracked.bool(True)
)
