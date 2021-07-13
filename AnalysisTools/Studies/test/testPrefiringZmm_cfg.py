import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

# process.load("Configuration.Geometry.GeometryIdeal_cff")
# process.load("Configuration.StandardSequences.MagneticField_cff")
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# process.load('Configuration.Geometry.GeometryRecoDB_cff')


# Need a proper geometry to do L1 muon object matching at low energy
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')

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

import commands
def add_files(path):
    files = commands.getoutput("find %s -type f -name '*root'"%path).split("\n")
    print "Number of files: %d" % len(files)
    for file in files:
        process.source.fileNames.extend(cms.untracked.vstring('file:%s'%file))

# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016B-07Aug17_ver2-v1_unprefirable/190714_073234')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016C-07Aug17-v1_unprefirable/190714_073319')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016D-07Aug17-v1_unprefirable/190714_073354')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016E-07Aug17-v1_unprefirable/190716_044107')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016F-07Aug17-v1_unprefirable/190714_073451')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016G-07Aug17-v1_unprefirable/190726_094125')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2016H-07Aug17-v1_unprefirable/190714_073545')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017B-17Nov2017-v1_unprefirable/190820_154737')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017C-17Nov2017-v1_unprefirable/190820_154812')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017D-17Nov2017-v1_unprefirable/190819_154116')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017E-17Nov2017-v1_unprefirable/190728_065825')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2017F-17Nov2017-v1_unprefirable/190820_154833')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018A-17Sep2018-v2_unprefirable/190725_130436')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018B-17Sep2018-v1_unprefirable/190715_070339')
# add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018C-17Sep2018-v1_unprefirable/190728_064656')
add_files('/eos/cms/store/group/phys_muon/dmytro/skims/SingleMuon/crab_skim_aod_Run2018D-22Jan2019-v2_unprefirable/190726_103909')

# process.prefireVetoFilter = cms.EDFilter("TriggerRulePrefireVetoFilter",
#     l1AcceptRecordLabel = cms.InputTag("scalersRawToDigi"),
#)

# Latest 2016 is cutBasedElectronID-Summer16-80X-V1-medium but close enough

process.ntuple = cms.EDAnalyzer("PrefiringMuAna",
    triggerRule = cms.InputTag("prefireVetoFilter:ruleIndex"),
    muonSrc = cms.InputTag("muons"),
    vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    l1Src = cms.InputTag("gmtStage2Digis:Muon"),
    l1GtSrc = cms.InputTag("gtStage2Digis"),
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    triggerEvents = cms.InputTag('hltTriggerSummaryAOD','','HLT'),
)

# process.skimPath = cms.Path(process.prefireVetoFilter+process.ntuple)
process.skimPath = cms.Path(process.ntuple)

process.TFileService = cms.Service("TFileService",
    # fileName = cms.string("Run2017F.root"),
    fileName = cms.string("Run2018D.root"),
    closeFileFast = cms.untracked.bool(True)
)
