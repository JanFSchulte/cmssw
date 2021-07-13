import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/dmytro/tmp/store+data+Run2018C+SingleMuon+AOD+17Sep2018-v1+120000+F12B0069-26A7-2E42-B573-412918602C2B.root'
    ),
)

process.output = cms.OutputModule("PoolOutputModule",
   SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('skimPath')),
    fileName = cms.untracked.string('output.root'),
   # dataset = cms.untracked.PSet(
   #     filterName = cms.untracked.string(''),
   #     dataTier = cms.untracked.string('')
   # )
)

# filter
process.prefireVetoFilter = cms.EDFilter(
    "TriggerRulePrefireVetoFilter",
    tcdsRecordLabel = cms.InputTag("tcdsDigis","tcdsRecord"),
)


process.skimSequence = cms.Sequence(
    process.prefireVetoFilter 
    )

process.skimPath = cms.Path(process.skimSequence)

# output
process.output_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.skimPath, process.output_step
)

# Spit out filter efficiency at the end.                                                                                                                                         
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
