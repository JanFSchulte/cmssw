# One-shot utility: dump the AK4PFHLT JEC payload (as used by
# patJetCorrFactors in scoutingToMiniAODDerivedCollections_cff.py) from the
# GlobalTag to standard JEC text files, for standalone (non-CMSSW) use in
# the offline PUPPI recalibration's FastJet-reclustered jets (analysis/
# puppi_refit.py), which otherwise has no JEC applied at all.
import FWCore.ParameterSet.Config as cms

process = cms.Process("JECDUMP")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "150X_mcRun3_2024_realistic_v2", "")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.dump = cms.EDAnalyzer("JetCorrectorDBReader",
    payloadName = cms.untracked.string("AK4PFHLT"),
    globalTag = cms.untracked.string("150X_mcRun3_2024_realistic_v2"),
    printScreen = cms.untracked.bool(False),
    createTextFile = cms.untracked.bool(True),
)
process.p = cms.Path(process.dump)
