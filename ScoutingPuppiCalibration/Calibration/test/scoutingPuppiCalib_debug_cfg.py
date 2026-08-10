# TEMPORARY debug driver: same customisation as scoutingPuppiCalib_mc_cfg.py,
# but writes a standard PoolOutputModule EDM file (keeps EventAux, so FWLite
# can read it) with the raw PuppiProducer diagnostic products kept
# (including the new PuppiDbgPValUsed/PuppiDbgChi2 debug instrumentation
# added to CommonTools/PileupAlgos for the id==0 weight-reproduction
# investigation), instead of the NanoAOD-table-only output the main cfg uses.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run3_2024_cff import Run3_2024

options = VarParsing.VarParsing('analysis')
options.maxEvents = 100
options.outputFile = 'scouting_puppi_calib_debug.root'
options.parseArguments()

process = cms.Process('NANO', Run3_2024)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.custom_run3scouting_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(2),
    numberOfStreams = cms.untracked.uint32(0),
    wantSummary = cms.untracked.bool(False),
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_mcRun3_2024_realistic_v2', '')

process.nanoAOD_step = cms.Path(process.scoutingNanoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.debugOutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_scoutingPuppiNominal_*_*',
        'keep *_packedPFCandidates_*_*',
        'keep *_puppiAlphaDiagnostics_*_*',
        'keep *_scoutingPuppiCalibCandTable_*_*',
    ),
)
process.debugOutput_step = cms.EndPath(process.debugOutput)

process.schedule = cms.Schedule(process.nanoAOD_step, process.endjob_step, process.debugOutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from ScoutingPuppiCalibration.Calibration.scoutingPuppiCalibration_cff import customiseScoutingPuppiCalibrationNano
process = customiseScoutingPuppiCalibrationNano(process, "NANO", variants={"nominal": {}}, referenceVariant="nominal")

process.source.delayReadingEventProducts = cms.untracked.bool(False)
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
