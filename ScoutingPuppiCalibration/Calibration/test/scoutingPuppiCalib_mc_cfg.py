# Standalone driver cfg for the scouting-PUPPI calibration study.
#
# Structurally the same as ScoutingNanoProduction/scoutingnano_mc_standalone2.py
# (same process/era/GlobalTag/output module), but parameterized via VarParsing
# and calling customiseScoutingPuppiCalibrationNano instead of just
# customiseScoutingNanoDerived, with a trimmed outputCommands list so files
# stay small when run over several QCD HT bins.
#
# Usage:
#   cmsRun scoutingPuppiCalib_mc_cfg.py maxEvents=300 \
#       inputFiles=root://xrootd-cms.infn.it//store/mc/.../someFile.root \
#       outputFile=calib_smoke.root
#
# NOTE: VarParsing('analysis')'s outputFile is a "tagString" that
# auto-appends "_numEvent<N>" whenever maxEvents>0 (standard CMSSW batch-job
# collision avoidance, not a bug here) -- e.g. outputFile=calib_smoke.root
# with maxEvents=300 actually writes calib_smoke_numEvent300.root.

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run3_2024_cff import Run3_2024

options = VarParsing.VarParsing('analysis')
options.maxEvents = 100
options.outputFile = 'scouting_puppi_calib.root'
options.parseArguments()

process = cms.Process('NANO', Run3_2024)

# import of standard configurations
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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(2),
    numberOfStreams = cms.untracked.uint32(0),
    wantSummary = cms.untracked.bool(False),
)

# Output definition: kept lean on purpose (see keepTables below) so running
# over several QCD HT bins doesn't produce bloated files.
process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_mcRun3_2024_realistic_v2', '')

process.nanoAOD_step = cms.Path(process.scoutingNanoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

process.schedule = cms.Schedule(process.nanoAOD_step, process.endjob_step, process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation: reuse the stock scouting NanoAOD customisation, then layer
# the PUPPI-variant calibration collections + tables on top.
from ScoutingPuppiCalibration.Calibration.scoutingPuppiCalibration_cff import customiseScoutingPuppiCalibrationNano
from ScoutingPuppiCalibration.Calibration.variants import DEFAULT_VARIANTS
process = customiseScoutingPuppiCalibrationNano(process, "NANO")

# Trim outputCommands to just what the calibration analysis needs: the
# baseline plain/CHS jet tables, every PUPPI variant's jet table, the shared
# candidate table, and pileup/gen-weight bookkeeping. Everything else from
# the stock NANOAODSIMEventContent (muons, electrons, tracks, V0s, SVs, ...)
# is dropped to keep files small across many HT bins.
#
# NOTE: 'keep'/'drop' patterns match on the EDM *module label*
# (nanoaodFlatTable_<moduleLabel>_<instance>_<process>), not the table's
# name= branch-name parameter -- these module labels must match the
# attribute names assigned in scoutingToMiniAODDerivedCollections_cff.py
# and scoutingPuppiCalibration_cff.py exactly.
keepModuleLabels = [
    "scoutingPFJetRecluster2Table", "scoutingPFJetRecluster2MCTable",
    "scoutingPFJetReclusterCHS2Table", "scoutingPFJetReclusterCHS2MCTable",
    "scoutingPuppiCalibCandTable",
    "puTable",
    "genJetTable",
    "rhoTable",
    # AK8 validation (_addAK8Jets, scoutingPuppiCalibration_cff.py) -- plain
    # AK8 plus one PUPPI-weighted AK8 table per variant, below, and the AK8
    # gen-truth table (a DIFFERENT branch/cut than AK4's GenJet -- needed to
    # resolve genJetAK8Idx, see io.py's genjet_idx_branch parameter).
    "scoutingFatPFJetRecluster2Table", "scoutingFatPFJetRecluster2MCTable",
    "genJetAK8Table",
    # Soft-drop-groomed mass companions (_addSoftDropTable in _addAK8Jets),
    # for boosted-decay (e.g. WWW) softdrop-mass validation.
    "scoutingFatPFJetRecluster2SoftDropTable",
    # MET/hadronic-recoil validation (_addMETTables, scoutingPuppiCalibration_cff.py).
    "pfMetTable", "rawPFMetTable", "genMetTable",
    # Muons, for dimuon Z reconstruction / hadronic-recoil validation
    # (_addMuonTable, scoutingPuppiCalibration_cff.py). NOT scoutingMuonTable
    # -- that name collides with a pre-existing, unrelated stock scouting
    # table (raw HLT Run3ScoutingMuon-based) -- see _addMuonTable's docstring.
    "scoutingPuppiCalibMuonTable",
]
for v, params in DEFAULT_VARIANTS.items():
    # module labels are camelCase, no underscores -- see the comment on
    # vcap in customiseForScoutingPuppiCalibration (scoutingPuppiCalibration_cff.py)
    vcap = v[0].upper() + v[1:]
    keepModuleLabels.append("scoutingPFJetReclusterPUPPI%sTable" % vcap)
    keepModuleLabels.append("scoutingPFJetReclusterPUPPI%sMCTable" % vcap)
    keepModuleLabels.append("scoutingFatPFJetReclusterPUPPI%sTable" % vcap)
    keepModuleLabels.append("scoutingFatPFJetReclusterPUPPI%sMCTable" % vcap)
    keepModuleLabels.append("scoutingFatPFJetReclusterPUPPI%sSoftDropTable" % vcap)
    if params.get("trackMatchedVertexAssociation"):
        # its own candidate table, separate from scoutingPuppiCalibCandTable
        # -- see _addTrackMatchedCandidateTable in scoutingPuppiCalibration_cff.py
        keepModuleLabels.append("scoutingPuppiCalibTrackMatchedVtxAssocCandTable")
    elif params.get("improvedVertexAssociation"):
        # its own candidate table, separate from scoutingPuppiCalibCandTable
        # -- see _addImprovedCandidateTable in scoutingPuppiCalibration_cff.py
        keepModuleLabels.append("scoutingPuppiCalibImprovedVtxAssocCandTable")

outputCommands = cms.untracked.vstring('drop *')
for m in keepModuleLabels:
    outputCommands.append('keep nanoaodFlatTable_%s_*_*' % m)
outputCommands += [
    'keep nanoaodFlatTable_genTable_*_*',
    'keep nanoaodMergeableCounterTable_*_*_*',
    'keep nanoaodUniqueString_nanoMetadata_*_*',
]
process.NANOAODSIMoutput.outputCommands = outputCommands

process.source.delayReadingEventProducts = cms.untracked.bool(False)
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
