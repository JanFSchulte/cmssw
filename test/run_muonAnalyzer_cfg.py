'''Author: g. karathanasis. georgios.karathanasis@cern.ch
cfg to run tag and probe ntuple for muon POG. It runs both on AOD and miniAOD
Modified by Andre Frankenthal (a.franken@cern.ch) -- September 2020
usage: cmsRun run_muonAnalyzer_cfg.py option1=value1 option2=value2
'''

# cmsRun run_muonAnalyzer_cfg.py resonance=Z isFullAOD=False isMC=True globalTag=124X_mcRun3_2022_realistic_v5 era=Run2022 includeJets=False maxEvents=1000
# cmsRun run_muonAnalyzer_cfg.py resonance=Z isFullAOD=True isMC=True globalTag=124X_mcRun3_2022_realistic_v5 era=Run2022 includeJets=False maxEvents=1000
# cmsRun run_muonAnalyzer_cfg.py isFullAOD=False isMC=True globalTag=122X_mcRun3_2021_realistic_v9 era=Run2022 maxEvents=102

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('resonance', 'Z',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set resonance ('Z'/'JPsi')"
)

options.register('isFullAOD', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to False for MiniAOD datatier"
)

options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to True for MC"
)

options.register('globalTag', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)

options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report frequency"
)

options.register('numThreads', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Number of CMSSW threads" 
)

# this parameter is added for Jet Branches (ID varies for different era)
options.register('era', 'Run2018',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "era"
)

options.register('includeJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to False to exclude jets information in output ntuples"
)

options.register('fromCRAB', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is config run from CRAB"
)

options.parseArguments()

# defaults

if options._beenSet['globalTag'] and options.globalTag != '':
    globaltag = options.globalTag
else:
    globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'

# Run local test if no input files provided
if len(options.inputFiles) == 0:
    if options.resonance == 'Z':
        if options.isFullAOD:
            if options.isMC:
                options.inputFiles.append('/store/mc/Run3Winter22DRPremix/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/AODSIM/122X_mcRun3_2021_realistic_v9_ext2-v2/40000/0027e0ac-873f-4f49-884b-4c5b67c85724.root')
            else:
                options.inputFiles.append('/store/data/Run2022D/Muon/AOD/PromptReco-v2/000/357/734/00000/011a59b3-42bc-493c-bc2b-73cbddcb0e79.root')
        else:
            if options.isMC:
                options.inputFiles.append('/store/mc/Run3Winter22MiniAOD/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/MINIAODSIM/122X_mcRun3_2021_realistic_v9_ext2-v2/40000/ff8f64b6-9dd6-4080-8bad-20c78f1c9e39.root')
            else:
                options.inputFiles.append('/store/data/Run2022D/Muon/MINIAOD/PromptReco-v2/000/357/734/00000/0326173e-e2c7-4efc-8fdd-466a81844cbd.root')
    elif options.resonance == 'JPsi':
        if options.isFullAOD:
            if options.isMC:
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/270001/FFF2FC1D-18CB-7244-9663-4E36963494B7.root')
            else:
              options.inputFiles.append('/store/data/Run2018D/SingleMuon/AOD/12Nov2019_UL2018-v8/120000/00E9106D-CE60-5D4D-805A-E086AD3F6EEA.root')
        else:
            if options.isMC:
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/270001/FFF2FC1D-18CB-7244-9663-4E36963494B7.root')
            else:
                options.inputFiles.append('/store/data/Run2018C/SingleMuon/MINIAOD/UL2018_MiniAODv2-v2/130000/644CB07D-0BBF-5E4B-A8CA-79FA2AA576D4.root')


if options.outputFile=="":
    options.outputFile="output"
    if options.isMC:
        options.outputFile+="_mc"
    else:
        options.outputFile+="_data" 
    if options.isFullAOD:
        options.outputFile+="_full"
    else:
        options.outputFile+="_mini"
    options.outputFile+=".root"


process = cms.Process("MuonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag, '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles),
        secondaryFileNames=cms.untracked.vstring(),
        inputCommands=cms.untracked.vstring(
            'keep *',
            'drop *_ctppsPixelClusters_*_*'
        )
)

if options.includeJets:
    # for b-tagging
    process.load("RecoBTag.ImpactParameter.impactParameter_cff")
    process.load("RecoBTag.SecondaryVertex.secondaryVertex_cff")
    process.load("RecoBTag.SoftLepton.softLepton_cff")
    process.load("RecoBTag.Combined.combinedMVA_cff")
    process.load("RecoBTag.CTagging.cTagging_cff")
    process.load("RecoBTag.Combined.deepFlavour_cff")
    process.load("JetMETCorrections.Configuration.JetCorrectors_cff")

# Include pat:packedCandidateCollection in AOD for miniPFIsolation
if options.isFullAOD:   
    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
    process.load("PhysicsTools.PatAlgos.slimming.offlineSlimmedPrimaryVertices_cfi")
    process.load("PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi")
    from PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi import packedPFCandidates
    process.packedCandsForMuons = packedPFCandidates.clone()
    process.packedCandsForMuons.PuppiSrc=cms.InputTag("")
    process.packedCandsForMuons.PuppiNoLepSrc=cms.InputTag("")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(options.numThreads)
)

from MuonAnalysis.MuonAnalyzer.tools.ntuple_tools import *
if options.isFullAOD:
    if options.resonance == 'Z':
        process = muonAnalysis_customizeFullAOD_Z(process)
    else:
        process = muonAnalysis_customizeFullAOD_JPsi(process)
    if not options.isMC:
        process.muon.jetCorrector = cms.InputTag(
            "ak4PFCHSL1FastL2L3ResidualCorrector")
else:
    if options.resonance == 'Z':
        process = muonAnalysis_customizeMiniAOD_Z(process)
    else:
        process = muonAnalysis_customizeMiniAOD(process)

process.muon.isMC = options.isMC
process.muon.includeJets = options.includeJets
process.muon.era = options.era

# Trigger matching
muonSrc = "muons" if options.isFullAOD else "slimmedMuons"
from MuonAnalysis.MuonAssociators.muonL1Match_cfi import muonL1Match as _muonL1Match
process.muonL1Info = _muonL1Match.clone(
    src = cms.InputTag(muonSrc),
    useMB2InOverlap = cms.bool(True),
    useStage2L1 = cms.bool(True),
    preselection = cms.string(""),
    matched = cms.InputTag("gmtStage2Digis:Muon:"),

    useStation2 = cms.bool(True),
    cosmicPropagationHypothesis = cms.bool(False),
    propagatorAlong = cms.ESInputTag("", "SteppingHelixPropagatorAlong"),
    propagatorAny = cms.ESInputTag("", "SteppingHelixPropagatorAny"),
    propagatorOpposite = cms.ESInputTag("", "SteppingHelixPropagatorOpposite"),
    fallbackToME1 = cms.bool(False)
)
process.muonL1InfoByQ = process.muonL1Info.clone(
    sortBy = cms.string("quality"),
    sortByQuality  = cms.bool(True),
    sortByDeltaPhi = cms.bool(False),
    sortByDeltaEta = cms.bool(False),
    sortByPt       = cms.bool(False)
)

process.muon.fallbackToME1 = cms.bool(False)
process.muon.cosmicPropagationHypothesis = cms.bool(False)
process.muon.useMB2InOverlap = cms.bool(True)
process.muon.propagatorAlong = cms.ESInputTag("", "SteppingHelixPropagatorAlong")
process.muon.propagatorAny = cms.ESInputTag("", "SteppingHelixPropagatorAny")
process.muon.propagatorOpposite = cms.ESInputTag("", "SteppingHelixPropagatorOpposite")

from MuonAnalysis.MuonAnalyzer.hltInfo_cff import getHLTInfo, selectTriggers
hltInfo = getHLTInfo(options.resonance, options.era)
excludeDSA = (not options.isFullAOD)
process.muon.triggerPaths = cms.vstring(selectTriggers(hltInfo['triggerPaths'], True, False, excludeDSA))
process.muon.tagFilters = cms.vstring(selectTriggers(hltInfo['tagFilters'], not options.isFullAOD, True, excludeDSA))
process.muon.probeFilters = cms.vstring(selectTriggers(hltInfo['probeFilters'], not options.isFullAOD, True, excludeDSA))

# Standard selectors
from MuonAnalysis.MuonAnalyzer.selectorInfo_cff import getSelectorNamesAndBits
selectorNames, selectorBits = getSelectorNamesAndBits(options.era, options.isFullAOD)
process.muon.probeSelectorNames = cms.vstring(selectorNames)
process.muon.probeSelectorBits = cms.vuint32(selectorBits)

if options.isFullAOD:
    if options.includeJets:
        if not options.isMC:
	        process.analysis_step = cms.Path(
                process.primaryVertexAssociation +
                process.offlineSlimmedPrimaryVertices +
                process.packedCandsForMuons +
                process.muonL1Info +
                process.muonL1InfoByQ +
                process.ak4PFCHSL1FastL2L3ResidualCorrectorChain +
                process.muSequence
            )
        else:
            process.analysis_step = cms.Path(
                process.primaryVertexAssociation +
                process.offlineSlimmedPrimaryVertices +
                process.packedCandsForMuons +
                process.muonL1Info +
                process.muonL1InfoByQ +
                process.ak4PFCHSL1FastL2L3CorrectorChain +
                process.muSequence
	    )
    else:
        process.analysis_step = cms.Path(
            process.primaryVertexAssociation +
            process.offlineSlimmedPrimaryVertices +
            process.packedCandsForMuons +
            process.muonL1Info +
	        process.muonL1InfoByQ +
            process.muSequence
        )
else:
    if options.includeJets:
        if not options.isMC:
            process.analysis_step = cms.Path(
                process.muonL1Info +
                process.muonL1InfoByQ +
                process.ak4PFCHSL1FastL2L3ResidualCorrectorChain +
                process.muSequence
            )
        else:
            process.analysis_step = cms.Path(
                process.muonL1Info +
                process.muonL1InfoByQ +
                process.ak4PFCHSL1FastL2L3CorrectorChain +
                process.muSequence
            )
    else:
        process.analysis_step = cms.Path(
            process.muonL1Info +
            process.muonL1InfoByQ +
            process.muSequence
        )



process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile)
)
process.endjob_step = cms.EndPath(process.endOfProcess)


process.schedule = cms.Schedule(process.analysis_step, process.endjob_step)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

