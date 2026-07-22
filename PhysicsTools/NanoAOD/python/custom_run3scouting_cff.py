import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.run3scouting_cff import *
from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits
from PhysicsTools.NanoAOD.globals_cff import puTable
from PhysicsTools.NanoAOD.genWeightsTable_cfi import genWeightsTable
from PhysicsTools.NanoAOD.jetMC_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles, genIso, genParticleTable, genParticleTask, genParticleTablesTask
#from PhysicsTools.NanoAOD.nanogen_cff import customizeNanoGENFromMini 
############################
### Sub Task Definitions ###
############################

# Task contains all dependent tasks
# ExtensionTask must be run on top of another Task

#############################
# Scouting Original Objects #
#############################

# Scouting Muon
scoutingMuonTableTask = cms.Task(scoutingMuonTable)
scoutingMuonDisplacedVertexTableTask = cms.Task(scoutingMuonDisplacedVertexTable)

# from 2024, there are two muon collections (https://its.cern.ch/jira/browse/CMSHLT-3089)
# only the vertex-constrained collection is kept, to slim the output
(run3_scouting_2024 | run3_scouting_2025).toReplaceWith(scoutingMuonTableTask, cms.Task(scoutingMuonVtxTable))\
    .toReplaceWith(scoutingMuonDisplacedVertexTableTask, cms.Task(scoutingMuonVtxDisplacedVertexTable))

# muon matched gen particle index (only for MC); era-switched like scoutingMuonTableTask above
scoutingMuonGenPartMatchTask = cms.Task(scoutingMuonGenPartMatch)
(run3_scouting_2024 | run3_scouting_2025).toReplaceWith(scoutingMuonGenPartMatchTask, cms.Task(scoutingMuonVtxGenPartMatch))

# Scouting Electron
scoutingElectronTableTask = cms.Task(scoutingElectronTable)

# from 2023, scouting electron's tracks are added as std::vector since multiple tracks can be associated to a scouting electron
# plugin to select the best track to reduce to a single track per scouting electron is added
(run3_scouting_2023 | run3_scouting_2024 | run3_scouting_2025).toReplaceWith(
     scoutingElectronTableTask, cms.Task(scoutingElectronBestTrack, scoutingElectronTable)
)

# other collections are directly from original Run3Scouting objects, so unnessary to define tasks

############################
# Scouting Derived Objects #
############################

scoutingPFCandidateTask = cms.Task(scoutingPFCandidate, scoutingPFCandidateTable)

scoutingFatPFJetReclusterTask = cms.Task(
    scoutingPFCandidate, # translate to reco::PFCandidate, used as input
    scoutingFatPFJetRecluster, # jet clustering
    scoutingFatPFJetReclusterGlobalParticleTransformerJetTagInfos, scoutingFatPFJetReclusterGlobalParticleTransformerJetTags, # jet tagging with Global Particle Transformer
    scoutingFatPFJetReclusterSoftDrop, scoutingFatPFJetReclusterSoftDropMass, # softdrop mass
    scoutingFatPFJetReclusterEcfNbeta1, scoutingFatPFJetReclusterNjettiness, # substructure variables
    scoutingFatPFJetReclusterTable
)
scoutingFatPFJetReclusterMatchGenExtensionTask = cms.Task(
    genJetsAK8ForMatch, # pt>100 filtered slimmedGenJetsAK8, shared with GenJetAK8 table so indices line up
    scoutingFatPFJetReclusterMatchGen, # gen jet matching
    scoutingFatPFJetReclusterMatchFlavourAssociation, scoutingFatPFJetReclusterFlavourOnlyPATJets, scoutingFatPFJetReclusterMatchFlavour, # hadron/parton flavour
    scoutingFatPFJetReclusterTopWCategory, # top/W merging category
    scoutingFatPFJetReclusterGloParTCategory, # GloParT tagger truth category
    scoutingFatPFJetReclusterGenParticleMatch, scoutingFatPFJetReclusterGenPartIdxTable, # matched gen particle index
    scoutingFatPFJetReclusterMatchGenExtensionTable
)

scoutingFatPFJetReclusterCHSTask = cms.Task(
    scoutingPFCandidateCHS, # CHS translation to reco::PFCandidate
    scoutingFatPFJetReclusterCHS, # jet clustering
    scoutingFatPFJetReclusterCHSSoftDrop, scoutingFatPFJetReclusterCHSSoftDropMass, # softdrop mass
    scoutingFatPFJetReclusterCHSEcfNbeta1, scoutingFatPFJetReclusterCHSNjettiness, # substructure variables
    scoutingFatPFJetReclusterCHSTable
)
scoutingFatPFJetReclusterCHSMatchGenExtensionTask = cms.Task(
    scoutingFatPFJetReclusterCHSMatchGen, # gen jet matching
    scoutingFatPFJetReclusterCHSMatchFlavourAssociation, scoutingFatPFJetReclusterCHSFlavourOnlyPATJets, scoutingFatPFJetReclusterCHSMatchFlavour, # hadron/parton flavour
    scoutingFatPFJetReclusterCHSTopWCategory, # top/W merging category
    scoutingFatPFJetReclusterCHSGenParticleMatch, scoutingFatPFJetReclusterCHSGenPartIdxTable, # matched gen particle index
    scoutingFatPFJetReclusterCHSMatchGenExtensionTable
)

############################
# Trigger Bits and Objects #
############################

## L1 decisions
gtStage2DigisScouting = gtStage2Digis.clone(InputLabel="hltFEDSelectorL1")
l1bitsScouting = l1bits.clone(src="gtStage2DigisScouting") 

## L1 objects
from PhysicsTools.NanoAOD.l1trig_cff import *
l1MuScoutingTable = l1MuTable.clone(src=cms.InputTag("gtStage2DigisScouting", "Muon"))
l1EGScoutingTable = l1EGTable.clone(src=cms.InputTag("gtStage2DigisScouting", "EGamma"))
l1TauScoutingTable = l1TauTable.clone(src=cms.InputTag("gtStage2DigisScouting", "Tau"))
l1JetScoutingTable = l1JetTable.clone(src=cms.InputTag("gtStage2DigisScouting", "Jet"))
l1EtSumScoutingTable = l1EtSumTable.clone(src=cms.InputTag("gtStage2DigisScouting", "EtSum"))

# reduce the variables to the core variables as only these are available in gtStage2Digis
l1MuScoutingTable.variables = cms.PSet(l1MuonReducedVars)
l1EGScoutingTable.variables = cms.PSet(l1EGReducedVars)
l1TauScoutingTable.variables = cms.PSet(l1TauReducedVars)
l1JetScoutingTable.variables = cms.PSet(l1JetReducedVars)
l1EtSumScoutingTable.variables = cms.PSet(l1EtSumReducedVars)

##############################
### Main Tasks Definitions ###
##############################

# default configuration for ScoutingNano common for both data and MC
def prepareScoutingNanoTaskCommon():
    # Scouting original objects
    # all scouting objects are saved except PF Candidate and Track
    scoutingNanoTaskCommon = cms.Task()
    scoutingNanoTaskCommon.add(scoutingMuonTableTask, scoutingMuonDisplacedVertexTableTask)
    scoutingNanoTaskCommon.add(scoutingElectronTableTask)
    scoutingNanoTaskCommon.add(scoutingPrimaryVertexTable)
    scoutingNanoTaskCommon.add(scoutingMETTable)

    # Scouting derived objects
    scoutingNanoTaskCommon.add(scoutingFatPFJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingFatPFJetReclusterCHSTask)

    return scoutingNanoTaskCommon

# tasks related to trigger bits and objects
# L1 object tables (l1Mu/l1EG/l1Tau/l1Jet/l1EtSumScoutingTable) are defined
# above but deliberately not added here, to slim the output; the definitions
# are kept because customiseScoutingNanoForScoutingPFMonitor/FromMini below
# still reference them (e.g. to repoint their src for non-scouting L1 input).
def prepareScoutingTriggerTask():
    scoutingTriggerTask = cms.Task(gtStage2DigisScouting, l1bitsScouting)
    return scoutingTriggerTask

# additional tasks for running on MC
def prepareScoutingNanoTaskMC():
    scoutingNanoTaskMC = cms.Task()
    scoutingNanoTaskMC.add(scoutingFatPFJetReclusterMatchGenExtensionTask)
    scoutingNanoTaskMC.add(scoutingFatPFJetReclusterCHSMatchGenExtensionTask)

    scoutingNanoTaskMC.add(puTable)
    scoutingNanoTaskMC.add(genWeightsTable)
    scoutingNanoTaskMC.add(genJetTable)
    scoutingNanoTaskMC.add(patJetPartonsNano)
    scoutingNanoTaskMC.add(genJetFlavourAssociation)
    scoutingNanoTaskMC.add(genJetFlavourTable)

    # GenJetAK8 table, sourced from the same pt>100 filtered collection
    # (genJetsAK8ForMatch) that scoutingFatPFJetReclusterMatchGen/CHS match
    # against, so genJetAK8Idx lines up exactly with this table's row numbers.
    # cut is cleared since genJetsAK8ForMatch already applied it upstream.
    genJetAK8Table.src = cms.InputTag("genJetsAK8ForMatch")
    genJetAK8Table.cut = cms.string("")
    genJetAK8FlavourTable.src = genJetAK8Table.src
    genJetAK8FlavourTable.cut = genJetAK8Table.cut
    scoutingNanoTaskMC.add(genJetAK8Table)
    scoutingNanoTaskMC.add(genJetAK8FlavourAssociation)
    scoutingNanoTaskMC.add(genJetAK8FlavourTable)

    # GenPart table
    scoutingNanoTaskMC.add(genParticleTask)
    scoutingNanoTaskMC.add(genParticleTablesTask)

    # lepton matched gen particle index
    scoutingNanoTaskMC.add(scoutingMuonGenPartMatchTask)
    scoutingNanoTaskMC.add(scoutingElectronGenPartMatch)

    return scoutingNanoTaskMC

# Common tasks added to main scoutingNanoSequence
scoutingNanoTaskCommon = prepareScoutingNanoTaskCommon()
scoutingNanoSequence = cms.Sequence(scoutingNanoTaskCommon)

# Specific tasks which will be added to sequence during customization
scoutingTriggerTask = prepareScoutingTriggerTask()
scoutingTriggerSequence = cms.Sequence(scoutingTriggerTask)
scoutingNanoTaskMC = prepareScoutingNanoTaskMC()

def customiseScoutingNano(process):
    # if running with standard NanoAOD, triggerSequence is already added
    # if running standalone, triggerSequence need to be added
    if not ((hasattr(process, "nanoSequence") and process.schedule.contains(process.nanoSequence))
            or hasattr(process, "nanoSequenceMC") and process.schedule.contains(process.nanoSequenceMC)):
        process.trigger_step = cms.Path(process.scoutingTriggerSequence)
        process.schedule.extend([process.trigger_step])

    # specific tasks when running on MC
    runOnMC = hasattr(process,"NANOEDMAODSIMoutput") or hasattr(process,"NANOAODSIMoutput")
    if runOnMC:
        #process.load('PhysicsTools.NanoAOD.nanogen_cff')
        #process = customizeNanoGENFromMini(process)
        process.scoutingNanoSequence.associate(scoutingNanoTaskMC)
    else:
        # scoutingMuonTable/scoutingMuonVtxTable/scoutingElectronTable carry a
        # genPartIdx externalVariable pointing at the *GenPartMatch producers,
        # but those producers are only scheduled by scoutingNanoTaskMC above;
        # on data the getByToken for that ValueMap<int> would otherwise fail
        # with ProductNotFound, so drop the branch here instead.
        for tableName in ("scoutingMuonTable", "scoutingMuonVtxTable", "scoutingElectronTable"):
            table = getattr(process, tableName, None)
            if table is not None and hasattr(table.externalVariables, "genPartIdx"):
                del table.externalVariables.genPartIdx

    return process

##############
### Filter ###
##############

# this filter selects only events triggered by scouting paths by checking scouting primary dataset bit(s)
# if scouting paths are triggered, scouting objects will be reconstructed, so this gurantees that scouting objects  exist
import HLTrigger.HLTfilters.hltHighLevel_cfi
scoutingTriggerPathFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = cms.vstring("Dataset_ScoutingPFRun3", "Dataset_ScoutingPF0", "Dataset_ScoutingPF1"),
    throw = cms.bool(False)
)

#####################
### Customisation ###
#####################
# these function are designed to be used with --customise flag in cmsDriver.py
# e.g. --customise PhysicsTools/NanoAOD/python/custom_run3scouting_cff.addScoutingPFCandidate

# additional customisation for running with ScoutingPFMonitor/RAW inputs
# should be used with default customiseScoutingNano
# this is suitable when ScoutingPFMonitor/RAW is involved, e.g. RAW, RAW-MiniAOD two-file solution, full chain RAW-MiniAOD-NanoAOD
# when running full chain RAW-MiniAOD-NanoAOD, this ensures that gtStage2Digis, gmtStage2Digis, and caloStage2Digis are run
def customiseScoutingNanoForScoutingPFMonitor(process):
    process = skipEventsWithoutScoutingByEra(process)

    # replace gtStage2DigisScouting with standard gtStage2Digis
    process.scoutingTriggerTask.remove(process.gtStage2DigisScouting)
    process.scoutingTriggerTask.add(process.gtStage2Digis)

    # add gmtStage2Digis
    process.load("EventFilter.L1TRawToDigi.gmtStage2Digis_cfi")
    process.scoutingTriggerTask.add(process.gmtStage2Digis)

    # add caloStage2Digis
    process.load("EventFilter.L1TRawToDigi.caloStage2Digis_cfi")
    process.scoutingTriggerTask.add(process.caloStage2Digis)

    # replace l1bitsScouting with standard l1bits
    process.scoutingTriggerTask.remove(process.l1bitsScouting)
    process.scoutingTriggerTask.add(process.l1bits)

    # change src for l1 objects
    process.l1MuScoutingTable.src = cms.InputTag("gmtStage2Digis", "Muon")
    process.l1EGScoutingTable.src = cms.InputTag("caloStage2Digis", "EGamma")
    process.l1TauScoutingTable.src = cms.InputTag("caloStage2Digis", "Tau")
    process.l1JetScoutingTable.src = cms.InputTag("caloStage2Digis", "Jet")
    process.l1EtSumScoutingTable.src = cms.InputTag("caloStage2Digis", "EtSum")

    return process

# additional customisation for running with ScoutingPFMonitor/MiniAOD inputs alone
# can also be used on MC input
# should be used with default customiseScoutingNano and NOT with customiseScoutingNanoForScoutingPFMonitor
def customiseScoutingNanoFromMini(process):
    # when running on data, assume ScoutingPFMonitor/MiniAOD dataset as inputs
    runOnData = hasattr(process,"NANOAODSIMoutput") or hasattr(process,"NANOAODoutput")
    if runOnData:
        process = skipEventsWithoutScoutingByEra(process)

    # remove gtStage2Digis since they are already run for Mini
    process.scoutingTriggerTask.remove(process.gtStage2DigisScouting)

    # replace l1bitsScouting with standard l1bits
    process.scoutingTriggerTask.remove(process.l1bitsScouting)
    process.scoutingTriggerTask.add(process.l1bits)

    # change src for l1 objects
    process.l1MuScoutingTable.src = cms.InputTag("gmtStage2Digis", "Muon")
    process.l1EGScoutingTable.src = cms.InputTag("caloStage2Digis", "EGamma")
    process.l1TauScoutingTable.src = cms.InputTag("caloStage2Digis", "Tau")
    process.l1JetScoutingTable.src = cms.InputTag("caloStage2Digis", "Jet")
    process.l1EtSumScoutingTable.src = cms.InputTag("caloStage2Digis", "EtSum")

    return process

# skip events without scouting object products
# this may be useful since ScoutingPFMonitor dataset contains some events which do not contain scouting object products in 2022-24
# this is fixed for 2025: https://its.cern.ch/jira/browse/CMSHLT-3331
def skipEventsWithoutScouting(process):
    # add filter/skim path to the process
    process.scoutingNanoSkim_step = cms.Path(process.scoutingTriggerPathFilter)

    # schedule filter/skim path to run
    process.schedule.extend([process.scoutingNanoSkim_step])

    if hasattr(process, "NANOAODoutput"):
        process.NANOAODoutput.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("scoutingNanoSkim_step"))

    if hasattr(process, "NANOEDMAODoutput"):
        process.NANOEDMAODoutput.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("scoutingNanoSkim_step"))

    if hasattr(process, "write_NANOAOD"): # PromptReco
        process.write_NANOAOD.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("scoutingNanoSkim_step")) 

    return process

# skip events without scouting object products by era
# this may be useful since ScoutingPFMonitor dataset contains some events which do not contain scouting object products in 2022-24
# this is fixed for 2025: https://its.cern.ch/jira/browse/CMSHLT-3331
def skipEventsWithoutScoutingByEra(process):
    # add filter/skim path to the process
    process.scoutingNanoSkim_step = cms.Path(process.scoutingTriggerPathFilter)

    if hasattr(process, "NANOAODoutput"):
        (~run3_scouting_2025).toModify(process.NANOAODoutput, SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("scoutingNanoSkim_step")))
        if hasattr(process.NANOAODoutput, "SelectEvents") and "scoutingNanoSkim_step" in process.NANOAODoutput.SelectEvents.SelectEvents:
            # schedule filter/skim path to run
            process.schedule.extend([process.scoutingNanoSkim_step])

    if hasattr(process, "NANOEDMAODoutput"):
        (~run3_scouting_2025).toModify(process.NANOEDMAODoutput, SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("scoutingNanoSkim_step")))
        if hasattr(process.NANOEDMAODoutput, "SelectEvents") and "scoutingNanoSkim_step" in process.NANOEDMAODoutput.SelectEvents.SelectEvents:
            # schedule filter/skim path to run
            process.schedule.extend([process.scoutingNanoSkim_step])

    if hasattr(process, "write_NANOAOD"): # PromptReco
        (~run3_scouting_2025).toModify(process.write_NANOAOD, SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("scoutingNanoSkim_step")))
        if hasattr(process.write_NANOAOD, "SelectEvents") and "scoutingNanoSkim_step" in process.write_NANOAOD.SelectEvents.SelectEvents:
            # schedule filter/skim path to run
            process.schedule.extend([process.scoutingNanoSkim_step])

    return process

def addScoutingTrack(process):
    process.scoutingNanoSequence.associate(cms.Task(scoutingTrackTable))
    return process

def addScoutingParticle(process):
    # original PF candidate without post-processing
    process.scoutingNanoSequence.associate(cms.Task(scoutingParticleTable))
    return process

def addScoutingPFCandidate(process):
    # PF candidate after translation to reco::PFCandidate
    process.scoutingNanoSequence.associate(scoutingPFCandidateTask)
    return process

# this adds all electron tracks in addition to best track selected
# this should be only used with ScoutingElectron format from 2023
def addScoutingElectronTrack(process):
    process.scoutingElectronTable.externalVariables.bestTrack_index\
            = ExtVar(cms.InputTag("scoutingElectronBestTrack", "Run3ScoutingElectronBestTrackIndex"), int, doc="best track index")

    process.scoutingElectronTable.collectionVariables = cms.PSet(
        ScoutingElectronTrack = cms.PSet(
            name = cms.string("ScoutingElectronTrack"),
            doc = cms.string("Scouting Electron Track"),
            useCount = cms.bool(True),
            useOffset = cms.bool(True),
            variables = cms.PSet(
                d0 = Var("trkd0", "float", doc="track d0"),
                dz = Var("trkdz", "float", doc="track dz"),
                pt = Var("trkpt", "float", doc="track pt"),
                eta = Var("trketa", "float", doc="track eta"),
                phi = Var("trkphi", "float", doc="track phi"),
                chi2overndf = Var("trkchi2overndf", "float", doc="track normalized chi squared"),
                charge = Var("trkcharge", "int", doc="track charge"),
            ),
        ),
    )
    
    # additional electron track variables added in 2024 in https://github.com/cms-sw/cmssw/pull/43744
    (run3_scouting_2024 | run3_scouting_2025).toModify(
        process.scoutingElectronTable.collectionVariables.variables,
        pMode = Var("trkpMode", "float", doc="track pMode"),
        etaMode = Var("trketaMode", "float", doc="track etaMode"),
        phiMode = Var("trkphiMode", "float", doc="track phiMode"),
        qoverpModeError = Var("trkqoverpModeError", "float", doc="track qoverpModeError"),
    )
    return process

# use for samples with no relevant top/W/Z gen truth chain (pure multijet
# QCD): every ScoutingFatPFJetRecluster jet's topWCategory is then
# unconditionally "Others" (see scoutingFatPFJetReclusterTopWCategory in
# run3scouting_cff.py). Should NOT be applied for ttbar/single-top/ttV/VV/
# Z+jets samples, which use the default (applyTopWMerging=True) to get the
# actual Top-merged/W-merged/Z-merged/Non-merged categorization.
def setScoutingFatJetTopWCategoryAsBackground(process):
    process.scoutingFatPFJetReclusterTopWCategory.applyTopWMerging = cms.bool(False)
    return process

# use for samples with no relevant resonance gen truth chain (pure multijet
# QCD): every ScoutingFatPFJetRecluster jet's gloParTCategory is then
# unconditionally "QCD" (see scoutingFatPFJetReclusterGloParTCategory in
# run3scouting_cff.py). Should NOT be applied for samples with a genuine
# resonance decay (Higgs, Z', or any other model), which use the default
# (applyGloParTMatching=True) to get the actual per-class categorization.
def setScoutingFatJetGloParTCategoryAsBackground(process):
    process.scoutingFatPFJetReclusterGloParTCategory.applyGloParTMatching = cms.bool(False)
    return process

