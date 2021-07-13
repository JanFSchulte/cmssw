import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # "file:/eos/cms/store/user/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+BuToKJpsi_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+PUPoissonAve20_102X_upgrade2018_realistic_v15-v2+60000+71247BCA-D0BB-3947-A624-CECE01B2BC19.root"
        # "file:/afs/cern.ch/work/d/dmytro/projects/bparking-reco/src/calo-miniaod.root"
        "file:/eos/cms/store/group/phys_muon/dmytro/tmp/bparking-reco-test/calo-miniaod.root"
    ),
)

from AnalysisTools.TinyNtupleMaker.trigger_names import *

trigger_info = get_triggers()
print "Trigger names:"
print trigger_info[0]
print "Trigger bits:"
print trigger_info[1]


process.tree = tree = cms.EDAnalyzer('HeavyFlavDileptonNtupleMakerMiniAOD',
      electrons     = cms.InputTag("slimmedElectrons"),
      muons         = cms.InputTag("slimmedMuons"),
      vertices      = cms.InputTag("offlineSlimmedPrimaryVertices"),
      beamSpot      = cms.InputTag('offlineBeamSpot'),                                     
      packed        = cms.InputTag("packedGenParticles"),
      pruned        = cms.InputTag("prunedGenParticles"),
      pfCands       = cms.InputTag("packedPFCandidates"),                      
      trigger       = cms.InputTag("TriggerResults","","HLT"),
      prescales     = cms.InputTag('patTrigger'),
      triggerNames  = trigger_info[0],
      triggerBits   = trigger_info[1],
      ### Object selection requirements ###
      electronMinPt = cms.double(10),
      electronId    = cms.uint32( 0),
      #
      muon1MinPt     = cms.double(3.5),
      muon2MinPt     = cms.double(1.5),
      hadMinPt       = cms.double(1.5),
      llMaxMass      = cms.double(4),
      muon1Id        = cms.uint32(0), 
      muon2Id        = cms.uint32(0), 
      minLLvtxProb   = cms.double(0.001),                               
      minBvtxProb    = cms.double(0.001),                               
      minTagPt       = cms.double(7),
      requireTag     = cms.bool(False),                               
      ### Event selection requirements ### 
      # Skim_MuonFakes            = 1UL<<0, // pt>20, LooseID, VLooseIso, prescaled triggers
      # Skim_ElectronFakes        = 1UL<<1, // pt>20, VetoID, VetoIso, prescaled triggers
      # Skim_SingleMuon           = 1UL<<2, // pt>30, LooseID, LooseIso
      # Skim_SingleElectron       = 1UL<<3, // pt>30, LooseID, LooseIso
      # Skim_Dilepton             = 1UL<<4, // pt 20/20
      # Skim_DileptonMET          = 1UL<<5, // Dilep + MET>30
      # Skim_MET                  = 1UL<<6, // MET>200
      # Skim_Recoil               = 1UL<<7, // Recoil (MET build ignoring, els, mus and photons)>200
      # Skim_Photon               = 1UL<<8, // pt>175 with loose selection requirements
      # requirement is treated as OR of different skims
      skimming      = cms.uint32( 0b111111111 )
)


# disable skimming
process.tree.skimming = cms.uint32(0)

process.TFileService = cms.Service("TFileService",
    # fileName = cms.string("HeavyFlavDileptonNtupleMakerMiniAOD.root"),
    fileName = cms.string("calomuon.root"),
    closeFileFast = cms.untracked.bool(True)
)


process.p = cms.Path(process.tree)
