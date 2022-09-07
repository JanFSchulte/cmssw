import FWCore.ParameterSet.Config as cms

hltTau3MuSonicFilter = cms.EDFilter('HLTTau3MuSonicFilter',
  Client = cms.PSet(
    mode = cms.string('Async'),
    allowedTries = cms.untracked.uint32(0),
    verbose = cms.untracked.bool(False),
    modelName = cms.string("GNN_full_dR_1_ts"),
    modelVersion = cms.string('5'),
    modelConfigPath = cms.FileInPath("HeterogeneousCore/SonicTriton/data/models/GNN_full_dR_1_ts/config.pbtxt"),
    preferredServer = cms.untracked.string(''),
    timeout = cms.untracked.uint32(300),
    useSharedMemory = cms.untracked.bool(True),
    compression = cms.untracked.string(''),
    outputs = cms.untracked.vstring()
  ),
  L1EMTFHitInputTag = cms.InputTag('simEmtfDigis'),
  zeropad = cms.bool(False),
  gnn_threshold = cms.double(0.5),
  max_n_hits = cms.uint32(999999),
  max_n_edges = cms.uint32(999999),
  mightGet = cms.optional.untracked.vstring
)

