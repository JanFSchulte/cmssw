import FWCore.ParameterSet.Config as cms

from ..tasks.HLTL3MuonLowPtTask_cfi import *

HLTL3MuonLowPtSequence = cms.Sequence(
    HLTL3MuonLowPtTask
)
