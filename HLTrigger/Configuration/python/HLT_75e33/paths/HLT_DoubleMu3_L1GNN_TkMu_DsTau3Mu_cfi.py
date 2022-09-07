import FWCore.ParameterSet.Config as cms

#from ..modules.hltL1TkMuons_cfi import *
#from ..modules.hltL1SingleMuFiltered5_cfi import *
#from ..modules.hltL1DoubleMuFiltered2_cfi import *
#from ..modules.hltL3fL1TkTripleMu533L31055DZFiltered0p2_cfi import *
#from ..modules.hltL3fL1TkTripleMu533L3Filtered1055_cfi import *
from ..modules.hltL3fL1TkDoubleMu3TkMuDsTau3MuSonicPreFiltered_cfi import *
#from ..modules.hltTripleMuon3DR0_cfi import *
#from ..modules.hltTripleMuon3DZ1p0_cfi import *
from ..modules.hltCsc2DRecHits_cfi import *
from ..modules.hltCscSegments_cfi import *
from ..modules.hltDimuonMergergingTkMu_cfi import * 
from ..modules.hltDimuonLinksTkMuMerge_cfi import *
from ..modules.hltdstau3mumuontrkFltr_cfi import *
from ..modules.hltdstau3muDisplaced3muVtxProducer_cfi import *
from ..modules.hltdstau3muDisplaced3muFltr_cfi import *
from ..modules.hltDt1DRecHits_cfi import *
from ..modules.hltDt4DSegments_cfi import *
from ..modules.hltGemRecHits_cfi import *
from ..modules.hltGemSegments_cfi import *
from ..modules.hltPhase2L3MuonCandidatesLowPt_cfi import *
from ..modules.hltGlbTrkMuonsLowPt_cfi import *
from ..modules.hltGlbTrkMuonsNoIDLowPt_cfi import *
from ..modules.hltL2MuonFromL1MuonCandidates_cfi import *
from ..modules.hltL2MuonSeedsFromL1Muon_cfi import *
from ..modules.hltL2MuonsFromL1Muon_cfi import *
from ..modules.hltL2OfflineMuonSeeds_cfi import *
from ..modules.hltMe0RecHits_cfi import *
from ..modules.hltMe0Segments_cfi import *
from ..modules.hltRpcRecHits_cfi import *
from ..modules.MeasurementTrackerEvent_cfi import *
from ..modules.siPhase2Clusters_cfi import *
from ..modules.siPixelClusters_cfi import *
from ..modules.siPixelClusterShapeCache_cfi import *
from ..modules.siPixelRecHits_cfi import *
from ..modules.hltTau3MuSonicFilter_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTTrackingV61Sequence_cfi import *
from ..sequences.HLTL3MuonLowPtSequence_cfi import *

HLT_DoubleMu3_L1GNN_TkMu_DsTau3Mu = cms.Path(
    HLTBeginSequence +
    hltTau3MuSonicFilter + 
    hltL3fL1TkDoubleMu3TkMuDsTau3MuSonicPreFiltered + 
    #hltL1TripleMuFiltered3 +
    HLTTrackingV61Sequence +
    HLTL3MuonLowPtSequence +
    #hltL1SingleMuFiltered5 +
    #hltTripleMuon3DZ1p0 +
    #hltTripleMuon3DR0 +
    #hltL3fL1TkTripleMu533PreFiltered555 +
    #hltL3fL1TkTripleMu533L3Filtered1055 +
    #hltL3fL1TkTripleMu533L31055DZFiltered0p2 +
    hltdstau3mumuontrkFltr + 
    hltdstau3muDisplaced3muFltr + 
    HLTEndSequence,
    cms.Task(
        MeasurementTrackerEvent,
        hltCsc2DRecHits,
        hltCscSegments,
        hltDt1DRecHits,
        hltDt4DSegments,
        hltDimuonMergergingTkMu,
        hltDimuonLinksTkMuMerge,
        hltdstau3muDisplaced3muVtxProducer,
        hltGemRecHits,
        hltGemSegments,
        hltGlbTrkMuonsNoIDLowPt,
        hltGlbTrkMuonsLowPt,
        hltPhase2L3MuonCandidatesLowPt,
        hltL2MuonFromL1MuonCandidates,
        hltL2MuonSeedsFromL1Muon,
        hltL2MuonsFromL1Muon,
        hltL2OfflineMuonSeeds,
        hltMe0RecHits,
        hltMe0Segments,
        hltRpcRecHits,
        siPhase2Clusters,
        siPixelClusterShapeCache,
        siPixelClusters,
        siPixelRecHits
    )
)
