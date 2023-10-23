#include "NtupleContent.h"

NtupleContent::NtupleContent() :
  branches{
    {"run", {-1}},
    {"event", {-1}},
    {"lumi", {-1}},
    {"fromFullAOD", {false}},
    {"genWeight", {-99.f}},
    {"BSpot_x", {-99.f}},
    {"BSpot_y", {-99.f}},
    {"BSpot_z", {-99.f}},
    {"pv_x", {-99.f}},
    {"pv_y", {-99.f}},
    {"pv_z", {-99.f}},
    {"nVertices", {0}},
    {"nTrueInteractions", {-99.f}},
    {"nPUInteractions", {-99.f}},
    {"rho", {-1.d}},
    {"nmuons", {0}},
    {"ntag", {0}},
    {"npairs", {0}},
    {"genmu1_pt", {-99.f}},
    {"genmu1_eta", {-99.f}},
    {"genmu1_phi", {-99.f}},
    {"genmu1_charge", {0}},
    {"genmu2_pt", {-99.f}},
    {"genmu2_eta", {-99.f}},
    {"genmu2_phi", {-99.f}},
    {"genmu2_charge", {0}},
    {"genMass", {-99.f}},
    {"genmuFSfromHP1_pt", {-99.f}},
    {"genmuFSfromHP1_eta", {-99.f}},
    {"genmuFSfromHP1_phi", {-99.f}},
    {"genmuFSfromHP1_charge", {0}},
    {"genmuFSfromHP2_pt", {-99.f}},
    {"genmuFSfromHP2_eta", {-99.f}},
    {"genmuFSfromHP2_phi", {-99.f}},
    {"genmuFSfromHP2_charge", {0}},
    {"genMassFSfromHP", {-99.f}},
    {"genJets_pt", {-99.f}},
    {"genJets_eta", {-99.f}},
    {"genJets_phi", {-99.f}},
    {"genJets_mass", {-99.f}},
    {"jets_pt", {-99.f}},
    {"jets_eta", {-99.f}},
    {"jets_phi", {-99.f}},
    {"jets_mass", {-99.f}},
    {"jets_isTight", {-99}},
    {"jets_isTightLepVeto", {-99}},
    {"jets_bTag_deepCSV", {-99.f}},
    {"jets_bTag_deepFlav", {-99.f}},
    {"nTightJets", {0}},
    {"nTightLepVetoJets", {0}},
    {"tag_pt", {-99.f}},
    {"tag_eta", {-99.f}},
    {"tag_phi", {-99.f}},
    {"tag_charge", {0}},
    {"tag_inner_pterr", {-99.f}},
    {"tag_dxy", {-99.f}},
    {"tag_dz", {-99.f}},
    {"tag_isPF", {false}},
    {"tag_isSA", {false}},
    {"tag_isdSA", {false}},
    {"tag_isTracker", {false}},
    {"tag_isGlobal", {false}},
    {"tag_isRPC", {false}},
    {"tag_isLoose", {false}},
    {"tag_isMedium", {false}},
    {"tag_isTight", {false}},
    {"tag_isSoft", {false}},
    {"tag_isHighPt", {false}},
    {"tag_tuneP_ExistingRefit", {false}},
    {"tag_ZprimeMatchedStations", {false}},
    {"tag_SIP3D", {-99.f}},
    {"tag_SIP3D_err", {-99.f}},
    {"tag_absTrkIso04", {-99.f}},
    {"tag_absTrkIso03", {-99.f}},
    {"tag_miniIso", {-1.f}},
    {"tag_miniIsoCharged", {0.f}},
    {"tag_miniIsoPhotons", {0.f}},
    {"tag_miniIsoNeutrals", {0.f}},
    {"tag_isMatchedGen", {false}},
    {"tag_minDR", {0.f}},
    {"tag_ptRel_minDR", {0.f}},
    {"tag_iso03_sumPt", {-99.f}},
    {"tag_pfIso03_charged", {-99.f}},
    {"tag_pfIso03_neutral", {-99.f}},
    {"tag_pfIso03_photon", {-99.f}},
    {"tag_pfIso03_sumPU", {-99.f}},
    {"tag_combRelIsoPF03dBeta", {-99.f}},
    {"tag_pfIso04_charged", {-99.f}},
    {"tag_pfIso04_neutral", {-99.f}},
    {"tag_pfIso04_photon", {-99.f}},
    {"tag_pfIso04_sumPU", {-99.f}},
    {"tag_combRelIsoPF04dBeta", {-99.f}},
    {"tag_tuneP_pt", {-99.f}},
    {"tag_tuneP_eta", {-99.f}},
    {"tag_tuneP_phi", {-99.f}},
    {"tag_tuneP_muonHits", {-99}},
    {"tag_tuneP_charge", {-99}},
    {"tag_tuneP_pterr", {-99.f}},
    {"tag_nsegments", {-99}},
    {"tag_inner_validFraction", {-99.f}},
    {"tag_inner_trackerLayers", {-99}},
    {"tag_inner_pixelLayers", {-99}},
    {"tag_inner_pixelHits", {-99}},
    {"tag_inner_pt", {-99.f}},
    {"tag_inner_eta", {-99.f}},
    {"tag_inner_phi", {-99.f}},
    {"tag_inner_charge", {0}},
    {"tag_GlobalValidHits", {-99}},
    {"tag_RPCLayers", {-99}},
    {"tag_tpfms_pt", {-99.f}},
    {"tag_tpfms_eta", {-99.f}},
    {"tag_tpfms_phi", {-99.f}},
    {"tag_tpfms_muonHits", {-99}},
    {"tag_tpfms_charge", {-99}},
    {"tag_tpfms_pterr", {-99.f}},
    {"tag_picky_pt", {-99.f}},
    {"tag_picky_eta", {-99.f}},
    {"tag_picky_phi", {-99.f}},
    {"tag_picky_muonHits", {-99}},
    {"tag_picky_charge", {-99}},
    {"tag_picky_pterr", {-99.f}},
    {"tag_dyt_pt", {-99.f}},
    {"tag_dyt_eta", {-99.f}},
    {"tag_dyt_phi", {-99.f}},
    {"tag_dyt_muonHits", {-99}},
    {"tag_dyt_charge", {-99}},
    {"tag_dyt_pterr", {-99.f}},
    {"iprobe", {0}},
    {"probe_pt", {-99.f}},
    {"probe_eta", {-99.f}},
    {"probe_phi", {-99.f}},
    {"probe_charge", {0}},
    {"probe_inner_pt", {-99.f}},
    {"probe_inner_eta", {-99.f}},
    {"probe_inner_phi", {-99.f}},
    {"probe_inner_charge", {0}},
    {"probe_outer_pt", {-99.f}},
    {"probe_outer_eta", {-99.f}},
    {"probe_outer_phi", {-99.f}},
    {"probe_outer_charge", {0}},
    {"probe_global_pt", {-99.f}},
    {"probe_global_eta", {-99.f}},
    {"probe_global_phi", {-99.f}},
    {"probe_global_charge", {0}},
    {"probe_best_pt", {-99.f}},
    {"probe_best_eta", {-99.f}},
    {"probe_best_phi", {-99.f}},
    {"probe_best_charge", {0}},
    {"probe_inner_pterr", {-99.f}},
    {"probe_dxy", {-99.f}},
    {"probe_dz", {-99.f}},
    {"probe_isPF", {false}},
    {"probe_isSA", {false}},
    {"probe_isTracker", {false}},
    {"probe_isGlobal", {false}},
    {"probe_isLoose", {false}},
    {"probe_isMedium", {false}},
    {"probe_isTight", {false}},
    {"probe_isSoft", {false}},
    {"probe_isHighPt", {false}},
    {"probe_isRPC", {false}},
    {"probe_isArbitratedTracker", {false}},
    {"probe_isMuMatched", {false}},
    {"probe_tuneP_ExistingRefit", {false}},
    {"probe_isdSA", {false}},
    {"probe_isdGlobal", {false}},
    {"probe_isCosmic", {false}},
    {"probe_ncosmic", {-99}},
    {"probe_cosmic_minDR", {+99.f}},
    {"probe_isHighPurity", {false}},
    {"probe_isDuplicated", {false}},
    {"probe_isBestPair", {false}},
    {"probe_inner_validFraction", {-99.f}},
    {"probe_trkChi2", {-99.f}},
    {"probe_positionChi2", {-99.f}},
    {"probe_trkKink", {-99.f}},
    {"probe_inner_trackerLayers", {-99}},
    {"probe_inner_pixelLayers", {-99}},
    {"probe_muonStations", {-99}},
    {"probe_muonHits", {-99}},
    {"probe_DTHits", {-99}},
    {"probe_CSCHits", {-99}},
    {"probe_SIP3D", {-99.f}},
    {"probe_SIP3D_err", {-99.f}},
    {"probe_absTrkIso04", {-99.f}},
    {"probe_absTrkIso03", {-99.f}},
    {"probe_miniIso", {-1.f}},
    {"probe_miniIsoCharged", {0.f}},
    {"probe_miniIsoPhotons", {0.f}},
    {"probe_miniIsoNeutrals", {0.f}},
    {"probe_isMatchedGen", {false}},
    {"probe_minDR", {0.f}},
    {"probe_ptRel_minDR", {0.f}},
    {"probe_iso03_sumPt", {-99.f}},
    {"probe_pfIso03_charged", {-99.f}},
    {"probe_pfIso03_neutral", {-99.f}},
    {"probe_pfIso03_photon", {-99.f}},
    {"probe_pfIso03_sumPU", {-99.f}},
    {"probe_pfIso04_charged", {-99.f}},
    {"probe_pfIso04_neutral", {-99.f}},
    {"probe_pfIso04_photon", {-99.f}},
    {"probe_pfIso04_sumPU", {-99.f}},
    {"probe_inner_pixelHits", {-99}},
    {"probe_matchedStations", {-99}},
    {"probe_expectedMatchedStations", {-99}},
    {"probe_RPCLayers", {-99}},
    {"probe_stationMask", {-99}},
    {"probe_nShowers", {-99}},
    {"probe_tuneP_pt", {-99.f}},
    {"probe_tuneP_eta", {-99.f}},
    {"probe_tuneP_phi", {-99.f}},
    {"probe_tuneP_charge", {-99}},
    {"probe_tuneP_pterr", {-99.f}},
    {"probe_tuneP_muonHits", {-99}},
    {"probe_nsegments", {-99}},
    {"probe_tpfms_pt", {-99.f}},
    {"probe_tpfms_eta", {-99.f}},
    {"probe_tpfms_phi", {-99.f}},
    {"probe_tpfms_charge", {-99.f}},
    {"probe_tpfms_pterr", {-99.f}},
    {"probe_tpfms_muonHits", {-99}},
    {"probe_picky_pt", {-99.f}},
    {"probe_picky_eta", {-99.f}},
    {"probe_picky_phi", {-99.f}},
    {"probe_picky_charge", {-99.f}},
    {"probe_picky_pterr", {-99.f}},
    {"probe_picky_muonHits", {-99}},
    {"probe_dyt_pt", {-99.f}},
    {"probe_dyt_eta", {-99.f}},
    {"probe_dyt_phi", {-99.f}},
    {"probe_dyt_charge", {-99.f}},
    {"probe_dyt_pterr", {-99.f}},
    {"probe_dyt_muonHits", {-99}},
    {"l1pt", {-99.f}},
    {"l1q", {-99.f}},
    {"l1dr", {-99.f}},
    {"l1ptByQ", {-99.f}},
    {"l1qByQ", {-99.f}},
    {"l1drByQ", {-99.f}},
    {"tag_l1pt", {-99.f}},
    {"tag_l1q", {-99.f}},
    {"tag_l1dr", {-99.f}},
    {"tag_l1ptByQ", {-99.f}},
    {"tag_l1qByQ", {-99.f}},
    {"tag_l1drByQ", {-99.f}},
    {"probe_dsa_segmentMatches", {-99}},
    {"probe_dsa_nsegments", {-99}},
    {"probe_dsa_muonStations", {-99}},
    {"probe_dsa_muonHits", {-99}},
    {"probe_dsa_DTHits", {-99}},
    {"probe_dsa_CSCHits", {-99}},
    {"probe_dsa_pterr", {0.f}},
    {"probe_dsa_dxy", {-99.f}},
    {"probe_dsa_dz", {-99.f}},
    {"probe_dsa_trkChi2", {-99.f}},
    {"probe_dsa_pt", {0.f}},
    {"probe_dsa_eta", {-99.f}},
    {"probe_dsa_phi", {-99.f}},
    {"probe_dsa_outerEta", {-99.f}},
    {"probe_dsa_outerPhi", {-99.f}},
    {"probe_dsa_minDR", {+99.f}},
    {"probe_dsa_charge", {-99}},
    {"tag_dsa_segmentMatches", {-99}},
    {"tag_dsa_nsegments", {-99}},
    {"tag_dsa_muonStations", {-99}},
    {"tag_dsa_muonHits", {-99}},
    {"tag_dsa_DTHits", {-99}},
    {"tag_dsa_CSCHits", {-99}},
    {"tag_dsa_pterr", {0.f}},
    {"tag_dsa_dxy", {-99.f}},
    {"tag_dsa_dz", {-99.f}},
    {"tag_dsa_trkChi2", {-99.f}},
    {"tag_dsa_pt", {0.f}},
    {"tag_dsa_eta", {-99.f}},
    {"tag_dsa_phi", {-99.f}},
    {"tag_dsa_outerEta", {-99.f}},
    {"tag_dsa_outerPhi", {-99.f}},
    {"tag_dsa_minDR", {+99.f}},
    {"tag_dsa_charge", {-99}},
    {"probe_dgl_segmentMatches", {-99.f}},
    {"probe_dgl_nsegments", {-99}},
    {"probe_dgl_muonStations", {-99}},
    {"probe_dgl_muonHits", {-99}},
    {"probe_dgl_totalHits", {-99}},
    {"probe_dgl_outerTrackerHits", {-99}},
    {"probe_dgl_trackerHits", {-99}},
    {"probe_dgl_DTHits", {-99}},
    {"probe_dgl_CSCHits", {-99}},
    {"probe_dgl_pterr", {0.f}},
    {"probe_dgl_dxy", {-99.f}},
    {"probe_dgl_dz", {-99.f}},
    {"probe_dgl_trkChi2", {-99.f}},
    {"probe_dgl_pt", {0.f}},
    {"probe_dgl_eta", {-99.f}},
    {"probe_dgl_phi", {-99.f}},
    {"probe_dgl_charge", {-99.f}},
    {"probe_dgl_minDR", {+99.f}},
    {"pair_pt", {-99.f}},
    {"pair_eta", {-99.f}},
    {"pair_phi", {-99.f}},
    {"pair_mass", {-99.f}},
    {"pair_fit_mass", {0.f}},
    {"pair_svprob", {0.f}},
    {"pair_normalchi2", {0.f}},
    {"pair_dz", {-99.f}},
    {"pair_dR", {-99.f}},
    {"pair_drM1", {-99.f}},
    {"pair_rank_vtx_prob", {-1}},
    {"pair_rank_dz_PV_SV", {-1.f}},
    {"pair_rank_dPhi_muons", {-1.f}},
    {"pair_rank_dM_Z_Mmumu", {-1.f}},
    {"pair_tuneP_pt", {-99.f}},
    {"pair_tuneP_eta", {-99.f}},
    {"pair_tuneP_phi", {-99.f}},
    {"pair_tuneP_mass", {-99.f}},
    {"pair_tuneP_fit_mass", {-99.f}},
    {"pair_tuneP_svprob", {-99.f}},
    {"pair_tuneP_normalchi2", {-99.f}},
    {"pair_tuneP_dz", {-99.f}},
    {"pair_tuneP_dR", {-99.f}},
    {"tag_simType", {-99}},
    {"tag_simExtType", {-99}},
    {"tag_simFlavour", {-99}},
    {"tag_simHeaviestMotherFlavour", {-99}},
    {"tag_simPdgId", {-99}},
    {"tag_simMotherPdgId", {-99}},
    {"tag_simBX", {-99}},
    {"tag_simProdRho", {-99.f}},
    {"tag_simProdZ", {-99.f}},
    {"tag_simPt", {-99.f}},
    {"tag_simEta", {-99.f}},
    {"tag_simPhi", {-99.f}},
    {"probe_simType", {-99}},
    {"probe_simExtType", {-99}},
    {"probe_simFlavour", {-99}},
    {"probe_simHeaviestMotherFlavour", {-99}},
    {"probe_simPdgId", {-99}},
    {"probe_simMotherPdgId", {-99}},
    {"probe_simBX", {-99}},
    {"probe_simProdRho", {-99.f}},
    {"probe_simProdZ", {-99.f}},
    {"probe_simPt", {-99.f}},
    {"probe_simEta", {-99.f}},
    {"probe_simPhi", {-99.f}},
}
{}

NtupleContent::~NtupleContent() {}

void NtupleContent::SetTree(TTree *mytree) { t1 = mytree; }

void NtupleContent::CreateBranches(const std::vector<std::string> &HLTs,
                                   const std::vector<std::string> &selectorNames) {
  // General
  for (auto & [name, branch] : branches) {
    if (std::holds_alternative<bool>(branch.value))
      t1->Branch(name, & std::get<bool>(branch.value));
    else if (std::holds_alternative<int>(branch.value))
      t1->Branch(name, & std::get<int>(branch.value));
    else if (std::holds_alternative<float>(branch.value))
      t1->Branch(name, & std::get<float>(branch.value));
    else if (std::holds_alternative<double>(branch.value))
      t1->Branch(name, & std::get<double>(branch.value));
    else if (std::holds_alternative<unsigned>(branch.value))
      t1->Branch(name, & std::get<unsigned>(branch.value));
    else if (std::holds_alternative<long unsigned>(branch.value))
      t1->Branch(name, & std::get<long unsigned>(branch.value));
    else if (std::holds_alternative<long long unsigned>(branch.value))
      t1->Branch(name, & std::get<long long unsigned>(branch.value));
    else if (std::holds_alternative<std::vector<int>>(branch.value))
      t1->Branch(name, & std::get<std::vector<int>>(branch.value));
    else if (std::holds_alternative<std::vector<float>>(branch.value))
      t1->Branch(name, & std::get<std::vector<float>>(branch.value));
    else if (std::holds_alternative<std::vector<bool>>(branch.value))
      t1->Branch(name, & std::get<std::vector<bool>>(branch.value));
  }

  // Trigger info
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++)
    t1->Branch(TString(HLTs[ihlt]), &trigger[ihlt]);

  // selectors for probe
  for (unsigned int isel = 0; isel < selectorNames.size(); ++isel) {
    t1->Branch(TString("probe_" + selectorNames[isel]), &probe_selectors[isel]);
  }

  

}

void NtupleContent::CreateExtraTrgBranches(const std::vector<std::string> &HLTs, bool isTag = false) {
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++) {
    if (isTag) {
      t1->Branch(TString("tag_" + HLTs[ihlt]), &tag_trg[ihlt]);
      t1->Branch(TString("tag_" + HLTs[ihlt] + "_pt"), &tag_trg_pt[ihlt]);
      t1->Branch(TString("tag_" + HLTs[ihlt] + "_eta"), &tag_trg_eta[ihlt]);
      t1->Branch(TString("tag_" + HLTs[ihlt] + "_phi"), &tag_trg_phi[ihlt]);
      t1->Branch(TString("tag_" + HLTs[ihlt] + "_dr"), &tag_trg_dr[ihlt]);
    } else {
      t1->Branch(TString("probe_" + HLTs[ihlt]), &probe_trg[ihlt]);
      t1->Branch(TString("probe_" + HLTs[ihlt] + "_pt"), &probe_trg_pt[ihlt]);
      t1->Branch(TString("probe_" + HLTs[ihlt] + "_eta"), &probe_trg_eta[ihlt]);
      t1->Branch(TString("probe_" + HLTs[ihlt] + "_phi"), &probe_trg_phi[ihlt]);
      t1->Branch(TString("probe_" + HLTs[ihlt] + "_dr"), &probe_trg_dr[ihlt]);
    }
  }
}

void NtupleContent::ClearBranches() {
 
  for (unsigned int itrg = 0; itrg < NTRIGGERMAX; itrg++) {
    trigger[itrg] = false;
    tag_trg[itrg] = false;
    tag_trg_pt[itrg] = -99;
    tag_trg_eta[itrg] = -99;
    tag_trg_phi[itrg] = -99;
    tag_trg_dr[itrg] = 99;
    probe_trg[itrg] = false;
    probe_trg_pt[itrg] = -99;
    probe_trg_eta[itrg] = -99;
    probe_trg_phi[itrg] = -99;
    probe_trg_dr[itrg] = 99;
  }

  trg_filter.clear();
  trg_pt.clear();
  trg_eta.clear();
  trg_phi.clear();
  prb_filter.clear();
  prb_pt.clear();
  prb_eta.clear();
  prb_phi.clear();

  for (unsigned int isel = 0; isel < 100; isel++) {
    probe_selectors[isel] = false;
  }

 
}
