//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// filling functions for aod and miniaod tag/probe

#ifndef MuonAnalysis_MuonAnalyzer_MuonBranches
#define MuonAnalysis_MuonAnalyzer_MuonBranches

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

#include <type_traits>
#include "NtupleContent.h"
#include "helper.h"

template <typename MUON, typename TRK>
inline void FillTagBranches(const MUON &muon,
                            const std::vector<TRK> &tracks,
                            NtupleContent &nt,
                            const reco::Vertex &vertex) {

  nt.branches["tag_pt"] = (float)muon.pt();
  nt.branches["tag_eta"] = (float)muon.eta();
  nt.branches["tag_phi"] = (float)muon.phi();
  nt.branches["tag_charge"] = (int)muon.charge();
  nt.branches["tag_dxy"] = (float)muon.innerTrack()->dxy(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
  nt.branches["tag_dz"] = (float)muon.innerTrack()->dz(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
  nt.branches["tag_isPF"] = (bool)muon.isPFMuon();
  nt.branches["tag_isSA"] = (bool)muon.isStandAloneMuon();
  nt.branches["tag_isTracker"] = (bool)muon.isTrackerMuon();
  nt.branches["tag_isGlobal"] = (bool)muon.isGlobalMuon();
  nt.branches["tag_isRPC"] = (bool)muon.isRPCMuon();
  
  // Use selectors instead of 'muon.passed' method which is only introduced in CMSSW_9_4_X
  nt.branches["tag_isLoose"] = (bool)muon::isLooseMuon(muon);
  nt.branches["tag_isMedium"] = (bool)muon::isMediumMuon(muon);
  nt.branches["tag_isTight"] = (bool)muon::isTightMuon(muon, vertex);
  nt.branches["tag_isSoft"] = (bool)muon::isSoftMuon(muon, vertex, false);
  nt.branches["tag_isHighPt"] = (bool)muon::isHighPtMuon(muon, vertex);
  float Trkiso04 = (TrackerEnergy04<TRK>(muon.eta(), muon.phi(), tracks) - muon.pt());
  nt.branches["tag_absTrkIso04"] = (float)(Trkiso04 > 0) ? Trkiso04 : 0;
  float Trkiso03 = (TrackerEnergy03<TRK>(muon.eta(), muon.phi(), tracks) - muon.pt());
  nt.branches["tag_absTrkIso03"] = (Trkiso03 > 0) ? Trkiso03 : 0;
  nt.branches["tag_iso03_sumPt"] = (float)muon.isolationR03().sumPt;
  nt.branches["tag_pfIso03_charged"] = (float)muon.pfIsolationR03().sumChargedHadronPt;
  nt.branches["tag_pfIso03_neutral"] = (float)muon.pfIsolationR03().sumNeutralHadronEt;
  nt.branches["tag_pfIso03_photon"] = (float)muon.pfIsolationR03().sumPhotonEt;
  nt.branches["tag_pfIso03_sumPU"] = (float)muon.pfIsolationR03().sumPUPt;
  nt.branches["tag_combRelIsoPF03dBeta"] = (float)((muon.pfIsolationR03().sumChargedHadronPt + TMath::Max(muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - muon.pfIsolationR03().sumPUPt/2.0,0.0))/muon.pt());
  nt.branches["tag_pfIso04_charged"] = (float)muon.pfIsolationR04().sumChargedHadronPt;
  nt.branches["tag_pfIso04_neutral"] = (float)muon.pfIsolationR04().sumNeutralHadronEt;
  nt.branches["tag_pfIso04_photon"] = (float)muon.pfIsolationR04().sumPhotonEt;
  nt.branches["tag_pfIso04_sumPU"] = (float)muon.pfIsolationR04().sumPUPt;
  nt.branches["tag_combRelIsoPF04dBeta"] = (float)((muon.pfIsolationR04().sumChargedHadronPt + TMath::Max(muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - muon.pfIsolationR04().sumPUPt/2.0,0.0))/muon.pt());
  if (muon.tunePMuonBestTrack().isNonnull()) {
    nt.branches["tag_tuneP_ExistingRefit"] = (bool)true;
    nt.branches["tag_tuneP_charge"] = (int)muon.tunePMuonBestTrack()->charge();
    nt.branches["tag_tuneP_pt"] = (float)muon.tunePMuonBestTrack()->pt();
    nt.branches["tag_tuneP_pterr"] = (float)muon.tunePMuonBestTrack()->ptError();
    nt.branches["tag_tuneP_eta"] = (float)muon.tunePMuonBestTrack()->eta();
    nt.branches["tag_tuneP_phi"] = (float)muon.tunePMuonBestTrack()->phi();
    nt.branches["tag_tuneP_muonHits"] = (int)muon.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits();
  } else {
    nt.branches["tag_tuneP_ExistingRefit"] = (bool)false;
    nt.branches["tag_tuneP_charge"] = (int)-99;
    nt.branches["tag_tuneP_pt"] = (float)-99.;
    nt.branches["tag_tuneP_pterr"] = (float)-99.;
    nt.branches["tag_tuneP_eta"] = (float)-99.;
    nt.branches["tag_tuneP_phi"] = (float)-99.;
    nt.branches["tag_tuneP_pterr"] = (float)-99.;
    nt.branches["tag_tuneP_muonHits"] = (int)-99;
  }
  int nsegments = 0;
  for (auto &chamber : muon.matches()) {
    if (chamber.id.det() != DetId::Muon)
      continue;
    if (chamber.id.subdetId() != MuonSubdetId::DT && chamber.id.subdetId() != MuonSubdetId::CSC)
      continue;
    nsegments += chamber.segmentMatches.size();
  }
  nt.branches["tag_nsegments"] = (int)nsegments;

  // for high-pt
  if (muon.innerTrack().isNonnull() && muon.innerTrack().isAvailable()) {
    nt.branches["tag_inner_validFraction"] = (float)muon.innerTrack()->validFraction();
    nt.branches["tag_inner_trackerLayers"] = (int)muon.innerTrack()->hitPattern().trackerLayersWithMeasurement();
    nt.branches["tag_inner_pixelLayers"] = (int)muon.innerTrack()->hitPattern().pixelLayersWithMeasurement();
    nt.branches["tag_inner_pterr"] = (float)muon.innerTrack()->ptError();
    nt.branches["tag_inner_pixelHits"] = (int)muon.innerTrack()->hitPattern().numberOfValidPixelHits();
    nt.branches["tag_inner_pt"] = (float)muon.innerTrack()->pt();
    nt.branches["tag_inner_eta"] = (float)muon.innerTrack()->eta();
    nt.branches["tag_inner_phi"] = (float)muon.innerTrack()->phi();
    nt.branches["tag_inner_charge"] = (int)muon.innerTrack()->charge();
  } else {
    nt.branches["tag_inner_validFraction"] = (float)-99;
    nt.branches["tag_inner_trackerLayers"] = (int)-99;
    nt.branches["tag_inner_pixelLayers"] = (int)-99;
    nt.branches["tag_inner_pterr"] = (float)-99;
    nt.branches["tag_inner_pixelHits"] = (int)-99;
    nt.branches["tag_inner_pt"] = (float)-99;
    nt.branches["tag_inner_eta"] = (float)-99;
    nt.branches["tag_inner_phi"] = (float)-99;
    nt.branches["tag_inner_charge"] = (int)0;
  }

  if (muon.tpfmsTrack().isNonnull() && muon.tpfmsTrack().isAvailable()) {
    nt.branches["tag_tpfms_charge"] = (int)muon.tpfmsTrack()->charge();
    nt.branches["tag_tpfms_pt"] = (float)muon.tpfmsTrack()->pt();
    nt.branches["tag_tpfms_pterr"] = (float)muon.tpfmsTrack()->ptError();
    nt.branches["tag_tpfms_eta"] = (float)muon.tpfmsTrack()->eta();
    nt.branches["tag_tpfms_phi"] = (float)muon.tpfmsTrack()->phi();
    nt.branches["tag_tpfms_muonHits"] = (int)muon.tpfmsTrack()->hitPattern().numberOfValidMuonHits();
  } else {
    nt.branches["tag_tpfms_charge"] = (int)-99;
    nt.branches["tag_tpfms_pt"] = (float)-99.;
    nt.branches["tag_tpfms_pterr"] = (float)-99.;
    nt.branches["tag_tpfms_eta"] = (float)-99.;
    nt.branches["tag_tpfms_phi"] = (float)-99.;
    nt.branches["tag_tpfms_pterr"] = (float)-99.;
    nt.branches["tag_tpfms_muonHits"] = (int)-99.;
  }

  if (muon.pickyTrack().isNonnull() && muon.pickyTrack().isAvailable()) {
    nt.branches["tag_picky_charge"] = (int)muon.pickyTrack()->charge();
    nt.branches["tag_picky_pt"] = (float)muon.pickyTrack()->pt();
    nt.branches["tag_picky_pterr"] = (float)muon.pickyTrack()->ptError();
    nt.branches["tag_picky_eta"] = (float)muon.pickyTrack()->eta();
    nt.branches["tag_picky_phi"] = (float)muon.pickyTrack()->phi();
    nt.branches["tag_picky_muonHits"] = (int)muon.pickyTrack()->hitPattern().numberOfValidMuonHits();
  } else {
    nt.branches["tag_picky_charge"] = (int)-99;
    nt.branches["tag_picky_pt"] = (float)-99.;
    nt.branches["tag_picky_pterr"] = (float)-99.;
    nt.branches["tag_picky_eta"] = (float)-99.;
    nt.branches["tag_picky_phi"] = (float)-99.;
    nt.branches["tag_picky_pterr"] = (float)-99.;
    nt.branches["tag_picky_muonHits"] = (int)-99.;
  }

  if (muon.dytTrack().isNonnull() && muon.dytTrack().isAvailable()) {
    nt.branches["tag_dyt_charge"] = (int)muon.dytTrack()->charge();
    nt.branches["tag_dyt_pt"] = (float)muon.dytTrack()->pt();
    nt.branches["tag_dyt_pterr"] = (float)muon.dytTrack()->ptError();
    nt.branches["tag_dyt_eta"] = (float)muon.dytTrack()->eta();
    nt.branches["tag_dyt_phi"] = (float)muon.dytTrack()->phi();
    nt.branches["tag_dyt_muonHits"] = (int)muon.dytTrack()->hitPattern().numberOfValidMuonHits();
  } else {
    nt.branches["tag_dyt_charge"] = (int)-99;
    nt.branches["tag_dyt_pt"] = (float)-99.;
    nt.branches["tag_dyt_pterr"] = (float)-99.;
    nt.branches["tag_dyt_eta"] = (float)-99.;
    nt.branches["tag_dyt_phi"] = (float)-99.;
    nt.branches["tag_dyt_pterr"] = (float)-99.;
    nt.branches["tag_dyt_muonHits"] = (int)-99.;
  }


  if (muon.globalTrack().isNonnull()) {
    nt.branches["tag_GlobalValidHits"] = (int)muon.globalTrack()->hitPattern().numberOfValidMuonHits();
  }
  nt.branches["tag_ZprimeMatchedStations"] = (bool)(muon.numberOfMatchedStations() > 1 || (muon.numberOfMatchedStations() == 1 && !(muon.stationMask() == 1 || muon.stationMask() == 16)) || (muon.numberOfMatchedStations() == 1 && (muon.stationMask() == 1 || muon.stationMask() == 16) && muon.numberOfMatchedRPCLayers() > 2));
    
  nt.branches["tag_RPCLayers"] = (int)muon.numberOfMatchedRPCLayers();
}

template <typename MUON, typename TRK>
inline void FillProbeBranches(
    const MUON &mu, const std::vector<TRK> &tracks, NtupleContent &nt, bool success, const reco::Vertex &vertex) {
  nt.branches["probe_pt"] = (float)mu.pt();
  nt.branches["probe_eta"] = (float)mu.eta();
  nt.branches["probe_phi"] = (float)mu.phi();
  nt.branches["probe_charge"] = (int)mu.charge();
  float Trkiso04 = (TrackerEnergy04<TRK>(mu.eta(), mu.phi(), tracks) - mu.pt());
  nt.branches["probe_absTrkIso04"] = (float)(Trkiso04 > 0) ? Trkiso04 : 0;
  float Trkiso03 = (TrackerEnergy03<TRK>(mu.eta(), mu.phi(), tracks) - mu.pt());
  nt.branches["probe_absTrkIso03"] = (float)(Trkiso03 > 0) ? Trkiso03 : 0;
  // success --> muon obj and track match in dR
  if (success) {
    nt.branches["probe_isLoose"] = (bool)muon::isLooseMuon(mu);
    nt.branches["probe_isMedium"] = (bool)muon::isMediumMuon(mu);
    nt.branches["probe_isTight"] = (bool)muon::isTightMuon(mu, vertex);
    nt.branches["probe_isSoft"] = (bool)muon::isSoftMuon(mu, vertex, false);
    nt.branches["probe_isHighPt"] = (bool)muon::isHighPtMuon(mu, vertex);
    nt.branches["probe_isArbitratedTracker"] = (bool)muon::isGoodMuon(mu, muon::TrackerMuonArbitrated);
    nt.branches["probe_isPF"] = (bool)mu.isPFMuon();
    nt.branches["probe_isSA"] = (bool)mu.isStandAloneMuon();
    nt.branches["probe_isTracker"] = (bool)mu.isTrackerMuon();
    nt.branches["probe_isGlobal"] = (bool)mu.isGlobalMuon();
    nt.branches["probe_isRPC"] = (bool)mu.isRPCMuon();
    nt.branches["probe_iso03_sumPt"] = (float)mu.isolationR03().sumPt;
    nt.branches["probe_pfIso03_charged"] = (float)mu.pfIsolationR03().sumChargedHadronPt;
    nt.branches["probe_pfIso03_neutral"] = (float)mu.pfIsolationR03().sumNeutralHadronEt;
    nt.branches["probe_pfIso03_photon"] = (float)mu.pfIsolationR03().sumPhotonEt;
    nt.branches["probe_pfIso03_sumPU"] = (float)mu.pfIsolationR03().sumPUPt;
    nt.branches["probe_pfIso04_charged"] = (float)mu.pfIsolationR04().sumChargedHadronPt;
    nt.branches["probe_pfIso04_neutral"] = (float)mu.pfIsolationR04().sumNeutralHadronEt;
    nt.branches["probe_pfIso04_photon"] = (float)mu.pfIsolationR04().sumPhotonEt;
    nt.branches["probe_pfIso04_sumPU"] = (float)mu.pfIsolationR04().sumPUPt;
    nt.branches["probe_matchedStations"] = (int)mu.numberOfMatchedStations();
    nt.branches["probe_expectedMatchedStations"] = (int)mu.expectedNnumberOfMatchedStations();
    nt.branches["probe_RPCLayers"] = (int)mu.numberOfMatchedRPCLayers();
    nt.branches["probe_stationMask"] = (int)mu.stationMask();
    nt.branches["probe_nShowers"] = (int)mu.numberOfShowers();
    if (mu.globalTrack().isNonnull()) {
      nt.branches["probe_muonHits"] = (int)mu.globalTrack()->hitPattern().numberOfValidMuonHits();
      nt.branches["probe_trkChi2"] = (float)mu.globalTrack()->normalizedChi2();
    } else if (mu.innerTrack().isNonnull() && mu.innerTrack().isAvailable()) {
      nt.branches["probe_trkChi2"] = (float)mu.innerTrack()->normalizedChi2();
      nt.branches["probe_muonHits"] = (int)mu.innerTrack()->hitPattern().numberOfValidMuonHits();
    } else {    
      nt.branches["probe_muonHits"] = (int)-99;
      nt.branches["probe_trkChi2"] = (float)-99;
    }
    if (mu.innerTrack().isNonnull() && mu.innerTrack().isAvailable()) {
      nt.branches["probe_inner_validFraction"] = (float)mu.innerTrack()->validFraction();
      nt.branches["probe_inner_trackerLayers"] = (int)mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
      nt.branches["probe_inner_pixelLayers"] = (int)mu.innerTrack()->hitPattern().pixelLayersWithMeasurement();
      nt.branches["probe_inner_pterr"] = (float)mu.innerTrack()->ptError();
      nt.branches["probe_dxy"] = (float)mu.innerTrack()->dxy(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
      nt.branches["probe_dz"] = (float)mu.innerTrack()->dz(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
      nt.branches["probe_inner_pixelHits"] = (int)mu.innerTrack()->hitPattern().numberOfValidPixelHits();
      nt.branches["probe_inner_pt"] = (float)mu.innerTrack()->pt();
      nt.branches["probe_inner_eta"] = (float)mu.innerTrack()->eta();
      nt.branches["probe_inner_phi"] = (float)mu.innerTrack()->phi();
      nt.branches["probe_inner_charge"] = (int)mu.innerTrack()->charge();
    } else {
      nt.branches["probe_inner_validFraction"] = (int)-99;
      nt.branches["probe_inner_trackerLayers"] = (int)-99;
      nt.branches["probe_inner_pixelLayers"] = (int)-99;
      nt.branches["probe_inner_pterr"] = (float)-99;
      nt.branches["probe_dxy"] = (float)-99;
      nt.branches["probe_dz"] = (float)-99;
      nt.branches["probe_inner_pixelHits"] = (int)-99;
      nt.branches["probe_inner_pt"] = (float)-99;
      nt.branches["probe_inner_eta"] = (float)-99;
      nt.branches["probe_inner_phi"] = (float)-99;
      nt.branches["probe_inner_charge"] = (int)-99;
    }
    if (mu.outerTrack().isNonnull() && mu.outerTrack().isAvailable()) {
      nt.branches["probe_muonStations"] = (int)mu.outerTrack()->hitPattern().muonStationsWithValidHits();
      nt.branches["probe_DTHits"] = (int)mu.outerTrack()->hitPattern().numberOfValidMuonDTHits();
      nt.branches["probe_CSCHits"] = (int)mu.outerTrack()->hitPattern().numberOfValidMuonCSCHits();
      nt.branches["probe_outer_pt"] = (float)mu.outerTrack()->pt();
      nt.branches["probe_outer_eta"] = (float)mu.outerTrack()->eta();
      nt.branches["probe_outer_phi"] = (float)mu.outerTrack()->phi();
      nt.branches["probe_outer_charge"] = (int)mu.outerTrack()->charge();
    }else{
      nt.branches["probe_muonStations"] = (int)-99;
      nt.branches["probe_DTHits"] = (int)-99;
      nt.branches["probe_CSCHits"] = (int)-99;
      nt.branches["probe_outer_pt"] = (float)-99;
      nt.branches["probe_outer_eta"] = (float)-99;
      nt.branches["probe_outer_phi"] = (float)-99;
      nt.branches["probe_outer_charge"] = (int)-99;
    }
    if (mu.globalTrack().isNonnull() && mu.globalTrack().isAvailable()) {
      nt.branches["probe_global_pt"] = (float)mu.globalTrack()->pt();
      nt.branches["probe_global_eta"] = (float)mu.globalTrack()->eta();
      nt.branches["probe_global_phi"] = (float)mu.globalTrack()->phi();
      nt.branches["probe_global_charge"] = (int)mu.globalTrack()->charge();
    }else{
      nt.branches["probe_global_pt"] = (float)-99;
      nt.branches["probe_global_eta"] = (float)-99;
      nt.branches["probe_global_phi"] = (float)-99;
      nt.branches["probe_global_charge"] = (int)-99;
    }
    if (mu.muonBestTrack().isNonnull() && mu.muonBestTrack().isAvailable()) {
      nt.branches["probe_best_pt"] = (float)mu.muonBestTrack()->pt();
      nt.branches["probe_best_eta"] = (float)mu.muonBestTrack()->eta();
      nt.branches["probe_best_phi"] = (float)mu.muonBestTrack()->phi();
      nt.branches["probe_best_charge"] = (int)mu.muonBestTrack()->charge();
    }else{
      nt.branches["probe_best_pt"] = (float)-99;
      nt.branches["probe_best_eta"] = (float)-99;
      nt.branches["probe_best_phi"] = (float)-99;
      nt.branches["probe_best_charge"]   = (int)-99;
    }
    
    if (mu.tunePMuonBestTrack().isNonnull() && mu.tunePMuonBestTrack().isAvailable()) {
      nt.branches["probe_tuneP_ExistingRefit"] = true;
      nt.branches["probe_tuneP_charge"] = (int)mu.tunePMuonBestTrack()->charge();
      nt.branches["probe_tuneP_pt"] = (float)mu.tunePMuonBestTrack()->pt();
      nt.branches["probe_tuneP_eta"] = (float)mu.tunePMuonBestTrack()->eta();
      nt.branches["probe_tuneP_phi"] = (float)mu.tunePMuonBestTrack()->phi();
      nt.branches["probe_tuneP_pterr"] = (float)mu.tunePMuonBestTrack()->ptError();
      nt.branches["probe_tuneP_muonHits"] = (int)mu.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits();
    }else{
      nt.branches["probe_tuneP_ExistingRefit"] = false;
      nt.branches["probe_tuneP_charge"] = (int)-99;
      nt.branches["probe_tuneP_pt"] = (float)-99;
      nt.branches["probe_tuneP_eta"] = (float)-99;
      nt.branches["probe_tuneP_phi"] = (float)-99;
      nt.branches["probe_tuneP_pterr"] = (float)-99;
      nt.branches["probe_tuneP_muonHits"] = (int)-99;
    }
    if (mu.tpfmsTrack().isNonnull() && mu.tpfmsTrack().isAvailable()) {
      nt.branches["probe_tpfms_charge"] = (int)mu.tpfmsTrack()->charge();
      nt.branches["probe_tpfms_pt"] = (float)mu.tpfmsTrack()->pt();
      nt.branches["probe_tpfms_pterr"] = (float)mu.tpfmsTrack()->ptError();
      nt.branches["probe_tpfms_eta"] = (float)mu.tpfmsTrack()->eta();
      nt.branches["probe_tpfms_phi"] = (float)mu.tpfmsTrack()->phi();
      nt.branches["probe_tpfms_muonHits"] = (int)mu.tpfmsTrack()->hitPattern().numberOfValidMuonHits();
    }else{
      nt.branches["probe_tpfms_charge"] = (int)-99.;
      nt.branches["probe_tpfms_pt"] = (float)-99.;
      nt.branches["probe_tpfms_pterr"] = (float)-99.;
      nt.branches["probe_tpfms_eta"] = (float)-99.;
      nt.branches["probe_tpfms_phi"] = (float)-99.;
      nt.branches["probe_tpfms_muonHits"]  = (int)-99.;
    }
    if (mu.pickyTrack().isNonnull() && mu.pickyTrack().isAvailable()) {
      nt.branches["probe_picky_charge"] = (int)mu.pickyTrack()->charge();
      nt.branches["probe_picky_pt"] = (float)mu.pickyTrack()->pt();
      nt.branches["probe_picky_pterr"] = (float)mu.pickyTrack()->ptError();
      nt.branches["probe_picky_eta"] = (float)mu.pickyTrack()->eta();
      nt.branches["probe_picky_phi"] = (float)mu.pickyTrack()->phi();
      nt.branches["probe_picky_muonHits"] = (int)mu.pickyTrack()->hitPattern().numberOfValidMuonHits();
    }else{
      nt.branches["probe_picky_charge"] = (int)-99.;
      nt.branches["probe_picky_pt"] = (float)-99.;
      nt.branches["probe_picky_pterr"] = (float)-99.;
      nt.branches["probe_picky_eta"] = (float)-99.;
      nt.branches["probe_picky_phi"] = (float)-99.;
      nt.branches["probe_picky_muonHits"] = (int)-99.;
    }
    if (mu.dytTrack().isNonnull() && mu.dytTrack().isAvailable()) {
      nt.branches["probe_dyt_charge"] = (int)mu.dytTrack()->charge();
      nt.branches["probe_dyt_pt"] = (float)mu.dytTrack()->pt();
      nt.branches["probe_dyt_pterr"] = (float)mu.dytTrack()->ptError();
      nt.branches["probe_dyt_eta"] = (float)mu.dytTrack()->eta();
      nt.branches["probe_dyt_phi"] = (float)mu.dytTrack()->phi();
      nt.branches["probe_dyt_muonHits"] = (int)mu.dytTrack()->hitPattern().numberOfValidMuonHits();
    }else{
      nt.branches["probe_dyt_charge"] = (int)-99.;
      nt.branches["probe_dyt_pt"] = (float)-99.;
      nt.branches["probe_dyt_pterr"] = (float)-99.;
      nt.branches["probe_dyt_eta"] = (float)-99.;
      nt.branches["probe_dyt_phi"] = (float)-99.;
      nt.branches["probe_dyt_muonHits"] = (int)-99;
    }
    nt.branches["probe_positionChi2"] = (float)mu.combinedQuality().chi2LocalPosition;
    nt.branches["probe_trkKink"] = (float)mu.combinedQuality().trkKink;
    nt.branches["probe_segmentCompatibility"] = (float)muon::segmentCompatibility(mu);
    nt.branches["probe_isMuMatched"] = (bool)true;
    int nsegments = 0;
    for (auto &chamber : mu.matches()) {
      if (chamber.id.det() != DetId::Muon)
        continue;
      if (chamber.id.subdetId() != MuonSubdetId::DT && chamber.id.subdetId() != MuonSubdetId::CSC)
        continue;
      nsegments += chamber.segmentMatches.size();
    }
    nt.branches["probe_nsegments"] = (int)nsegments;
  }
  // no successs (no match)
  else {
    nt.branches["probe_isLoose"] = (bool)false;
    nt.branches["probe_isMedium"] = (bool)false;
    nt.branches["probe_isTight"] = (bool)false;
    nt.branches["probe_isSoft"] = (bool)false;
    nt.branches["probe_isHighPt"] = (bool)false;
    nt.branches["probe_isMuMatched"] = (bool)false;
    nt.branches["probe_isPF"] = (bool)false;
    nt.branches["probe_isSA"] = (bool)false;
    nt.branches["probe_isTracker"] = (bool)false;
    nt.branches["probe_isGlobal"] = (bool)false;
    nt.branches["probe_inner_validFraction"] = (int)-99;
    nt.branches["probe_trkChi2"] = (float)-99;
    nt.branches["probe_positionChi2"] = (float)-99;
    nt.branches["probe_trkKink"] = (float)-99;
    nt.branches["probe_inner_trackerLayers"] = (int)-99;
    nt.branches["probe_inner_pixelLayers"] = (int)-99;
    nt.branches["probe_dxy"] = (float)-99;
    nt.branches["probe_dz"] = (float)-99;
    nt.branches["probe_muonStations"] = (int)-99;
    nt.branches["probe_muonHits"] = (int)-99;
    nt.branches["probe_DTHits"] = (int)-99;
    nt.branches["probe_CSCHits"] = (int)-99;
    nt.branches["probe_inner_pterr"] = (float)-99;
    nt.branches["probe_iso03_sumPt"] = (float)-99;
    nt.branches["probe_pfIso03_charged"] = (float)-99;
    nt.branches["probe_pfIso03_neutral"] = (float)-99;
    nt.branches["probe_pfIso03_photon"] = (float)-99;
    nt.branches["probe_pfIso03_sumPU"] = (float)-99;
    nt.branches["probe_pfIso04_charged"] = (float)-99;
    nt.branches["probe_pfIso04_neutral"] = (float)-99;
    nt.branches["probe_pfIso04_photon"] = (float)-99;
    nt.branches["probe_pfIso04_sumPU"] = (float)-99;
    nt.branches["probe_inner_pixelHits"] = (int)-99;
    nt.branches["probe_matchedStations"] = (int)-99;
    nt.branches["probe_expectedMatchedStations"] = (int)-99;
    nt.branches["probe_RPCLayers"] = (int)-99;
    nt.branches["probe_stationMask"] = (int)0;
    nt.branches["probe_nShowers"] = (int)-99;
    nt.branches["probe_tuneP_pt"] = (float)-99;
    nt.branches["probe_tuneP_pterr"] = (float)-99;
    nt.branches["probe_tuneP_muonHits"] = (int)-99;
    nt.branches["probe_nsegments"] = (int)-99;
    nt.branches["l1pt"] = (float)-99;
    nt.branches["l1q"] = (int)-99;
    nt.branches["l1dr"] = (float)-99;
    nt.branches["l1ptByQ"] = (float)-99;
    nt.branches["l1qByQ"] = (float)-99;
    nt.branches["l1drByQ"] = (float)-99;
    nt.branches["probe_inner_pt"] = (float)-99;
    nt.branches["probe_inner_eta"] = (float)-99;
    nt.branches["probe_inner_phi"] = (float)-99;
    nt.branches["probe_inner_charge"] = (int)-99;
    nt.branches["probe_outer_pt"] = (float)-99;
    nt.branches["probe_outer_eta"] = (float)-99;
    nt.branches["probe_outer_phi"] = (float)-99;
    nt.branches["probe_outer_charge"] = (int)-99;
    nt.branches["probe_global_pt"] = (float)-99;
    nt.branches["probe_global_eta"] = (float)-99;
    nt.branches["probe_global_phi"] = (float)-99;
    nt.branches["probe_global_charge"] = (int)-99;
    nt.branches["probe_best_pt"] = (float)-99;
    nt.branches["probe_best_eta"] = (float)-99;
    nt.branches["probe_best_phi"] = (float)-99;
    nt.branches["probe_best_charge"] = (int)-99;


    nt.branches["probe_SIP3D"] = (float)-99.0;
    nt.branches["probe_SIP3D_err"] = (float)-99.0;
    nt.branches["probe_dyt_eta"] = (float)-99.0;
    nt.branches["probe_dyt_phi"] = (float)-99.0;
    nt.branches["probe_dyt_pt"] = (float)-99.0;
    nt.branches["probe_dyt_pterr"] = (float)-99.0;
    nt.branches["probe_dyt_muonHits"] = (int)-99;
    nt.branches["probe_picky_eta"] = (float)-99.0;
    nt.branches["probe_picky_phi"] = (float)-99.0;
    nt.branches["probe_picky_pt"] = (float)-99.0;
    nt.branches["probe_picky_pterr"] = (float)-99.0;
    nt.branches["probe_tpfms_eta"] = (float)-99.0;
    nt.branches["probe_tpfms_phi"] = (float)-99.0;
    nt.branches["probe_tpfms_pt"] = (float)-99.0;
    nt.branches["probe_tpfms_pterr"] = (float)-99.0;
    nt.branches["probe_tuneP_ExistingRefit"] = (float)-99.0;
    nt.branches["probe_tuneP_charge"] = (int)-99;
    nt.branches["probe_tuneP_eta"] = (float)-99.0;
    nt.branches["probe_tuneP_phi"] = (float)-99.0;
    nt.branches["probe_cosmic_minDR"] = (float)+99.0;
    nt.branches["probe_dsa_minDR"] = (float)+99.0;
    nt.branches["probe_isArbitratedTracker"] = (bool)false;
    nt.branches["probe_isRPC"] = (bool)false;
    nt.branches["probe_muonHits"] = (int)-99;
    nt.branches["probe_picky_muonHits"] = (int)-99;
    nt.branches["probe_tpfms_muonHits"] = (int)-99;
    nt.branches["probe_dgl_outerTrackerHits"] = (int)-99;
    nt.branches["probe_dgl_pterr"] = (float)0.0;
    nt.branches["probe_dgl_segmentMatches"] = (float)-99.0;
    nt.branches["probe_dgl_totalHits"] = (int)-99;
    nt.branches["probe_dgl_trackerHits"] = (int)-99;
    nt.branches["probe_dgl_minDR"] = (float)+99.0;
  }
}

template <typename MUON>
inline void FillProbeBranchesSelector(const MUON &mu,
                                      NtupleContent &nt,
                                      const std::vector<unsigned> selectorBits,
                                      bool success) {
  for (unsigned int ibit = 0; ibit < selectorBits.size(); ++ibit) {
    if (success && mu.selectors() != 0) {
      unsigned bit = selectorBits.at(ibit);
      nt.probe_selectors[ibit] = mu.passed(1UL << bit);
    } else {
      nt.probe_selectors[ibit] = false;
    }
  }
}

template <typename TRK>
inline void FillProbeBranchesdSA(const TRK &trk, NtupleContent &nt, bool passdSA) {
  nt.branches["probe_isdSA"] = (bool)passdSA;
  nt.branches["probe_dsa_pt"] = (float)trk.pt();
  nt.branches["probe_dsa_eta"] = (float)trk.eta();
  nt.branches["probe_dsa_phi"] = (float)trk.phi();
  nt.branches["probe_dsa_charge"] = (int)trk.charge();
  if (passdSA) {
    nt.branches["probe_dsa_outerEta"] = (float)trk.outerEta();
    nt.branches["probe_dsa_outerPhi"] = (float)trk.outerPhi();
    nt.branches["probe_dsa_dxy"] = (float)trk.dxy(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
    nt.branches["probe_dsa_dz"] = (float)trk.dz(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
    nt.branches["probe_dsa_muonStations"] = (int)trk.hitPattern().muonStationsWithValidHits();
    nt.branches["probe_dsa_muonHits"] = (int)trk.hitPattern().numberOfValidMuonHits();
    nt.branches["probe_dsa_DTHits"] = (int)trk.hitPattern().numberOfValidMuonDTHits();
    nt.branches["probe_dsa_CSCHits"] = (int)trk.hitPattern().numberOfValidMuonCSCHits();
    nt.branches["probe_dsa_pterr"] = (float)(trk.ptError() / trk.pt());
    nt.branches["probe_dsa_trkChi2"] = (float)trk.normalizedChi2();
    // [Adapted from displaced dimuon analysis]
    // Number of DT+CSC segments
    unsigned int nsegments = 0;
    for (auto &hit : trk.recHits()) {
      if (!hit->isValid())
        continue;
      DetId id = hit->geographicalId();
      if (id.det() != DetId::Muon)
        continue;
      if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC)
        nsegments++;
    }
    nt.branches["probe_dsa_nsegments"] = (int)nsegments;
  } else {   
    nt.branches["probe_dsa_outerEta"] = (float)-99;
    nt.branches["probe_dsa_outerPhi"] = (float)-99;
    nt.branches["probe_dsa_dxy"] = (float)-99;
    nt.branches["probe_dsa_dz"] = (float)-99;
    nt.branches["probe_dsa_muonStations"] = (int)-99;
    nt.branches["probe_dsa_muonHits"] = (int)-99;
    nt.branches["probe_dsa_DTHits"] = (int)-99;
    nt.branches["probe_dsa_CSCHits"] = (int)-99;
    nt.branches["probe_dsa_pterr"] = (float)-99;
    nt.branches["probe_dsa_trkChi2"] = (float)-99;
    nt.branches["probe_dsa_nsegments"] = (int)-99;
  }
}



template <typename TRK>
inline void FillTagBranchesdSA(const TRK &trk, NtupleContent &nt, bool passdSA) {
  nt.branches["tag_isdSA"] = (bool)passdSA;
  nt.branches["tag_dsa_pt"] = (float)trk.pt();
  nt.branches["tag_dsa_eta"] = (float)trk.eta();
  nt.branches["tag_dsa_phi"] = (float)trk.phi();
  nt.branches["tag_dsa_charge"] = (int)trk.charge();
  if (passdSA) {
    nt.branches["tag_dsa_outerEta"] = (float)trk.outerEta();
    nt.branches["tag_dsa_outerPhi"] = (float)trk.outerPhi();
    nt.branches["tag_dsa_dxy"] = (float)trk.dxy(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
    nt.branches["tag_dsa_dz"] = (float)trk.dz(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
    nt.branches["tag_dsa_muonStations"] = (int)trk.hitPattern().muonStationsWithValidHits();
    nt.branches["tag_dsa_muonHits"] = (int)trk.hitPattern().numberOfValidMuonHits();
    nt.branches["tag_dsa_DTHits"] = (int)trk.hitPattern().numberOfValidMuonDTHits();
    nt.branches["tag_dsa_CSCHits"] = (int)trk.hitPattern().numberOfValidMuonCSCHits();
    nt.branches["tag_dsa_pterr"] = (float)(trk.ptError() / trk.pt());
    nt.branches["tag_dsa_trkChi2"] = (float)trk.normalizedChi2();
    // [Adapted from displaced dimuon analysis]
    // Number of DT+CSC segments
    unsigned int nsegments = 0;
    for (auto &hit : trk.recHits()) {
      if (!hit->isValid())
        continue;
      DetId id = hit->geographicalId();
      if (id.det() != DetId::Muon)
        continue;
      if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC)
        nsegments++;
    }
    nt.branches["tag_dsa_nsegments"] = (int)nsegments;
  } else {
    nt.branches["tag_dsa_outerEta"] = (float)-99;
    nt.branches["tag_dsa_outerPhi"] = (float)-99;
    nt.branches["tag_dsa_dxy"] = (float)-99;
    nt.branches["tag_dsa_dz"] = (float)-99;
    nt.branches["tag_dsa_muonStations"] = (int)-99;
    nt.branches["tag_dsa_muonHits"] = (int)-99;
    nt.branches["tag_dsa_DTHits"] = (int)-99;
    nt.branches["tag_dsa_CSCHits"] = (int)-99;
    nt.branches["tag_dsa_pterr"] = (float)-99;
    nt.branches["tag_dsa_trkChi2"] = (float)-99;
    nt.branches["tag_dsa_nsegments"] = (int)-99;
  }
}

template <typename TRK>
inline void FillProbeBranchesdgl(const TRK &trk, NtupleContent &nt, bool passdgl) {
  nt.branches["probe_isdGlobal"] = (bool)passdgl;
  nt.branches["probe_dgl_pt"] = (float)trk.pt();
  nt.branches["probe_dgl_eta"] = (float)trk.eta();
  nt.branches["probe_dgl_phi"] = (float)trk.phi();
  nt.branches["probe_dgl_charge"] = (int)trk.charge();
  if (passdgl) {
    nt.branches["probe_dgl_dxy"] = (float)trk.dxy(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
    nt.branches["probe_dgl_dz"] = (float)trk.dz(reco::TrackBase::Point(std::get<float>(nt.branches["pv_x"].value), std::get<float>(nt.branches["pv_y"].value), std::get<float>(nt.branches["pv_z"].value)));
    nt.branches["probe_dgl_muonStations"] = (int)trk.hitPattern().muonStationsWithValidHits();
    nt.branches["probe_dgl_muonHits"] = (int)trk.hitPattern().numberOfValidMuonHits();
    nt.branches["probe_dgl_outerTrackerHits"] = (int)trk.hitPattern().numberOfValidStripHits();
    nt.branches["probe_dgl_trackerHits"] = (int)trk.hitPattern().numberOfValidTrackerHits();
    nt.branches["probe_dgl_totalHits"] = (int)trk.hitPattern().numberOfValidHits();
    nt.branches["probe_dgl_DTHits"] = (int)trk.hitPattern().numberOfValidMuonDTHits();
    nt.branches["probe_dgl_CSCHits"] = (int)trk.hitPattern().numberOfValidMuonCSCHits();
    nt.branches["probe_dgl_pterr"] = (float)(trk.ptError() / trk.pt());
    nt.branches["probe_dgl_trkChi2"] = (float)trk.normalizedChi2();
    // [Adapted from displaced dimuon analysis]
    // Number of DT+CSC segments
    unsigned int nsegments = 0;
    for (auto &hit : trk.recHits()) {
      if (!hit->isValid())
        continue;
      DetId id = hit->geographicalId();
      if (id.det() != DetId::Muon)
        continue;
      if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC)
        nsegments++;
    }
    nt.branches["probe_dgl_nsegments"] = (int)nsegments;
  } else {
    nt.branches["probe_dgl_dxy"] = (float)-99;
    nt.branches["probe_dgl_dz"] = (float)-99;
    nt.branches["probe_dgl_muonStations"] = (int)-99;
    nt.branches["probe_dgl_muonHits"] = (int)-99;
    nt.branches["probe_dgl_DTHits"] = (int)-99;
    nt.branches["probe_dgl_CSCHits"] = (int)-99;
    nt.branches["probe_dgl_pterr"] = (float)-99.0;
    nt.branches["probe_dgl_trkChi2"] = (float)-99;
    nt.branches["probe_dgl_nsegments"] = (int)-99;
  }
}
  

template <typename TRK>
inline void FillProbeBranchesCosmic(const TRK &trk, NtupleContent &nt, bool passcosmic) {
  nt.branches["probe_isCosmic"] = (bool)passcosmic;
}
template <typename MUOTT, typename TRKTT>
inline void FillPairBranches(const MUOTT &muon, const TRKTT &trk, NtupleContent &nt, PropagateToMuon &prop1_) {
  math::PtEtaPhiMLorentzVector mu1(muon.first.pt(), muon.first.eta(), muon.first.phi(), MU_MASS);
  math::PtEtaPhiMLorentzVector mu2(trk.first.pt(), trk.first.eta(), trk.first.phi(), MU_MASS);
  nt.branches["pair_pt"] = (float)(mu1 + mu2).pt();
  nt.branches["pair_mass"] = (float)(mu1 + mu2).mass();
  nt.branches["pair_eta"] = (float)(mu1 + mu2).eta();
  nt.branches["pair_phi"] = (float)(mu1 + mu2).phi();
  nt.branches["pair_dz"] = (float)(muon.first.vz() - trk.first.vz());
  nt.branches["pair_dR"] = (float)deltaR(muon.first.eta(), muon.first.phi(), trk.first.eta(), trk.first.phi());
  if (muon.second.impactPointStateAvailable() && trk.second.impactPointStateAvailable()) {
    FreeTrajectoryState trajectory_state_muon = muon.second.impactPointTSCP().theState();
    FreeTrajectoryState trajectory_state_trk = trk.second.impactPointTSCP().theState();

    TrajectoryStateOnSurface prop1_M1 = prop1_.extrapolate(trajectory_state_muon);
    TrajectoryStateOnSurface prop2_M1 = prop1_.extrapolate(trajectory_state_trk);

    if (prop1_M1.isValid() && prop2_M1.isValid()) {
      float dphiM1 = deltaPhi<float>(prop1_M1.globalPosition().phi(), prop2_M1.globalPosition().phi());
      nt.branches["pair_drM1"] = (float)hypot(dphiM1, std::abs(prop1_M1.globalPosition().eta() - prop2_M1.globalPosition().eta()));
    } else
      nt.branches["pair_drM1"] = (float)1000;
  } else
    nt.branches["pair_drM1"] = (float)1000;
}

template <typename MUO, typename TRK>
inline void FillTunePPairBranches(const MUO &muon, const TRK &trk, NtupleContent &nt) {
  math::PtEtaPhiMLorentzVector mu1(muon.pt(), muon.eta(), muon.phi(), MU_MASS);
  math::PtEtaPhiMLorentzVector mu2(trk.pt(), trk.eta(), trk.phi(), MU_MASS);
  nt.branches["pair_tuneP_pt"] = (float)(mu1 + mu2).pt();
  nt.branches["pair_tuneP_mass"] = (float)(mu1 + mu2).mass();
  nt.branches["pair_tuneP_eta"] = (float)(mu1 + mu2).eta();
  nt.branches["pair_tuneP_phi"] = (float)(mu1 + mu2).phi();
  nt.branches["pair_tuneP_dz"] = (float)(muon.vz() - trk.vz());
  nt.branches["pair_tuneP_dR"] = (float)deltaR(muon.eta(), muon.phi(), trk.eta(), trk.phi());
}

inline void FillTunePPairBranchesDummy(NtupleContent &nt) {
  nt.branches["pair_tuneP_pt"] = (float)-99;
  nt.branches["pair_tuneP_mass"] = (float)-99;
  nt.branches["pair_tuneP_eta"] = (float)-99;
  nt.branches["pair_tuneP_phi"] = (float)-99;
  nt.branches["pair_tuneP_dz"] = (float)-99;
  nt.branches["pair_tuneP_dR"] = (float)-99;
  nt.branches["pair_tuneP_fit_mass"] = (float)-99;
  nt.branches["pair_tuneP_svprob"] = (float)-99;
  nt.branches["pair_tuneP_normalchi2"] = (float)-99;
}

inline void FillSimMatchingBranches(const pat::Muon &mu, NtupleContent &nt, bool isTag) {
  if (isTag) {
    nt.branches["tag_simType"] = (int)mu.simType();
    nt.branches["tag_simExtType"] = (int)mu.simExtType();
    nt.branches["tag_simFlavour"] = (int)mu.simFlavour();
    nt.branches["tag_simHeaviestMotherFlavour"] = (int)mu.simHeaviestMotherFlavour();
    nt.branches["tag_simPdgId"] = (int)mu.simPdgId();
    nt.branches["tag_simMotherPdgId"] = (int)mu.simMotherPdgId();
    nt.branches["tag_simBX"] = (int)mu.simBX();
    nt.branches["tag_simProdRho"] = (float)mu.simProdRho();
    nt.branches["tag_simProdZ"] = (float)mu.simProdZ();
    nt.branches["tag_simPt"] = (float)mu.simPt();
    nt.branches["tag_simEta"] = (float)mu.simEta();
    nt.branches["tag_simPhi"] = (float)mu.simPhi();
  } else {
    nt.branches["probe_simType"] = (int)mu.simType();
    nt.branches["probe_simExtType"] = (int)mu.simExtType();
    nt.branches["probe_simFlavour"] = (int)mu.simFlavour();
    nt.branches["probe_simHeaviestMotherFlavour"] = (int)mu.simHeaviestMotherFlavour();
    nt.branches["probe_simPdgId"] = (int)mu.simPdgId();
    nt.branches["probe_simMotherPdgId"] = (int)mu.simMotherPdgId();
    nt.branches["probe_simBX"] = (int)mu.simBX();
    nt.branches["probe_simProdRho"] = (float)mu.simProdRho();
    nt.branches["probe_simProdZ"] = (float)mu.simProdZ();
    nt.branches["probe_simPt"] = (float)mu.simPt();
    nt.branches["probe_simEta"] = (float)mu.simEta();
    nt.branches["probe_simPhi"] = (float)mu.simPhi();
  }
}

inline void FillSimMatchingBranchesAOD(const reco::MuonSimInfo &msi, NtupleContent &nt, bool isTag) {
  if (isTag) {
    nt.branches["tag_simType"] = (int)msi.primaryClass;
    nt.branches["tag_simExtType"] = (int)msi.extendedClass;
    nt.branches["tag_simFlavour"] = (int)msi.flavour;
    nt.branches["tag_simHeaviestMotherFlavour"] = (int)msi.heaviestMotherFlavour;
    nt.branches["tag_simPdgId"] = (int)msi.pdgId;
    nt.branches["tag_simMotherPdgId"] = (int)msi.motherPdgId;
    nt.branches["tag_simBX"] = (int)msi.tpBX;
    nt.branches["tag_simProdRho"] = (float)msi.vertex.Rho();
    nt.branches["tag_simProdZ"] = (float)msi.vertex.Z();
    nt.branches["tag_simPt"] = (float)msi.p4.pt();
    nt.branches["tag_simEta"] = (float)msi.p4.eta();
    nt.branches["tag_simPhi"] = (float)msi.p4.phi();
  } else {
    nt.branches["probe_simType"] = (int)msi.primaryClass;
    nt.branches["probe_simExtType"] = (int)msi.extendedClass;
    nt.branches["probe_simFlavour"] = (int)msi.flavour;
    nt.branches["probe_simHeaviestMotherFlavour"] = (int)msi.heaviestMotherFlavour;
    nt.branches["probe_simPdgId"] = (int)msi.pdgId;
    nt.branches["probe_simMotherPdgId"] = (int)msi.motherPdgId;
    nt.branches["probe_simBX"] = (int)msi.tpBX;
    nt.branches["probe_simProdRho"] = (float)msi.vertex.Rho();
    nt.branches["probe_simProdZ"] = (float)msi.vertex.Z();
    nt.branches["probe_simPt"] = (float)msi.p4.pt();
    nt.branches["probe_simEta"] = (float)msi.p4.eta();
    nt.branches["probe_simPhi"] = (float)msi.p4.phi();
  }
}

inline void FillSimMatchingBranchesDummy(NtupleContent &nt, bool isTag) {
  if (isTag) {
    nt.branches["tag_simType"] = (int)-99;
    nt.branches["tag_simExtType"] = (int)-99;
    nt.branches["tag_simFlavour"] = (int)-99;
    nt.branches["tag_simHeaviestMotherFlavour"] = (int)-99;
    nt.branches["tag_simPdgId"] = (int)-99;
    nt.branches["tag_simMotherPdgId"] = (int)-99;
    nt.branches["tag_simBX"] = (int)-99;
    nt.branches["tag_simProdRho"] = (float)-99;
    nt.branches["tag_simProdZ"] = (float)-99;
    nt.branches["tag_simPt"] = (float)-99;
    nt.branches["tag_simEta"] = (float)-99;
    nt.branches["tag_simPhi"] = (float)-99;
  } else {
    nt.branches["probe_simType"] = (int)-99;
    nt.branches["probe_simExtType"] = (int)-99;
    nt.branches["probe_simFlavour"] = (int)-99;
    nt.branches["probe_simHeaviestMotherFlavour"] = (int)-99;
    nt.branches["probe_simPdgId"] = (int)-99;
    nt.branches["probe_simMotherPdgId"] = (int)-99;
    nt.branches["probe_simBX"] = (int)-99;
    nt.branches["probe_simProdRho"] = (float)-99;
    nt.branches["probe_simProdZ"] = (float)-99;
    nt.branches["probe_simPt"] = (float)-99;
    nt.branches["probe_simEta"] = (float)-99;
    nt.branches["probe_simPhi"] = (float)-99;
  }
}

#endif
