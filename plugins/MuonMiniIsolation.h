//
// Original Author:
//                wing yan wong (jess)
//         Created:  Wed, 14 Oct 2020 19:40:23 GMT
//
// functions for getting miniIsolation of muons

#ifndef MuonAnalysis_MuonAnalyzer_MuonMiniIsolation
#define MuonAnalysis_MuonAnalyzer_MuonMiniIsolation

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <type_traits>
#include "NtupleContent.h"
#include "helper.h"

template <typename PFCANDS>
inline bool is_fromPV(PFCANDS &pfc) {
  return false;
}

template <>
inline bool is_fromPV<reco::PFCandidate>(reco::PFCandidate &pfc) {
  if (pfc.bestTrack())
    return (pfc.bestTrack()->dz() < 0.0 ? true : false);
  return false;
}

template <>
inline bool is_fromPV<pat::PackedCandidate>(pat::PackedCandidate &pfc) {
  return (pfc.fromPV() > 1 ? true : false);
}

inline float miniIsoDr(const reco::Candidate::PolarLorentzVector &p4, float mindr, float maxdr, float kt_scale) {
  return std::max(mindr, std::min(maxdr, float(kt_scale / p4.pt())));
}

// Adapted from PhysicsTools/PatUtils/src/MiniIsolation.cc
template <typename PFCANDS>
pat::PFIsolation getMiniPFIsolation(const std::vector<PFCANDS> & pfcands,
                                const reco::Candidate::PolarLorentzVector &p4,
                                float mindr = 0.05,
                                float maxdr = 0.2,
                                float kt_scale = 10.0,
                                float ptthresh = 0.5,
                                float deadcone_ch = 0.0001,
                                float deadcone_pu = 0.01,
                                float deadcone_ph = 0.01,
                                float deadcone_nh = 0.01,
                                float dZ_cut = 0.0) {
  float chiso = 0, nhiso = 0, phiso = 0, puiso = 0;
  float drcut = miniIsoDr(p4, mindr, maxdr, kt_scale);
  for (auto const pc : pfcands) {
    float dr2 = deltaR2(p4, pc);
    if (dr2 > drcut * drcut)
      continue;
    float pt = pc.p4().pt();
    int id = pc.pdgId();
    if (std::abs(id) == 211) {
      // bool fromPV = (pc.fromPV() > 1 || fabs(pc.dz()) < dZ_cut);
      bool fromPV = is_fromPV(pc);
      if (fromPV && dr2 > deadcone_ch * deadcone_ch) {
        // if charged hadron and from primary vertex, add to charged hadron isolation
        chiso += pt;
      } else if (!fromPV && pt > ptthresh && dr2 > deadcone_pu * deadcone_pu) {
        // if charged hadron and NOT from primary vertex, add to pileup isolation
        puiso += pt;
      }
    }
    // if neutral hadron, add to neutral hadron isolation
    if (std::abs(id) == 130 && pt > ptthresh && dr2 > deadcone_nh * deadcone_nh)
      nhiso += pt;
    // if photon, add to photon isolation
    if (std::abs(id) == 22 && pt > ptthresh && dr2 > deadcone_ph * deadcone_ph)
      phiso += pt;
  }
  return pat::PFIsolation(chiso, nhiso, phiso, puiso);
}

// Adapted from implementation from Daniel Li (Brown)
template <typename MUON, typename PFCANDS>
inline void FillMiniIso(
  const std::vector<PFCANDS> & pfcands, 
  const MUON &mu,
  const double rho,
  NtupleContent &nt,
  bool isTag  
) {
  
  double Aeff_Fall17[5] = { 0.0566, 0.0562, 0.0363, 0.0119, 0.0064 };
  double EA;

  auto iso = getMiniPFIsolation<PFCANDS>(pfcands, mu.polarP4());
  // auto iso = mu.miniPFIsolation();
  auto chg = iso.chargedHadronIso();
  auto neu = iso.neutralHadronIso();
  auto pho = iso.photonIso();
    
  if( TMath::Abs(mu.eta()) < 0.8 ) EA = Aeff_Fall17[0];
  else if( TMath::Abs(mu.eta()) < 1.3 ) EA = Aeff_Fall17[1];
  else if( TMath::Abs(mu.eta()) < 2.0 ) EA = Aeff_Fall17[2];
  else if( TMath::Abs(mu.eta()) < 2.2 ) EA = Aeff_Fall17[3];
  else EA = Aeff_Fall17[4];

  float R = 10.0 / std::min( std::max( mu.pt(), 50.0 ), 200.0 );
  EA *= std::pow( R / 0.3, 2 );
    
  float miniIso = ( chg + TMath::Max( 0.0, neu + pho - (rho) * EA ) ) / mu.pt();

  if( isTag ){
    nt.tag_miniIso = miniIso;
    nt.tag_miniIsoCharged = chg / mu.pt();
    nt.tag_miniIsoPhotons = pho / mu.pt();
    nt.tag_miniIsoNeutrals = neu / mu.pt();
  }
  else {
    nt.probe_miniIso = miniIso;
    nt.probe_miniIsoCharged = chg / mu.pt();
    nt.probe_miniIsoPhotons = pho / mu.pt();
    nt.probe_miniIsoNeutrals = neu / mu.pt();
  }
}

#endif
