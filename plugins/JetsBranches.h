//
// Original Author:
//                Jess Wing Yan Wong
//         Created:  Nov 2020
//
// filling functions for aod and miniaod jets

#ifndef MuonAnalysis_MuonAnalyzer_JetBranches
#define MuonAnalysis_MuonAnalyzer_JetBranches

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <type_traits>
#include "TLorentzVector.h"
#include "NtupleContent.h"
#include "helper.h"

template <typename JET, typename MU>
inline bool CrossClean(const JET &jet, std::vector<MU> &muForJetCleaning) {
  bool muFakeJet = false;
  for (const auto &mu : muForJetCleaning) {
    double dR = deltaR(jet, mu);
    if (dR < 0.4)
      muFakeJet = true;
  }
  return muFakeJet;
}

template <typename JET>
inline void FillJetBranches(const JET &jet, const JET &corrJet, NtupleContent &nt, std::string era) {
  // Jet ID
  bool jetTightID = true, jetTightLepVeto = true;
  if (era.find("2016") != std::string::npos) {
    if (abs(jet.eta()) <= 2.7) {
      if (jet.neutralHadronEnergyFraction() >= 0.9 || jet.neutralEmEnergyFraction() >= 0.9 ||
          jet.numberOfDaughters() <= 1) {
        jetTightID = false;
        jetTightLepVeto = false;
      }
      if (jet.muonEnergyFraction() >= 0.8)
        jetTightLepVeto = false;

      if (abs(jet.eta()) < 2.4) {
        if (jet.chargedHadronEnergyFraction() <= 0 || jet.chargedMultiplicity() <= 0) {
          jetTightID = false;
          jetTightLepVeto = false;
        }
        if (jet.chargedEmEnergyFraction() >= 0.9)
          jetTightLepVeto = false;
        if (jet.chargedEmEnergyFraction() >= 0.99)
          jetTightID = false;
      }
    }
  } else if (era.find("2017") != std::string::npos) {
    if (abs(jet.eta()) <= 2.7) {
      if (jet.neutralHadronEnergyFraction() >= 0.9 || jet.neutralEmEnergyFraction() >= 0.9 ||
          jet.numberOfDaughters() <= 1) {
        jetTightID = false;
        jetTightLepVeto = false;
      }
      if (jet.muonEnergyFraction() >= 0.8)
        jetTightLepVeto = false;

      if (abs(jet.eta()) < 2.4) {
        if (jet.chargedHadronEnergyFraction() <= 0 || jet.chargedMultiplicity() <= 0) {
          jetTightID = false;
          jetTightLepVeto = false;
        }
        if (jet.chargedEmEnergyFraction() >= 0.8)
          jetTightLepVeto = false;
      }
    }
  } else if (era.find("2018") != std::string::npos) {
    if (abs(jet.eta()) <= 2.6) {
      if (jet.neutralHadronEnergyFraction() >= 0.9 || jet.neutralEmEnergyFraction() >= 0.9 ||
          jet.numberOfDaughters() <= 1 || jet.chargedHadronEnergyFraction() <= 0 || jet.chargedMultiplicity() <= 0) {
        jetTightID = false;
        jetTightLepVeto = false;
      }
      if (jet.muonEnergyFraction() >= 0.8 || jet.chargedEmEnergyFraction() >= 0.8) {
        jetTightLepVeto = false;
      }
    } else if (abs(jet.eta()) <= 2.7) {
      if (jet.neutralHadronEnergyFraction() >= 0.9 || jet.neutralEmEnergyFraction() >= 0.99 ||
          jet.chargedMultiplicity() <= 0) {
        jetTightID = false;
        jetTightLepVeto = false;
      }
      if (jet.muonEnergyFraction() >= 0.8 || jet.chargedEmEnergyFraction() >= 0.8) {
        jetTightLepVeto = false;
      }
    }
  }
  if (abs(jet.eta()) > 2.7 && abs(jet.eta()) <= 3.0) {
    if (jet.neutralEmEnergyFraction() <= 0.02 || jet.neutralEmEnergyFraction() >= 0.99 ||
        jet.neutralMultiplicity() <= 2) {
      jetTightID = false;
      jetTightLepVeto = false;
    }
  } else if (abs(jet.eta()) > 3.0) {
    if (jet.neutralHadronEnergyFraction() <= 0.2 || jet.neutralEmEnergyFraction() >= 0.9 ||
        jet.neutralMultiplicity() <= 10) {
      jetTightID = false;
      jetTightLepVeto = false;
    }
  }

  if (jetTightID)
    nt.nTightJets++;
  if (jetTightLepVeto)
    nt.nTightLepVetoJets++;

  // Store Jet Information
  nt.branches["jets_isTight"].push_back((int)jetTightID);
  nt.branches["jets_isTightLepVeto"].push_back((bool)jetTightLepVeto);
  nt.branches["jets_pt"].push_back((float)corrJet.pt());
  nt.branches["jets_eta"].push_back((float)corrJet.eta());
  nt.branches["jets_phi"].push_back((float)corrJet.phi());
  nt.branches["jets_mass"].push_back((float)corrJet.mass());
}

template <typename JET, typename MUON>
inline void FindJetProbePair(const std::vector<JET> &jets, const MUON &mu, NtupleContent &nt) {
  TLorentzVector muonP4;
  muonP4.SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu.mass());
  float minDR = 999;
  JET closestJet;
  TLorentzVector closestJetP4;
  for (const auto &jet : jets) {
    TLorentzVector jetP4;
    jetP4.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.mass());
    // find jet with minDR from mu
    if (jetP4.DeltaR(muonP4) < minDR) {
      minDR = jetP4.DeltaR(muonP4);
      closestJet = jet;
      closestJetP4 = jetP4;
    }
  }
  nt.branches["probe_minDR"] = (float)minDR;
  nt.branches["probe_ptRel_minDR"] =
      (float)(muonP4.P() * (closestJetP4.Vect().Cross(muonP4.Vect()).Mag() / closestJetP4.P() / muonP4.P()));
}
#endif
