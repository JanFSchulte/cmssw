#include "MuonGenAnalyzer.h"

MuonGenAnalyzer::MuonGenAnalyzer(){};

MuonGenAnalyzer::~MuonGenAnalyzer(){};

void MuonGenAnalyzer::SetInputs(const edm::Event& iEvent,
                                const edm::EDGetTokenT<edm::View<reco::GenParticle>>& gens_,
                                const int& momPdg_) {
  iEvent.getByToken(gens_, gens);
  const reco::Candidate* mother1st = nullptr;
  bool found1 = false;
  bool found2 = false;
  gmuon1.SetPtEtaPhiM(0, 0, 0, 0);
  gmuon2.SetPtEtaPhiM(0, 0, 0, 0);
  bool foundFSfromHP1 = false;
  bool foundFSfromHP2 = false;
  gmuonFSfromHP1.SetPtEtaPhiM(0, 0, 0, 0);
  gmuonFSfromHP2.SetPtEtaPhiM(0, 0, 0, 0);
  for (const auto& gen : *gens) {
    if (fabs(gen.pdgId()) != 13)
      continue;

    // look for final state muons from hard process
    if (gen.fromHardProcessFinalState()) {
      if (gen.charge() < 0) {
        if (foundFSfromHP1) {
          std::cout
              << "Warning|MuonGenAnalyzer::SetInputs| there are more than 2 final state mu- from hard process -> skip"
              << std::endl;
        } else {
          gmuonFSfromHP1.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
          foundFSfromHP1 = true;
        }
      } else {
        if (foundFSfromHP2) {
          std::cout
              << "Warning|MuonGenAnalyzer::SetInputs| there are more than 2 final state mu+ from hard process -> skip"
              << std::endl;
        } else {
          gmuonFSfromHP2.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
          foundFSfromHP2 = true;
        }
      }
    }

    // look for final state muons
    if (gen.status() != 1)
      continue;

    bool fromSameMother = false;
    const reco::Candidate* mother_tmp = gen.mother();
    while (mother_tmp) {
      if (mother_tmp->pdgId() == momPdg_) {
        // if there are more than 2 resonances, take the first one for now...
        if (mother1st == nullptr) {
          mother1st = mother_tmp;
        }
        fromSameMother = (mother1st == mother_tmp);
        break;
      }
      mother_tmp = mother_tmp->mother();
    }
    if (!fromSameMother)
      continue;

    if (gen.charge() < 0.) {
      if (found1) {
        std::cout << "Warning|MuonGenAnalyzer::SetInputs| there are more than 2 final state mu- from the same "
                     "resonance -> skip"
                  << std::endl;
      } else {
        gmuon1.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
        found1 = true;
      }
    } else {
      if (found2) {
        std::cout << "Warning|MuonGenAnalyzer::SetInputs| there are more than 2 final state mu+ from the same "
                     "resonance -> skip"
                  << std::endl;
      } else {
        gmuon2.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
        found2 = true;
      }
    }
  }
}

void MuonGenAnalyzer::FillNtuple(NtupleContent& nt) {
  if (gmuon1.Pt() > 0. && gmuon2.Pt() > 0.) {
    nt.branches["genmu1_pt"] = (float)gmuon1.Pt();
    nt.branches["genmu1_eta"] = (float)gmuon1.Eta();
    nt.branches["genmu1_phi"] = (float)gmuon1.Phi();
    nt.branches["genmu1_charge"] = (int)-1;
    nt.branches["genmu2_pt"] = (float)gmuon2.Pt();
    nt.branches["genmu2_eta"] = (float)gmuon2.Eta();
    nt.branches["genmu2_phi"] = (float)gmuon2.Phi();
    nt.branches["genmu2_charge"] = (int)1;
    nt.branches["genMass"] = (float)(gmuon1 + gmuon2).M();
  } else {
    nt.branches["genmu1_pt"] = (float)-99;
    nt.branches["genmu1_eta"] = (float)-99;
    nt.branches["genmu1_phi"] = (float)-99;
    nt.branches["genmu1_charge"] = (int)-99;
    nt.branches["genmu2_pt"] = (float)-99;
    nt.branches["genmu2_eta"] = (float)-99;
    nt.branches["genmu2_phi"] = (float)-99;
    nt.branches["genmu2_charge"] = (int)-99;
    nt.branches["genMass"] = (float)-99;
  }

  // Fill final state muons from hard process
  if (gmuonFSfromHP1.Pt() > 0. && gmuonFSfromHP2.Pt() > 0.) {
    nt.branches["genmuFSfromHP1_pt"] = (float)gmuonFSfromHP1.Pt();
    nt.branches["genmuFSfromHP1_eta"] = (float)gmuonFSfromHP1.Eta();
    nt.branches["genmuFSfromHP1_phi"] = (float)gmuonFSfromHP1.Phi();
    nt.branches["genmuFSfromHP1_charge"] = (int)-1;
    nt.branches["genmuFSfromHP2_pt"] = (float)gmuonFSfromHP2.Pt();
    nt.branches["genmuFSfromHP2_eta"] = (float)gmuonFSfromHP2.Eta();
    nt.branches["genmuFSfromHP2_phi"] = (float)gmuonFSfromHP2.Phi();
    nt.branches["genmuFSfromHP2_charge"] = (int)1;
    nt.branches["genMassFSfromHP"] = (float)(gmuonFSfromHP1 + gmuonFSfromHP2).M();
  } else {
    nt.branches["genmuFSfromHP1_pt"] = (float)-99;
    nt.branches["genmuFSfromHP1_eta"] = (float)-99;
    nt.branches["genmuFSfromHP1_phi"] = (float)-99;
    nt.branches["genmuFSfromHP1_charge"] = (int)-99;
    nt.branches["genmuFSfromHP2_pt"] = (float)-99;
    nt.branches["genmuFSfromHP2_eta"] = (float)-99;
    nt.branches["genmuFSfromHP2_phi"] = (float)-99;
    nt.branches["genmuFSfromHP2_charge"] = (int)-99;
    nt.branches["genMassFSfromHP"] = (float)-99;
  }
}
