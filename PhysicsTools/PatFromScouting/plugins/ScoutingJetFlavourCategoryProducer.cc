// Categorizes jets into the heavy/light-flavour training categories used for
// b-tagging, following the scheme implemented in
// DeepNTuples/DeepNtuplizer/src/helpers.cc (deep_ntuples::jet_flavour).
//
// The muon/electron/tau sub-categories (MU, ELE, TAU, TAUP1H0P, ...) of that
// scheme require gen lepton/tau matching that DeepNtuplizer itself does not
// exercise (ntuple_JetInfo.cc hard-codes the tau collection to stay empty and
// never runs the muon/electron matching for this call), so they are omitted
// here; the enum values below are otherwise numbered identically to
// DeepNTuples/DeepNtuplizer/interface/helpers.h::JetFlavor so that PU keeps
// its original value of 26.

#include <cstdlib>
#include <vector>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace {

  enum JetFlavorCategory {
    UNDEFINED = 0,
    G = 1,
    U = 2,
    D = 3,
    S = 4,
    C = 5,
    GCC = 6,
    CC = 7,
    B = 8,
    GBB = 9,
    BB = 10,
    LeptonicB = 11,
    LeptonicB_C = 12,
    PU = 26,
  };

  int jetFlavourCategory(const pat::Jet& jet,
                          const std::vector<reco::GenParticle>& gToBB,
                          const std::vector<reco::GenParticle>& gToCC,
                          const std::vector<reco::GenParticle>& neutrinosLepB,
                          const std::vector<reco::GenParticle>& neutrinosLepB_C,
                          bool usePhysForLightAndUndefined) {
    int hflav = std::abs(jet.hadronFlavour());
    int pflav = std::abs(jet.partonFlavour());
    int physflav = 0;

    if (!jet.genJet())
      return pflav == 0 ? static_cast<int>(PU) : static_cast<int>(UNDEFINED);

    if (jet.genParton())
      physflav = std::abs(jet.genParton()->pdgId());

    std::size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
    std::size_t ncs = jet.jetFlavourInfo().getcHadrons().size();

    unsigned int nbFromGSP = 0;
    for (const auto& p : gToBB) {
      if (reco::deltaR(jet, p) < 0.4)
        ++nbFromGSP;
    }
    unsigned int ncFromGSP = 0;
    for (const auto& p : gToCC) {
      if (reco::deltaR(jet, p) < 0.4)
        ++ncFromGSP;
    }

    auto lightOrUndefined = [&]() -> int {
      if (usePhysForLightAndUndefined) {
        if (physflav == 21)
          return G;
        if (physflav == 3)
          return S;
        if (physflav == 2)
          return U;
        if (physflav == 1)
          return D;
      }
      return UNDEFINED;
    };

    if (hflav == 5) {  // B jet
      if (nbs > 1)
        return nbFromGSP > 0 ? GBB : BB;
      if (nbs == 1) {
        for (const auto& n : neutrinosLepB)
          if (reco::deltaR(n.eta(), n.phi(), jet.eta(), jet.phi()) < 0.4)
            return LeptonicB;
        for (const auto& n : neutrinosLepB_C)
          if (reco::deltaR(n.eta(), n.phi(), jet.eta(), jet.phi()) < 0.4)
            return LeptonicB_C;
        return B;
      }
      return lightOrUndefined();
    } else if (hflav == 4) {  // C jet
      if (ncs > 1)
        return ncFromGSP > 0 ? GCC : CC;
      return C;
    } else {  // not a heavy jet
      if (std::abs(pflav) == 4 || std::abs(pflav) == 5 || nbs || ncs)
        return lightOrUndefined();
      if (usePhysForLightAndUndefined)
        return lightOrUndefined();
      if (pflav == 21)
        return G;
      if (pflav == 3)
        return S;
      if (pflav == 2)
        return U;
      if (pflav == 1)
        return D;
      return UNDEFINED;
    }
  }

}  // namespace

class ScoutingJetFlavourCategoryProducer : public edm::stream::EDProducer<> {
public:
  explicit ScoutingJetFlavourCategoryProducer(const edm::ParameterSet& iConfig)
      : jetsToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
        genParticlesToken_(
            consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
        usePhysForLightAndUndefined_(iConfig.getParameter<bool>("usePhysForLightAndUndefined")) {
    produces<edm::ValueMap<int>>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets")->setComment(
        "input pat::Jet collection; must have hadronFlavour/partonFlavour/jetFlavourInfo/genJet/genParton "
        "available");
    desc.add<edm::InputTag>("genParticles")->setComment(
        "prunedGenParticles-like collection used to build the gluon-splitting and leptonic-B categories");
    desc.add<bool>("usePhysForLightAndUndefined", false)
        ->setComment("fall back to the matched parton pdgId for light/gluon/undefined jets");
    descriptions.add("scoutingJetFlavourCategoryProducer", desc);
  }

private:
  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(jetsToken_, jets);

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);

    std::vector<reco::GenParticle> gToBB, gToCC, neutrinosLepB, neutrinosLepB_C;
    for (const reco::GenParticle& gen : *genParticles) {
      int id = std::abs(gen.pdgId());
      if ((id == 12 || id == 14 || id == 16) && gen.mother() != nullptr) {
        int momId = std::abs(gen.mother()->pdgId());
        if ((momId > 500 && momId < 600) || (momId > 5000 && momId < 6000))
          neutrinosLepB.push_back(gen);
        if ((momId > 400 && momId < 500) || (momId > 4000 && momId < 5000))
          neutrinosLepB_C.push_back(gen);
      }
      if (id == 21 && gen.status() >= 21 && gen.status() <= 59 && gen.numberOfDaughters() == 2) {
        const reco::Candidate* d0 = gen.daughter(0);
        const reco::Candidate* d1 = gen.daughter(1);
        if (std::abs(d0->pdgId()) == 5 && std::abs(d1->pdgId()) == 5 && d0->pdgId() * d1->pdgId() < 0 &&
            reco::deltaR(*d0, *d1) < 0.4)
          gToBB.push_back(gen);
        if (std::abs(d0->pdgId()) == 4 && std::abs(d1->pdgId()) == 4 && d0->pdgId() * d1->pdgId() < 0 &&
            reco::deltaR(*d0, *d1) < 0.4)
          gToCC.push_back(gen);
      }
    }

    std::vector<int> categories;
    categories.reserve(jets->size());
    for (const pat::Jet& jet : *jets)
      categories.push_back(
          jetFlavourCategory(jet, gToBB, gToCC, neutrinosLepB, neutrinosLepB_C, usePhysForLightAndUndefined_));

    auto out = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler filler(*out);
    filler.insert(jets, categories.begin(), categories.end());
    filler.fill();
    iEvent.put(std::move(out));
  }

  const edm::EDGetTokenT<edm::View<pat::Jet>> jetsToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  const bool usePhysForLightAndUndefined_;
};

DEFINE_FWK_MODULE(ScoutingJetFlavourCategoryProducer);
