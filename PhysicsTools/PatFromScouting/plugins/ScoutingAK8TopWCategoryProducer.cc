// Categorizes AK8 jets into gen-level top/W/Z merging categories, following
// slide 5 of https://cds.cern.ch/record/2941747/files/DP2025_052.pdf (categories
// 1-4) plus a Z-merged extension (category 5, not from that note):
//
//   1 Top-merged: the b quark and both W-decay quarks from a hadronically
//     decaying top are all within dR < 0.8 of the fatjet.
//   2 W-merged: both W-decay quarks are within dR < 0.8 of the fatjet, and
//     either (a) the b quark from the same top is outside that cone, or
//     (b) the merged W does not come from a hadronic top decay at all (e.g.
//     the associated on-shell W in single-top tW production).
//   3 Non-merged: none of the above (top/W/Z background present, but not
//     merged into this jet).
//   4 Others: processes with no relevant top/W/Z truth chain at all (pure
//     multijet QCD). This is a process-level flag set by the user per sample
//     (applyTopWMerging=False), not something derived from an absence of
//     gen tops/Zs in a given event.
//   5 Z-merged: both quarks from a hadronically-decaying Z boson are within
//     dR < 0.8 of the fatjet. Checked after Top-merged/W-merged, so a jet
//     that (very unusually) also happens to overlap an unrelated Z is still
//     reported as the higher-priority top/W category.
//
// The gen matching walks prunedGenParticles: for every top quark, its last
// same-flavour copy is found, then its W and b daughters; the W is followed
// to its own last copy and, if it decays hadronically (into two quarks), its
// two daughters and the top's b become one (w1, w2, b) candidate. Separately,
// every hadronically-decaying W boson not descending from a top quark is
// recorded as a (w1, w2) candidate with no associated b (case 2b above), and
// every hadronically-decaying Z boson is recorded as a (q1, q2) candidate.

#include <array>
#include <cmath>
#include <utility>
#include <vector>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace {

  enum TopWCategory {
    TopMerged = 1,
    WMerged = 2,
    NonMerged = 3,
    Others = 4,
    ZMerged = 5,
  };

  // follows same-pdgId daughters down to the final copy before it actually decays
  const reco::Candidate* lastCopy(const reco::Candidate* p) {
    for (unsigned i = 0; i < p->numberOfDaughters(); ++i) {
      if (p->daughter(i)->pdgId() == p->pdgId())
        return lastCopy(p->daughter(i));
    }
    return p;
  }

  bool isLastCopy(const reco::Candidate& p) {
    for (unsigned i = 0; i < p.numberOfDaughters(); ++i)
      if (p.daughter(i)->pdgId() == p.pdgId())
        return false;
    return true;
  }

  bool isDescendantOfTop(const reco::Candidate& p) {
    const reco::Candidate* anc = &p;
    for (int depth = 0; depth < 20 && anc->numberOfMothers() > 0; ++depth) {
      anc = anc->mother(0);
      if (std::abs(anc->pdgId()) == 6)
        return true;
    }
    return false;
  }

  struct TopWB {
    const reco::Candidate* w1;
    const reco::Candidate* w2;
    const reco::Candidate* b;
  };

  int categorize(const reco::Jet& jet, const std::vector<TopWB>& topHadronicWs,
                 const std::vector<std::pair<const reco::Candidate*, const reco::Candidate*>>& standaloneHadronicWs,
                 const std::vector<std::pair<const reco::Candidate*, const reco::Candidate*>>& hadronicZs) {
    for (const auto& t : topHadronicWs) {
      if (reco::deltaR(jet, *t.w1) < 0.8 && reco::deltaR(jet, *t.w2) < 0.8 && reco::deltaR(jet, *t.b) < 0.8)
        return TopMerged;
    }
    for (const auto& t : topHadronicWs) {
      if (reco::deltaR(jet, *t.w1) < 0.8 && reco::deltaR(jet, *t.w2) < 0.8 && reco::deltaR(jet, *t.b) >= 0.8)
        return WMerged;
    }
    for (const auto& w : standaloneHadronicWs) {
      if (reco::deltaR(jet, *w.first) < 0.8 && reco::deltaR(jet, *w.second) < 0.8)
        return WMerged;
    }
    for (const auto& z : hadronicZs) {
      if (reco::deltaR(jet, *z.first) < 0.8 && reco::deltaR(jet, *z.second) < 0.8)
        return ZMerged;
    }
    return NonMerged;
  }

}  // namespace

class ScoutingAK8TopWCategoryProducer : public edm::stream::EDProducer<> {
public:
  explicit ScoutingAK8TopWCategoryProducer(const edm::ParameterSet& iConfig)
      : jetsToken_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
        genParticlesToken_(
            consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
        applyTopWMerging_(iConfig.getParameter<bool>("applyTopWMerging")) {
    produces<edm::ValueMap<int>>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets")->setComment("input AK8 jet collection (any reco::Jet-derived collection)");
    desc.add<edm::InputTag>("genParticles")->setComment("prunedGenParticles-like collection with the full top/W decay chain");
    desc.add<bool>("applyTopWMerging", true)
        ->setComment(
            "set to False for processes with no relevant top/W/Z truth chain (pure multijet QCD); every jet is "
            "then categorized as Others without touching genParticles");
    descriptions.add("scoutingAK8TopWCategoryProducer", desc);
  }

private:
  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<edm::View<reco::Jet>> jets;
    iEvent.getByToken(jetsToken_, jets);

    std::vector<int> categories(jets->size(), Others);

    if (applyTopWMerging_) {
      edm::Handle<reco::GenParticleCollection> genParticles;
      iEvent.getByToken(genParticlesToken_, genParticles);

      std::vector<TopWB> topHadronicWs;
      std::vector<std::pair<const reco::Candidate*, const reco::Candidate*>> standaloneHadronicWs;

      for (const reco::GenParticle& gen : *genParticles) {
        if (std::abs(gen.pdgId()) != 6 || !isLastCopy(gen))
          continue;
        const reco::Candidate* wDau = nullptr;
        const reco::Candidate* bDau = nullptr;
        for (unsigned i = 0; i < gen.numberOfDaughters(); ++i) {
          const reco::Candidate* d = gen.daughter(i);
          if (std::abs(d->pdgId()) == 24)
            wDau = d;
          else if (std::abs(d->pdgId()) == 5)
            bDau = d;
        }
        if (!wDau || !bDau)
          continue;
        wDau = lastCopy(wDau);
        bDau = lastCopy(bDau);
        if (wDau->numberOfDaughters() != 2)
          continue;
        const reco::Candidate* q0 = wDau->daughter(0);
        const reco::Candidate* q1 = wDau->daughter(1);
        if (std::abs(q0->pdgId()) > 5 || std::abs(q1->pdgId()) > 5)
          continue;  // leptonic W decay
        topHadronicWs.push_back({q0, q1, bDau});
      }

      for (const reco::GenParticle& gen : *genParticles) {
        if (std::abs(gen.pdgId()) != 24 || !isLastCopy(gen) || gen.numberOfDaughters() != 2)
          continue;
        const reco::Candidate* q0 = gen.daughter(0);
        const reco::Candidate* q1 = gen.daughter(1);
        if (std::abs(q0->pdgId()) > 5 || std::abs(q1->pdgId()) > 5)
          continue;  // leptonic W decay
        if (!isDescendantOfTop(gen))
          standaloneHadronicWs.emplace_back(q0, q1);
      }

      std::vector<std::pair<const reco::Candidate*, const reco::Candidate*>> hadronicZs;
      for (const reco::GenParticle& gen : *genParticles) {
        if (gen.pdgId() != 23 || !isLastCopy(gen) || gen.numberOfDaughters() != 2)
          continue;
        const reco::Candidate* q0 = gen.daughter(0);
        const reco::Candidate* q1 = gen.daughter(1);
        if (std::abs(q0->pdgId()) > 5 || std::abs(q1->pdgId()) > 5)
          continue;  // leptonic/invisible Z decay
        hadronicZs.emplace_back(q0, q1);
      }

      for (std::size_t i = 0; i < jets->size(); ++i)
        categories[i] = categorize(jets->at(i), topHadronicWs, standaloneHadronicWs, hadronicZs);
    }

    auto out = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler filler(*out);
    filler.insert(jets, categories.begin(), categories.end());
    filler.fill();
    iEvent.put(std::move(out));
  }

  const edm::EDGetTokenT<edm::View<reco::Jet>> jetsToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  const bool applyTopWMerging_;
};

DEFINE_FWK_MODULE(ScoutingAK8TopWCategoryProducer);
