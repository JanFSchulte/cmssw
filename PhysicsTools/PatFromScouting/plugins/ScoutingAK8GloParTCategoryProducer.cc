// Categorizes AK8 jets by which GlobalParticleTransformer (GloParT) truth
// class they best match. Category values are numbered to match the index of
// the corresponding score in scoutingFatPFJetReclusterGlobalParticleTransformerJetTags's
// flav_names list (run3scouting_cff.py): probQCD=0, probXbb=1, probXcc=2,
// probXss=3, probXqq=4, probXbs=5, probXgg=6, probXtauhtaue=9,
// probXtauhtaum=10, probXtauhtauh=11, probXbc=12, probXcs=13, probXud=14 --
// so a category value can be used directly as an index into the score array.
// Indices 7 (Xee) and 8 (Xmm) are reserved but never produced by this
// producer (no truth support for X->ee/mumu yet; the corresponding score
// branches aren't currently exposed in NanoAOD either). NonMerged=15 is a
// sentinel appended after the last real class.
//
// "X" is any generic (i.e. not restricted to a specific PDG ID, since
// GloParT is a mass-decorrelated, model-agnostic tagger) hard-process,
// last-copy particle with exactly two direct decay daughters -- including
// SM W/Z/top. This is intentionally independent of and NOT coordinated with
// ScoutingAK8TopWCategoryProducer's Top-merged/W-merged/Z-merged
// categories: a jet from a hadronic W or Z decay is categorized here purely
// by its two daughters' flavour (e.g. Xud/Xqq/Xbb/Xtauhtauh/...), same as
// any other resonance would be. (A top's own decay, top -> W b, never
// matches any category here since W isn't a quark/gluon/lepton, so top
// itself has no direct effect on this categorization either way.) The
// hard-process requirement excludes ordinary hadron decays (e.g. J/psi/
// Upsilon -> mumu) from being mistaken for a genuine X -> mumu/tautau/etc.
//
// Category assignment from the flavour of X's two (last-copy) daughters:
//   Xbb: b bbar        Xcc: c cbar       Xss: s sbar
//   Xgg: g g           Xbc: b/c mixed    Xbs: b/s mixed
//   Xcs: c/s mixed     Xud: u/d mixed    Xqq: any other quark pair
//     (e.g. u ubar, d dbar, u+s, u+c, ...)
//   Xtauhtauh/Xtauhtaue/Xtauhtaum: tau tau, split by each tau's decay mode
//     (hadronic if neither direct daughter of the tau's last copy is an
//     electron or muon, else the corresponding lepton flavour). Other
//     tau-tau decay-mode combinations (e.g. both leptonic) aren't covered by
//     a GloParT class and fall through to NonMerged.
//   NonMerged: an X candidate exists somewhere in the event but its decay
//     products aren't both within the jet cone, or its decay doesn't match
//     any class above.
//   QCD: process-level flag (applyGloParTMatching=False) for samples with no
//     relevant resonance truth chain at all (pure multijet QCD), mirroring
//     ScoutingAK8TopWCategoryProducer's Others/applyTopWMerging.

#include <cmath>
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

  enum GloParTCategory {
    QCD = 0,
    Xbb = 1,
    Xcc = 2,
    Xss = 3,
    Xqq = 4,
    Xbs = 5,
    Xgg = 6,
    // 7 = Xee, 8 = Xmm: reserved, not produced
    Xtauhtaue = 9,
    Xtauhtaum = 10,
    Xtauhtauh = 11,
    Xbc = 12,
    Xcs = 13,
    Xud = 14,
    NonMerged = 15,
  };

  enum class TauDecay { Hadronic, Electron, Muon, Other };

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

  TauDecay tauDecayMode(const reco::Candidate* tau) {
    const reco::Candidate* last = lastCopy(tau);
    for (unsigned i = 0; i < last->numberOfDaughters(); ++i) {
      int id = std::abs(last->daughter(i)->pdgId());
      if (id == 11)
        return TauDecay::Electron;
      if (id == 13)
        return TauDecay::Muon;
    }
    return TauDecay::Hadronic;
  }

  struct XCandidate {
    int category;
    const reco::Candidate* d1;
    const reco::Candidate* d2;
  };

  std::vector<XCandidate> findXCandidates(const reco::GenParticleCollection& genParticles) {
    std::vector<XCandidate> result;
    for (const reco::GenParticle& gen : genParticles) {
      if (!gen.statusFlags().fromHardProcess())
        continue;
      if (!isLastCopy(gen) || gen.numberOfDaughters() != 2)
        continue;

      const reco::Candidate* d0 = lastCopy(gen.daughter(0));
      const reco::Candidate* d1 = lastCopy(gen.daughter(1));
      int id0 = std::abs(d0->pdgId());
      int id1 = std::abs(d1->pdgId());
      auto isPair = [&](int a, int b) { return (id0 == a && id1 == b) || (id0 == b && id1 == a); };

      int category = -1;
      if (id0 >= 1 && id0 <= 5 && id1 >= 1 && id1 <= 5) {
        if (id0 == id1) {
          if (id0 == 5)
            category = Xbb;
          else if (id0 == 4)
            category = Xcc;
          else if (id0 == 3)
            category = Xss;
          else
            category = Xqq;  // u ubar or d dbar
        } else if (isPair(5, 4)) {
          category = Xbc;
        } else if (isPair(5, 3)) {
          category = Xbs;
        } else if (isPair(4, 3)) {
          category = Xcs;
        } else if (isPair(1, 2)) {
          category = Xud;
        } else {
          category = Xqq;
        }
      } else if (id0 == 21 && id1 == 21) {
        category = Xgg;
      } else if (id0 == 15 && id1 == 15) {
        TauDecay t0 = tauDecayMode(d0);
        TauDecay t1 = tauDecayMode(d1);
        bool h0 = t0 == TauDecay::Hadronic;
        bool h1 = t1 == TauDecay::Hadronic;
        if (h0 && h1)
          category = Xtauhtauh;
        else if ((h0 && t1 == TauDecay::Electron) || (h1 && t0 == TauDecay::Electron))
          category = Xtauhtaue;
        else if ((h0 && t1 == TauDecay::Muon) || (h1 && t0 == TauDecay::Muon))
          category = Xtauhtaum;
        // other tau-tau decay-mode combinations (both leptonic, etc.) aren't
        // covered by a GloParT class; category stays -1 and this X is skipped
      }
      // Xee/Xmm (id0==id1==11 or 13) intentionally not categorized yet

      if (category >= 0)
        result.push_back({category, d0, d1});
    }
    return result;
  }

  int categorize(const reco::Jet& jet, const std::vector<XCandidate>& candidates) {
    for (const auto& c : candidates) {
      if (reco::deltaR(jet, *c.d1) < 0.8 && reco::deltaR(jet, *c.d2) < 0.8)
        return c.category;
    }
    return NonMerged;
  }

}  // namespace

class ScoutingAK8GloParTCategoryProducer : public edm::stream::EDProducer<> {
public:
  explicit ScoutingAK8GloParTCategoryProducer(const edm::ParameterSet& iConfig)
      : jetsToken_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
        genParticlesToken_(
            consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
        applyGloParTMatching_(iConfig.getParameter<bool>("applyGloParTMatching")) {
    produces<edm::ValueMap<int>>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets")->setComment("input AK8 jet collection (any reco::Jet-derived collection)");
    desc.add<edm::InputTag>("genParticles")->setComment("prunedGenParticles-like collection with the full decay chain");
    desc.add<bool>("applyGloParTMatching", true)
        ->setComment(
            "set to False for processes with no relevant resonance truth chain (pure multijet QCD); every jet "
            "is then categorized as QCD without touching genParticles");
    descriptions.add("scoutingAK8GloParTCategoryProducer", desc);
  }

private:
  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<edm::View<reco::Jet>> jets;
    iEvent.getByToken(jetsToken_, jets);

    std::vector<int> categories(jets->size(), QCD);

    if (applyGloParTMatching_) {
      edm::Handle<reco::GenParticleCollection> genParticles;
      iEvent.getByToken(genParticlesToken_, genParticles);

      std::vector<XCandidate> candidates = findXCandidates(*genParticles);

      for (std::size_t i = 0; i < jets->size(); ++i)
        categories[i] = categorize(jets->at(i), candidates);
    }

    auto out = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler filler(*out);
    filler.insert(jets, categories.begin(), categories.end());
    filler.fill();
    iEvent.put(std::move(out));
  }

  const edm::EDGetTokenT<edm::View<reco::Jet>> jetsToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  const bool applyGloParTMatching_;
};

DEFINE_FWK_MODULE(ScoutingAK8GloParTCategoryProducer);
