// Nearest-deltaR match of a bare (non-reco::Candidate) scouting lepton
// collection (Run3ScoutingMuon, Run3ScoutingElectron) to reco::GenParticle,
// producing an edm::ValueMap<int> of "index into the matched collection"
// (intended to be finalGenParticles, i.e. the GenPart NanoAOD table's
// source), or -1 if unmatched.
//
// Modeled on CommonTools/RecoAlgos/plugins/JetDeltaRValueMapProducer.cc's
// greedy one-to-one (per-reco-object lock) matching, and on
// PhysicsTools/HepMCCandAlgos/plugins/MCTruthMatchers.cc's MCMatcher
// pdgId/status filtering -- but templated so it also works on the
// non-reco::Candidate scouting lepton formats (only .pt()/.eta()/.phi() are
// required of T).
//
// genPartFlav (the flavour byte CandMCMatchTableProducer also emits for
// standard NanoAOD leptons) is deliberately not produced here.

#include <vector>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Math/interface/deltaR.h"

template <typename T>
class ScoutingLeptonGenParticleMatcher : public edm::stream::EDProducer<> {
public:
  explicit ScoutingLeptonGenParticleMatcher(edm::ParameterSet const& params)
      : srcToken_(consumes<edm::View<T>>(params.getParameter<edm::InputTag>("src"))),
        matchedToken_(consumes<edm::View<reco::GenParticle>>(params.getParameter<edm::InputTag>("matched"))),
        mcPdgId_(params.getParameter<std::vector<int>>("mcPdgId")),
        mcStatus_(params.getParameter<std::vector<int>>("mcStatus")),
        maxDeltaR_(params.getParameter<double>("maxDeltaR")) {
    produces<edm::ValueMap<int>>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("input scouting lepton collection");
    desc.add<edm::InputTag>("matched")->setComment("gen particle collection to match against, e.g. finalGenParticles");
    desc.add<std::vector<int>>("mcPdgId")->setComment("absolute-value pdgId accept-list");
    desc.add<std::vector<int>>("mcStatus")->setComment("status accept-list");
    desc.add<double>("maxDeltaR", 0.3);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  bool passesSelection(const reco::GenParticle& gp) const {
    if (!mcPdgId_.empty()) {
      bool ok = false;
      for (int id : mcPdgId_) {
        if (std::abs(gp.pdgId()) == std::abs(id)) {
          ok = true;
          break;
        }
      }
      if (!ok)
        return false;
    }
    if (!mcStatus_.empty()) {
      bool ok = false;
      for (int st : mcStatus_) {
        if (gp.status() == st) {
          ok = true;
          break;
        }
      }
      if (!ok)
        return false;
    }
    return true;
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<edm::View<T>> src;
    iEvent.getByToken(srcToken_, src);
    edm::Handle<edm::View<reco::GenParticle>> matched;
    iEvent.getByToken(matchedToken_, matched);

    std::vector<int> result(src->size(), -1);
    std::vector<bool> srcLocked(src->size(), false);

    // Outer loop over the FULL matched (finalGenParticles) collection, so
    // that the recorded index always corresponds to a row in that
    // collection, never to a filtered subset.
    for (unsigned int gi = 0; gi < matched->size(); ++gi) {
      const reco::GenParticle& gp = matched->at(gi);
      if (!passesSelection(gp))
        continue;

      float bestDR2 = maxDeltaR_ * maxDeltaR_;
      int bestSrc = -1;
      for (unsigned int si = 0; si < src->size(); ++si) {
        if (srcLocked[si])
          continue;
        const T& obj = src->at(si);
        float dr2 = reco::deltaR2(obj.eta(), obj.phi(), gp.eta(), gp.phi());
        if (dr2 < bestDR2) {
          bestDR2 = dr2;
          bestSrc = static_cast<int>(si);
        }
      }
      if (bestSrc >= 0) {
        srcLocked[bestSrc] = true;
        result[bestSrc] = static_cast<int>(gi);
      }
    }

    auto vm = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler filler(*vm);
    filler.insert(src, result.begin(), result.end());
    filler.fill();
    iEvent.put(std::move(vm));
  }

  const edm::EDGetTokenT<edm::View<T>> srcToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> matchedToken_;
  const std::vector<int> mcPdgId_;
  const std::vector<int> mcStatus_;
  const double maxDeltaR_;
};

typedef ScoutingLeptonGenParticleMatcher<Run3ScoutingMuon> ScoutingMuonGenParticleMatcher;
typedef ScoutingLeptonGenParticleMatcher<Run3ScoutingElectron> ScoutingElectronGenParticleMatcher;

DEFINE_FWK_MODULE(ScoutingMuonGenParticleMatcher);
DEFINE_FWK_MODULE(ScoutingElectronGenParticleMatcher);
