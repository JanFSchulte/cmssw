/**
 * PuppiRawAlphaToValueMapProducer
 *
 * PuppiProducer's puppiDiagnostics=True output (PuppiRawAlphas/PuppiAlphasMed/
 * PuppiAlphasRms) is a set of bare event-level std::vector<double> products,
 * not ValueMaps keyed to the candidate collection -- so they can't be
 * attached directly as SimplePFCandidateFlatTableProducer externalVariables.
 * This producer republishes them as three parallel ValueMap<float> keyed to
 * the same candidate collection PuppiProducer consumed, so they line up as
 * sibling columns next to a variant's puppiWeight_<variant> column.
 *
 * IMPORTANT, and easy to get wrong (caught by cross-checking a reproduced
 * weight against the stored puppiWeight_nominal -- see analysis/
 * validate_recompute.py): PuppiRawAlphas and PuppiAlphasMed/PuppiAlphasRms
 * do NOT share a layout.
 *  - PuppiRawAlphas (PuppiContainer::getRawAlphas) is laid out as
 *    [algoBlock][candidate], one block per top-level PuppiProducer "algos"
 *    VPSet entry (PuppiContainer::PuppiContainer does one
 *    fPuppiAlgo.emplace_back(algos) per entry -- confirmed via direct
 *    source read), each block the same length as the candidate collection.
 *    Puppi_cff.py's production config has TWO such entries -- algos[0]
 *    ("central", |eta|<2.5) and algos[1] ("forward", a SINGLE PuppiAlgo
 *    object internally covering BOTH 2.5<=|eta|<3.0 and |eta|>=3.0 via its
 *    own etaMin/etaMax sub-bins, since cone/rmsPtMin -- the only things that
 *    affect the raw alpha value itself, see goodVar()/coneSize() -- are
 *    shared by both forward sub-bins in Puppi_cff.py's puppiForward PSet --
 *    so raw alpha only ever needs a 2-way (not 3-way) block choice for this
 *    config; MedEtaSF/RMSEtaSF do differ per forward sub-bin, but those are
 *    already correctly captured per-candidate in alphaMed/alphaRms below,
 *    not in rawAlpha. A single fixed algoBlock_ (as this producer used to
 *    have, hardcoded to 0/central) is therefore only correct for candidates
 *    actually in the central region -- every forward candidate got the
 *    central block's cone=0.4/rmsPtMin=0.1 raw alpha instead of the forward
 *    block's cone=0.4/rmsPtMin=0.5 value, silently wrong for ~14% of
 *    candidates (the |eta|>=2.5 population, measured on real v5 data).
 *    Fixed below: each candidate's own eta picks its block via
 *    etaBoundaries (ascending upper edges of every block but the last),
 *    mirroring PuppiContainer::getPuppiId's per-candidate eta lookup
 *    (simplified vs. the real getPuppiId: this assumes every "algos" ptMin
 *    is 0, true for every block in Puppi_cff.py's actual central/forward
 *    definitions, so the ptMin-based fallback getPuppiId also implements
 *    never triggers and doesn't need reproducing here).
 *  - PuppiAlphasMed/PuppiAlphasRms (PuppiContainer::calculatePuppiWeights,
 *    filled via fPuppiAlgo[pPupId].median()/.rms() inside the per-candidate
 *    loop) are ALREADY per-candidate: entry i is candidate i's own eta-region
 *    algo's median/RMS, one value per candidate, length == nCandidates. They
 *    must be passed straight through, NOT block-indexed -- indexing into
 *    them by a block offset (as an earlier version of this file did) reads
 *    one arbitrary candidate's value and broadcasts it to the whole event,
 *    which silently corrupts every reproduced weight (and produces
 *    NaN/negative-rms garbage whenever that arbitrary candidate happens to
 *    be one PuppiContainer flagged outside its fiducial region, sentinel
 *    value -10).
 *
 * With these three columns correctly filled, PuppiAlgo::compute()'s
 * chi2 = (alpha-med)*|alpha-med|/rms^2, weight = chisquared_cdf(chi2, ndof)
 * can be recomputed offline for any hypothetical MedEtaSF/RMSEtaSF (which
 * only rescale med/rms) -- see analysis/puppi_recompute.py for the full
 * reproduction including the MinNeutralPt/MinNeutralPtSlope floor and the
 * CHS id==1/2 overrides. Only a cone/rmsPtMin change requires a new
 * production pass, since those affect the alpha computation itself.
 */

#include <cmath>
#include <memory>
#include <limits>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

class PuppiRawAlphaToValueMapProducer : public edm::stream::EDProducer<> {
public:
  explicit PuppiRawAlphaToValueMapProducer(const edm::ParameterSet&);
  ~PuppiRawAlphaToValueMapProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<reco::PFCandidateCollection> candToken_;
  const edm::EDGetTokenT<std::vector<double>> rawAlphasToken_;
  const edm::EDGetTokenT<std::vector<double>> alphasMedToken_;
  const edm::EDGetTokenT<std::vector<double>> alphasRmsToken_;
  const std::vector<double> etaBoundaries_;
};

PuppiRawAlphaToValueMapProducer::PuppiRawAlphaToValueMapProducer(const edm::ParameterSet& iConfig)
    : candToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("candidates"))),
      rawAlphasToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("rawAlphas"))),
      alphasMedToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("alphasMed"))),
      alphasRmsToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("alphasRms"))),
      etaBoundaries_(iConfig.getParameter<std::vector<double>>("etaBoundaries")) {
  produces<edm::ValueMap<float>>("rawAlpha");
  produces<edm::ValueMap<float>>("alphaMed");
  produces<edm::ValueMap<float>>("alphaRms");
}

void PuppiRawAlphaToValueMapProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  auto candHandle = iEvent.getHandle(candToken_);
  const size_t nCand = candHandle->size();

  const auto& rawAlphas = iEvent.get(rawAlphasToken_);
  const auto& alphasMed = iEvent.get(alphasMedToken_);
  const auto& alphasRms = iEvent.get(alphasRmsToken_);

  constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
  const size_t nBlocks = etaBoundaries_.size() + 1;

  // PuppiRawAlphas: [algoBlock][candidate] -- each candidate's own |eta|
  // picks its block (etaBoundaries_ = ascending upper edges of every block
  // but the last), mirroring which top-level "algos" PuppiAlgo object
  // PuppiContainer::getPuppiId would have assigned it -- see class doc
  // comment for why a single fixed block was wrong for forward candidates.
  std::vector<float> rawAlphaOut(nCand, kNaN);
  for (size_t i = 0; i < nCand; ++i) {
    const double absEta = std::abs((*candHandle)[i].eta());
    size_t block = nBlocks - 1;
    for (size_t b = 0; b < etaBoundaries_.size(); ++b) {
      if (absEta < etaBoundaries_[b]) {
        block = b;
        break;
      }
    }
    const size_t offset = block * nCand;
    if (offset + i < rawAlphas.size()) {
      rawAlphaOut[i] = static_cast<float>(rawAlphas[offset + i]);
    }
  }

  // PuppiAlphasMed/PuppiAlphasRms: already one entry per candidate -- direct
  // passthrough, no algoBlock_ indexing.
  std::vector<float> medOut(nCand, kNaN);
  std::vector<float> rmsOut(nCand, kNaN);
  if (alphasMed.size() == nCand && alphasRms.size() == nCand) {
    for (size_t i = 0; i < nCand; ++i) {
      medOut[i] = static_cast<float>(alphasMed[i]);
      rmsOut[i] = static_cast<float>(alphasRms[i]);
    }
  }

  auto fill = [&](const char* label, std::vector<float>& values) {
    auto vm = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler filler(*vm);
    filler.insert(candHandle, values.begin(), values.end());
    filler.fill();
    iEvent.put(std::move(vm), label);
  };
  fill("rawAlpha", rawAlphaOut);
  fill("alphaMed", medOut);
  fill("alphaRms", rmsOut);
}

void PuppiRawAlphaToValueMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("candidates")
      ->setComment("The same candidate collection/instance the source PuppiProducer ran over "
                   "(candName), e.g. packedPFCandidates:recoCands -- ValueMap order must match.");
  desc.add<edm::InputTag>("rawAlphas")->setComment("PuppiProducer's <module>:PuppiRawAlphas output.");
  desc.add<edm::InputTag>("alphasMed")->setComment("PuppiProducer's <module>:PuppiAlphasMed output.");
  desc.add<edm::InputTag>("alphasRms")->setComment("PuppiProducer's <module>:PuppiAlphasRms output.");
  desc.add<std::vector<double>>("etaBoundaries", {2.5})
      ->setComment("Ascending |eta| upper edges of every PuppiRawAlphas block but the last -- "
                   "each candidate's own |eta| picks its block from these (that product only, see "
                   "class doc comment; alphaMed/alphaRms are always per-candidate already). Default "
                   "{2.5} gives 2 blocks (0: |eta|<2.5 'central', 1: |eta|>=2.5 'forward'), matching "
                   "Puppi_cff.py's production 'algos' VPSet (algos[0]=central, algos[1]=forward -- "
                   "the forward entry internally covers both its 2.5-3.0 and 3.0-10.0 sub-bins with "
                   "one PuppiAlgo object/raw-alpha block, since they share cone/rmsPtMin).");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(PuppiRawAlphaToValueMapProducer);
