#ifndef PhysicsTools_PatAlgos_SoftMuonMvaEstimatorRun3_h
#define PhysicsTools_PatAlgos_SoftMuonMvaEstimatorRun3_h

#include <memory>
#include <string>
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

//

namespace pat {
  class Muon;
}

namespace edm {
  class FileInPath;
}

namespace pat {

  typedef std::pair<const reco::MuonChamberMatch*, const reco::MuonSegmentMatch*> MatchPair;

  class SoftMuonMvaEstimatorRun3 {
  public:
    SoftMuonMvaEstimatorRun3(const edm::FileInPath &weightsfile);
    ~SoftMuonMvaEstimatorRun3() = default;

    static void fillDescriptions(edm::ConfigurationDescriptions &);
    static void globalEndJob(const cms::Ort::ONNXRuntime *);
    const MatchPair& getBetterMatch(const MatchPair&, const MatchPair&) const;
    float dX(const MatchPair&) const;
    float pullX(const MatchPair&) const;
    float pullDxDz(const MatchPair&) const;
    float dY(const MatchPair&) const;
    float pullY(const MatchPair&) const;
    float pullDyDz(const MatchPair&) const;
    std::vector<float> computeMVAID(const pat::Muon &imuon) const;


  private:
    std::vector<std::string> flav_names_;   // names of the output scores
    std::vector<std::string> input_names_;  // names of each input group - the ordering is important!
    std::unique_ptr<const cms::Ort::ONNXRuntime> hgb_;
  };
};  // namespace pat
#endif
