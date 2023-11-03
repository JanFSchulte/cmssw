#include "PhysicsTools/PatAlgos/interface/SoftMuonMvaEstimatorRun3.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

using namespace pat;
using namespace cms::Ort;

SoftMuonMvaEstimatorRun3::SoftMuonMvaEstimatorRun3(const edm::FileInPath &weightsfile) {
  hgb_ = std::make_unique<ONNXRuntime>(weightsfile.fullPath());
  LogDebug("MuonMvaIDEstimator") << hgb_.get();
}

void SoftMuonMvaEstimatorRun3::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::FileInPath>("mvaIDTrainingFile", edm::FileInPath("RecoMuon/MuonIdentification/data/rf_HGB_weighted.onnx"));
  desc.add<std::vector<std::string>>("flav_names",
                                     std::vector<std::string>{
                                         "probGOOD",
                                     });

  descriptions.addWithDefaultLabel(desc);
}

void SoftMuonMvaEstimatorRun3::globalEndJob(const cms::Ort::ONNXRuntime *cache) {}


const MatchPair& SoftMuonMvaEstimatorRun3::getBetterMatch(const MatchPair& match1, const MatchPair& match2) const{

  if (match2.first->detector() == MuonSubdetId::DT and
      match1.first->detector() != MuonSubdetId::DT)
    return match2;
 
  if ( abs(match1.first->x - match1.second->x) >
       abs(match2.first->x - match2.second->x) )
    return match2;
    
  return match1;
}


float SoftMuonMvaEstimatorRun3::dX(const MatchPair& match) const{
  if (match.first and match.second->hasPhi())
    return (match.first->x - match.second->x);
  else
    return 9999.;
}

float SoftMuonMvaEstimatorRun3::pullX(const MatchPair& match) const{
  if (match.first and match.second->hasPhi())
    return dX(match) /
      sqrt(pow(match.first->xErr, 2) + pow(match.second->xErr, 2));
  else
    return 9999.;
}

float SoftMuonMvaEstimatorRun3::pullDxDz(const MatchPair& match) const{
  if (match.first and match.second->hasPhi())
    return (match.first->dXdZ - match.second->dXdZ) /
           sqrt(pow(match.first->dXdZErr, 2) + pow(match.second->dXdZErr, 2));
  else
    return 9999.;
}

float SoftMuonMvaEstimatorRun3::dY(const MatchPair& match) const{
  if (match.first and match.second->hasZed())
    return (match.first->y - match.second->y);
  else
    return 9999.;
}

float SoftMuonMvaEstimatorRun3::pullY(const MatchPair& match) const{
  if (match.first and match.second->hasZed())
    return dY(match) /
      sqrt(pow(match.first->yErr, 2) + pow(match.second->yErr, 2));
  else
    return 9999.;
}

float SoftMuonMvaEstimatorRun3::pullDyDz(const MatchPair& match) const{
  if (match.first and match.second->hasZed())
    return (match.first->dYdZ - match.second->dYdZ) /
           sqrt(pow(match.first->dYdZErr, 2) + pow(match.second->dYdZErr, 2));
  else
    return 9999.;
}

const reco::Muon::ArbitrationType arbitrationType = reco::Muon::SegmentAndTrackArbitration;
std::vector<float> SoftMuonMvaEstimatorRun3::computeMVAID(const pat::Muon &muon) const {


  const float eta = muon.eta();
  const float charge = muon.charge();

  float chargeProduct = 0;
  float staValidHits  = 0;
  float glbNormChi2 = 9999.;
  float staNormChi2 = 9999.;
   if (muon.isGlobalMuon()) {
       chargeProduct = muon.innerTrack()->charge()*muon.outerTrack()->charge();
       staValidHits = muon.outerTrack()->hitPattern().muonStationsWithValidHits();
       glbNormChi2 = muon.globalTrack()->normalizedChi2();
       staNormChi2 = muon.outerTrack()->normalizedChi2();
  }

  const float isGlobal = muon.isGlobalMuon();
  const float isTracker = muon.isTrackerMuon();
  const float isStandalone = muon.isStandAloneMuon();
  const float isPF = muon.isPFMuon();

  const float trkKink = muon.combinedQuality().trkKink;
  const float glbTrackProbability = muon.combinedQuality().glbTrackProbability;
  const float chi2LocalPosition = muon.combinedQuality().chi2LocalPosition;
  const float chi2LocalMomentum = muon.combinedQuality().chi2LocalMomentum;
  const float trkRelChi2 = muon.combinedQuality().trkRelChi2;
  const float staRelChi2 = muon.combinedQuality().staRelChi2;
  const float nStations = muon.numberOfMatchedStations();
  const float segmentComp = muon.segmentCompatibility(arbitrationType);

  float trkValidFrac = 0.;
  float trkNormChi2 = 9999;
  float nPixels = 0;
  float nValidHits = 0;
  float nLostHitsInner = 0;
  float nLostHitsOn = 0;
  float nLostHitsOuter = 0;
  float trkLayers = 0;
  float trkLostLayersInner = 0;
  float trkLostLayersOn = 0;
  float trkLostLayersOuter = 0;
  float highPurity = 0;

  if (muon.isTrackerMuon() or muon.isGlobalMuon()){

	trkValidFrac = muon.innerTrack()->validFraction();
	trkNormChi2 = muon.innerTrack()->normalizedChi2();
	
	nPixels = muon.innerTrack()->hitPattern().numberOfValidPixelHits();
	nValidHits = muon.innerTrack()->hitPattern().numberOfValidTrackerHits();
	nLostHitsInner = muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS);
	nLostHitsOn = muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::TRACK_HITS);
	nLostHitsOuter = muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_OUTER_HITS);
	trkLostLayersInner = muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS);
	trkLostLayersOn = muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
	trkLostLayersOuter = muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_OUTER_HITS);

	highPurity = muon.innerTrack()->quality(reco::Track::highPurity);
  }

  // do matching magic
  const int n_stations = 2;
  std::vector<MatchPair> matches;
  for (unsigned int i=0; i < n_stations; ++i)
    matches.push_back(std::pair(nullptr, nullptr));

  for (auto& chamberMatch : muon.matches()){
    unsigned int station = chamberMatch.station() - 1;
    if (station >= n_stations) continue;

    for (auto& segmentMatch : chamberMatch.segmentMatches){
      if ( not segmentMatch.isMask(reco::MuonSegmentMatch::BestInStationByDR) ||
	   not segmentMatch.isMask(reco::MuonSegmentMatch::BelongsToTrackByDR) )
	continue;


      auto match_pair = MatchPair(&chamberMatch, &segmentMatch);
      
      if (matches[station].first)
	matches[station] = getBetterMatch(matches[station], match_pair);
      else
	matches[station] = match_pair;
    }
  }


  float match1_dX = dX(matches[0]);
  float match1_pullX = pullX(matches[0]);
  float match1_pullDxDz = pullDxDz(matches[0]);
  float match1_dY = dY(matches[0]);
  float match1_pullY = pullY(matches[0]);
  float match1_pullDyDz = pullDyDz(matches[0]);
  float match2_dX = dX(matches[0]);
  float match2_pullX = pullX(matches[0]);
  float match2_pullDxDz = pullDxDz(matches[0]);
  float match2_dY = dY(matches[0]);
  float match2_pullY = pullY(matches[0]);
  float match2_pullDyDz = pullDyDz(matches[0]);


  const std::vector<std::string> input_names_{"X"};
  std::vector<float> vars = {eta,
                             charge,
                             chargeProduct,
                             isGlobal,
                             isTracker,
                             isStandalone,
                             isPF,
                             trkKink,
                             glbTrackProbability,
                             chi2LocalPosition,
                             glbNormChi2,
                             trkValidFrac,
			     chi2LocalMomentum,
                             trkRelChi2,
                             staRelChi2,
                             trkNormChi2,
                             staNormChi2,
                             nStations,
                             match1_dX,
                             match1_pullX,
                             match1_pullDxDz,
                             match1_dY,
                             match1_pullY,
                             match1_pullDyDz,
                             match2_dX,
                             match2_pullX,
                             match2_pullDxDz,
                             match2_dY,
                             match2_pullY,
                             match2_pullDyDz,
                             nPixels,
                             segmentComp,
                             nValidHits,
                             nLostHitsInner,
                             nLostHitsOn,
                             nLostHitsOuter,
                             trkLayers,
                             trkLostLayersInner,
                             trkLostLayersOn,
                             trkLostLayersOuter,
                             staValidHits,
                             highPurity};


  const std::vector<std::string> flav_names_{"probBAD", "probGOOD"};
  cms::Ort::FloatArrays input_values_;
  input_values_.emplace_back(vars);
  std::vector<float> outputs;
  LogDebug("MuonMvaIDEstimator") << hgb_.get();
  outputs = hgb_->run(input_names_, input_values_, {}, {"probabilities"})[0];
  assert(outputs.size() == flav_names_.size());
  return outputs;
}
