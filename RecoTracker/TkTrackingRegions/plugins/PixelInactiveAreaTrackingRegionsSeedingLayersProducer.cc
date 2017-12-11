#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TrackingTools/TransientTrackingRecHit/interface/SeedingLayerSetsHits.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionsSeedingLayerSets.h"

#include "VertexBeamspotOrigins.h"
#include "Candidates.h"
#include "AreaSeededTrackingRegionsBuilder.h"
#include "PixelInactiveAreaFinder.h"

#include <vector>
#include <utility>
#include <memory>

class PixelInactiveAreaTrackingRegionsSeedingLayersProducer: public edm::stream::EDProducer<> {
public:
  PixelInactiveAreaTrackingRegionsSeedingLayersProducer(const edm::ParameterSet& iConfig);
  ~PixelInactiveAreaTrackingRegionsSeedingLayersProducer() override = default;

  enum class MODE {GLOBAL, CANDIDATE_SEEDED};
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

private:
  SeedingLayerSetsBuilder seedingLayerSetsBuilder_;
  VertexBeamspotOrigins origins_;
  Candidates candidates_;
  PixelInactiveAreaFinder inactiveAreaFinder_;
  AreaSeededTrackingRegionsBuilder trackingRegionsBuilder_;
 };

PixelInactiveAreaTrackingRegionsSeedingLayersProducer::PixelInactiveAreaTrackingRegionsSeedingLayersProducer(const edm::ParameterSet& iConfig):
  seedingLayerSetsBuilder_(iConfig, consumesCollector()),
  origins_(iConfig.getParameter<edm::ParameterSet>("RegionPSet"), consumesCollector()),
  candidates_(iConfig.getParameter<edm::ParameterSet>("RegionPSet"), consumesCollector()),
  inactiveAreaFinder_(iConfig, seedingLayerSetsBuilder_.layers(), seedingLayerSetsBuilder_.seedingLayerSetsLooper(), consumesCollector()),
  trackingRegionsBuilder_(iConfig.getParameter<edm::ParameterSet>("RegionPSet"), consumesCollector())

{
  
  produces<SeedingLayerSetsHits>();
  produces<TrackingRegionsSeedingLayerSets>();
}

void PixelInactiveAreaTrackingRegionsSeedingLayersProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  edm::ParameterSetDescription descRegion;
  VertexBeamspotOrigins::fillDescriptions(descRegion);
  Candidates::fillDescriptions(descRegion);
  AreaSeededTrackingRegionsBuilder::fillDescriptions(descRegion);
  desc.add<edm::ParameterSetDescription>("RegionPSet", descRegion);

  PixelInactiveAreaFinder::fillDescriptions(desc);
  SeedingLayerSetsBuilder::fillDescriptions(desc);

  descriptions.add("pixelInactiveAreaTrackingRegionsAndSeedingLayers", desc);
}

void PixelInactiveAreaTrackingRegionsSeedingLayersProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto orphanHandle = iEvent.put(seedingLayerSetsBuilder_.hits(iEvent, iSetup));
  const SeedingLayerSetsHits *seedingLayers = orphanHandle.product();

  auto regions = std::make_unique<TrackingRegionsSeedingLayerSets>(seedingLayers);

  const auto origins = origins_.origins(iEvent);
  const auto candidates = candidates_.objects(iEvent);
  int n_objects = 0;
  if( candidates.first.isValid() ){
    n_objects = candidates.first->size();
  }
  const auto builder = trackingRegionsBuilder_.beginEvent(iEvent);

  const auto allAreas = inactiveAreaFinder_.inactiveAreas(iEvent, iSetup);
  for(const auto& origin: origins) {
    auto areasLayerSets = allAreas.areasAndLayerSets(origin.first, origin.second); // point, half length in z
    LogTrace("PixelInactiveAreaTrackingRegionsSeedingLayersProducer") << "Origin " << origin.first.x() << "," << origin.first.y() << "," << origin.first.z() << " z half lengh " << origin.second;
    for(auto& areasLayerSet: areasLayerSets) {
      auto region = builder.region(origin, areasLayerSet.first);

      auto etaPhiRegion = dynamic_cast<const RectangularEtaPhiTrackingRegion *>(region.get());
#ifdef EDM_ML_DEBUG
      std::stringstream ss;
      for(const auto& ind: areasLayerSet.second) {
        ss << ind << ",";
      
      LogTrace("PixelInactiveAreaTrackingRegionsSeedingLayersProducer") << " region eta,phi " << region->direction().eta() << "," << region->direction().phi() << " eta range " << etaPhiRegion->etaRange().min() << "," << etaPhiRegion->etaRange().max() << " phi range " << (region->direction().phi()-etaPhiRegion->phiMargin().left()) << "," << (region->direction().phi()+etaPhiRegion->phiMargin().right()) << " layer sets " << ss.str();
#endif
      bool store = true;
      if (n_objects > 0){
		store = false;
		float dEta_Cand = candidates.second.first;
		float dPhi_Cand = candidates.second.second;
		double eta_Point = region->direction().eta();
          	double phi_Point = region->direction().phi();
		double m_deltaEta_Point = std::max(abs(region->direction().eta()-etaPhiRegion->etaRange().min()),abs(abs(region->direction().eta()-etaPhiRegion->etaRange().max())));
 		double m_deltaPhi_Point = std::max(etaPhiRegion->phiMargin().left(),etaPhiRegion->phiMargin().right());
	      	for(const auto& object : *candidates.first) {
			 
       			 double eta_Cand = object.eta();
                	 double phi_Cand = object.phi();
          		 double dEta_Cand_Point = std::abs(eta_Cand-eta_Point);
          		 double dPhi_Cand_Point = std::abs(deltaPhi(phi_Cand,phi_Point));
			 if(dEta_Cand_Point > (dEta_Cand + m_deltaEta_Point) || dPhi_Cand_Point > (dPhi_Cand + m_deltaPhi_Point)) continue;
			 store = true;
		}
      }
      if (store) regions->emplace_back(std::move(region), std::move(areasLayerSet.second));
    }
  }

  iEvent.put(std::move(regions));
}

DEFINE_FWK_MODULE(PixelInactiveAreaTrackingRegionsSeedingLayersProducer);
