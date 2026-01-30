#include <memory>
#include <set>

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/TrackVertexArbitration.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/TrackVertexArbitrationNoBeamSpot.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/TTHelpers.h"

//#define VTXDEBUG

inline const unsigned int nTracks(const reco::Vertex &sv) { return sv.nTracks(); }
inline const unsigned int nTracks(const reco::VertexCompositePtrCandidate &sv) {
  return sv.numberOfSourceCandidatePtrs();
}

template <class InputContainer, class VTX>
class TemplatedVertexArbitratorScouting : public edm::stream::EDProducer<> {
public:
  typedef std::vector<VTX> Product;
  TemplatedVertexArbitratorScouting(const edm::ParameterSet &params);

  static void fillDescriptions(edm::ConfigurationDescriptions &cdesc) {
    edm::ParameterSetDescription pdesc;
    pdesc.add<edm::InputTag>("primaryVertices", edm::InputTag("scoutingPrimaryVertexReco"));
    if (std::is_same<VTX, reco::Vertex>::value) {
      pdesc.add<edm::InputTag>("tracks", edm::InputTag("scoutingTrackReco"));
      pdesc.add<edm::InputTag>("secondaryVertices", edm::InputTag("vertexMergerScouting"));
    } else if (std::is_same<VTX, reco::VertexCompositePtrCandidate>::value) {
      pdesc.add<edm::InputTag>("tracks", edm::InputTag("scoutingPFCandidateReco"));
      pdesc.add<edm::InputTag>("secondaryVertices", edm::InputTag("candidateVertexMergerScouting"));
    } else {
      pdesc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
      pdesc.add<edm::InputTag>("secondaryVertices", edm::InputTag("vertexMerger"));
    }
    pdesc.add<edm::InputTag>("valueMapNValidPixelHits", edm::InputTag("scoutingTrackReco", "nValidPixelHits"));
    pdesc.add<edm::InputTag>("valueMapNTrackerLayersWithMeasurements", edm::InputTag("scoutingTrackReco", "nTrackerLayersWithMeasurement"));

    pdesc.add<double>("dLenFraction", 0.333);
    pdesc.add<double>("dRCut", 0.4);
    pdesc.add<double>("distCut", 0.04);
    pdesc.add<double>("sigCut", 5.0);
    pdesc.add<double>("fitterSigmacut", 3.0);
    pdesc.add<double>("fitterTini", 256);
    pdesc.add<double>("fitterRatio", 0.25);
    pdesc.add<int>("trackMinLayers", 4);
    pdesc.add<double>("trackMinPt", 0.4);
    pdesc.add<int>("trackMinPixels", 1);
    pdesc.add<double>("maxTimeSignificance", 3.5);
    pdesc.add<bool>("isScouting", false);
    if (std::is_same<VTX, reco::Vertex>::value) {
      cdesc.add("trackVertexArbitratorScoutingDefault", pdesc);
    } else if (std::is_same<VTX, reco::VertexCompositePtrCandidate>::value) {
      cdesc.add("candidateVertexArbitratorScoutingDefault", pdesc);
    } else {
      cdesc.addDefault(pdesc);
    }
  }

  void produce(edm::Event &event, const edm::EventSetup &es) override;

private:
  bool trackFilter(const reco::TrackRef &track) const;

  edm::EDGetTokenT<reco::VertexCollection> token_primaryVertex;
  edm::EDGetTokenT<Product> token_secondaryVertex;
  edm::EDGetTokenT<InputContainer> token_tracks;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> token_trackBuilder;
  edm::EDGetTokenT<edm::ValueMap<int> > nValidPixelHitsValueMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > nTrackerLayersWithMeasurementsValueMapToken_;

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> esTokenMF_;
  std::unique_ptr<TrackVertexArbitrationNoBeamSpot<VTX> > theArbitratorNoBS;

  bool isScouting;
};

template <class InputContainer, class VTX>
TemplatedVertexArbitratorScouting<InputContainer, VTX>::TemplatedVertexArbitratorScouting(const edm::ParameterSet &params) : esTokenMF_(esConsumes()) {
  token_primaryVertex = consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("primaryVertices"));
  token_secondaryVertex = consumes<Product>(params.getParameter<edm::InputTag>("secondaryVertices"));
  token_tracks = consumes<InputContainer>(params.getParameter<edm::InputTag>("tracks"));
  nValidPixelHitsValueMapToken_ = consumes<edm::ValueMap<int> >(params.getParameter<edm::InputTag>("valueMapNValidPixelHits"));
  nTrackerLayersWithMeasurementsValueMapToken_ = consumes<edm::ValueMap<int> >(params.getParameter<edm::InputTag>("valueMapNTrackerLayersWithMeasurements"));
  token_trackBuilder =
      esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  produces<Product>();
  theArbitratorNoBS.reset(new TrackVertexArbitrationNoBeamSpot<VTX>(params));
}

template <class InputContainer, class VTX>
void TemplatedVertexArbitratorScouting<InputContainer, VTX>::produce(edm::Event &event, const edm::EventSetup &es) {
  using namespace reco;

  edm::Handle<Product> secondaryVertices;
  event.getByToken(token_secondaryVertex, secondaryVertices);
  Product theSecVertexColl = *(secondaryVertices.product());

  edm::Handle<VertexCollection> primaryVertices;
  event.getByToken(token_primaryVertex, primaryVertices);

  const MagneticField* theMagneticField = &es.getData(esTokenMF_);


  auto recoVertices = std::make_unique<Product>();
  if (!primaryVertices->empty()) {
    const reco::Vertex &pv = (*primaryVertices)[0];

    edm::Handle<InputContainer> tracks;
    event.getByToken(token_tracks, tracks);

    edm::ESHandle<TransientTrackBuilder> trackBuilder = es.getHandle(token_trackBuilder);


    //        const edm::RefVector< TrackCollection > tracksForArbitration= selectedTracks;:/
    //
    Product theRecoVertices;

    edm::Handle<edm::ValueMap<int> > nValidPixelHitsValueMap;
    edm::Handle<edm::ValueMap<int> > nTrackerLayersWithMeasurementsValueMap;

    event.getByToken(nValidPixelHitsValueMapToken_, nValidPixelHitsValueMap);
    event.getByToken(nTrackerLayersWithMeasurementsValueMapToken_, nTrackerLayersWithMeasurementsValueMap);
    const edm::ValueMap<int>& nValidPixelHitsMap = *nValidPixelHitsValueMap;
    const edm::ValueMap<int>& nTrackerLayersWithMeasurementsMap = *nTrackerLayersWithMeasurementsValueMap;


    std::vector<int> n_pixel_hits;
    std::vector<int> n_tracker_layers;
    std::vector<TransientTrack> selectedTracks;
    for (typename InputContainer::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {

        reco::TrackRef ref;
        if constexpr (std::is_same_v<InputContainer, reco::TrackCollection>) {
            ref = reco::TrackRef(tracks, track - tracks->begin());
        }
        else if constexpr (std::is_same_v<InputContainer, edm::View<reco::Candidate>>) {
	      const reco::Candidate& cand = *track;
	      const reco::PFCandidate* tmpCand = dynamic_cast<const reco::PFCandidate*>(&cand);
	      ref = tmpCand->trackRef();
        }

        //else if constexpr (std::is_same_v<InputContainer, edm::View<reco::PFCandidate>>) {
	//   ref = track->trackRef();
       // }
       //
        int n_ph = 99;
        int n_tl = 99;
        if (ref.isNonnull()) {
	    n_ph = nValidPixelHitsMap[ref];
	    n_tl = nTrackerLayersWithMeasurementsMap[ref];
	}
	n_pixel_hits.push_back(n_ph);
	n_tracker_layers.push_back(n_tl);

        TransientTrack tt(tthelpers::buildTT(tracks, trackBuilder, track - tracks->begin()));
       // reco::TransientTrack tt(*ref, theMagneticField);
        selectedTracks.push_back(tt);
    }
    
    theRecoVertices = theArbitratorNoBS->trackVertexArbitratorNoBeamSpot(pv, selectedTracks, n_pixel_hits, n_tracker_layers, theSecVertexColl);
	

    for (unsigned int ivtx = 0; ivtx < theRecoVertices.size(); ivtx++) {
      if (!(nTracks(theRecoVertices[ivtx]) > 1))
        continue;
      recoVertices->push_back(theRecoVertices[ivtx]);
    }
  }
  event.put(std::move(recoVertices));
}

typedef TemplatedVertexArbitratorScouting<reco::TrackCollection, reco::Vertex> TrackVertexArbitratorScouting;
typedef TemplatedVertexArbitratorScouting<edm::View<reco::Candidate>, reco::VertexCompositePtrCandidate>
    CandidateVertexArbitratorScouting;

DEFINE_FWK_MODULE(TrackVertexArbitratorScouting);
DEFINE_FWK_MODULE(CandidateVertexArbitratorScouting);
