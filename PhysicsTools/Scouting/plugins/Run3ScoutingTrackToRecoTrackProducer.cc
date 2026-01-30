// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "fastjet/contrib/SoftKiller.hh"

class Run3ScoutingTrackToRecoTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit Run3ScoutingTrackToRecoTrackProducer(const edm::ParameterSet &);
  ~Run3ScoutingTrackToRecoTrackProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
  void endStream() override {}

  void createTracks(edm::Handle<std::vector<Run3ScoutingTrack>> scoutingtrackHandle,
                          std::unique_ptr<reco::TrackCollection> &tracks);
  reco::Track createTrack(Run3ScoutingTrack scoutingtrack);

  void clearVars();

private:
  const edm::EDGetTokenT<std::vector<Run3ScoutingTrack>> input_scoutingtrack_token_;

  std::vector<int> nValidPixelHits_;
  std::vector<int> nTrackerLayersWithMeasurement_;
  std::vector<int> nValidStripHits_;

};

//
// constructors and destructor
//
Run3ScoutingTrackToRecoTrackProducer::Run3ScoutingTrackToRecoTrackProducer(
    const edm::ParameterSet &iConfig)
    : input_scoutingtrack_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingtrack"))) {
  //register products
  produces<reco::TrackCollection>();
  produces<edm::ValueMap<int>>("nValidPixelHits");
  produces<edm::ValueMap<int>>("nTrackerLayersWithMeasurement");
  produces<edm::ValueMap<int>>("nValidStripHits");
}

Run3ScoutingTrackToRecoTrackProducer::~Run3ScoutingTrackToRecoTrackProducer() = default;

reco::Track Run3ScoutingTrackToRecoTrackProducer::createTrack(Run3ScoutingTrack scoutingtrack) {


  reco::Track::Point v(scoutingtrack.tk_vx(), scoutingtrack.tk_vy(), scoutingtrack.tk_vz());
  reco::Track::Vector p(math::RhoEtaPhiVector(scoutingtrack.tk_pt(), scoutingtrack.tk_eta(), scoutingtrack.tk_phi()));

  reco::TrackBase::CovarianceMatrix cov;
  cov(0, 0) = pow(scoutingtrack.tk_qoverp_Error(), 2);
  cov(0, 1) = scoutingtrack.tk_qoverp_lambda_cov();
  cov(0, 2) = scoutingtrack.tk_qoverp_phi_cov();
  cov(0, 3) = scoutingtrack.tk_qoverp_dxy_cov();
  cov(0, 4) = scoutingtrack.tk_qoverp_dsz_cov();
  cov(1, 1) = pow(scoutingtrack.tk_lambda_Error(), 2);
  cov(1, 2) = scoutingtrack.tk_lambda_phi_cov();
  cov(1, 3) = scoutingtrack.tk_lambda_dxy_cov();
  cov(1, 4) = scoutingtrack.tk_lambda_dsz_cov();
  cov(2, 2) = pow(scoutingtrack.tk_phi_Error(), 2);
  cov(2, 3) = scoutingtrack.tk_phi_dxy_cov();
  cov(2, 4) = scoutingtrack.tk_phi_dsz_cov();
  cov(3, 3) = pow(scoutingtrack.tk_dxy_Error(), 2);
  cov(3, 4) = scoutingtrack.tk_dxy_dsz_cov();
  cov(4, 4) = pow(scoutingtrack.tk_dsz_Error(), 2);

  nValidPixelHits_.push_back(scoutingtrack.tk_nValidPixelHits());
  nTrackerLayersWithMeasurement_.push_back(scoutingtrack.tk_nTrackerLayersWithMeasurement());
  nValidStripHits_.push_back(scoutingtrack.tk_nValidStripHits());

  return reco::Track(scoutingtrack.tk_chi2(), scoutingtrack.tk_ndof(), v, p, scoutingtrack.tk_charge(), cov);

}

void Run3ScoutingTrackToRecoTrackProducer::createTracks(
    edm::Handle<std::vector<Run3ScoutingTrack>> scoutingtrackHandle,
    std::unique_ptr<reco::TrackCollection> &tracks) {
  for (unsigned int itrack = 0; itrack < scoutingtrackHandle->size(); ++itrack) {
    auto &scoutingtrack = (*scoutingtrackHandle)[itrack];


    auto track = createTrack(scoutingtrack);
    //if (track.p() != 0)
      tracks->push_back(track);
  }
}


// ------------ method called to produce the data  ------------
void Run3ScoutingTrackToRecoTrackProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup) {
  using namespace edm;


  Handle<std::vector<Run3ScoutingTrack>> scoutingtrackHandle;
  iEvent.getByToken(input_scoutingtrack_token_, scoutingtrackHandle);

  auto tracks = std::make_unique<reco::TrackCollection>();

  createTracks(scoutingtrackHandle, tracks);
  std::cout << tracks->size() << " " << nValidPixelHits_.size() << std::endl;
  edm::OrphanHandle<reco::TrackCollection> oh = iEvent.put(std::move(tracks));

  std::unique_ptr<edm::ValueMap<int>> nValidPixelHits_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_nValidPixelHits(*nValidPixelHits_VM);
  filler_nValidPixelHits.insert(oh, nValidPixelHits_.begin(), nValidPixelHits_.end());
  filler_nValidPixelHits.fill();
  iEvent.put(std::move(nValidPixelHits_VM), "nValidPixelHits");

  std::unique_ptr<edm::ValueMap<int>> nTrackerLayersWithMeasurement_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_nTrackerLayersWithMeasurement(*nTrackerLayersWithMeasurement_VM);
  filler_nTrackerLayersWithMeasurement.insert(oh, nTrackerLayersWithMeasurement_.begin(), nTrackerLayersWithMeasurement_.end());
  filler_nTrackerLayersWithMeasurement.fill();
  iEvent.put(std::move(nTrackerLayersWithMeasurement_VM), "nTrackerLayersWithMeasurement");

  std::unique_ptr<edm::ValueMap<int>> nValidStripHits_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_nValidStripHits(*nValidStripHits_VM);
  filler_nValidStripHits.insert(oh, nValidStripHits_.begin(), nValidStripHits_.end());
  filler_nValidStripHits.fill();
  iEvent.put(std::move(nValidStripHits_VM), "nValidStripHits");

  const edm::ValueMap<int>& nValidStripHitsMap = *nValidStripHits_VM;
  clearVars();
}


void Run3ScoutingTrackToRecoTrackProducer::clearVars() {
  nValidPixelHits_.clear();
  nTrackerLayersWithMeasurement_.clear();
  nValidStripHits_.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingTrackToRecoTrackProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingtrack", edm::InputTag("hltScoutingTrackPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingTrackToRecoTrackProducer);
