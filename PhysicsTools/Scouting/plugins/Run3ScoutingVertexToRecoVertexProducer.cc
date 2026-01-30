// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
//#include "fastjet/contrib/SoftKiller.hh"

class Run3ScoutingVertexToRecoVertexProducer : public edm::stream::EDProducer<> {
public:
    explicit Run3ScoutingVertexToRecoVertexProducer(const edm::ParameterSet &);
    ~Run3ScoutingVertexToRecoVertexProducer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
    void endStream() override {}

    reco::Vertex createVertex(Run3ScoutingVertex scoutingvertex);
    void createVertices(edm::Handle<std::vector<Run3ScoutingVertex>> scoutingvertexHandle,
                        std::unique_ptr<reco::VertexCollection> &vertices);

private:
    const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> input_scoutingvertex_token_; 
};

Run3ScoutingVertexToRecoVertexProducer::Run3ScoutingVertexToRecoVertexProducer(
    const edm::ParameterSet &iConfig)
    : input_scoutingvertex_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))){
    
    //register products
    produces<reco::VertexCollection>();
}

Run3ScoutingVertexToRecoVertexProducer::~Run3ScoutingVertexToRecoVertexProducer() = default;

reco::Vertex Run3ScoutingVertexToRecoVertexProducer::createVertex(Run3ScoutingVertex scoutingVertex){
    // fill point coordinate
    reco::Vertex::Point point(scoutingVertex.x(), scoutingVertex.y(), scoutingVertex.z());

    // fill error
    std::vector<float> error_vec(6);
    error_vec[0] = scoutingVertex.xError() * scoutingVertex.xError(); // cov(0, 0)
    error_vec[1] = 0; // cov(0, 1)
    error_vec[2] = 0; // cov(0, 2)
    error_vec[3] = scoutingVertex.yError() * scoutingVertex.yError(); // cov(1, 1)
    error_vec[4] = 0; // cov(1, 2)
    error_vec[5] = scoutingVertex.zError() * scoutingVertex.zError(); // cov(2, 2)

    // off-diagonal errors are added in the begining of 2024
    // see https://github.com/cms-sw/cmssw/pull/43758
    try {
      error_vec[1] = scoutingVertex.xyCov();
      error_vec[2] = scoutingVertex.xzCov();
      error_vec[4] = scoutingVertex.yzCov();
    }
    catch (...) { // do nothing
    }

    reco::Vertex::Error error(error_vec.begin(), error_vec.end());

    return scoutingVertex.isValidVtx() ? reco::Vertex(point, error, scoutingVertex.chi2(), scoutingVertex.ndof(), scoutingVertex.tracksSize()): reco::Vertex(point, error);
}

void Run3ScoutingVertexToRecoVertexProducer::createVertices(
    edm::Handle<std::vector<Run3ScoutingVertex>> scoutingvertexHandle,
    std::unique_ptr<reco::VertexCollection> &vertices) {

    for (unsigned int ivertex = 0; ivertex < scoutingvertexHandle->size(); ++ivertex) {
        auto &scoutingvertex = (*scoutingvertexHandle)[ivertex];
        vertices->push_back(createVertex(scoutingvertex));
    } 
}

void Run3ScoutingVertexToRecoVertexProducer::produce(edm::Event &iEvent, edm::EventSetup const &setep){
    using namespace edm;
    Handle<std::vector<Run3ScoutingVertex>> scoutingvertexHandle;
    iEvent.getByToken(input_scoutingvertex_token_, scoutingvertexHandle);

    auto vertices = std::make_unique<reco::VertexCollection>();
    createVertices(scoutingvertexHandle, vertices);
    iEvent.put(std::move(vertices));
}

void Run3ScoutingVertexToRecoVertexProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingVertexToRecoVertexProducer);
