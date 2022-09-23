#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonEDFilter.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"

#include "DataFormats/Math/interface/deltaR.h"

using namespace l1t;

class HLTTau3MuSonicFilter : public TritonEDFilter<> {
public:
  explicit HLTTau3MuSonicFilter(const edm::ParameterSet&);
  void acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) override;
  bool filter(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const edm::EDGetTokenT< EMTFHitCollection >   mu_hitToken;
  const double gnn_threshold_;
  const unsigned int max_n_hits_;
  const unsigned int max_n_edges_;
  const bool zeropad_;
};

HLTTau3MuSonicFilter::HLTTau3MuSonicFilter(const edm::ParameterSet& cfg)
    //: TritonEDFilter<>(cfg, "HLTTau3MuSonicFilter"),
    : TritonEDFilter<>(cfg),
      mu_hitToken (consumes< EMTFHitCollection > (cfg.getParameter<edm::InputTag>("L1EMTFHitInputTag"))),
      gnn_threshold_(cfg.getParameter<double>("gnn_threshold")),
      max_n_hits_(cfg.getParameter<unsigned int>("max_n_hits")),
      max_n_edges_(cfg.getParameter<unsigned int>("max_n_edges")),
      zeropad_(cfg.getParameter<bool>("zeropad")) {}


void HLTTau3MuSonicFilter::acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) {

    edm::Handle<EMTFHitCollection> l1muhitsH;
    iEvent.getByToken(mu_hitToken, l1muhitsH);  
    const EMTFHitCollection& l1muhits = (*l1muhitsH.product());

    auto& nodes = iInput.at("x__0");
    auto nodedata = nodes.allocate<float>();
    auto& vnodedata = (*nodedata)[0];

    auto& edgeIndices = iInput.at("edge_index__1");
    auto edgeIndexData = edgeIndices.allocate<int64_t>();
    auto& vEdgeIndexData = (*edgeIndexData)[0];

    auto& edgeFeatures = iInput.at("edge_attr__2");
    auto edgeFeaturesData = edgeFeatures.allocate<float>();
    auto& vEdgeFeaturesData = (*edgeFeaturesData)[0];




    std::vector<float> hit_z;
    std::vector<float> hit_eta;
    std::vector<float> hit_phi;
    std::vector<float> hit_bend;

    size_t i_hits = 1;

    for (const auto& hit : l1muhits)
    {    
        if ((hit.Station() != 1) || (hit.Neighbor() != 0) || (hit.Subsystem() == 0)) continue;
        vnodedata.push_back(hit.Z_sim());
        vnodedata.push_back(hit.Eta_sim());
        vnodedata.push_back(hit.Bend());

        hit_z.push_back(hit.Z_sim());
        hit_eta.push_back(hit.Eta_sim());
        hit_phi.push_back(hit.Phi_sim());
        hit_bend.push_back(hit.Bend());

        ++i_hits;
        if (i_hits == max_n_hits_) {
            break;  // output a warning?
        }
    }
    //add virtual global node
    vnodedata.push_back(0);
    vnodedata.push_back(0);
    vnodedata.push_back(0);

    hit_z.push_back(0);
    hit_eta.push_back(0);
    hit_phi.push_back(0);
    hit_bend.push_back(0);



    size_t i_edges = 0;
    // starting both loops at zero ensures that the graph is undirected and has self-loops
    for (size_t index1 = 0; index1 < i_hits; index1++){
        for (size_t index2 = 0; index2 < i_hits; index2++){

            if (deltaR(hit_eta.at(index1),hit_phi.at(index1),hit_eta.at(index2),hit_phi.at(index2)) > 1.0 && index1 !=0 && index2 !=0) continue;
            vEdgeIndexData.push_back(index1);            
            vEdgeIndexData.push_back(index2);            

            vEdgeFeaturesData.push_back(hit_z.at(index1)-hit_z.at(index2));
            vEdgeFeaturesData.push_back(hit_eta.at(index1)-hit_eta.at(index2));
            vEdgeFeaturesData.push_back(hit_bend.at(index1)-hit_bend.at(index2));
            vEdgeFeaturesData.push_back(hit_phi.at(index1)-hit_phi.at(index2));

            ++i_edges;
            if (i_edges == max_n_edges_) {
                break;  // output a warning?
            }
        } 
    }

    nodes.setShape(0,i_hits);
    edgeIndices.setShape(1,i_edges);
    edgeFeatures.setShape(0,i_edges);

    //std::cout << "found " << i_hits << " hits and " << i_edges << " edges" << std::endl;
    if (zeropad_){
        vnodedata.resize(3 * max_n_hits_);
        vnodedata.resize(2 * max_n_edges_);
        vnodedata.resize(4 * max_n_edges_);
    }
    nodes.toServer(nodedata);
    //std::cout << "nodes to server" << std::endl;
    edgeIndices.toServer(edgeIndexData);
    //std::cout << "edge indices to server" << std::endl;
    edgeFeatures.toServer(edgeFeaturesData);
    //std::cout << "edge features to server, running inference now" << std::endl;
}

bool HLTTau3MuSonicFilter::filter(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) {

    const auto& output1 = iOutput.begin()->second;
    const auto& outputs = output1.fromServer<float>();
    //std::cout << "GNN score: "  << outputs[0][0] << std::endl;
    return (outputs[0][0] > gnn_threshold_);

}

void HLTTau3MuSonicFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<edm::InputTag>("L1EMTFHitInputTag",edm::InputTag("simEmtfDigis"));
  desc.add<bool>("zeropad", false);
  desc.add<double>("gnn_threshold", 0.5);
  desc.add<unsigned int>("max_n_hits", 999999);
  desc.add<unsigned int>("max_n_edges", 999999);
  descriptions.add("hltTau3MuSonicFilter", desc);
}
// register as framework plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTTau3MuSonicFilter);
