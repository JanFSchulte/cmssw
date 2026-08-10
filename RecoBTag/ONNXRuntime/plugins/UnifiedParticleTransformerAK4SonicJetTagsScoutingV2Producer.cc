#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/BTauReco/interface/UnifiedParticleTransformerAK4TagInfo.h"
#include "DataFormats/BTauReco/interface/UnifiedParticleTransformerAK4Features.h"

#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonData.h"

#include "RecoBTag/ONNXRuntime/interface/tensor_fillers.h"
#include "RecoBTag/ONNXRuntime/interface/tensor_configs.h"

class UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer : public TritonEDProducer<> {
public:
  explicit UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer(const edm::ParameterSet &);
  ~UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer() override;

  void acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput) override;

  void produce(edm::Event &iEvent, edm::EventSetup const &iSetup, Output const &iOutput) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &);

private:
  typedef std::vector<reco::UnifiedParticleTransformerAK4TagInfo> TagInfoCollection;
  typedef reco::JetTagCollection JetTagCollection;

  const edm::EDGetTokenT<TagInfoCollection> src_;
  std::vector<std::string> flav_names_;
  std::vector<std::string> input_names_;
  std::vector<std::string> output_names_;

  bool skippedInference_ = false;

  // the scouting v2 model only takes 4 real inputs (charged candidates, neutral candidates,
  // vertices, lost tracks) -- the *4Vec groups defined in ScoutingUparT::InputFeatures are unused,
  // matching UnifiedParticleTransformerAK4ONNXJetTagsScoutingv2Producer.cc
  static const std::vector<ScoutingUparT::InputFeatures> active_features_;
};

const std::vector<ScoutingUparT::InputFeatures> UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer::active_features_{
    ScoutingUparT::kChargedCandidates,
    ScoutingUparT::kNeutralCandidates,
    ScoutingUparT::kVertices,
    ScoutingUparT::kLostTracks};

UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer::UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer(
    const edm::ParameterSet &iConfig)
    : TritonEDProducer<>(iConfig),
      src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      flav_names_(iConfig.getParameter<std::vector<std::string>>("flav_names")),
      input_names_(iConfig.getParameter<std::vector<std::string>>("input_names")),
      output_names_(iConfig.getParameter<std::vector<std::string>>("output_names")) {
  // get output names from flav_names
  for (const auto &flav_name : flav_names_) {
    produces<JetTagCollection>(flav_name);
  }
}

UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer::~UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer() {}

void UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer::fillDescriptions(
    edm::ConfigurationDescriptions &descriptions) {
  // pfUnifiedParticleTransformerAK4SonicJetTagsScoutingV2
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<edm::InputTag>("src", edm::InputTag("pfUnifiedParticleTransformerAK4TagInfos"));
  desc.add<std::vector<std::string>>("input_names", {"input_1", "input_2", "input_3", "input_4"});
  desc.add<std::vector<std::string>>("output_names", {"softmax"});
  desc.add<std::vector<std::string>>(
      "flav_names",
      std::vector<std::string>{"probb",
                               "probbb",
                               "problepb",
                               "probc",
                               "probs",
                               "probu",
                               "probd",
                               "probg",
                               "ptcorr",
                               "ptreshigh",
                               "ptreslow"});

  descriptions.add("pfUnifiedParticleTransformerAK4SonicJetTagsScoutingV2", desc);
}

void UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer::acquire(edm::Event const &iEvent,
                                                                          edm::EventSetup const &iSetup,
                                                                          Input &iInput) {
  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(src_, tag_infos);
  client_->setBatchSize(tag_infos->size());
  skippedInference_ = false;
  if (tag_infos->empty())
    return;

  // Unlike the non-scouting (V01) UParT model, the scouting v2 ONNX models were exported with a fixed
  // (non-dynamic) candidate-count axis per input -- see get_input_sizes()'s use_dynamic_axes_==false
  // branch in UnifiedParticleTransformerAK4ONNXJetTagsScoutingv2Producer.cc, and the Triton config.pbtxt
  // in RecoBTag/CombinedScouting/data/models/, which declare fixed dims (32/32/4/8) rather than -1. So,
  // unlike the non-scouting Sonic producer, the tensors sent to the server must always be padded up to
  // the full accepted size -- sizing per-event to the actual max candidate count (as the non-scouting
  // producer does) would send a shape the server rejects whenever an event's jets have fewer candidates
  // than the accepted maximum.
  bool any_candidates = false;
  for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
    const auto &features = ((*tag_infos)[jet_n]).features();
    if (!features.c_pf_features.empty() || !features.n_pf_features.empty() || !features.sv_features.empty() ||
        !features.lt_features.empty()) {
      any_candidates = true;
      break;
    }
  }

  // If an event has no jet, or all jets have zero n_cpf, n_npf, n_vtx and n_lt, the inference is skipped.
  if (!any_candidates) {
    client_->setBatchSize(0);
    skippedInference_ = true;
    return;
  }

  const unsigned int target_n_cpf = ScoutingUparT::n_cpf_accept;
  const unsigned int target_n_npf = ScoutingUparT::n_npf_accept;
  const unsigned int target_n_vtx = ScoutingUparT::n_sv_accept;
  const unsigned int target_n_lt = ScoutingUparT::n_lt_accept;
  const std::map<ScoutingUparT::InputFeatures, unsigned int> target_n{{ScoutingUparT::kChargedCandidates, target_n_cpf},
                                                                      {ScoutingUparT::kNeutralCandidates, target_n_npf},
                                                                      {ScoutingUparT::kVertices, target_n_vtx},
                                                                      {ScoutingUparT::kLostTracks, target_n_lt}};

  for (const auto ifeature : active_features_) {
    const auto &group_name = input_names_[ifeature];
    auto &input = iInput.at(group_name);
    input.setShape(0, target_n.at(ifeature));
    auto tdata = input.allocate<float>(true);
    for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
      const auto &taginfo = (*tag_infos)[jet_n];
      const auto &features = taginfo.features();
      auto &vdata = (*tdata)[jet_n];

      if (ifeature == ScoutingUparT::kChargedCandidates)
        ScoutingUparT_tensor_filler(vdata, ifeature, features.c_pf_features, target_n_cpf);
      else if (ifeature == ScoutingUparT::kNeutralCandidates)
        ScoutingUparT_tensor_filler(vdata, ifeature, features.n_pf_features, target_n_npf);
      else if (ifeature == ScoutingUparT::kVertices)
        ScoutingUparT_tensor_filler(vdata, ifeature, features.sv_features, target_n_vtx);
      else if (ifeature == ScoutingUparT::kLostTracks)
        ScoutingUparT_tensor_filler(vdata, ifeature, features.lt_features, target_n_lt);
    }
    input.toServer(tdata);
  }
}

void UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer::produce(edm::Event &iEvent,
                                                                          const edm::EventSetup &iSetup,
                                                                          Output const &iOutput) {
  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(src_, tag_infos);

  // initialize output collection
  std::vector<std::unique_ptr<JetTagCollection>> output_tags;
  if (!tag_infos->empty()) {
    auto jet_ref = tag_infos->begin()->jet();
    auto ref2prod = edm::makeRefToBaseProdFrom(jet_ref, iEvent);
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>(ref2prod));
    }
  } else {
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>());
    }
  }
  if (!tag_infos->empty()) {
    if (!skippedInference_) {
      const auto &output1 = iOutput.begin()->second;
      const auto &outputs_from_server = output1.fromServer<float>();

      for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
        const auto &taginfo = (*tag_infos)[jet_n];
        const auto &jet_ref = tag_infos->at(jet_n).jet();

        if (taginfo.features().is_filled) {
          for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
            (*(output_tags[flav_n]))[jet_ref] = outputs_from_server[jet_n][flav_n];
          }
        } else {
          for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
            (*(output_tags[flav_n]))[jet_ref] = -1.0;
          }
        }
      }
    } else {
      for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
        const auto &jet_ref = tag_infos->at(jet_n).jet();
        for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
          (*(output_tags[flav_n]))[jet_ref] = -1.0;
        }
      }
    }
  }
  // put into the event
  for (std::size_t flav_n = 0; flav_n < flav_names_.size(); ++flav_n) {
    iEvent.put(std::move(output_tags[flav_n]), flav_names_[flav_n]);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(UnifiedParticleTransformerAK4SonicJetTagsScoutingV2Producer);
