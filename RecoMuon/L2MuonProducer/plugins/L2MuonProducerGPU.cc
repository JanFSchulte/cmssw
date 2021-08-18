#include <cuda_runtime.h>

#include "CUDADataFormats/Common/interface/Product.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/RunningAverage.h"
#include "HeterogeneousCore/CUDACore/interface/ScopedContext.h"


#include "CUDADataFormats/Track/interface/L2MuonTrackHeterogeneous.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDA.h"
#include "RecoMuon/L2MuonProducer/plugins/L2MuonGeneratorOnGPU.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"

class L2MuonProducerGPU : public edm::global::EDProducer<> {
public:
  explicit L2MuonProducerGPU(const edm::ParameterSet& iConfig);
  ~L2MuonProducerGPU() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

  edm::EDGetTokenT<cms::cuda::Product<MuonSegmentsCUDA>> tokenSegmentsGPU_;
  edm::EDPutTokenT<cms::cuda::Product<L2MuonTrackHeterogeneous>> tokenTrackGPU_;

  L2MuonGeneratorOnGPU gpuAlgo_;

};

L2MuonProducerGPU::L2MuonProducerGPU(const edm::ParameterSet& iConfig):
    gpuAlgo_(iConfig, consumesCollector()) {
    tokenSegmentsGPU_ =
        consumes<cms::cuda::Product<MuonSegmentsCUDA>>(iConfig.getParameter<edm::InputTag>("muonSegmentsSource"));
    tokenTrackGPU_ = produces<cms::cuda::Product<L2MuonTrackHeterogeneous>>();
}

void L2MuonProducerGPU::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("muonSegmentsSource", edm::InputTag("hltMuonSegmentsToCUDA"));

  L2MuonGeneratorOnGPU::fillDescriptions(desc);

  descriptions.add("L2MuonProducerGPU", desc);
}

void L2MuonProducerGPU::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& es) const {

    auto bf = 1. / PixelRecoUtilities::fieldInInvGev(es);
     
    edm::Handle<cms::cuda::Product<MuonSegmentsCUDA>> muonSegments;
    iEvent.getByToken(tokenSegmentsGPU_, muonSegments);

    cms::cuda::ScopedContextProduce ctx{*muonSegments};
    auto const& muonSegments_h = ctx.get(*muonSegments);


    ctx.emplace(iEvent, tokenTrackGPU_, gpuAlgo_.makeTuplesAsync(muonSegments_h, bf, ctx.stream()));

}

DEFINE_FWK_MODULE(L2MuonProducerGPU);   

