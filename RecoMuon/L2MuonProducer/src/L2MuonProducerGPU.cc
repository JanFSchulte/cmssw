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
#include "CUDADataFormats/CSCRecHit/interface/CSCSegmentCUDA.h"
#include "CUDADataFormats/DTRecHit/interface/DTRecSegment4DCUDA.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentContainerCUDA.h"
#include "RecoMuon/L2MuonProducer/src/L2MuonGeneratorOnGPU.h"
class L2MuonProducerGPU : public edm::global::EDProducer<> {
public:
  explicit L2MuonProducerGPU(const edm::ParameterSet& iConfig);
  ~L2MuonProducerGPU() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

  edm::EDGetTokenT<cms::cuda::Product<CSCSegmentCUDA>> tokenCSCSegmentsGPU_;
  edm::EDGetTokenT<cms::cuda::Product<DTRecSegment4DCUDA>>  tokenDTSegmentsGPU_;
  edm::EDPutTokenT<cms::cuda::Product<L2MuonTrackHeterogeneous>> tokenTrackGPU_;

  L2MuonGeneratorOnGPU gpuAlgo_;

};

L2MuonProducerGPU::L2MuonProducerGPU(const edm::ParameterSet& iConfig):
    gpuAlgo_(iConfig, consumesCollector()) {
    tokenCSCSegmentsGPU_ =
        consumes<cms::cuda::Product<CSCSegmentCUDA>>(iConfig.getParameter<edm::InputTag>("cscSegmentsSource"));
    tokenDTSegmentsGPU_ =
        consumes<cms::cuda::Product<DTRecSegment4DCUDA>>(iConfig.getParameter<edm::InputTag>("dtSegmentsSource"));
    tokenTrackGPU_ = produces<cms::cuda::Product<L2MuonTrackHeterogeneous>>();
}

void L2MuonProducerGPU::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("cscSegmentsSource", edm::InputTag("hltCSCSegmentsToCUDA"));
  desc.add<edm::InputTag>("dtSegmentsSource", edm::InputTag("hltDTSegmentsToCUDA"));

  L2MuonGeneratorOnGPU::fillDescriptions(desc);

  descriptions.add("L2MuonProducerGPU", desc);
}

void L2MuonProducerGPU::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& es) const {

    edm::Handle<cms::cuda::Product<DTRecSegment4DCUDA>> dtSegments;
    iEvent.getByToken(tokenDTSegmentsGPU_, dtSegments);

    cms::cuda::ScopedContextProduce ctx{*dtSegments};
    auto const& dtSegments_h = ctx.get(*dtSegments);

    edm::Handle<cms::cuda::Product<CSCSegmentCUDA>> cscSegments;
    iEvent.getByToken(tokenCSCSegmentsGPU_, cscSegments);

    auto const& cscSegments_h = ctx.get(*cscSegments);

    ctx.emplace(iEvent, tokenTrackGPU_, gpuAlgo_.makeTuplesAsync(dtSegments_h, cscSegments_h, ctx.stream()));

}

DEFINE_FWK_MODULE(L2MuonProducerGPU);   

