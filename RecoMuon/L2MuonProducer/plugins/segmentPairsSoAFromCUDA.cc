#include <cuda_runtime.h>

#include "CUDADataFormats/Common/interface/Product.h"
#include "CUDADataFormats/Common/interface/HostProduct.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentPairsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentPairsHeterogeneous.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/CUDACore/interface/ScopedContext.h"

class segmentPairsSoAFromCUDA : public edm::stream::EDProducer<edm::ExternalWork> {
public:
  explicit segmentPairsSoAFromCUDA(const edm::ParameterSet& iConfig);
  ~segmentPairsSoAFromCUDA() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void acquire(edm::Event const& iEvent,
               edm::EventSetup const& iSetup,
               edm::WaitingTaskWithArenaHolder waitingTaskHolder) override;
  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override;

  edm::EDGetTokenT<cms::cuda::Product<MuonSegmentPairsHeterogeneous>> tokenCUDA_;
  edm::EDPutTokenT<MuonSegmentPairsHeterogeneous> tokenSOA_;

  cms::cuda::host::unique_ptr<MuonSegmentPairsCUDA> soa_;
};

segmentPairsSoAFromCUDA::segmentPairsSoAFromCUDA(const edm::ParameterSet& iConfig)
    : tokenCUDA_(consumes<cms::cuda::Product<MuonSegmentPairsHeterogeneous>>(iConfig.getParameter<edm::InputTag>("src"))),
      tokenSOA_(produces<MuonSegmentPairsHeterogeneous>()) {}

void segmentPairsSoAFromCUDA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src", edm::InputTag("pixelTracksCUDA"));
  descriptions.add("segmentPairsSoA", desc);
}


void segmentPairsSoAFromCUDA::acquire(edm::Event const& iEvent,
                                    edm::EventSetup const& iSetup,
                                    edm::WaitingTaskWithArenaHolder waitingTaskHolder) {
  cms::cuda::Product<MuonSegmentPairsHeterogeneous> const& inputDataWrapped = iEvent.get(tokenCUDA_);
  cms::cuda::ScopedContextAcquire ctx{inputDataWrapped, std::move(waitingTaskHolder)};
  auto const& inputData = ctx.get(inputDataWrapped);

  soa_ = inputData.toHostAsync(ctx.stream());
}

void segmentPairsSoAFromCUDA::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {
#ifdef PIXEL_DEBUG_PRODUCE
  auto const& tsoa = *soa_;
  auto maxTracks = tsoa.stride();
  std::cout << "size of SoA" << sizeof(tsoa) << " stride " << maxTracks << std::endl;

  int32_t nt = 0;
  for (int32_t it = 0; it < maxTracks; ++it) {
    auto nHits = tsoa.nHits(it);
    assert(nHits == int(tsoa.hitIndices.size(it)));
    if (nHits == 0)
      break;  // this is a guard: maybe we need to move to nTracks...
    nt++;
  }
  std::cout << "found " << nt << " tracks in cpu SoA at " << &tsoa << std::endl;
#endif

  iEvent.emplace(tokenSOA_, MuonSegmentPairsHeterogeneous(std::move(soa_)));

  assert(!soa_);
}

DEFINE_FWK_MODULE(segmentPairsSoAFromCUDA);


