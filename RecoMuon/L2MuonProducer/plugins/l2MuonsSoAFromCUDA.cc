#include <cuda_runtime.h>

#include "CUDADataFormats/Common/interface/Product.h"
#include "CUDADataFormats/Common/interface/HostProduct.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentNtupletsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentNtupletsHeterogeneous.h"
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

class l2MuonsSoAFromCUDA : public edm::stream::EDProducer<edm::ExternalWork> {
public:
  explicit l2MuonsSoAFromCUDA(const edm::ParameterSet& iConfig);
  ~l2MuonsSoAFromCUDA() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void acquire(edm::Event const& iEvent,
               edm::EventSetup const& iSetup,
               edm::WaitingTaskWithArenaHolder waitingTaskHolder) override;
  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override;

  edm::EDGetTokenT<cms::cuda::Product<MuonSegmentNtupletsHeterogeneous>> tokenCUDA_;
  edm::EDPutTokenT<MuonSegmentNtupletsHeterogeneous> tokenSOA_;

  cms::cuda::host::unique_ptr<MuonSegmentNtupletsCUDA> soa_;
};

l2MuonsSoAFromCUDA::l2MuonsSoAFromCUDA(const edm::ParameterSet& iConfig)
    : tokenCUDA_(consumes<cms::cuda::Product<MuonSegmentNtupletsHeterogeneous>>(iConfig.getParameter<edm::InputTag>("src"))),
      tokenSOA_(produces<MuonSegmentNtupletsHeterogeneous>()) {}

void l2MuonsSoAFromCUDA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src", edm::InputTag("pixelTracksCUDA"));
  descriptions.add("l2MuonsSoA", desc);
}


void l2MuonsSoAFromCUDA::acquire(edm::Event const& iEvent,
                                    edm::EventSetup const& iSetup,
                                    edm::WaitingTaskWithArenaHolder waitingTaskHolder) {
  cms::cuda::Product<MuonSegmentNtupletsHeterogeneous> const& inputDataWrapped = iEvent.get(tokenCUDA_);
  cms::cuda::ScopedContextAcquire ctx{inputDataWrapped, std::move(waitingTaskHolder)};
  auto const& inputData = ctx.get(inputDataWrapped);

  soa_ = inputData.toHostAsync(ctx.stream());
}

void l2MuonsSoAFromCUDA::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {
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

  iEvent.emplace(tokenSOA_, MuonSegmentNtupletsHeterogeneous(std::move(soa_)));

  assert(!soa_);
}

DEFINE_FWK_MODULE(l2MuonsSoAFromCUDA);


