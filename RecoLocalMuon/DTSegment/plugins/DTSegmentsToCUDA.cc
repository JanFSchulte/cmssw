// C++ includes
#include <memory>
#include <string>
#include <vector>

// CMSSW includes
#include "CUDADataFormats/Common/interface/Product.h"
#include "CUDADataFormats/DTRecHit/interface/DTRecSegment4DCUDA.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DContainerCUDA.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HeterogeneousCore/CUDACore/interface/ScopedContext.h"
#include "HeterogeneousCore/CUDAServices/interface/CUDAService.h"
#include "HeterogeneousCore/CUDAUtilities/interface/copyAsync.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_noncached_unique_ptr.h"

namespace {

  class DTSegmentHost {
  public:
    DTSegmentHost() : data_h_{cms::cuda::make_host_noncached_unique<DTRecSegment4DContainerCUDA>(cudaHostAllocWriteCombined)} {}

    DTSegmentHost(DTSegmentHost const&) = delete;
    DTSegmentHost(DTSegmentHost&&) = default;

    DTSegmentHost& operator=(DTSegmentHost const&) = delete;
    DTSegmentHost& operator=(DTSegmentHost&&) = default;

    DTRecSegment4DContainerCUDA* data() { return data_h_.get(); }
    DTRecSegment4DContainerCUDA const* data() const { return data_h_.get(); }

    cms::cuda::host::noncached::unique_ptr<DTRecSegment4DContainerCUDA>& ptr() { return data_h_; }
    cms::cuda::host::noncached::unique_ptr<DTRecSegment4DContainerCUDA> const& ptr() const { return data_h_; }

  private:
    cms::cuda::host::noncached::unique_ptr<DTRecSegment4DContainerCUDA> data_h_;
  };

}  // namespace

class DTSegmentsToCUDA : public edm::global::EDProducer<edm::StreamCache<DTSegmentHost>> {
public:
  explicit DTSegmentsToCUDA(const edm::ParameterSet& iConfig);
  ~DTSegmentsToCUDA() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::unique_ptr<DTSegmentHost> beginStream(edm::StreamID) const override {
    edm::Service<CUDAService> cs;
    if (cs->enabled()) {
      return std::make_unique<DTSegmentHost>();
    } else {
      return nullptr;
    }
  }
  void produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;
private:

  edm::EDGetTokenT<DTRecSegment4DCollection> segmentsGetToken_;
  edm::EDPutTokenT<cms::cuda::Product<DTRecSegment4DCUDA>> segmentsPutToken_;

};

DTSegmentsToCUDA::DTSegmentsToCUDA(const edm::ParameterSet& iConfig)
    : segmentsGetToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      segmentsPutToken_(produces<cms::cuda::Product<DTRecSegment4DCUDA>>()) {

}

void DTSegmentsToCUDA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltDTSegments"));
  descriptions.add("DTSegmentsToCUDA",desc);
}

void DTSegmentsToCUDA::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  cms::cuda::ScopedContextProduce ctx{streamID};

  const DTRecSegment4DCollection& segments = iEvent.get(segmentsGetToken_);

  auto& segmentHost = streamCache(streamID)->ptr();

  int index = 0;
  for (DTRecSegment4DCollection::const_iterator it = segments.begin(); it != segments.end(); it++) {
	  segmentHost->x[index] = (*it).parameters()[2];
	  segmentHost->y[index] = (*it).parameters()[3];
	  segmentHost->dxdz[index] = (*it).parameters()[0];
	  segmentHost->dydz[index] = (*it).parameters()[1];

	  segmentHost->sigmaX[index] = (*it).parametersError()[2][2];
	  segmentHost->sigmaY[index] = (*it).parametersError()[3][3];
	  segmentHost->sigmaDXDZ[index] = (*it).parametersError()[0][0];
	  segmentHost->sigmaDYDZ[index] = (*it).parametersError()[0][0];
	  index++;
  }
  segmentHost->nSegments = index;
  DTRecSegment4DCUDA segmentsDevice(ctx.stream());
  cms::cuda::copyAsync(segmentsDevice.ptr(), segmentHost, ctx.stream());

  ctx.emplace(iEvent, segmentsPutToken_, std::move(segmentsDevice));
} 

DEFINE_FWK_MODULE(DTSegmentsToCUDA);
