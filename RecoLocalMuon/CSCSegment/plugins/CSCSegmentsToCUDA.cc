// C++ includes
#include <memory>
#include <string>
#include <vector>

// CMSSW includes
#include "CUDADataFormats/Common/interface/Product.h"
#include "CUDADataFormats/CSCRecHit/interface/CSCSegmentCUDA.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentContainerCUDA.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
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

  class CSCSegmentHost {
  public:
    CSCSegmentHost() : data_h_{cms::cuda::make_host_noncached_unique<CSCSegmentContainerCUDA>(cudaHostAllocWriteCombined)} {}

    CSCSegmentHost(CSCSegmentHost const&) = delete;
    CSCSegmentHost(CSCSegmentHost&&) = default;

    CSCSegmentHost& operator=(CSCSegmentHost const&) = delete;
    CSCSegmentHost& operator=(CSCSegmentHost&&) = default;

    CSCSegmentContainerCUDA* data() { return data_h_.get(); }
    CSCSegmentContainerCUDA const* data() const { return data_h_.get(); }

    cms::cuda::host::noncached::unique_ptr<CSCSegmentContainerCUDA>& ptr() { return data_h_; }
    cms::cuda::host::noncached::unique_ptr<CSCSegmentContainerCUDA> const& ptr() const { return data_h_; }

  private:
    cms::cuda::host::noncached::unique_ptr<CSCSegmentContainerCUDA> data_h_;
  };

}  // namespace

class CSCSegmentsToCUDA : public edm::global::EDProducer<edm::StreamCache<CSCSegmentHost>> {
public:
  explicit CSCSegmentsToCUDA(const edm::ParameterSet& iConfig);
  ~CSCSegmentsToCUDA() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::unique_ptr<CSCSegmentHost> beginStream(edm::StreamID) const override {
    edm::Service<CUDAService> cs;
    if (cs->enabled()) {
      return std::make_unique<CSCSegmentHost>();
    } else {
      return nullptr;
    }
  }
  void produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;
private:

  edm::EDGetTokenT<CSCSegmentCollection> segmentsGetToken_;
  edm::EDPutTokenT<cms::cuda::Product<CSCSegmentCUDA>> segmentsPutToken_;

};

CSCSegmentsToCUDA::CSCSegmentsToCUDA(const edm::ParameterSet& iConfig)
    : segmentsGetToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      segmentsPutToken_(produces<cms::cuda::Product<CSCSegmentCUDA>>()) {

}

void CSCSegmentsToCUDA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltCSCSegments"));
  descriptions.add("CSCSegmentsToCUDA",desc);
}

void CSCSegmentsToCUDA::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  cms::cuda::ScopedContextProduce ctx{streamID};

  const CSCSegmentCollection& segments = iEvent.get(segmentsGetToken_);

  auto& segmentHost = streamCache(streamID)->ptr();

  int index = 0;
  for (CSCSegmentCollection::const_iterator it = segments.begin(); it != segments.end(); it++) {
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
  CSCSegmentCUDA segmentsDevice(ctx.stream());
  cms::cuda::copyAsync(segmentsDevice.ptr(), segmentHost, ctx.stream());

  ctx.emplace(iEvent, segmentsPutToken_, std::move(segmentsDevice));
} 

DEFINE_FWK_MODULE(CSCSegmentsToCUDA);
