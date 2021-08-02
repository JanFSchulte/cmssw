#ifndef RecoMuon_L2MuonProducer_src_L2MuonGeneratorOnGPU_h
#define RecoMuon_L2MuonProducer_src_L2MuonGeneratorOnGPU_h

#include <cuda_runtime.h>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

//#include "CUDADataFormats/TrackingRecHit/interface/TrackingRecHit2DHeterogeneous.h"
#include "CUDADataFormats/CSCRecHit/interface/CSCSegmentCUDA.h"
#include "CUDADataFormats/DTRecHit/interface/DTRecSegment4DCUDA.h"
#include "CUDADataFormats/Track/interface/L2MuonTrackHeterogeneous.h"

#include "RecoMuon/L2MuonProducer/src/L2MuonGeneratorKernels.h"

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSetDescription;
}  // namespace edm




class L2MuonGeneratorOnGPU {

public:

  L2MuonGeneratorOnGPU(const edm::ParameterSet& cfg, edm::ConsumesCollector&& iC)
      : L2MuonGeneratorOnGPU(cfg, iC) {}
  L2MuonGeneratorOnGPU(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC);

  ~L2MuonGeneratorOnGPU();

  static void fillDescriptions(edm::ParameterSetDescription& desc);
  static const char* fillDescriptionsLabel() { return "l2MuonGeneratorOnGPU"; }

  L2MuonTrackHeterogeneous makeTuplesAsync(DTRecSegment4DCUDA const& dtSegments_d, CSCSegmentCUDA const& cscSegments_d, cudaStream_t stream) const;

public:

  using Params = l2MuonGenerator::Params;
  using Counters = l2MuonGenerator::Counters;

private:
  //void buildDoublets(HitsOnCPU const& hh, cudaStream_t stream) const;

  //void hitNtuplets(HitsOnCPU const& hh, const edm::EventSetup& es, bool useRiemannFit, cudaStream_t cudaStream);

  //void launchKernels(HitsOnCPU const& hh, bool useRiemannFit, cudaStream_t cudaStream) const;

  Params m_params;

  Counters* m_counters = nullptr;
};

#endif  // RecoMuon_L2MuonProducer_src_L2MuonGeneratorOnGPU_h
