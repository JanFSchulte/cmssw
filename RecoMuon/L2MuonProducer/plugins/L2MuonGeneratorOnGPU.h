#ifndef RecoMuon_L2MuonProducer_plugins_L2MuonGeneratorOnGPU_h
#define RecoMuon_L2MuonProducer_plugins_L2MuonGeneratorOnGPU_h

#include <cuda_runtime.h>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"


#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentPairsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentNtupletsCUDA.h"
#include "CUDADataFormats/Track/interface/L2MuonTrackHeterogeneous.h"

#include "CUDADataFormats/Muon/interface/MuonSegmentPairsHeterogeneous.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentNtupletsHeterogeneous.h"
#include "RecoMuon/L2MuonProducer/plugins/L2MuonGeneratorKernels.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/HelixFitOnGPU.h"

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

  L2MuonTrackHeterogeneous makeTuplesAsync(MuonSegmentsCUDA const& muonSegments_h, float bfield, cudaStream_t stream) const;
  MuonSegmentNtupletsHeterogeneous makeTuplesAsyncForReturn(MuonSegmentsCUDA const& muonSegments_h, float bfield, cudaStream_t stream) const;
  MuonSegmentPairsHeterogeneous makeDoubletsAsync(MuonSegmentsCUDA const& muonSegments_h, float bfield, cudaStream_t stream) const;

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

#endif  // RecoMuon_L2MuonProducer_plugins_L2MuonGeneratorOnGPU_h
