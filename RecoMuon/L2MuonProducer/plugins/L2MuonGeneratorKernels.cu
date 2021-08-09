// CUDA runtime
#include <cuda_runtime.h>

// CMSSW headers
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "RecoMuon/L2MuonProducer/plugins/L2MuonGeneratorKernels.h"
#include "RecoMuon/L2MuonProducer/plugins/L2MuonGeneratorKernelsImpl.h"

template <>
void L2MuonGeneratorKernelsGPU::buildL2Muons(MuonSegmentsCUDA const& muonSegments_h,
						          L2MuonTrack::TrackSoA* l2Muons_d,
                                                          cudaStream_t stream) const {
    auto nSegments = muonSegments_h.nSegments();

    int threadsPerBlock = 128;
    int blocks = 60;  // number of sectors in DT (???)

#ifdef GPU_DEBUG
    std::cout << "launching createL2Muons kernel for " << blocks << " blocks" << std::endl;
#endif
    // protect from empty events
    if (blocks) {
      makeL2Muon<<<blocks, threadsPerBlock, 0, stream>>>(
          l2Muons_d, muonSegments_h.view());
      cudaCheck(cudaGetLastError());
#ifdef GPU_DEBUG
      cudaCheck(cudaDeviceSynchronize());
#endif
    }


  }
