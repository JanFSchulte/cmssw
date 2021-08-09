#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "L2MuonGeneratorKernels.h"

template<>
void L2MuonGeneratorKernelsGPU::allocateOnGPU(int32_t nSegments, cudaStream_t stream) {

  nSegments++;  // storage requires one more counter;
  assert(nSegments > 0);


  device_storage_ = Traits::template make_unique<cms::cuda::AtomicPairCounter::c_type[]>(3, stream);

  device_segmentTuple_apc_ = (cms::cuda::AtomicPairCounter*)device_storage_.get();

}
