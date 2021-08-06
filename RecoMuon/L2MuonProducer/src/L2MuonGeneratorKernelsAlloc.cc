#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "L2MuonGeneratorKernels.h"

template<>
void L2MuonGeneratorKernelsGPU::allocateOnGPU(int32_t nSegmentsDT, int32_t nSegmentsCSC, cudaStream_t stream) {

  nSegmentsDT++;  // storage requires one more counter;
  nSegmentsCSC++;  // storage requires one more counter;
  assert(nSegmentsDT > 0);
  assert(nSegmentsCSC > 0);


  device_storage_ = Traits::template make_unique<cms::cuda::AtomicPairCounter::c_type[]>(3, stream);

  device_segmentTupleDT_apc_ = (cms::cuda::AtomicPairCounter*)device_storage_.get();
  device_segmentTupleCSC_apc_ = (cms::cuda::AtomicPairCounter*)device_storage_.get();

}
