#ifndef CUDADataFormats_CSCRecHit_interface_CSCSegmentCUDA_h
#define CUDADataFormats_CSCRecHit_interface_CSCSegmentCUDA_h



#include <cuda_runtime.h>

#include "DataFormats/CSCRecHit/interface/CSCSegmentContainerCUDA.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"


class CSCSegmentCUDA {
public:
  // default constructor, required by cms::cuda::Product<CSCSegmentCUDA>
  CSCSegmentCUDA() = default;
  // constructor that allocates cached device memory on the given CUDA stream
   CSCSegmentCUDA(cudaStream_t stream) { data_d_ = cms::cuda::make_device_unique<CSCSegmentContainerCUDA>(stream); };

  // movable, non-copiable
  CSCSegmentCUDA(CSCSegmentCUDA const&) = delete;
  CSCSegmentCUDA(CSCSegmentCUDA&&) = default;
  CSCSegmentCUDA& operator=(CSCSegmentCUDA const&) = delete;
  CSCSegmentCUDA& operator=(CSCSegmentCUDA&&) = default;


  CSCSegmentContainerCUDA* data() { return data_d_.get(); }
  CSCSegmentContainerCUDA const* data() const { return data_d_.get(); }

  cms::cuda::device::unique_ptr<CSCSegmentContainerCUDA>& ptr() { return data_d_; }
  cms::cuda::device::unique_ptr<CSCSegmentContainerCUDA> const& ptr() const { return data_d_; } 


private:

  cms::cuda::device::unique_ptr<CSCSegmentContainerCUDA> data_d_;

};

#endif  // CUDADataFormats_CSCRecHit_interface_CSCSegmentCUDA_h
