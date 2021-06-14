#ifndef CUDADataFormats_DTRecHit_interface_DTRecSegment4DCUDA_h
#define CUDADataFormats_DTRecHit_interface_DTRecSegment4DCUDA_h



#include <cuda_runtime.h>

#include "DataFormats/DTRecHit/interface/DTRecSegment4DContainerCUDA.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"


class DTRecSegment4DCUDA {
public:
  // default constructor, required by cms::cuda::Product<DTRecSegment4DCUDA>
  DTRecSegment4DCUDA() = default;
  // constructor that allocates cached device memory on the given CUDA stream
   DTRecSegment4DCUDA(cudaStream_t stream) { data_d_ = cms::cuda::make_device_unique<DTRecSegment4DContainerCUDA>(stream); };

  // movable, non-copiable
  DTRecSegment4DCUDA(DTRecSegment4DCUDA const&) = delete;
  DTRecSegment4DCUDA(DTRecSegment4DCUDA&&) = default;
  DTRecSegment4DCUDA& operator=(DTRecSegment4DCUDA const&) = delete;
  DTRecSegment4DCUDA& operator=(DTRecSegment4DCUDA&&) = default;


  DTRecSegment4DContainerCUDA* data() { return data_d_.get(); }
  DTRecSegment4DContainerCUDA const* data() const { return data_d_.get(); }

  cms::cuda::device::unique_ptr<DTRecSegment4DContainerCUDA>& ptr() { return data_d_; }
  cms::cuda::device::unique_ptr<DTRecSegment4DContainerCUDA> const& ptr() const { return data_d_; } 


private:

  cms::cuda::device::unique_ptr<DTRecSegment4DContainerCUDA> data_d_;

};

#endif  // CUDADataFormats_DTRecHit_interface_DTRecSegment4DCUDA_h
