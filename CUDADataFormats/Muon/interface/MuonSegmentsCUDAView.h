#ifndef CUDADataFormats_Muon_interface_MuonSegmentsCUDAView_h
#define CUDADataFormats_Muon_interface_MuonSegmentsCUDAView_h

#include <cuda_runtime.h>

class MuonSegmentsCUDAView {
  public:

    friend class MuonSegmentsCUDA;

    __device__ __forceinline__ float lx(int i) const { return __ldg(lx_d_ + i); }
    __device__ __forceinline__ float ly(int i) const { return __ldg(ly_d_ + i); }
    __device__ __forceinline__ float ldxdz(int i) const { return __ldg(ldxdz_d_ + i); }
    __device__ __forceinline__ float ldydz(int i) const { return __ldg(ldydz_d_ + i); }
    __device__ __forceinline__ float lSigmaX(int i) const { return __ldg(lSigmaX_d_ + i); }
    __device__ __forceinline__ float lSigmaY(int i) const { return __ldg(lSigmaY_d_ + i); }
    __device__ __forceinline__ float lSigmaDXDZ(int i) const { return __ldg(lSigmaDXDZ_d_ + i); }
    __device__ __forceinline__ float lSigmaDYDZ(int i) const { return __ldg(lSigmaDYDZ_d_ + i); }

    __device__ __forceinline__ float gx(int i) const { return __ldg(gx_d_ + i); }
    __device__ __forceinline__ float gy(int i) const { return __ldg(gy_d_ + i); }
 
    __device__ __forceinline__ uint32_t layerID(int i) const { return __ldg(layerID_d_ + i); }

    __device__ __forceinline__ int nSemgents() const { return nSegments_d_; }


    float *lx_d_;
    float *ly_d_;

    float *ldxdz_d_;
    float *ldydz_d_;

    float *lSigmaX_d_;
    float *lSigmaY_d_;
    float *lSigmaDXDZ_d_;
    float *lSigmaDYDZ_d_;

    float *gx_d_;
    float *gy_d_;

    uint32_t *layerID_d_;

    int nSegments_d_;
};
 

#endif  // CUDADataFormats_Muon_interface_MuonSegmentsCUDAView_h
