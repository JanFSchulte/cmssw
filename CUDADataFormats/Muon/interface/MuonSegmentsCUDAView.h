#ifndef CUDADataFormats_Muon_interface_MuonSegmentsCUDAView_h
#define CUDADataFormats_Muon_interface_MuonSegmentsCUDAView_h

#include <cuda_runtime.h>

class MuonSegmentsCUDAView {
  public:

    friend class MuonSegmentsCUDA;

    using hindex_type = uint32_t;

    __device__ __forceinline__ float lx(uint32_t i) const { return __ldg(lx_d_ + i); }
    __device__ __forceinline__ float ly(uint32_t i) const { return __ldg(ly_d_ + i); }
    __device__ __forceinline__ float ldxdz(uint32_t i) const { return __ldg(ldxdz_d_ + i); }
    __device__ __forceinline__ float ldydz(uint32_t i) const { return __ldg(ldydz_d_ + i); }
    __device__ __forceinline__ float lSigmaX(uint32_t i) const { return __ldg(lSigmaX_d_ + i); }
    __device__ __forceinline__ float lSigmaY(uint32_t i) const { return __ldg(lSigmaY_d_ + i); }
    __device__ __forceinline__ float lSigmaDXDZ(uint32_t i) const { return __ldg(lSigmaDXDZ_d_ + i); }
    __device__ __forceinline__ float lSigmaDYDZ(uint32_t i) const { return __ldg(lSigmaDYDZ_d_ + i); }

    __device__ __forceinline__ float gx(uint32_t i) const { return __ldg(gx_d_ + i); }
    __device__ __forceinline__ float gy(uint32_t i) const { return __ldg(gy_d_ + i); }
    __device__ __forceinline__ float gz(uint32_t i) const { return __ldg(gz_d_ + i); }
    __device__ __forceinline__ float gr(uint32_t i) const { return __ldg(gr_d_ + i); }
    __device__ __forceinline__ float gdx(uint32_t i) const { return __ldg(gdx_d_ + i); }
    __device__ __forceinline__ float gdy(uint32_t i) const { return __ldg(gdy_d_ + i); }
    __device__ __forceinline__ float gdz(uint32_t i) const { return __ldg(gdz_d_ + i); }

    __device__ __forceinline__ float phi(uint32_t i) const { return __ldg(phi_d_ + i); }
 
    __device__ __forceinline__ uint32_t layerID(uint32_t i) const { return __ldg(layerID_d_ + i); }
    __device__ __forceinline__ uint32_t offset(uint32_t i) const { return offsets_d_[i]; }
    __device__ __forceinline__ uint32_t* offsets() const { return offsets_d_; }

    __device__ __forceinline__ uint32_t nSegments() const { return nSegments_d_; }


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
    float *gz_d_;
    float *gr_d_;
    float *gdx_d_;
    float *gdy_d_;
    float *gdz_d_;

    float *phi_d_;

    uint32_t *layerID_d_;
    uint32_t *offsets_d_;

    uint32_t nSegments_d_;
};
 

#endif  // CUDADataFormats_Muon_interface_MuonSegmentsCUDAView_h
