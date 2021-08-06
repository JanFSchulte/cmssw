#ifndef CUDADataFormats_Muon_interface_MuonSegmentsCUDA_h
#define CUDADataFormats_Muon_interface_MuonSegmentsCUDA_h

#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCompat.h"

#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDAView.h"

#include <cuda_runtime.h>

class MuonSegmentsCUDA {
public:
  MuonSegmentsCUDA() = default;
  explicit MuonSegmentsCUDA(size_t maxSegments, cudaStream_t stream);
  ~MuonSegmentsCUDA() = default;

  MuonSegmentsCUDA(const MuonSegmentsCUDA &) = delete;
  MuonSegmentsCUDA &operator=(const MuonSegmentsCUDA &) = delete;
  MuonSegmentsCUDA(MuonSegmentsCUDA &&) = default;
  MuonSegmentsCUDA &operator=(MuonSegmentsCUDA &&) = default;

  MuonSegmentsCUDAView* view() { return view_d_.get(); }
  MuonSegmentsCUDAView const* view() const { return view_d_.get(); }

  void fillViewAndCopy(cudaStream_t stream);

  void setNSegents(uint32_t nSegments) { nSegments_ = nSegments; }
  void fillOffsets(uint32_t i, uint32_t offset) { offsets_[i] = offset; }
  void fillLocalX(int i, float x) { lx_h_[i] = x;}
  void fillLocalY(int i, float y) { ly_h_[i] = y;}
  void fillLocalDXDZ(int i, float dxdz) { ldxdz_h_[i] = dxdz;}
  void fillLocalDYDZ(int i, float dydz) { ldxdz_h_[i] = dydz;}
  void fillLocalSigmaX(int i, float sigmaX) { lSigmaX_h_[i] = sigmaX;}
  void fillLocalSigmaY(int i, float sigmaY) { lSigmaY_h_[i] = sigmaY;}
  void fillLocalSigmaDXDZ(int i, float sigmaDXDZ) { lSigmaDXDZ_h_[i] = sigmaDXDZ;}
  void fillLocalSigmaDYDZ(int i, float sigmaDYDZ) { lSigmaDYDZ_h_[i] = sigmaDYDZ;}
  void fillGlobalX(int i, float x) { gx_h_[i] = x;}
  void fillGlobalY(int i, float y) { gy_h_[i] = y;}
  void fillLayerID(int i, uint32_t layerID) { layerID_[i] = layerID;}
 
  uint32_t nSegments() const { return nSegments_; }
  uint32_t getOffset(int i) const { return offsets_[i]; }



private:

  //local position of the segments
  cms::cuda::host::unique_ptr<float[]> lx_h_; 
  cms::cuda::host::unique_ptr<float[]> ly_h_;
  //local direction of the segments
  cms::cuda::host::unique_ptr<float[]> ldxdz_h_; 
  cms::cuda::host::unique_ptr<float[]> ldydz_h_;
  //parameter uncertainties
  cms::cuda::host::unique_ptr<float[]> lSigmaX_h_; 
  cms::cuda::host::unique_ptr<float[]> lSigmaY_h_;
  cms::cuda::host::unique_ptr<float[]> lSigmaDXDZ_h_; 
  cms::cuda::host::unique_ptr<float[]> lSigmaDYDZ_h_;

  //global position of the segments
  cms::cuda::host::unique_ptr<float[]> gx_h_; 
  cms::cuda::host::unique_ptr<float[]> gy_h_;
 
  cms::cuda::host::unique_ptr<uint32_t[]> layerID_;

  cms::cuda::host::unique_ptr<uint32_t[]> offsets_;

  cms::cuda::host::unique_ptr<MuonSegmentsCUDAView> view_h_;
  cms::cuda::device::unique_ptr<MuonSegmentsCUDAView> view_d_;

  size_t nSegments_; 
};

#endif  // CUDADataFormats_Muon_interface_MuonSegmentsCUDA_h
