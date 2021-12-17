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

  void fillLocalX(uint32_t i, float x) { lx_h_[i] = x;}
  void fillLocalY(uint32_t i, float y) { ly_h_[i] = y;}
  void fillLocalDXDZ(uint32_t i, float dxdz) { ldxdz_h_[i] = dxdz;}
  void fillLocalDYDZ(uint32_t i, float dydz) { ldxdz_h_[i] = dydz;}
  void fillLocalSigmaX(uint32_t i, float sigmaX) { lSigmaX_h_[i] = sigmaX;}
  void fillLocalSigmaY(uint32_t i, float sigmaY) { lSigmaY_h_[i] = sigmaY;}
  void fillLocalSigmaDXDZ(uint32_t i, float sigmaDXDZ) { lSigmaDXDZ_h_[i] = sigmaDXDZ;}
  void fillLocalSigmaDYDZ(uint32_t i, float sigmaDYDZ) { lSigmaDYDZ_h_[i] = sigmaDYDZ;}
  void fillGlobalX(uint32_t i, float x) { gx_h_[i] = x;}
  void fillGlobalY(uint32_t i, float y) { gy_h_[i] = y;}
  void fillGlobalZ(uint32_t i, float z) { gz_h_[i] = z;}
  void fillGlobalR(uint32_t i, float r) { gr_h_[i] = r;}
  void fillGlobalDX(uint32_t i, float dx) { gdx_h_[i] = dx;}
  void fillGlobalDY(uint32_t i, float dy) { gdy_h_[i] = dy;}  
  void fillGlobalDZ(uint32_t i, float dz) { gdz_h_[i] = dz;}  
  void fillPhi(uint32_t i, float phi) { phi_h_[i] = phi;}
  void fillLayerID(uint32_t i, uint32_t layerID) { layerID_[i] = layerID;}
 
  uint32_t nSegments() const { return nSegments_; }
  uint32_t getOffset(uint32_t i) const { return offsets_[i]; }

  float localX(uint32_t i) const { return lx_h_[i];}
  float localY(uint32_t i) const { return ly_h_[i];}
  float localDXDZ(uint32_t i) const { return ldxdz_h_[i];}
  float localDYDZ(uint32_t i) const { return ldydz_h_[i];}
  float localSigmaX(uint32_t i) const { return lSigmaX_h_[i];}
  float localSigmaY(uint32_t i) const { return lSigmaY_h_[i];}
  float localSigmaDXDZ(uint32_t i) const { return lSigmaDXDZ_h_[i];}
  float localSigmaDYDZ(uint32_t i) const { return lSigmaDYDZ_h_[i];}  
  float globalX(uint32_t i) const { return gx_h_[i];}
  float globalY(uint32_t i) const { return gy_h_[i];}
  float globalZ(uint32_t i) const { return gz_h_[i];}
  float globalR(uint32_t i) const { return gr_h_[i];}
  float globalDX(uint32_t i) const { return gdx_h_[i];}
  float globalDY(uint32_t i) const { return gdy_h_[i];}
  float globalDZ(uint32_t i) const { return gdz_h_[i];}
  float phi(uint32_t i) const { return phi_h_[i];}
  uint32_t layerID(uint32_t i) const { return layerID_[i];} 

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
  cms::cuda::host::unique_ptr<float[]> gz_h_;
  cms::cuda::host::unique_ptr<float[]> gr_h_;
  //global direction of the segments
  cms::cuda::host::unique_ptr<float[]> gdx_h_; 
  cms::cuda::host::unique_ptr<float[]> gdy_h_;
  cms::cuda::host::unique_ptr<float[]> gdz_h_;

  cms::cuda::host::unique_ptr<float[]> phi_h_;
 
  cms::cuda::host::unique_ptr<uint32_t[]> layerID_;

  cms::cuda::host::unique_ptr<uint32_t[]> offsets_;

  cms::cuda::host::unique_ptr<MuonSegmentsCUDAView> view_h_;
  cms::cuda::device::unique_ptr<MuonSegmentsCUDAView> view_d_;

  size_t nSegments_; 
};

#endif  // CUDADataFormats_Muon_interface_MuonSegmentsCUDA_h
