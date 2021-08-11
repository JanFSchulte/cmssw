#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDA.h"
#include "HeterogeneousCore/CUDAUtilities/interface/copyAsync.h"


MuonSegmentsCUDA::MuonSegmentsCUDA(size_t maxSegments, cudaStream_t stream) {

  lx_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  ly_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  ldxdz_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  ldydz_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  lSigmaX_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  lSigmaY_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  lSigmaDXDZ_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  lSigmaDYDZ_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  gx_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  gy_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  gz_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  gr_h_ = cms::cuda::make_host_unique<float[]>(maxSegments, stream);;
  layerID_ = cms::cuda::make_host_unique<uint32_t[]>(maxSegments, stream);;
  offsets_ = cms::cuda::make_host_unique<uint32_t[]>(12, stream);;

  view_h_ = cms::cuda::make_host_unique<MuonSegmentsCUDAView>(stream);;


}

void MuonSegmentsCUDA::fillViewAndCopy(cudaStream_t stream){

  view_h_->lx_d_ = lx_h_.get();
  view_h_->ly_d_ = ly_h_.get();
  view_h_->ldxdz_d_ = ldxdz_h_.get();
  view_h_->ldydz_d_ = ldydz_h_.get();
  view_h_->lSigmaX_d_ = lSigmaX_h_.get();
  view_h_->lSigmaY_d_ = lSigmaY_h_.get();
  view_h_->lSigmaDXDZ_d_ = lSigmaDXDZ_h_.get();
  view_h_->lSigmaDYDZ_d_ = lSigmaDYDZ_h_.get();
  view_h_->gx_d_ = gx_h_.get();
  view_h_->gy_d_ = gy_h_.get();
  view_h_->gz_d_ = gz_h_.get();
  view_h_->gr_d_ = gr_h_.get();
  view_h_->layerID_d_ = layerID_.get();

  view_h_->offsets_d_ = offsets_.get();

  view_d_ = cms::cuda::make_device_unique<MuonSegmentsCUDAView>(stream);
  cms::cuda::copyAsync(view_d_, view_h_, stream);
}
