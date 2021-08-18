// CUDA runtime
#include <cuda_runtime.h>

// CMSSW headers
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "RecoMuon/L2MuonProducer/plugins/L2MuonGeneratorKernels.h"
#include "RecoMuon/L2MuonProducer/plugins/L2MuonGeneratorKernelsImpl.h"
#include "RecoMuon/L2MuonProducer/plugins/gpuMuonDoublets.h"

template <>
void L2MuonGeneratorKernelsGPU::buildL2Muons(MuonSegmentsCUDA const& muonSegments_h,
						          L2MuonTrack::TrackSoA* l2Muons_d,
                                                          cudaStream_t stream) const {

  // these are pointer on GPU!
  auto *tuples_d = &l2Muons_d->hitIndices;
  auto *quality_d = l2Muons_d->qualityData();

  // zero tuples
  cms::cuda::launchZero(tuples_d, stream);

  auto nSegments = muonSegments_h.nSegments();


  auto nthTot = 64;
  auto stride = 4;
  auto blockSize = nthTot / stride;
  auto numberOfBlocks = nDoubletBlocks(blockSize);
  auto rescale = numberOfBlocks / 65536;
  blockSize *= (rescale + 1);
  numberOfBlocks = nDoubletBlocks(blockSize);
  assert(numberOfBlocks < 65536);
  assert(blockSize > 0 && 0 == blockSize % 16);
  dim3 blks(1, numberOfBlocks, 1);
  dim3 thrs(stride, blockSize, 1);

  kernel_connect<<<blks, thrs, 0, stream>>>(
      device_hitTuple_apc_,
      device_hitToTuple_apc_,  // needed only to be reset, ready for next kernel
      muonSegments_h.view(),
      device_theCells_.get(),
      device_nCells_,
      device_theCellNeighbors_.get(),
      device_isOuterHitOfCell_.get(),
      params_.hardCurvCut_,
      params_.ptmin_,
      params_.CAThetaCutBarrel_,
      params_.CAThetaCutForward_,
      params_.dcaCutInnerTriplet_,
      params_.dcaCutOuterTriplet_);
  cudaCheck(cudaGetLastError());

  blockSize = 64;
  numberOfBlocks = (3 * params_.maxNumberOfDoublets_ / 4 + blockSize - 1) / blockSize;
  kernel_find_ntuplets<<<numberOfBlocks, blockSize, 0, stream>>>(muonSegments_h.view(),
                                                                     device_theCells_.get(),
                                                                     device_nCells_,
                                                                     device_theCellTracks_.get(),
                                                                     tuples_d,
                                                                     device_hitTuple_apc_,
                                                                     quality_d,
                                                                     params_.minHitsPerNtuplet_);
  cudaCheck(cudaGetLastError());


  if (params_.doStats_)
    kernel_mark_used<<<numberOfBlocks, blockSize, 0, stream>>>(muonSegments_h.view(), device_theCells_.get(), device_nCells_);
  cudaCheck(cudaGetLastError());

#ifdef GPU_DEBUG
  cudaDeviceSynchronize();
  cudaCheck(cudaGetLastError());
#endif

  blockSize = 128;
  numberOfBlocks = (HitContainer::ctNOnes() + blockSize - 1) / blockSize;
  cms::cuda::finalizeBulk<<<numberOfBlocks, blockSize, 0, stream>>>(device_hitTuple_apc_, tuples_d);
  // remove duplicates (tracks that share a doublet)
  numberOfBlocks = nDoubletBlocks(blockSize);
  kernel_earlyDuplicateRemover<<<numberOfBlocks, blockSize, 0, stream>>>(
      device_theCells_.get(), device_nCells_, tuples_d, quality_d, params_.dupPassThrough_);
  cudaCheck(cudaGetLastError());

  blockSize = 128;
  numberOfBlocks = (3 * caConstants::maxTuples / 4 + blockSize - 1) / blockSize;
  kernel_countMultiplicity<<<numberOfBlocks, blockSize, 0, stream>>>(
      tuples_d, quality_d, device_tupleMultiplicity_.get());
  cms::cuda::launchFinalize(device_tupleMultiplicity_.get(), stream);
  kernel_fillMultiplicity<<<numberOfBlocks, blockSize, 0, stream>>>(
      tuples_d, quality_d, device_tupleMultiplicity_.get());
  cudaCheck(cudaGetLastError());

#ifdef GPU_DEBUG
  cudaDeviceSynchronize();
  cudaCheck(cudaGetLastError());
#endif

}



template <>
void L2MuonGeneratorKernelsGPU::buildDoublets(MuonSegmentsCUDA const &muonSegments_h, cudaStream_t stream) {
  int32_t nSegments = muonSegments_h.nSegments();


#ifdef GPU_DEBUG
  cudaDeviceSynchronize();
  cudaCheck(cudaGetLastError());
#endif

  device_isOuterHitOfCell_ = cms::cuda::make_device_unique<GPUCACellMuon::OuterHitOfCell[]>(std::max(1, nSegments), stream);
  assert(device_isOuterHitOfCell_.get());

  cellStorage_ = cms::cuda::make_device_unique<unsigned char[]>(
      caConstants::maxNumOfActiveDoublets * sizeof(GPUCACellMuon::CellNeighbors) +
          caConstants::maxNumOfActiveDoublets * sizeof(GPUCACellMuon::CellTracks),
      stream);
  device_theCellNeighborsContainer_ = (GPUCACellMuon::CellNeighbors *)cellStorage_.get();
  device_theCellTracksContainer_ = (GPUCACellMuon::CellTracks *)(cellStorage_.get() + caConstants::maxNumOfActiveDoublets *
                                                                                      sizeof(GPUCACellMuon::CellNeighbors));



  {
    int threadsPerBlock = 128;
    // at least one block!
    int blocks = (std::max(1, nSegments) + threadsPerBlock - 1) / threadsPerBlock;
    gpuMuonDoublets::initDoublets<<<blocks, threadsPerBlock, 0, stream>>>(device_isOuterHitOfCell_.get(),
                                                                           nSegments,
                                                                           device_theCellNeighbors_.get(),
                                                                           device_theCellNeighborsContainer_,
                                                                           device_theCellTracks_.get(),
                                                                           device_theCellTracksContainer_);
    cudaCheck(cudaGetLastError());
  }
  device_theCells_ = cms::cuda::make_device_unique<GPUCACellMuon[]>(params_.maxNumberOfDoublets_, stream);

#ifdef GPU_DEBUG
  cudaDeviceSynchronize();
  cudaCheck(cudaGetLastError());
#endif

  if (0 == nSegments)
    return;  // protect against empty events

  // take all layer pairs into account
  auto nActualPairs = gpuMuonDoublets::nPairs;
  if (not params_.includeJumpingForwardDoublets_) {
    // exclude forward "jumping" layer pairs
    nActualPairs = gpuMuonDoublets::nPairsForTriplets;
  }
  if (params_.minHitsPerNtuplet_ > 3) {
    // for quadruplets, exclude all "jumping" layer pairs
    nActualPairs = gpuMuonDoublets::nPairsForQuadruplets;
  }

  assert(nActualPairs <= gpuMuonDoublets::nPairs);
  int stride = 4;
  int threadsPerBlock = gpuMuonDoublets::getDoubletsFromHistoMaxBlockSize / stride;
  int blocks = (4 * nSegments + threadsPerBlock - 1) / threadsPerBlock;
  dim3 blks(1, blocks, 1);
  dim3 thrs(stride, threadsPerBlock, 1);
  gpuMuonDoublets::getDoubletsFromHisto<<<blks, thrs, 0, stream>>>(device_theCells_.get(),
                                                                    device_nCells_,
                                                                    device_theCellNeighbors_.get(),
                                                                    device_theCellTracks_.get(),
                                                                    muonSegments_h.view(),
                                                                    device_isOuterHitOfCell_.get(),
                                                                    nActualPairs,
                                                                    params_.doZ0Cut_,
                                                                    params_.doPtCut_,
                                                                    params_.maxNumberOfDoublets_);
  cudaCheck(cudaGetLastError());

#ifdef GPU_DEBUG
  cudaDeviceSynchronize();
  cudaCheck(cudaGetLastError());
#endif
}

template <>
void L2MuonGeneratorKernelsGPU::fillHitDetIndices(MuonSegmentsCUDAView const *hv, L2MuonTrack::TrackSoA *tracks_d, cudaStream_t cudaStream) {
  auto blockSize = 128;
  auto numberOfBlocks = (HitContainer::ctCapacity() + blockSize - 1) / blockSize;

  kernel_fillHitDetIndices<<<numberOfBlocks, blockSize, 0, cudaStream>>>(
      &tracks_d->hitIndices, hv, &tracks_d->detIndices);
  cudaCheck(cudaGetLastError());
#ifdef GPU_DEBUG
  cudaDeviceSynchronize();
  cudaCheck(cudaGetLastError());
#endif
}

