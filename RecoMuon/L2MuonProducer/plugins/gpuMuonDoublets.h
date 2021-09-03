#ifndef RecoMuon_L2MuonProducer_plugins_gpuMuonDoublets_h
#define RecoMuon_L2MuonProducer_plugins_gpuMuonDoublets_h

#include "RecoMuon/L2MuonProducer/plugins/gpuMuonDoubletsAlgos.h"

#define CONSTANT_VAR __constant__

namespace gpuMuonDoublets {

  constexpr int nPairsForQuadruplets = 15;                     // quadruplets require segments in all layers
  constexpr int nPairsForTriplets = nPairsForQuadruplets + 2;  // include barrel "jumping" layer pairs
  constexpr int nPairs = nPairsForTriplets + 6;                // include forward "jumping" layer pairs
  static_assert(nPairs <= caConstants::maxNumberOfLayerPairs);

  // start constants
  // clang-format off

  CONSTANT_VAR const uint8_t layerPairs[2 * nPairs] = {
      0, 1, 0, 4, 0, 8,              // MB1 (3)
      1, 2, 1, 4, 1, 8,              // MB2 (6)
      4, 5, 8, 9,                    // ME1 (8)
      2, 3, 2, 4, 2, 8, 5, 6, 9, 10, // MB3 & ME2 (13)
      6, 7, 10, 11,                  // ME3 (15)       
      0, 2, 1, 3,                    // Jumping Barrel (17)
      0, 5, 0, 9,                    // Jumping Forward (MB1,ME2) (19)
      4, 6, 5, 7, 8, 10, 9, 11       // Jumping Forward (23)
  };


  CONSTANT_VAR const float phicuts[nPairs]{0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.7,
                                             0.7,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,
                                             0.5,};
  //   phi0p07, phi0p07, phi0p06,phi0p06, phi0p06,phi0p06};  // relaxed cuts

  CONSTANT_VAR float const minz[nPairs] = {
      -750., 0.,   -750., -750., 0.,  -750., 500., -800., -750., -750., 0, 800., -900., 900., -1000., -750, -750., 0., -750.,500., 750., -800.,-800};
  CONSTANT_VAR float const maxz[nPairs] = {
       750., 750., 0.,     750., 750., 0.,   800., -500.,  750., 0.,  750., 900.,-800.,  1000., -900., 750., 750., 750., 0., 800., 800., -500.,-750 };
  CONSTANT_VAR float const maxr[nPairs] = {
      150,250.,300.,150.,200., 200., 200., 200., 200., 100 , 100., 150., 150., 150., 150., 300., 300., 250., 250., 200., 200., 200., 200.};
  CONSTANT_VAR float const minr[nPairs] = {
      0.,-200.,-200.,0.,0., 0., -100., -100., 0. , -100., -100., -100., -150., -150., -300., 0., 0., 0., 0., 0., -100., 0.,-100.};


  // end constants
  // clang-format on

  using CellNeighbors = caConstants::CellNeighbors;
  using CellTracks = caConstants::CellTracks;
  using CellNeighborsVector = caConstants::CellNeighborsVector;
  using CellTracksVector = caConstants::CellTracksVector;

  __global__ void initDoublets(GPUCACellMuon::OuterHitOfCell* isOuterHitOfCell,
                               int nHits,
                               CellNeighborsVector* cellNeighbors,
                               CellNeighbors* cellNeighborsContainer,
                               CellTracksVector* cellTracks,
                               CellTracks* cellTracksContainer) {
    assert(isOuterHitOfCell);
    int first = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = first; i < nHits; i += gridDim.x * blockDim.x)
      isOuterHitOfCell[i].reset();

    if (0 == first) {
      cellNeighbors->construct(caConstants::maxNumOfActiveDoublets, cellNeighborsContainer);
      cellTracks->construct(caConstants::maxNumOfActiveDoublets, cellTracksContainer);
      auto i = cellNeighbors->extend();
      assert(0 == i);
      (*cellNeighbors)[0].reset();
      i = cellTracks->extend();
      assert(0 == i);
      (*cellTracks)[0].reset();
    }
  }

  constexpr auto getDoubletsFromHistoMaxBlockSize = 64;  // for both x and y
  constexpr auto getDoubletsFromHistoMinBlocksPerMP = 16;

  __global__
#ifdef __CUDACC__
  __launch_bounds__(getDoubletsFromHistoMaxBlockSize, getDoubletsFromHistoMinBlocksPerMP)
#endif
      void getDoubletsFromHisto(GPUCACellMuon* cells,
                                uint32_t* nCells,
                                CellNeighborsVector* cellNeighbors,
                                CellTracksVector* cellTracks,
                                MuonSegmentsCUDAView const* __restrict__ hhp,
                                GPUCACellMuon::OuterHitOfCell* isOuterHitOfCell,
                                int nActualPairs,
                                bool doZ0Cut,
                                bool doPtCut,
                                uint32_t maxNumOfDoublets) {
    auto const& __restrict__ hh = *hhp;
    doubletsFromHisto(layerPairs,
                      nActualPairs,
                      cells,
                      nCells,
                      cellNeighbors,
                      cellTracks,
                      hh,
                      isOuterHitOfCell,
                      phicuts,
                      minz,
                      maxz,
                      maxr,
                      minr,
                      doZ0Cut,
                      doPtCut,
                      maxNumOfDoublets);
  }

}  // namespace gpuMuonDoublets

#endif  // RecoMuon_L2MuonProducer_plugins_gpuMuonDoublets_h
