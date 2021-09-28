#ifndef RecoMuon_L2MuonProducer_plugins_gpuMuonDoubletsAlgos_h
#define RecoMuon_L2MuonProducer_plugins_gpuMuonDoubletsAlgos_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>

#include "CUDADataFormats/TrackingRecHit/interface/TrackingRecHit2DHeterogeneous.h"
#include "DataFormats/Math/interface/approx_atan2.h"
#include "HeterogeneousCore/CUDAUtilities/interface/VecArray.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cuda_assert.h"

#include "CAConstantsMuon.h"
#include "GPUCACellMuon.h"

namespace gpuMuonDoublets {

  using CellNeighbors = caConstants::CellNeighbors;
  using CellTracks = caConstants::CellTracks;
  using CellNeighborsVector = caConstants::CellNeighborsVector;
  using CellTracksVector = caConstants::CellTracksVector;

  __device__ __forceinline__ void doubletsFromHisto(uint8_t const* __restrict__ layerPairs,
                                                    uint32_t nPairs,
                                                    GPUCACellMuon* cells,
                                                    uint32_t* nCells,
                                                    CellNeighborsVector* cellNeighbors,
                                                    CellTracksVector* cellTracks,
                                                    MuonSegmentsCUDAView const& __restrict__ hh,
                                                    GPUCACellMuon::OuterHitOfCell* isOuterHitOfCell,
                                                    float const* __restrict__ phicuts,
                                                    float const* __restrict__ minz,
                                                    float const* __restrict__ maxz,
                                                    float const* __restrict__ maxr,
                                                    float const* __restrict__ minr,
                                                    float const* __restrict__ mindz,
                                                    float const* __restrict__ maxdz,
                                                    float const* __restrict__ maxdist,
                                                    bool doZ0Cut,
                                                    bool doPtCut,
                                                    uint32_t maxNumOfDoublets) {

    uint32_t const* __restrict__ offsets = hh.offsets();
    assert(offsets);
    auto layerSize = [=](uint8_t li) { return offsets[li+1] - offsets[li]; };

    // nPairsMax to be optimized later (originally was 64).
    // If it should be much bigger, consider using a block-wide parallel prefix scan,
    // e.g. see  https://nvlabs.github.io/cub/classcub_1_1_warp_scan.html
    const int nPairsMax = caConstants::maxNumberOfLayerPairs;
    assert(nPairs <= nPairsMax);
    __shared__ uint32_t innerLayerCumulativeSize[nPairsMax];
    __shared__ uint32_t ntot;
    if (threadIdx.y == 0 && threadIdx.x == 0) {
      innerLayerCumulativeSize[0] = layerSize(layerPairs[0]);
      for (uint32_t i = 1; i < nPairs; ++i) {
        innerLayerCumulativeSize[i] = innerLayerCumulativeSize[i - 1] + layerSize(layerPairs[2 * i]);
      }
      ntot = innerLayerCumulativeSize[nPairs - 1];
    }
    __syncthreads();

    // x runs faster
    auto idy = blockIdx.y * blockDim.y + threadIdx.y;
    auto first = threadIdx.x;
    auto stride = blockDim.x;
    //printf("going in\n");
    uint32_t pairLayerId = 0;  // cannot go backward
    for (auto j = idy; j < ntot; j += blockDim.y * gridDim.y) {
      while (j >= innerLayerCumulativeSize[pairLayerId++])
        ;
      --pairLayerId;  // move to lower_bound ??
      assert(pairLayerId < nPairs);
      assert(j < innerLayerCumulativeSize[pairLayerId]);
      assert(0 == pairLayerId || j >= innerLayerCumulativeSize[pairLayerId - 1]);

      uint8_t inner = layerPairs[2 * pairLayerId];
      uint8_t outer = layerPairs[2 * pairLayerId + 1];
      assert(outer > inner);
      auto i = (0 == pairLayerId) ? j : j - innerLayerCumulativeSize[pairLayerId - 1];
      i += offsets[inner];
      //if (!(inner == 2 && outer ==4)) continue;
      //printf("inner: %d outer %d\n",inner, outer);
      assert(i >= offsets[inner]);
      assert(i < offsets[inner+1]);

      // found hit corresponding to our cuda thread, now do the job
      auto mez = hh.gz(i);
      //if (mez < minz[pairLayerId] || mez > maxz[pairLayerId])
      //  continue;
      auto mep = hh.phi(i);
      auto mer = hh.gr(i);
      // all cuts: true if fails
      auto z0cut = maxdist[pairLayerId];      // cm
      auto z0cutoff = [&](int j) {
        auto zo = hh.gz(j);
        auto ro = hh.gr(j);
        auto dr = ro - mer;
        return (dr > maxr[pairLayerId] || dr < minr[pairLayerId] || (std::abs((mez * ro - mer * zo))/dr) > z0cut);
      };
      auto dzcutoff = [&](int j) {
        auto zo = hh.gz(j);
        auto dz = zo - mez;
        return (dz < mindz[pairLayerId] || dz > maxdz[pairLayerId]);
      };

      auto iphicut = phicuts[pairLayerId];

#ifdef GPU_DEBUG
      int tot = 0;
      int nmin = 0;
      int tooMany = 0;
#endif
      uint32_t p = offsets[outer];
      uint32_t e = offsets[outer+1];
      //printf("%d %d\n", p, e);
      p += first;
      //printf("before the loop\n");
      for (; p < e; p += stride) {
        auto oi = p;
        assert(oi >= offsets[outer]);
        assert(oi < offsets[outer+1]);
        if (doZ0Cut && (z0cutoff(oi) || dzcutoff(oi)))
          continue;

        //printf("after cut1\n");
        auto mop = hh.phi(oi);
        float dphi = std::min(std::abs(mop - mep), std::abs(mep - mop));
        if (dphi > iphicut)
          continue;
        //printf("after cut 2\n");
        auto ind = atomicAdd(nCells, 1);
        if (ind >= maxNumOfDoublets) {
          atomicSub(nCells, 1);
          break;
        }  // move to SimpleVector??
        // int layerPairId, int doubletId, int innerHitId, int outerHitId)
        cells[ind].init(*cellNeighbors, *cellTracks, hh, pairLayerId, ind, i, oi);
        isOuterHitOfCell[oi].push_back(ind);
#ifdef GPU_DEBUG
        if (isOuterHitOfCell[oi].full())
          ++tooMany;
        ++tot;
#endif
      }
#ifdef GPU_DEBUG
      if (tooMany > 0)
        printf("OuterHitOfCell full for %d in layer %d/%d, %d,%d %d\n", i, inner, outer, nmin, tot, tooMany);
#endif
    }  // loop in block...
  }

}  // namespace gpuMuonDoublets

#endif  // RecoMuon_L2MuonProducer_plugins_gpuMuonDoubletsAlgos_h
