#include <cmath>
#include <cstdint>
#include <limits>

#include <cuda_runtime.h>

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cuda_assert.h"

#include "L2MuonGeneratorKernels.h"
#include "CAConstantsMuon.h"
#include "GPUCACellMuon.h"
#include "gpuMuonDoublets.h"


using HitsOnGPU = MuonSegmentsCUDAView;
using HitsOnCPU = MuonSegmentsCUDA;

using HitToTuple = caConstants::HitToTuple;
using TupleMultiplicity = caConstants::TupleMultiplicity;

using Quality = L2MuonTrack::Quality;
using TkSoA = L2MuonTrack::TrackSoA;
using HitContainer = L2MuonTrack::HitContainer;

__global__ void kernel_connect(cms::cuda::AtomicPairCounter *apc1,
                               cms::cuda::AtomicPairCounter *apc2,  // just to zero them,
                               GPUCACellMuon::Hits const *__restrict__ hhp,
                               GPUCACellMuon *cells,
                               uint32_t const *__restrict__ nCells,
                               gpuMuonDoublets::CellNeighborsVector *cellNeighbors,
                               GPUCACellMuon::OuterHitOfCell const *__restrict__ isOuterHitOfCell,
                               float hardCurvCut,
                               float ptmin,
                               float CAThetaCutBarrel,
                               float CAThetaCutForward,
                               float dcaCutInnerTriplet,
                               float dcaCutOuterTriplet) {
  auto const &hh = *hhp;

  auto firstCellIndex = threadIdx.y + blockIdx.y * blockDim.y;
  auto first = threadIdx.x;
  auto stride = blockDim.x;

  if (0 == (firstCellIndex + first)) {
    (*apc1) = 0;
    (*apc2) = 0;
  }  // ready for next kernel

  for (int idx = firstCellIndex, nt = (*nCells); idx < nt; idx += gridDim.y * blockDim.y) {
    auto cellIndex = idx;
    auto &thisCell = cells[idx];
    auto innerHitId = thisCell.inner_hit_id();
    int numberOfPossibleNeighbors = isOuterHitOfCell[innerHitId].size();
    auto vi = isOuterHitOfCell[innerHitId].data();

    auto ri = thisCell.inner_r(hh);
    auto zi = thisCell.inner_z(hh);

    auto ro = thisCell.outer_r(hh);
    auto zo = thisCell.outer_z(hh);
    auto isBarrel = thisCell.inner_detIndex(hh) < caConstants::last_barrel_detIndex;

    for (int j = first; j < numberOfPossibleNeighbors; j += stride) {
      auto otherCell = __ldg(vi + j);
      auto &oc = cells[otherCell];
      auto r1 = oc.inner_r(hh);
      auto z1 = oc.inner_z(hh);
      bool aligned = GPUCACellMuon::areAlignedRZ(
          r1,
          z1,
          ri,
          zi,
          ro,
          zo,
          ptmin,
          isBarrel ? CAThetaCutBarrel : CAThetaCutForward);  // 2.f*thetaCut); // FIXME tune cuts

      auto dDirMax = 0.;
      auto innerID = thisCell.inner_detIndex(hh);
      auto outerID = thisCell.outer_detIndex(hh);

      if (innerID < 5 && outerID < 5) dDirMax = 0.6;
      else if (innerID < 5 && outerID >= 5) dDirMax = 0.3;
      else dDirMax = 0.75;
      aligned = true;
      if (aligned && thisCell.dDirCut(hh,
                                     oc,
                                     dDirMax)) {  // FIXME tune cuts
        oc.addOuterNeighbor(cellIndex, *cellNeighbors);
        thisCell.setUsedBit(1);
        oc.setUsedBit(1);
      }
    }  // loop on inner cells
  }    // loop on outer cells
}



__global__ void kernel_find_ntuplets(GPUCACellMuon::Hits const *__restrict__ hhp,
                                     GPUCACellMuon *__restrict__ cells,
                                     uint32_t const *nCells,
                                     gpuMuonDoublets::CellTracksVector *cellTracks,
                                     HitContainer *foundNtuplets,
                                     cms::cuda::AtomicPairCounter *apc,
                                     Quality *__restrict__ quality,
                                     unsigned int minHitsPerNtuplet) {

  // recursive: not obvious to widen
  auto const &hh = *hhp;
  auto first = threadIdx.x + blockIdx.x * blockDim.x;
  for (int idx = first, nt = (*nCells); idx < nt; idx += gridDim.x * blockDim.x) {
    auto const &thisCell = cells[idx];
    if (thisCell.isKilled())
      continue;  // cut by earlyFishbone
    // we require at least three hits...
    if (thisCell.outerNeighbors().empty())
      continue;
    auto pid = thisCell.layerPairId();
    auto doit = true;
    if (doit) {
      GPUCACellMuon::TmpTuple stack;
      stack.reset();
      thisCell.find_ntuplets(hh, cells, *cellTracks, *foundNtuplets, *apc, quality, stack, minHitsPerNtuplet, pid < 3);
      assert(stack.empty());
    }
  }
}


__global__ void kernel_mark_used(GPUCACellMuon::Hits const *__restrict__ hhp,
                                 GPUCACellMuon *__restrict__ cells,
                                 uint32_t const *nCells) {
  auto first = threadIdx.x + blockIdx.x * blockDim.x;
  for (int idx = first, nt = (*nCells); idx < nt; idx += gridDim.x * blockDim.x) {
    auto &thisCell = cells[idx];
    if (!thisCell.tracks().empty())
      thisCell.setUsedBit(2);
  }
}


__global__ void kernel_earlyDuplicateRemover(GPUCACellMuon const *cells,
                                             uint32_t const *__restrict__ nCells,
                                             HitContainer *foundNtuplets,
                                             Quality *quality,
                                             bool dupPassThrough) {
  // quality to mark rejected
  constexpr auto reject = L2MuonTrack::Quality::edup;  /// cannot be loose

  assert(nCells);
  auto first = threadIdx.x + blockIdx.x * blockDim.x;
  for (int idx = first, nt = (*nCells); idx < nt; idx += gridDim.x * blockDim.x) {
    auto const &thisCell = cells[idx];

    if (thisCell.tracks().size() < 2)
      continue;
    //if (0==thisCell.theUsed) continue;
    // if (thisCell.theDoubletId < 0) continue;

    uint32_t maxNh = 0;

    // find maxNh
    for (auto it : thisCell.tracks()) {
      auto nh = foundNtuplets->size(it);
      maxNh = std::max(nh, maxNh);
    }

    // quad pass through (leave it her for tests)
    //  maxNh = std::min(4U, maxNh);

    for (auto it : thisCell.tracks()) {
      if (foundNtuplets->size(it) < maxNh)
        quality[it] = reject;  //no race:  simple assignment of the same constant
    }
  }
}

__global__ void kernel_countMultiplicity(HitContainer const *__restrict__ foundNtuplets,
                                         Quality const *__restrict__ quality,
                                         caConstants::TupleMultiplicity *tupleMultiplicity) {
  auto first = blockIdx.x * blockDim.x + threadIdx.x;
  for (int it = first, nt = foundNtuplets->nOnes(); it < nt; it += gridDim.x * blockDim.x) {
    auto nhits = foundNtuplets->size(it);
    if (nhits < 3)
      continue;
    if (quality[it] == L2MuonTrack::Quality::edup)
      continue;
    assert(quality[it] == L2MuonTrack::Quality::bad);
    if (nhits > 5)
      printf("wrong mult %d %d\n", it, nhits);
    assert(nhits < 8);
    tupleMultiplicity->count(nhits);
  }
}

__global__ void kernel_fillMultiplicity(HitContainer const *__restrict__ foundNtuplets,
                                        Quality const *__restrict__ quality,
                                        caConstants::TupleMultiplicity *tupleMultiplicity) {
  auto first = blockIdx.x * blockDim.x + threadIdx.x;
  for (int it = first, nt = foundNtuplets->nOnes(); it < nt; it += gridDim.x * blockDim.x) {
    auto nhits = foundNtuplets->size(it);
    if (nhits < 3)
      continue;
    if (quality[it] == L2MuonTrack::Quality::edup)
      continue;
    assert(quality[it] == L2MuonTrack::Quality::bad);
    if (nhits > 5)
      printf("wrong mult %d %d\n", it, nhits);
    assert(nhits < 8);
    tupleMultiplicity->fill(nhits, it);
  }
}

__global__ void kernel_fillHitDetIndices(HitContainer const *__restrict__ tuples,
                                         MuonSegmentsCUDAView const *__restrict__ hhp,
                                         HitContainer *__restrict__ hitDetIndices) {
  int first = blockDim.x * blockIdx.x + threadIdx.x;
  // copy offsets
  for (int idx = first, ntot = tuples->totOnes(); idx < ntot; idx += gridDim.x * blockDim.x) {
    hitDetIndices->off[idx] = tuples->off[idx];
  }
  // fill hit indices
  auto const &hh = *hhp;
  auto nSegments = hh.nSegments();
  for (int idx = first, ntot = tuples->size(); idx < ntot; idx += gridDim.x * blockDim.x) {
    assert(tuples->content[idx] < nSegments);
    hitDetIndices->content[idx] = hh.layerID(tuples->content[idx]);
  }
}

__global__ void kernel_extractNtuplets(MuonSegmentsCUDAView const *__restrict__ hhp,
					 MuonSegmentNtupletsCUDA *ntuplets,
					 HitContainer const *__restrict__ tuples) {

  auto const& __restrict__ hh = *hhp;
  int first = blockDim.x * blockIdx.x + threadIdx.x;
  if (first == 0){
	for (int k = 0; k < 100; k++) {
		ntuplets->nNtuplets = 0;
		ntuplets->segmentsInNtuplet[k] = 0;
		ntuplets->gx1[k] = -999.;
		ntuplets->gy1[k] = -999.;
		ntuplets->gz1[k] = -999.;
		ntuplets->gphi1[k] = - 999.;
		ntuplets->gr1[k] = -999.;
		ntuplets->layerID1[k] = -999;

		ntuplets->gx2[k] = -999.;
		ntuplets->gy2[k] = -999.;
		ntuplets->gz2[k] = -999.;
		ntuplets->gphi2[k] = - 999.;
		ntuplets->gr2[k] = -999.;
		ntuplets->layerID2[k] = -999;

		ntuplets->gx3[k] = -999.;
		ntuplets->gy3[k] = -999.;
		ntuplets->gz3[k] = -999.;
		ntuplets->gphi3[k] = - 999.;
		ntuplets->gr3[k] = -999.;
		ntuplets->layerID3[k] = -999;

		ntuplets->gx4[k] = -999.;
		ntuplets->gy4[k] = -999.;
		ntuplets->gz4[k] = -999.;
		ntuplets->gphi4[k] = - 999.;
		ntuplets->gr4[k] = -999.;
		ntuplets->layerID4[k] = -999;

	}

  }
  for (int idx = first, ntot = tuples->size(); idx < ntot; idx += gridDim.x * blockDim.x) {
	ntuplets->nNtuplets = ntot;
	ntuplets->segmentsInNtuplet[idx] = tuples->size(idx);
        auto const *segmentID = tuples->begin(idx);

        for (unsigned int i = 0; i < tuples->size(idx); ++i) {
        	auto index = segmentID[i]; 	
		if (i == 0){	
			ntuplets->gx1[idx] = hh.gx(index) ;
	  		ntuplets->gy1[idx] = hh.gy(index);
	  		ntuplets->gz1[idx] = hh.gz(index);
	  		ntuplets->gphi1[idx] = hh.phi(index);
	  		ntuplets->gr1[idx] = hh.gr(index);
	  		ntuplets->layerID1[idx] = hh.layerID(index);
		}
		if (i == 1){	  
			ntuplets->gx2[idx] = hh.gx(index) ;
			ntuplets->gy2[idx] = hh.gy(index);
			ntuplets->gz2[idx] = hh.gz(index);
			ntuplets->gphi2[idx] = hh.phi(index);
			ntuplets->gr2[idx] = hh.gr(index);
			ntuplets->layerID2[idx] = hh.layerID(index);
		}
		if (i == 2){	  
			ntuplets->gx3[idx] = hh.gx(index) ;
			ntuplets->gy3[idx] = hh.gy(index);
			ntuplets->gz3[idx] = hh.gz(index);
			ntuplets->gphi3[idx] = hh.phi(index);
			ntuplets->gr3[idx] = hh.gr(index);
			ntuplets->layerID3[idx] = hh.layerID(index);
		}
		if (i == 3){	  
			ntuplets->gx4[idx] = hh.gx(index) ;
			ntuplets->gy4[idx] = hh.gy(index);
			ntuplets->gz4[idx] = hh.gz(index);
			ntuplets->gphi4[idx] = hh.phi(index);
			ntuplets->gr4[idx] = hh.gr(index);
			ntuplets->layerID4[idx] = hh.layerID(index);
		}


	}
  }
}
