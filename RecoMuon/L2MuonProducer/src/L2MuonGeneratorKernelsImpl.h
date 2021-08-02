#include <cmath>
#include <cstdint>
#include <limits>

#include <cuda_runtime.h>

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cuda_assert.h"

#include "L2MuonGeneratorKernels.h"
#include "CUDADataFormats/Track/interface/L2MuonTrackHeterogeneous.h"


__global__ void makeL2Muon(L2MuonTrack::TrackSoA* l2Muons_d, DTRecSegment4DCUDA const& dtSegments_d, CSCSegmentCUDA const& cscSegments_d) {




}
