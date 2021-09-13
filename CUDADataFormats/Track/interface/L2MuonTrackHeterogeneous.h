#ifndef CUDADataFormats_Track_L2MuonTrackHeterogeneous_h
#define CUDADataFormats_Track_L2MuonTrackHeterogeneous_h

#include "CUDADataFormats/Common/interface/HeterogeneousSoA.h"
#include "CUDADataFormats/Track/interface/L2MuonTrackSoAHeterogeneousT.h"

using L2MuonTrackHeterogeneous = HeterogeneousSoA<L2MuonTrack::TrackSoA>;

#endif  // #ifndef CUDADataFormats_Track_L2MuonTrackHeterogeneous_h
