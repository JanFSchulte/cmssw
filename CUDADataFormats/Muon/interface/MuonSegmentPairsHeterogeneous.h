#ifndef CUDADataFormats_Muon_MuonSegmentPairsHeterogeneous_h
#define CUDADataFormats_Muon_MuonSegmentPairsHeterogeneous_h

#include "CUDADataFormats/Common/interface/HeterogeneousSoA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentPairsCUDA.h"

using MuonSegmentPairsHeterogeneous = HeterogeneousSoA<MuonSegmentPairsCUDA>;

#endif  // #ifndef CUDADataFormats_Muon_MuonSegmentPairsHeterogeneous_h
