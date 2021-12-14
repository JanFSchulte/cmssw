#ifndef CUDADataFormats_Muon_MuonSegmentNtupletsHeterogeneous_h
#define CUDADataFormats_Muon_MuonSegmentNtupletsHeterogeneous_h

#include "CUDADataFormats/Common/interface/HeterogeneousSoA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentNtupletsCUDA.h"

using MuonSegmentNtupletsHeterogeneous = HeterogeneousSoA<MuonSegmentNtupletsCUDA>;

#endif  // #ifndef CUDADataFormats_Muon_MuonSegmentNtupletsHeterogeneous_h
