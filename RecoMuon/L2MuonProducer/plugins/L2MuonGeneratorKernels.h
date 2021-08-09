#ifndef RecoMuon_L2MuonProducer_plugins_L2MuonGeneratorKernels_h
#define RecoMuon_L2MuonProducer_plugins_L2MuonGeneratorKernels_h

#include "CUDADataFormats/Track/interface/L2MuonTrackHeterogeneous.h"
//#include "GPUCACell.h" // That's were the hard part will have to go
#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDA.h"

namespace l2MuonGenerator{

  struct Counters {
    unsigned long long nEvents;
    unsigned long long nHits;
    unsigned long long nCells;
    unsigned long long nTuples;
    unsigned long long nFitTracks;
    unsigned long long nLooseTracks;
    unsigned long long nGoodTracks;
    unsigned long long nUsedHits;
    unsigned long long nDupHits;
    unsigned long long nKilledCells;
    unsigned long long nEmptyCells;
    unsigned long long nZeroTrackCells;
  };

  struct Params {
    Params(bool onGPU,
           bool doStats)
        : onGPU_(onGPU), 
          doStats_(doStats){}

    const bool onGPU_;
    const bool doStats_;
                                                                                                                                                                                                                                             };  // Params

}
template <typename TTraits>
class L2MuonGeneratorKernels {
public:
  using Traits = TTraits;

  using Params = l2MuonGenerator::Params;
  using Counters = l2MuonGenerator::Counters;

  template <typename T>
  using unique_ptr = typename Traits::template unique_ptr<T>;


  L2MuonGeneratorKernels(Params const& params)
      : params_(params) {}
  ~L2MuonGeneratorKernels() = default;

  //void launchKernels(HitsOnCPU const& hh, TkSoA* tuples_d, cudaStream_t cudaStream);

  void buildL2Muons(MuonSegmentsCUDA const& hSegments, L2MuonTrack::TrackSoA* l2Muons_d, cudaStream_t stream) const;//
  void allocateOnGPU(int32_t nSegments, cudaStream_t stream);
  //void cleanup(cudaStream_t cudaStream);

  //static void printCounters(Counters const* counters);
  //void setCounters(Counters* counters) { counters_ = counters; }


private:
//  Counters* counters_ = nullptr;

  Params const& params_;
  
  cms::cuda::AtomicPairCounter* device_segmentTuple_apc_ = nullptr;

  unique_ptr<cms::cuda::AtomicPairCounter::c_type[]> device_storage_;
};

using L2MuonGeneratorKernelsGPU = L2MuonGeneratorKernels<cms::cudacompat::GPUTraits>;
#endif
#
