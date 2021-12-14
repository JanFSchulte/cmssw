#ifndef RecoMuon_L2MuonProducer_plugins_L2MuonGeneratorKernels_h
#define RecoMuon_L2MuonProducer_plugins_L2MuonGeneratorKernels_h

#include "CUDADataFormats/Track/interface/L2MuonTrackHeterogeneous.h"
#include "GPUCACellMuon.h" // That's were the hard part will have to go
#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentPairsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentNtupletsCUDA.h"

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

  using HitsView = MuonSegmentsCUDAView;
  using HitsOnGPU = MuonSegmentsCUDAView;

  using HitToTuple = caConstants::HitToTuple;
  using TupleMultiplicity = caConstants::TupleMultiplicity;

  using Quality = pixelTrack::Quality;
  using TkSoA = pixelTrack::TrackSoA;
  using HitContainer = pixelTrack::HitContainer;


  struct QualityCuts {
    // chi2 cut = chi2Scale * (chi2Coeff[0] + pT/GeV * (chi2Coeff[1] + pT/GeV * (chi2Coeff[2] + pT/GeV * chi2Coeff[3])))
    float chi2Coeff[4];
    float chi2MaxPt;  // GeV
    float chi2Scale;

    struct Region {
      float maxTip;  // cm
      float minPt;   // GeV
      float maxZip;  // cm
    };

    Region triplet;
    Region quadruplet;
  };

  struct Params {
    Params(bool onGPU,
           bool doStats,
           uint32_t minHitsPerNtuplet,
           uint32_t maxNumberOfDoublets,
           bool includeJumpingForwardDoublets,
           bool doZ0Cut,
           bool doPtCut,
           float ptmin,
           float CAThetaCutBarrel,
           float CAThetaCutForward,
           float hardCurvCut,
           float dcaCutInnerTriplet,
           float dcaCutOuterTriplet,
           bool dupPassThrough)
        : onGPU_(onGPU), 
          doStats_(doStats),
          minHitsPerNtuplet_(minHitsPerNtuplet),
          maxNumberOfDoublets_(maxNumberOfDoublets),
          includeJumpingForwardDoublets_(includeJumpingForwardDoublets),
          doZ0Cut_(doZ0Cut),
          doPtCut_(doPtCut),
          ptmin_(ptmin),
          CAThetaCutBarrel_(CAThetaCutBarrel),
          CAThetaCutForward_(CAThetaCutForward),
          hardCurvCut_(hardCurvCut),
          dcaCutInnerTriplet_(dcaCutInnerTriplet),
          dcaCutOuterTriplet_(dcaCutOuterTriplet),
          dupPassThrough_(dupPassThrough){}

    const bool onGPU_;
    const bool doStats_;
    const uint32_t minHitsPerNtuplet_;
    const uint32_t maxNumberOfDoublets_;
    const bool includeJumpingForwardDoublets_;
    const bool doZ0Cut_;
    const bool doPtCut_; 
    const float ptmin_;
    const float CAThetaCutBarrel_;
    const float CAThetaCutForward_;
    const float hardCurvCut_;
    const float dcaCutInnerTriplet_;
    const float dcaCutOuterTriplet_;
    const bool dupPassThrough_;
  };  // Params

}
template <typename TTraits>
class L2MuonGeneratorKernels {
public:
  using Traits = TTraits;

  using QualityCuts = l2MuonGenerator::QualityCuts;
  using Params = l2MuonGenerator::Params;
  using Counters = l2MuonGenerator::Counters;

  template <typename T>
  using unique_ptr = typename Traits::template unique_ptr<T>;

  using SegmentsView = MuonSegmentsCUDAView;
  using SegmentsOnGPU = MuonSegmentsCUDAView;
  using SegmentsOnCPU = MuonSegmentsCUDA;

  using HitToTuple = caConstants::HitToTuple;
  using TupleMultiplicity = caConstants::TupleMultiplicity;

  using Quality = pixelTrack::Quality;
  using TkSoA = L2MuonTrack::TrackSoA;
  using HitContainer = pixelTrack::HitContainer;

  L2MuonGeneratorKernels(Params const& params)
      : params_(params),  paramsMaxDoubletes3Quarters_(3 * params.maxNumberOfDoublets_ / 4) {}
  ~L2MuonGeneratorKernels() = default;

  void launchKernels(SegmentsOnCPU const& hh, TkSoA* tuples_d, cudaStream_t cudaStream);

  void buildDoublets(SegmentsOnCPU const& hSegments, cudaStream_t stream);
  void buildAndReturnDoublets(SegmentsOnCPU const& hSegments, MuonSegmentPairsCUDA* pairs_d, cudaStream_t stream);
  void extractNtuplets(SegmentsOnCPU const& hSegments, MuonSegmentNtupletsCUDA* ntuplets_d, L2MuonTrack::TrackSoA* l2Muons_d, cudaStream_t stream);
  void buildL2Muons(SegmentsOnCPU const& hSegments, L2MuonTrack::TrackSoA* l2Muons_d, cudaStream_t stream) const;//
  void allocateOnGPU(int32_t nSegments, cudaStream_t stream);
  void fillHitDetIndices(SegmentsView const* hv, TkSoA* tuples_d, cudaStream_t cudaStream);
  //void cleanup(cudaStream_t cudaStream);

  //static void printCounters(Counters const* counters);
  //void setCounters(Counters* counters) { counters_ = counters; }


private:
  Counters* counters_ = nullptr;

  unique_ptr<unsigned char[]> cellStorage_;
  unique_ptr<caConstants::CellNeighborsVector> device_theCellNeighbors_;
  caConstants::CellNeighbors* device_theCellNeighborsContainer_;
  unique_ptr<caConstants::CellTracksVector> device_theCellTracks_;
  caConstants::CellTracks* device_theCellTracksContainer_;

  unique_ptr<GPUCACellMuon[]> device_theCells_;
  unique_ptr<GPUCACellMuon::OuterHitOfCell[]> device_isOuterHitOfCell_;
  uint32_t* device_nCells_ = nullptr;

  unique_ptr<HitToTuple> device_hitToTuple_;
  unique_ptr<HitToTuple::Counter[]> device_hitToTupleStorage_;
  HitToTuple::View hitToTupleView_;

  cms::cuda::AtomicPairCounter* device_hitToTuple_apc_ = nullptr;

  cms::cuda::AtomicPairCounter* device_hitTuple_apc_ = nullptr;

  unique_ptr<TupleMultiplicity> device_tupleMultiplicity_;

  unique_ptr<cms::cuda::AtomicPairCounter::c_type[]> device_storage_;
  // params
  Params const& params_;

  const uint32_t paramsMaxDoubletes3Quarters_;
  /// Compute the number of doublet blocks for block size
  inline uint32_t nDoubletBlocks(uint32_t blockSize) const {
  // We want (3 * params_.maxNumberOfDoublets_ / 4 + blockSize - 1) / blockSize, but first part is pre-computed.
  	return (paramsMaxDoubletes3Quarters_ + blockSize - 1) / blockSize;
  }
  /// Compute the number of quadruplet blocks for block size
  inline uint32_t nQuadrupletBlocks(uint32_t blockSize) const {
  // caConstants::maxNumberOfQuadruplets is a constexpr, so the compiler will pre compute the 3*max/4
  	return (3 * caConstants::maxNumberOfQuadruplets / 4 + blockSize - 1) / blockSize;
  }

};

using L2MuonGeneratorKernelsGPU = L2MuonGeneratorKernels<cms::cudacompat::GPUTraits>;
#endif
#
