#include <array>
#include <cassert>
#include <functional>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "HeterogeneousCore/CUDAServices/interface/CUDAService.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"

#include "L2MuonGeneratorOnGPU.h"


L2MuonGeneratorOnGPU::L2MuonGeneratorOnGPU(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC)
    : m_params(cfg.getParameter<bool>("onGPU"),
               cfg.getParameter<bool>("doStats"),
               cfg.getParameter<int>("minHitsPerNtuplet"),
               cfg.getParameter<int>("maxNumberOfDoublets"),
               cfg.getParameter<bool>("includeJumpingForwardDoublets"),
               cfg.getParameter<bool>("doZ0Cut"),
               cfg.getParameter<bool>("doPtCut"),
               cfg.getParameter<double>("ptmin"),
               cfg.getParameter<double>("CAThetaCutBarrel"),
               cfg.getParameter<double>("CAThetaCutForward"),
               cfg.getParameter<double>("hardCurvCut"),
               cfg.getParameter<double>("dcaCutInnerTriplet"),
               cfg.getParameter<double>("dcaCutOuterTriplet"),
               cfg.getParameter<bool>("dupPassThrough")) {

  if (m_params.onGPU_) {
    // allocate pinned host memory only if CUDA is available
    edm::Service<CUDAService> cs;
    if (cs and cs->enabled()) {
      cudaCheck(cudaMalloc(&m_counters, sizeof(Counters)));
      cudaCheck(cudaMemset(m_counters, 0, sizeof(Counters)));
    }
  } else {
    m_counters = new Counters();
    memset(m_counters, 0, sizeof(Counters));
  }
}

L2MuonGeneratorOnGPU::~L2MuonGeneratorOnGPU() {
  if (m_params.onGPU_) {
    // print the gpu statistics and free pinned host memory only if CUDA is available
    edm::Service<CUDAService> cs;
    if (cs and cs->enabled()) {
      if (m_params.doStats_) {
        // crash on multi-gpu processes
        //CAHitNtupletGeneratorKernelsGPU::printCounters(m_counters);
      }
      cudaFree(m_counters);
    }
  } else {
    if (m_params.doStats_) {
      //CAHitNtupletGeneratorKernelsCPU::printCounters(m_counters);
    }
    delete m_counters;
  }
}


void L2MuonGeneratorOnGPU::fillDescriptions(edm::ParameterSetDescription& desc) {
     desc.add<bool>("onGPU", true);
     desc.add<bool>("doStats", true);
     desc.add<int>("minHitsPerNtuplet",3);
     desc.add<int>("maxNumberOfDoublets",1000);
     desc.add<bool>("includeJumpingForwardDoublets",true);
     desc.add<bool>("doZ0Cut",true);
     desc.add<bool>("doPtCut",true);
     desc.add<double>("ptmin", 0.9f)->setComment("Cut on minimum pt");
     desc.add<double>("CAThetaCutBarrel", 0.002f)->setComment("Cut on RZ alignement for Barrel");
     desc.add<double>("CAThetaCutForward", 0.003f)->setComment("Cut on RZ alignment for Forward");
     desc.add<double>("hardCurvCut", 1.f / (0.35 * 87.f))->setComment("Cut on minimum curvature");
     desc.add<double>("dcaCutInnerTriplet", 0.15f)->setComment("Cut on origin radius when the inner hit is on BPix1");
     desc.add<double>("dcaCutOuterTriplet", 0.25f)->setComment("Cut on origin radius when the outer hit is on BPix1");
     desc.add<bool>("dupPassThrough", false)->setComment("Do not reject duplicate");
}

L2MuonTrackHeterogeneous L2MuonGeneratorOnGPU::makeTuplesAsync(MuonSegmentsCUDA const& muonSegments_h,
                                                                    float bfield,
                                                                    cudaStream_t stream) const {
  L2MuonTrackHeterogeneous tracks(cms::cuda::make_device_unique<L2MuonTrack::TrackSoA>(stream));

  auto* soa = tracks.get();

  L2MuonGeneratorKernelsGPU kernels(m_params);
  //kernels.setCounters(m_counters);
  int32_t nSegments = muonSegments_h.nSegments();
  kernels.allocateOnGPU(nSegments,stream);
  kernels.buildDoublets(muonSegments_h, stream);
  kernels.buildL2Muons(muonSegments_h, soa, stream);
  kernels.fillHitDetIndices(muonSegments_h.view(), soa, stream);  // in principle needed only if Hits not "available"


 // bool fit5as4 = true;
 // HelixFitOnGPU fitter(bfield, fit5as4);
 // fitter.allocateOnGPU(&(soa->hitIndices), kernels.tupleMultiplicity(), soa);
  //if (m_params.useRiemannFit_) {
  //  fitter.launchRiemannKernels(hits_d.view(), hits_d.nHits(), caConstants::maxNumberOfQuadruplets, stream);
  //} else {
 // fitter.launchBrokenLineKernels(muonSegments_h.view(), muonSegments_h.nSegments(), caConstants::maxNumberOfQuadruplets, stream);
  //}
  //kernels.classifyTuples(hits_d, soa, stream);

  return tracks;
}
