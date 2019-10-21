//---------------------------------------------------------------------------
// class SeedingOTProducer
// author: ebrondol
// date: July, 2016
//---------------------------------------------------------------------------

#ifndef RecoTracker_TkSeedGenerator_SeedingOTProducer_h
#define RecoTracker_TkSeedGenerator_SeedingOTProducer_h

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackerRecHit2D/interface/VectorHit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include <ostream>
//#include <iostream>

class TrajectoryStateUpdator;

class SeedingOTProducer : public edm::stream::EDProducer<> {
public:
  explicit SeedingOTProducer(const edm::ParameterSet&);
  ~SeedingOTProducer() override;
  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

  TrajectorySeedCollection run(edm::Handle<VectorHitCollectionNew>);
  unsigned int checkLayer(unsigned int iidd);
  std::vector<VectorHit> collectVHsOnLayer(edm::Handle<VectorHitCollectionNew>, unsigned int);
  void printVHsOnLayer(edm::Handle<VectorHitCollectionNew>, unsigned int, std::ostream & stream);
  const TrajectoryStateOnSurface buildInitialTSOS(VectorHit&);
  AlgebraicSymMatrix assign44To55(AlgebraicSymMatrix);
  std::pair<bool, TrajectoryStateOnSurface> propagateAndUpdate(const TrajectoryStateOnSurface initialTSOS,
                                                               const Propagator&,
                                                               const TrackingRecHit& hit);
  float computeGlobalThetaError(const VectorHit& vh, const double sigmaZ_beamSpot);
  float computeInverseMomentumError(VectorHit& vh,
                                    const float globalTheta,
                                    const MagneticField* magField,
                                    const double sigmaZ_beamSpot);

  TrajectorySeed createSeed(const TrajectoryStateOnSurface& tsos,
                            const edm::OwnVector<TrackingRecHit>& container,
                            const DetId& id,
                            const Propagator& prop);

  struct isInvalid {
    bool operator()(const TrajectoryMeasurement& measurement) {
      return (((measurement).recHit() == nullptr) || !((measurement).recHit()->isValid()) ||
              !((measurement).updatedState().isValid()));
    }
  };

private:
  edm::EDGetTokenT<VectorHitCollectionNew> vhProducerToken_;
  const TrackerTopology* tkTopo_;
  const MeasurementTracker* measurementTracker_;
  std::unique_ptr<LayerMeasurements> layerMeasurements_;
  const MeasurementEstimator* estimator_;
  const Propagator* propagator_;
  const MagneticField* magField_;
  const TrajectoryStateUpdator* theUpdator_;
  const edm::EDGetTokenT<MeasurementTrackerEvent> tkMeasEventToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const reco::BeamSpot* beamSpot_;
  std::string updatorName_;
};

#endif
