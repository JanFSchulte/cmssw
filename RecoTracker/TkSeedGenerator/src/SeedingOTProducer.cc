#include "RecoTracker/TkSeedGenerator/interface/SeedingOTProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"

#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MeasurementDet/interface/TrajectoryMeasurementGroup.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

SeedingOTProducer::SeedingOTProducer(edm::ParameterSet const& conf)
    : theUpdator_(nullptr),
      tkMeasEventToken_(consumes<MeasurementTrackerEvent>(conf.getParameter<edm::InputTag>("trackerEvent"))) {
  vhProducerToken_ = consumes<VectorHitCollectionNew>(edm::InputTag(conf.getParameter<edm::InputTag>("src")));
  beamSpotToken_ = consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpotLabel"));
  updatorName_ = conf.getParameter<std::string>("updator");
  produces<TrajectorySeedCollection>();
}

SeedingOTProducer::~SeedingOTProducer() {}

void SeedingOTProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("siPhase2VectorHits", "vectorHitsAccepted"));
  desc.add<edm::InputTag>("trackerEvent", edm::InputTag("MeasurementTrackerEvent"));
  desc.add<edm::InputTag>("beamSpotLabel", edm::InputTag("offlineBeamSpot"));
  desc.add<std::string>("updator", std::string("KFUpdator"));
  descriptions.add("SeedingOTProducer", desc);
}

void SeedingOTProducer::produce(edm::Event& event, const edm::EventSetup& es) {
  std::unique_ptr<TrajectorySeedCollection> seedsWithVHs(new TrajectorySeedCollection());

  edm::ESHandle<TrackerTopology> tTopoHandle;
  es.get<TrackerTopologyRcd>().get(tTopoHandle);
  tkTopo_ = tTopoHandle.product();

  edm::ESHandle<MeasurementTracker> measurementTrackerHandle;
  es.get<CkfComponentsRecord>().get(measurementTrackerHandle);
  measurementTracker_ = measurementTrackerHandle.product();
  edm::Handle<MeasurementTrackerEvent> measurementTrackerEvent;
  event.getByToken(tkMeasEventToken_, measurementTrackerEvent);

  layerMeasurements_ = std::make_unique<LayerMeasurements>(*measurementTrackerHandle, *measurementTrackerEvent);

  edm::ESHandle<Chi2MeasurementEstimatorBase> est;
  es.get<TrackingComponentsRecord>().get("Chi2", est);
  estimator_ = est.product();

  edm::ESHandle<Propagator> prop;
  es.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", prop);
  propagator_ = prop.product();

  edm::ESHandle<MagneticField> magFieldHandle;
  es.get<IdealMagneticFieldRecord>().get(magFieldHandle);
  magField_ = magFieldHandle.product();

  edm::ESHandle<TrajectoryStateUpdator> updatorHandle;
  es.get<TrackingComponentsRecord>().get(updatorName_, updatorHandle);
  theUpdator_ = updatorHandle.product();

  edm::Handle<reco::BeamSpot> beamSpotH;
  event.getByToken(beamSpotToken_, beamSpotH);
  if (beamSpotH.isValid()) {
    beamSpot_ = beamSpotH.product();
  }

  // Get the vector hits
  edm::Handle<VectorHitCollectionNew> vhs;
  event.getByToken(vhProducerToken_, vhs);

  TrajectorySeedCollection const& tempSeeds = run(vhs);
  for (TrajectorySeedCollection::const_iterator qIt = tempSeeds.begin(); qIt < tempSeeds.end(); ++qIt) {
    seedsWithVHs->push_back(*qIt);
  }

  seedsWithVHs->shrink_to_fit();
  event.put(std::move(seedsWithVHs));
}

TrajectorySeedCollection SeedingOTProducer::run(edm::Handle<VectorHitCollectionNew> VHs) {
  TrajectorySeedCollection result;

  //check if all the first three layers have VHs
  std::vector<VectorHit> VHseedsL1 = collectVHsOnLayer(VHs, 1);
  std::vector<VectorHit> VHseedsL2 = collectVHsOnLayer(VHs, 2);
  std::vector<VectorHit> VHseedsL3 = collectVHsOnLayer(VHs, 3);
  if (VHseedsL1.empty() || VHseedsL2.empty() || VHseedsL3.empty()) {
    return result;
  }

  //seeds are built in the L3 of the OT
  const BarrelDetLayer* barrelOTLayer2 = measurementTracker_->geometricSearchTracker()->tobLayers().at(1);

  //the search propag directiondepend on the sign of signZ*signPz, while the building is always the contrary
  Propagator* searchingPropagator = &*propagator_->clone();
  Propagator* buildingPropagator = &*propagator_->clone();
  buildingPropagator->setPropagationDirection(alongMomentum);

  for (auto hitL3 : VHseedsL3) {
    //building a tsos out of a VectorHit
    const TrajectoryStateOnSurface initialTSOS = buildInitialTSOS(hitL3);
    float signZ = copysign(1.0, initialTSOS.globalPosition().z());
    float signPz = copysign(1.0, initialTSOS.globalMomentum().z());

    //set the direction of the propagator
    if (signZ * signPz > 0.0)
      searchingPropagator->setPropagationDirection(oppositeToMomentum);
    if (signZ * signPz < 0.0)
      searchingPropagator->setPropagationDirection(alongMomentum);

    //find vHits in layer 2
    std::vector<TrajectoryMeasurement> measurementsL2 =
        layerMeasurements_->measurements(*barrelOTLayer2, initialTSOS, *searchingPropagator, *estimator_);

    std::vector<TrajectoryMeasurement>::iterator measurementsL2end =
        std::remove_if(measurementsL2.begin(), measurementsL2.end(), isInvalid());
    measurementsL2.erase(measurementsL2end, measurementsL2.end());

    if (!measurementsL2.empty()) {
      //not sure if building it everytime takes time/memory
      const DetLayer* barrelOTLayer1 = measurementTracker_->geometricSearchTracker()->tobLayers().at(0);

      for (auto mL2 : measurementsL2) {
        const TrackingRecHit* hitL2 = mL2.recHit().get();

        //propagate to the L2 and update the TSOS
        std::pair<bool, TrajectoryStateOnSurface> updatedTSOS =
            propagateAndUpdate(initialTSOS, *searchingPropagator, *hitL2);
        if (!updatedTSOS.first)
          continue;

        //searching possible VHs in L1
        std::vector<TrajectoryMeasurement> measurementsL1 =
            layerMeasurements_->measurements(*barrelOTLayer1, updatedTSOS.second, *searchingPropagator, *estimator_);
        std::vector<TrajectoryMeasurement>::iterator measurementsL1end =
            std::remove_if(measurementsL1.begin(), measurementsL1.end(), isInvalid());
        measurementsL1.erase(measurementsL1end, measurementsL1.end());

        if (!measurementsL1.empty()) {
          for (auto mL1 : measurementsL1) {
            const TrackingRecHit* hitL1 = mL1.recHit().get();

            //propagate to the L1 and update the TSOS
            std::pair<bool, TrajectoryStateOnSurface> updatedTSOSL1 =
                propagateAndUpdate(updatedTSOS.second, *searchingPropagator, *hitL1);
            if (!updatedTSOSL1.first)
              continue;

            edm::OwnVector<TrackingRecHit> container;
            container.push_back(hitL1->clone());
            container.push_back(hitL2->clone());
            container.push_back(hitL3.clone());

            //building trajectory inside-out
            if (searchingPropagator->propagationDirection() == alongMomentum) {
              buildingPropagator->setPropagationDirection(oppositeToMomentum);
            } else if (searchingPropagator->propagationDirection() == oppositeToMomentum) {
              buildingPropagator->setPropagationDirection(alongMomentum);
            }

            updatedTSOSL1.second.rescaleError(100);

            TrajectoryStateOnSurface updatedTSOSL1_final = theUpdator_->update(updatedTSOSL1.second, *hitL1);
            if
              UNLIKELY(!updatedTSOSL1_final.isValid()) continue;
            std::pair<bool, TrajectoryStateOnSurface> updatedTSOSL2_final =
                propagateAndUpdate(updatedTSOSL1_final, *buildingPropagator, *hitL2);
            std::pair<bool, TrajectoryStateOnSurface> updatedTSOSL3_final =
                propagateAndUpdate(updatedTSOSL2_final.second, *buildingPropagator, hitL3);
            TrajectorySeed ts =
                createSeed(updatedTSOSL3_final.second, container, hitL3.geographicalId(), *buildingPropagator);
            result.push_back(ts);
          }
        }
      }
    }
  }

  return result;
}

unsigned int SeedingOTProducer::checkLayer(unsigned int iidd) {
  StripSubdetector strip = StripSubdetector(iidd);
  unsigned int subid = strip.subdetId();
  if (subid == StripSubdetector::TIB || subid == StripSubdetector::TOB) {
    return tkTopo_->layer(iidd);
  }
  return 0;
}

std::vector<VectorHit> SeedingOTProducer::collectVHsOnLayer(edm::Handle<VectorHitCollectionNew> VHs,
                                                              unsigned int layerNumber) {
  const VectorHitCollectionNew& input = *VHs;
  std::vector<VectorHit> VHsOnLayer;
  if (!input.empty()) {
    for (auto DSViter : input) {
      if (checkLayer(DSViter.id()) == layerNumber) {
        for (auto vh : DSViter) {
          VHsOnLayer.push_back(vh);
        }
      }
    }
  }

  return VHsOnLayer;
}

void SeedingOTProducer::printVHsOnLayer(edm::Handle<VectorHitCollectionNew> VHs, unsigned int layerNumber, std::ostream& stream) {
  const VectorHitCollectionNew& input = *VHs;
  if (!input.empty()) {
    for (auto DSViter : input) {
      for (auto vh : DSViter) {
        if (checkLayer(DSViter.id()) == layerNumber)
          stream << " VH in layer " << layerNumber << " >> " << vh;
      }
    }
  } else {
    stream << " No VHs in layer " << layerNumber << ".";
  }
}

const TrajectoryStateOnSurface SeedingOTProducer::buildInitialTSOS(VectorHit& vHit) {
  // having fun with theta
  Global3DVector gv(vHit.globalPosition().x(), vHit.globalPosition().y(), vHit.globalPosition().z());
  float theta = gv.theta();
  // gv transform to local (lv)
  const Local3DVector lv(vHit.det()->surface().toLocal(gv));

  //FIXME::charge is fine 1 every two times!!
  int charge = 1;
  float p = vHit.momentum(magField_);
  float x = vHit.localPosition().x();
  float y = vHit.localPosition().y();
  float dx = vHit.localDirection().x();
  // for dy use second component of the lv renormalized to the z component
  float dy = (lv.z() != 0 ? lv.y() / lv.z() : 0.);

  // Pz and Dz should have the same sign
  float signPz = copysign(1.0, vHit.globalPosition().z());

  LocalTrajectoryParameters ltpar2(charge / p, dx, dy, x, y, signPz);
  AlgebraicSymMatrix mat = assign44To55(vHit.parametersError());
  // set the error on 1/p
  mat[0][0] = pow(computeInverseMomentumError(vHit, theta, magField_, beamSpot_->sigmaZ()), 2);

  //building tsos
  LocalTrajectoryError lterr(asSMatrix<5>(mat));
  const TrajectoryStateOnSurface tsos(ltpar2, lterr, vHit.det()->surface(), magField_);

  return tsos;
}

AlgebraicSymMatrix SeedingOTProducer::assign44To55(AlgebraicSymMatrix mat44) {
  if (mat44.num_row() != 4 || mat44.num_col() != 4)
    assert("Wrong dimension! This should be a 4x4 matrix!");

  AlgebraicSymMatrix result(5, 0);
  for (int i = 1; i < 5; i++) {
    for (int j = 1; j < 5; j++) {
      result[i][j] = mat44[i - 1][j - 1];
    }
  }
  return result;
}

std::pair<bool, TrajectoryStateOnSurface> SeedingOTProducer::propagateAndUpdate(
    const TrajectoryStateOnSurface initialTSOS, const Propagator& prop, const TrackingRecHit& hit) {
  TrajectoryStateOnSurface propTSOS = prop.propagate(initialTSOS, hit.det()->surface());
  TrajectoryStateOnSurface updatedTSOS = theUpdator_->update(propTSOS, hit);
  if
    UNLIKELY(!updatedTSOS.isValid()) return std::make_pair(false, updatedTSOS);
  return std::make_pair(true, updatedTSOS);
}

float SeedingOTProducer::computeGlobalThetaError(const VectorHit& vh, const double sigmaZ_beamSpot) {
  double derivative =
      vh.globalPosition().perp() / (pow(vh.globalPosition().z(), 2) + pow(vh.globalPosition().perp(), 2));
  double derivative2 = pow(derivative, 2);
  return pow(derivative2 * vh.lowerGlobalPosErr().czz() + derivative2 * pow(sigmaZ_beamSpot, 2), 0.5);
}

float SeedingOTProducer::computeInverseMomentumError(VectorHit& vh,
                                                       const float globalTheta,
                                                       const MagneticField* magField,
                                                       const double sigmaZ_beamSpot) {
  //for pT > 2GeV, 1/pT has sigma = 1/sqrt(12)
  float varianceInverseTransvMomentum = 1. / 12;
  double derivativeTheta2 = pow(cos(globalTheta) / vh.transverseMomentum(magField), 2);
  double derivativeInverseTransvMomentum2 = pow(sin(globalTheta), 2);
  float thetaError = computeGlobalThetaError(vh, sigmaZ_beamSpot);
  return pow(derivativeTheta2 * pow(thetaError, 2) + derivativeInverseTransvMomentum2 * varianceInverseTransvMomentum,
             0.5);
}

TrajectorySeed SeedingOTProducer::createSeed(const TrajectoryStateOnSurface& tsos,
                                               const edm::OwnVector<TrackingRecHit>& container,
                                               const DetId& id,
                                               const Propagator& prop) {
  PTrajectoryStateOnDet seedTSOS = trajectoryStateTransform::persistentState(tsos, id.rawId());
  return TrajectorySeed(seedTSOS, container, prop.propagationDirection());
}
