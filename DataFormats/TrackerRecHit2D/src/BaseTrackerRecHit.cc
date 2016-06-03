//#define DO_THROW_UNINITIALIZED
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"
#include "DataFormats/Math/interface/ProjectMatrix.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace {
#if defined(DO_THROW_UNINITIALIZED) || defined(DO_INTERNAL_CHECKS_BTR)
  inline void throwExceptionUninitialized(const char *where) {
    throw cms::Exception("BaseTrackerRecHit")
        << "Trying to access " << where
        << " for a RecHit that was read from disk, but since CMSSW_2_1_X local positions are transient.\n"
        << "If you want to get coarse position/error estimation from disk, please set: "
           "ComputeCoarseLocalPositionFromDisk = True \n "
        << " to the TransientTrackingRecHitBuilder you are using from "
           "RecoTracker/TransientTrackingRecHit/python/TTRHBuilders_cff.py";
  }
#endif
  void obsolete() { throw cms::Exception("BaseTrackerRecHit") << "CLHEP is obsolete for Tracker Hits"; }
}  // namespace

#if !defined(VI_DEBUG) && defined(DO_INTERNAL_CHECKS_BTR)
void BaseTrackerRecHit::check() const {
  if (!hasPositionAndError())
    throwExceptionUninitialized("localPosition or Error");
}
#endif

bool BaseTrackerRecHit::hasPositionAndError() const {
  return det();

  //  return (err_.xx() != 0) || (err_.yy() != 0) || (err_.xy() != 0) ||
  //       (pos_.x()  != 0) || (pos_.y()  != 0) || (pos_.z()  != 0);
}

void BaseTrackerRecHit::getKfComponents1D(KfComponentsHolder &holder) const {
#if defined(DO_THROW_UNINITIALIZED)
  if (!hasPositionAndError())
    throwExceptionUninitialized("getKfComponents");
#endif
  AlgebraicVector1 &pars = holder.params<1>();
  pars[0] = pos_.x();

  AlgebraicSymMatrix11 &errs = holder.errors<1>();
  errs(0, 0) = err_.xx();

  ProjectMatrix<double, 5, 1> &pf = holder.projFunc<1>();
  pf.index[0] = 3;

  holder.measuredParams<1>() = AlgebraicVector1(holder.tsosLocalParameters().At(3));
  holder.measuredErrors<1>() = holder.tsosLocalErrors().Sub<AlgebraicSymMatrix11>(3, 3);
}

void BaseTrackerRecHit::getKfComponents2D(KfComponentsHolder &holder) const {
#if defined(DO_THROW_UNINITIALIZED)
  if (!hasPositionAndError())
    throwExceptionUninitialized("getKfComponents");
#endif
  AlgebraicVector2 &pars = holder.params<2>();
  pars[0] = pos_.x();
  pars[1] = pos_.y();

  AlgebraicSymMatrix22 &errs = holder.errors<2>();
  errs(0, 0) = err_.xx();
  errs(0, 1) = err_.xy();
  errs(1, 1) = err_.yy();

  ProjectMatrix<double, 5, 2> &pf = holder.projFunc<2>();
  pf.index[0] = 3;
  pf.index[1] = 4;

  holder.measuredParams<2>() = AlgebraicVector2(&holder.tsosLocalParameters().At(3), 2);
  holder.measuredErrors<2>() = holder.tsosLocalErrors().Sub<AlgebraicSymMatrix22>(3, 3);
}
/*
void 
BaseTrackerRecHit::getKfComponents4D( KfComponentsHolder & holder ) const {

  //if (!hasPositionAndError()) throwExceptionUninitialized("getKfComponents");
  AlgebraicVector4 & pars = holder.params<4>();
  pars[0] = theDirection.x();
  pars[1] = theDirection.y();
  pars[2] = thePosition.x();
  pars[3] = thePosition.y();

  AlgebraicSymMatrix44 & errs = holder.errors<4>();
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      errs(i,j) = theCovMatrix[i][j];
    }
  }

  AlgebraicMatrix45 & proj = holder.projection<4>();
  proj(0,1) = 1;
  proj(1,2) = 1;
  proj(2,3) = 1;
  proj(3,4) = 1;

  ProjectMatrix<double,5,4>  & pf = holder.projFunc<4>();
  pf.index[0] = 1;
  pf.index[1] = 2;
  pf.index[2] = 3;
  pf.index[3] = 4;
  holder.doUseProjFunc();

  holder.measuredParams<4>() = AlgebraicVector4( & holder.tsosLocalParameters().At(1), 4 );
  holder.measuredParams<4>() = AlgebraicVector4( & holder.tsosLocalParameters().At(1), 4 );

}
*/
// obsolete (for what tracker is concerned...) interface
AlgebraicVector BaseTrackerRecHit::parameters() const {
  obsolete();
  return AlgebraicVector();
}

AlgebraicSymMatrix BaseTrackerRecHit::parametersError() const {
  obsolete();
  return AlgebraicSymMatrix();
}

AlgebraicMatrix BaseTrackerRecHit::projectionMatrix() const {
  obsolete();
  return AlgebraicMatrix();
}
