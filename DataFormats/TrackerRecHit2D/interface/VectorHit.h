#ifndef TrackerRecHit2D_VectorHit_h
#define TrackerRecHit2D_VectorHit_h

#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/VectorHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"

#include "TkCloner.h"

class VectorHit final : public BaseTrackerRecHit {
public:
  typedef OmniClusterRef::Phase2Cluster1DRef ClusterRef;

  VectorHit() : thePosition(), theDirection(), theCovMatrix() { setType(bad); }

  VectorHit(const VectorHit& vh);

  VectorHit(const GeomDet& idet,
            const LocalPoint& posInner,
            const LocalVector& dir,
            const AlgebraicSymMatrix& covMatrix,
            const double& Chi2,
            OmniClusterRef const& lower,
            OmniClusterRef const& upper);

  VectorHit(const GeomDet& idet,
            const VectorHit2D& vh2Dzx,
            const VectorHit2D& vh2Dzy,
            OmniClusterRef const& lower,
            OmniClusterRef const& upper);

  ~VectorHit() override;

   enum curvPhiSwitch { CURV = 0, PHI = 1}; 

  VectorHit* clone() const override { return new VectorHit(*this); }
  RecHitPointer cloneSH() const override { return std::make_shared<VectorHit>(*this); }

  bool sharesInput(const TrackingRecHit* other, SharedInputType what) const override;
  bool sharesClusters(VectorHit const& h1, VectorHit const& h2, SharedInputType what) const;

  // Parameters of the segment, for the track fit
  // For a 4D segment: (dx/dz,dy/dz,x,y)
  bool hasPositionAndError() const override {
    return true;
  };

  AlgebraicVector parameters() const override;
  void getKfComponents(KfComponentsHolder& holder) const override { getKfComponents4D(holder); }
  void getKfComponents4D(KfComponentsHolder& holder) const;

  // returning methods
  LocalPoint localPosition() const override { return thePosition; }
  virtual LocalVector localDirection() const { return theDirection; }
  AlgebraicSymMatrix parametersError() const override;
  LocalError localPositionError() const override;
  virtual LocalError localDirectionError() const;
  Global3DVector globalDirection() const;

  virtual double chi2() const { return theChi2; }
  int dimension() const override { return theDimension; }

  std::pair<double, double> curvatureORphi(curvPhiSwitch curvORPhi ) const;
  const float transverseMomentum(const MagneticField* magField);
  const float momentum(const MagneticField* magField);

  ClusterRef lowerCluster() const { return theLowerCluster.cluster_phase2OT(); }
  ClusterRef upperCluster() const { return theUpperCluster.cluster_phase2OT(); }
  OmniClusterRef const lowerClusterRef() const { return theLowerCluster; }
  OmniClusterRef const upperClusterRef() const { return theUpperCluster; }

  //FIXME::to update with a proper CPE maybe...
  Global3DPoint lowerGlobalPos() const;
  Global3DPoint upperGlobalPos() const;
  Global3DPoint phase2clusterGlobalPos(const PixelGeomDetUnit* geomDet, ClusterRef cluster) const;
  GlobalError lowerGlobalPosErr() const;
  GlobalError upperGlobalPosErr() const;
  GlobalError phase2clusterGlobalPosErr(const PixelGeomDetUnit* geomDet) const;

  bool isPhase2() const override { return true; }

  //FIXME: I have always two clusters in a VH
  OmniClusterRef const& firstClusterRef() const final { return theLowerCluster; }
  ClusterRef cluster() const { return theLowerCluster.cluster_phase2OT(); }

  //This method returns the delta in global coordinates
  Global3DVector globalDelta() const;
  float theta();

  /// The projection matrix relates the trajectory state parameters to the segment parameters().
  AlgebraicMatrix projectionMatrix() const override;

  // Access to component RecHits (if any)
  std::vector<const TrackingRecHit*> recHits() const override;
  std::vector<TrackingRecHit*> recHits() override;

private:
  // double dispatch
  VectorHit* clone_(TkCloner const& cloner, TrajectoryStateOnSurface const& tsos) const override {
    return cloner(*this, tsos).release();
  }
  RecHitPointer cloneSH_(TkCloner const& cloner, TrajectoryStateOnSurface const& tsos) const override {
    return cloner.makeShared(*this, tsos);
  }

  LocalPoint thePosition;
  LocalVector theDirection;

  // the covariance matrix, has the following meaning
  // mat[0][0]=var(dx/dz)
  // mat[1][1]=var(dy/dz)
  // mat[2][2]=var(x)
  // mat[3][3]=var(y)
  // mat[0][2]=cov(dx/dz,x)
  // mat[1][3]=cov(dy/dz,y)
  AlgebraicSymMatrix theCovMatrix;
  double theChi2;
  static constexpr int theDimension = 4;
  OmniClusterRef theLowerCluster;
  OmniClusterRef theUpperCluster;
};

inline bool operator<(const VectorHit& one, const VectorHit& other) {
  return one.chi2() < other.chi2();
}

std::ostream& operator<<(std::ostream& os, const VectorHit& vh);

typedef edmNew::DetSetVector<VectorHit> VectorHitCollection;
typedef VectorHitCollection VectorHitCollectionNew;

#endif
