#include "DataFormats/TrackerRecHit2D/interface/VectorHit.h"
#include "Geometry/CommonDetUnit/interface/StackGeomDet.h"

VectorHit::VectorHit(const VectorHit& vh)
    : BaseTrackerRecHit(*vh.det(), trackerHitRTTI::vector),
      thePosition(vh.localPosition()),
      theDirection(vh.localDirection()),
      theCovMatrix(vh.parametersError()),
      theChi2(vh.chi2()),
      theLowerCluster(vh.lowerClusterRef()),
      theUpperCluster(vh.upperClusterRef()) {}

VectorHit::VectorHit(const GeomDet& idet,
                     const LocalPoint& posLower,
                     const LocalVector& dir,
                     const AlgebraicSymMatrix& covMatrix,
                     const double& Chi2,
                     OmniClusterRef const& lower,
                     OmniClusterRef const& upper)
    : BaseTrackerRecHit(idet, trackerHitRTTI::vector),
      thePosition(posLower),
      theDirection(dir),
      theCovMatrix(covMatrix),
      theChi2(Chi2),
      theLowerCluster(lower),
      theUpperCluster(upper) {}

VectorHit::VectorHit(const GeomDet& idet,
                     const VectorHit2D& vh2Dzx,
                     const VectorHit2D& vh2Dzy,
                     OmniClusterRef const& lower,
                     OmniClusterRef const& upper)
    : BaseTrackerRecHit(idet, trackerHitRTTI::vector), theLowerCluster(lower), theUpperCluster(upper) {
  thePosition = LocalPoint(vh2Dzx.localPosition().x(), vh2Dzy.localPosition().x(), 0.);

  theDirection = LocalVector(vh2Dzx.localDirection().x(), vh2Dzy.localDirection().x(), 1.);

  //building the cov matrix 4x4 starting from the 2x2
  AlgebraicSymMatrix22 covMatZX = vh2Dzx.covMatrix();
  AlgebraicSymMatrix22 covMatZY = vh2Dzy.covMatrix();

  theCovMatrix = AlgebraicSymMatrix(theDimension);
  theCovMatrix[0][0] = covMatZX[0][0];  // var(dx/dz)
  theCovMatrix[1][1] = covMatZY[0][0];  // var(dy/dz)
  theCovMatrix[2][2] = covMatZX[1][1];  // var(x)
  theCovMatrix[3][3] = covMatZY[1][1];  // var(y)
  theCovMatrix[0][2] = covMatZX[0][1];  // cov(dx/dz,x)
  theCovMatrix[1][3] = covMatZY[0][1];  // cov(dy/dz,y)

  theChi2 = vh2Dzx.chi2() + vh2Dzy.chi2();
}

bool VectorHit::sharesInput(const TrackingRecHit* other, SharedInputType what) const {
  if (what == all && (geographicalId() != other->geographicalId()))
    return false;

  if (!sameDetModule(*other))
    return false;

  if (trackerHitRTTI::isVector(*other)) {
    const VectorHit* otherVh = static_cast<const VectorHit*>(other);
    return sharesClusters(*this, *otherVh, what);
  }

  auto const& otherClus = reinterpret_cast<const BaseTrackerRecHit*>(other)->firstClusterRef();
  return (otherClus == lowerClusterRef()) || (otherClus == upperClusterRef());
}

bool VectorHit::sharesClusters(VectorHit const& h1, VectorHit const& h2, SharedInputType what) const {
  bool lower = h1.lowerClusterRef() == h2.lowerClusterRef();
  bool upper = h1.upperClusterRef() == h2.upperClusterRef();

  return (what == TrackingRecHit::all) ? (lower && upper) : (upper || lower);
}

void VectorHit::getKfComponents4D(KfComponentsHolder& holder) const {
  
  AlgebraicVector4& pars = holder.params<theDimension>();
  pars[0] = theDirection.x();
  pars[1] = theDirection.y();
  pars[2] = thePosition.x();
  pars[3] = thePosition.y();

  AlgebraicSymMatrix44& errs = holder.errors<theDimension>();
  for (int i = 0; i < theDimension; i++) {
    for (int j = 0; j < theDimension; j++) {
      errs(i, j) = theCovMatrix[i][j];
    }
  }

  ProjectMatrix<double, 5, theDimension>& pf = holder.projFunc<theDimension>();
  pf.index[0] = 1;
  pf.index[1] = 2;
  pf.index[2] = 3;
  pf.index[3] = 4;

  holder.measuredParams<theDimension>() = AlgebraicVector4(&holder.tsosLocalParameters().At(1), 4);
  holder.measuredErrors<theDimension>() = holder.tsosLocalErrors().Sub<AlgebraicSymMatrix44>(1, 1);
}

VectorHit::~VectorHit() {}

AlgebraicVector VectorHit::parameters() const {
  // (dx/dz,dy/dz,x,y)
  AlgebraicVector result(theDimension);

  result[0] = theDirection.x();
  result[1] = theDirection.y();
  result[2] = thePosition.x();
  result[3] = thePosition.y();
  return result;
}

Global3DPoint VectorHit::lowerGlobalPos() const {
  const PixelGeomDetUnit* geomDetLower = dynamic_cast<const PixelGeomDetUnit*>(dynamic_cast<const StackGeomDet*>(det())->lowerDet());
  return phase2clusterGlobalPos(geomDetLower, lowerCluster());
}

Global3DPoint VectorHit::upperGlobalPos() const {
  const PixelGeomDetUnit* geomDetUpper = dynamic_cast<const PixelGeomDetUnit*>(dynamic_cast<const StackGeomDet*>(det())->upperDet());
  return phase2clusterGlobalPos(geomDetUpper, upperCluster());
}

Global3DPoint VectorHit::phase2clusterGlobalPos(const PixelGeomDetUnit* geomDet, ClusterRef cluster) const {
  const PixelTopology* topo = &geomDet->specificTopology();
  float ix = cluster->center();
  float iy = cluster->column() + 0.5;                    // halfway the column
  LocalPoint lp(topo->localX(ix), topo->localY(iy), 0);  // x, y, z
  Global3DPoint gp = geomDet->surface().toGlobal(lp);
  return gp;
}

GlobalError VectorHit::lowerGlobalPosErr() const {
  const PixelGeomDetUnit* geomDetLower = dynamic_cast<const PixelGeomDetUnit*>(dynamic_cast<const StackGeomDet*>(det())->lowerDet());
  return phase2clusterGlobalPosErr(geomDetLower);
}

GlobalError VectorHit::upperGlobalPosErr() const {
  const PixelGeomDetUnit* geomDetUpper = dynamic_cast<const PixelGeomDetUnit*>(dynamic_cast<const StackGeomDet*>(det())->upperDet());
  return phase2clusterGlobalPosErr(geomDetUpper);
}

GlobalError VectorHit::phase2clusterGlobalPosErr(const PixelGeomDetUnit* geomDet) const {
  const PixelTopology* topo = &geomDet->specificTopology();
  float pitchX = topo->pitch().first;
  float pitchY = topo->pitch().second;
  LocalError le(pow(pitchX, 2) / 12., 0, pow(pitchY, 2) / 12.);  // e2_xx, e2_xy, e2_yy
  GlobalError ge(ErrorFrameTransformer().transform(le, geomDet->surface()));
  return ge;
}

Global3DVector VectorHit::globalDelta() const {
  Local3DVector theLocalDelta =
      LocalVector(theDirection.x() * theDirection.z(), theDirection.y() * theDirection.z(), theDirection.z());
  Global3DVector g = det()->surface().toGlobal(theLocalDelta);
  return g;
}

Global3DVector VectorHit::globalDirection() const { return (det()->surface().toGlobal(localDirection())); }

std::pair<double, double> VectorHit::curvatureORphi(curvPhiSwitch  curvORphi) const {
  double curvature = -999.;
  double errorCurvature = -999.;
  double phi = -999.;

  //global pos and errors
  Global3DPoint gPositionLower = lowerGlobalPos();
  Global3DPoint gPositionUpper = upperGlobalPos();

  GlobalError gErrorLower = lowerGlobalPosErr();
  GlobalError gErrorUpper = upperGlobalPosErr();

  //insert lower and upper in the global sor
  if (gPositionLower.perp() > gPositionUpper.perp()) {
    gPositionLower = upperGlobalPos();
    gPositionUpper = lowerGlobalPos();
    gErrorLower = upperGlobalPosErr();
    gErrorUpper = lowerGlobalPosErr();
  }

  double lowerX = gPositionLower.x();
  double lowerY = gPositionLower.y();
  double upperX = gPositionUpper.x();
  double upperY = gPositionUpper.y();

  double h1 = lowerX * upperY - upperX * lowerY;

  //determine sign of curvature
  AlgebraicVector2 n1;
  n1[0] = -lowerY;
  n1[1] = lowerX;
  AlgebraicVector2 n2;
  n2[0] = upperX - lowerX;
  n2[1] = upperY - lowerY;

  double n3 = n1[0] * n2[0] + n1[1] * n2[1];
  double signCurv = -copysign(1.0, n3);
  double phi1 = atan2(upperY - lowerY, upperX - lowerX);

  if (h1 != 0) {
    double h2 = 2 * h1;
    double r12 = pow(lowerX, 2) + pow(lowerY, 2);
    double r22 = pow(upperX, 2) + pow(upperY, 2);
    double h3 = pow(lowerX - upperX, 2) + pow(lowerY - upperY, 2);
    double h4 = -pow(lowerX, 2) * upperX + lowerX * pow(upperX, 2) +
                lowerX * pow(upperY, 2) - upperX * pow(lowerY, 2);
    double h5 = pow(lowerX, 2) * upperY - pow(upperX, 2) * lowerY +
                pow(lowerY, 2) * upperY - lowerY * pow(upperY, 2);

    //radius of circle
    double invRho2 = (4. * h1 * h1) / (r12 * r22 * h3);
    curvature = sqrt(invRho2);

    //center of circle
    double xcentre = h5 / h2;
    double ycentre = h4 / h2;

    //to compute phi at the cluster points
    double xtg = lowerY - ycentre;
    double ytg = -(lowerX - xcentre);

    //to compute phi at the origin
    phi = atan2(ytg, xtg);

    AlgebraicROOTObject<theDimension, theDimension>::Matrix jacobian;
    for (int i = 0; i < theDimension; i++) {
      for (int j = 0; j < theDimension; j++) {
        jacobian[i][j] = 0.0;
      }
    }
    double denom1 = 1. / sqrt(r12 * r22 * h3);
    double denom2 = 1. / (pow(r12 * r22 * h3, 1.5));
    jacobian[0][0] = 1.0;  // dx1/dx1 dx1/dy1 dx2/dx1 dy2/dx1
    jacobian[1][1] = 1.0;  //dy1/dx1 dy1/dy1 dy2/dx1 dy2/dx1
    // dkappa/dx1
    jacobian[2][0] = (h1 * (2. * lowerX * r22 * h3 + (2. * lowerX - 2. * upperX) * r12 * r22)) * denom2 - (2. * upperY) * denom1;
    // dkappa/dy1
    jacobian[2][1] = (2. * upperX) * denom1 + (h1 * (2. * lowerY * r22 * h3 + r12 * r22 * (2. * lowerY - 2. * upperY))) * denom2;
    // dkappa/dx2
    jacobian[2][2] = (2. * lowerY) * denom1 + (h1 * (2. * upperX * r12 * h3 - 2. * (lowerX - upperX) * r12 * r22)) * denom2;
    // dkappa/dy2
    jacobian[2][3] = (h1 * (2. * upperY * r12 * h3 - r12 * r22 * 2. * (lowerY - upperY))) * denom2 - (2. * lowerX) * denom1;

    for (int i = 0; i < theDimension; i++) {
      jacobian[2][i] = -jacobian[2][i];
    }

    AlgebraicVector2 vecM;
    //to compute phi at the cluster points
    vecM[0] = (lowerY - ycentre) * invRho2;   // dphi/dxcentre
    vecM[1] = -(lowerX - xcentre) * invRho2;  // dphi/dycentre
    //to compute phi at the origin

    AlgebraicROOTObject<2, theDimension>::Matrix matrixK;
    // dxm/dx1
    matrixK[0][0] = (2. * lowerX * upperY) / h2 - (2. * upperY * h5) / pow(h2, 2);
    // dxm/dy1
    matrixK[0][1] = (2. * upperX * h5) / pow(h2, 2) - (pow(upperX, 2) + pow(upperY, 2) - 2. * lowerY * upperY) / h2;
    // dxm/dx2
    matrixK[0][2] =  (2. * lowerY * h5) / pow(h2, 2) - (2. * upperX * lowerY) / h2;
    // dxm/dy2
    matrixK[0][3] = (pow(lowerX, 2) + pow(lowerY, 2) - 2. * upperY * lowerY) / h2 - (2. * lowerX * h5) / pow(h2, 2);
    // dym/dx1
    matrixK[1][0] = (pow(upperX, 2) - 2. * lowerX * upperX + pow(upperY, 2)) / h2 - (2. * upperY * h4) / pow(h2, 2);
    // dym/dy1
    matrixK[1][1] = (2. * upperX * h4) / pow(h2, 2) - (2. * upperX * lowerY) / h2;
    // dym/dx2
    matrixK[1][2] = (2. * lowerY * h4) / pow(h2, 2) - (pow(lowerX, 2) - 2. * upperX * lowerX + pow(lowerY, 2)) / h2;
    // dym/dy2
    matrixK[1][3] = (2. * lowerX * upperY) / h2 - (2. * lowerX * h4) / pow(h2, 2);

    AlgebraicVector4 vecN = vecM * matrixK;
    jacobian[3][0] = vecN[0];  // dphi/(dx1,dy1,dx2,dy2)
    jacobian[3][1] = vecN[1];  // dphi/(dx1,dy1,dx2,dy2)
    jacobian[3][2] = vecN[2];  // dphi/(dx1,dy1,dx2,dy2)
    jacobian[3][3] = vecN[3];  // dphi/(dx1,dy1,dx2,dy2)

    //assign correct sign to the curvature errors
    if ((signCurv < 0 && curvature > 0) || (signCurv > 0 && curvature < 0)) {
      curvature = -curvature;
      for (int i = 0; i < theDimension; i++) {
        jacobian[2][i] = -jacobian[2][i];
      }
    }

    // bring phi in the same quadrant as phi1
    if (abs(phi - phi1) > M_PI / 2.) {
      phi = phi + M_PI;
      if (phi > M_PI)
        phi = phi - 2. * M_PI;
    }

    //computing the curvature error
    AlgebraicVector4 curvatureJacobian;
    for (int i = 0; i < theDimension; i++) {
      curvatureJacobian[i] = jacobian[2][i];
    }

    AlgebraicROOTObject<theDimension, theDimension>::Matrix gErrors;
    for (int i = 0; i < theDimension; i++) {
      for (int j = 0; j < theDimension; j++) {
        gErrors[i][j] = 0.0;
      }
    }

    gErrors[0][0] = gErrorLower.cxx();
    gErrors[0][1] = gErrorLower.cyx();
    gErrors[1][0] = gErrorLower.cyx();
    gErrors[1][1] = gErrorLower.cyy();
    gErrors[2][2] = gErrorUpper.cxx();
    gErrors[2][3] = gErrorUpper.cyx();
    gErrors[3][2] = gErrorUpper.cyx();
    gErrors[3][3] = gErrorUpper.cyy();

    AlgebraicVector4 temp = curvatureJacobian;
    temp = temp * gErrors;
    errorCurvature = temp[0] * curvatureJacobian[0] + temp[1] * curvatureJacobian[1] + temp[2] * curvatureJacobian[2] +
                     temp[3] * curvatureJacobian[3];

  } else {
    return std::make_pair(0.0, 0.0);
  }

  if (curvORphi == CURV)
    return std::make_pair(curvature, errorCurvature);
  else if (curvORphi == PHI)
    return std::make_pair(phi, 0.0);
  else
    return std::make_pair(0.0, 0.0);
}

const float VectorHit::transverseMomentum(const MagneticField* magField){
  GlobalPoint center(0.0, 0.0, 0.0);
  float magnT = magField->inTesla(center).mag();
  double rho = 1.f / curvatureORphi(CURV).first;
  //0.003 is because the curvature (rho) is in cm and not in m
  return (0.003f * magnT * rho);
}

const float VectorHit::momentum(const MagneticField* magField){ return transverseMomentum(magField) / (1. * sin(theta())); }

float VectorHit::theta() { return globalDirection().theta(); }

AlgebraicMatrix VectorHit::projectionMatrix() const {
  // obsolete (for what tracker is concerned...) interface
  static const AlgebraicMatrix the4DProjectionMatrix(theDimension, 5, 0);
  return the4DProjectionMatrix;
}

LocalError VectorHit::localPositionError() const {
  return LocalError(theCovMatrix[2][2], theCovMatrix[2][3], theCovMatrix[3][3]);
}

LocalError VectorHit::localDirectionError() const {
  return LocalError(theCovMatrix[0][0], theCovMatrix[0][1], theCovMatrix[1][1]);
}

AlgebraicSymMatrix VectorHit::parametersError() const {
  //think about a more efficient method
  AlgebraicSymMatrix result(4);
  for (int i = 0; i < theDimension; i++) {
    for (int j = 0; j < theDimension; j++) {
      result[i][j] = theCovMatrix[i][j];
    }
  }
  return result;
}

std::ostream& operator<<(std::ostream& os, const VectorHit& vh) {
  os << " VectorHit create in the DetId#: " << vh.geographicalId() << "\n"
     << " Vectorhit local position      : " << vh.localPosition() << "\n"
     << " Vectorhit local direction     : " << vh.localDirection() << "\n"
     << " Vectorhit global direction    : " << vh.globalDirection() << "\n"
     << " Lower cluster global position : " << vh.lowerGlobalPos() << "\n"
     << " Upper cluster global position : " << vh.upperGlobalPos();

  return os;
}

/// Access to component RecHits (if any)
std::vector<const TrackingRecHit*> VectorHit::recHits() const {
  std::vector<const TrackingRecHit*> pointersOfRecHits;
  return pointersOfRecHits;
}

/// Non-const access to component RecHits (if any)
std::vector<TrackingRecHit*> VectorHit::recHits() {
  std::vector<TrackingRecHit*> pointersOfRecHits;
  return pointersOfRecHits;
}
