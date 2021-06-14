#ifndef CUDADataFormats_DTRecHit_interface_DTRecSegment4DContainerCUDA_h
#define CUDADataFormats_DTRecHit_interface_DTRecSegment4DContainerCUDA_h

class DTRecSegment4DContainerCUDA {
public:
  DTRecSegment4DContainerCUDA() = default;
  ~DTRecSegment4DContainerCUDA() = default;

  //local position of the segments
  float x[1000]; 
  float y[1000];
  //local direction of the segments
  float dxdz[1000]; 
  float dydz[100];
  //parameter uncertainties
  float sigmaX[1000]; 
  float sigmaY[1000];
  float sigmaDXDZ[1000]; 
  float sigmaDYDZ[1000];

  size_t nSegments; 
};

#endif  // CUDADataFormats_CSCRecHit_interface_DTRecSegment4DContainerCUDA_h
