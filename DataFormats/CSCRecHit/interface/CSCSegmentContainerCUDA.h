#ifndef CUDADataFormats_CSCRecHit_interface_CSCSegmentContainerCUDA_h
#define CUDADataFormats_CSCRecHit_interface_CSCSegmentContainerCUDA_h

class CSCSegmentContainerCUDA {
public:
  CSCSegmentContainerCUDA() = default;
  ~CSCSegmentContainerCUDA() = default;

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

#endif  // CUDADataFormats_CSCRecHit_interface_CSCSegmentContainerCUDA_h
