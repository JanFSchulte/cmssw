#ifndef CUDADataFormats_Muon_interface_MuonSegmentPairs_h
#define CUDADataFormats_Muon_interface_MuonSegmentPairs_h


class MuonSegmentPairsCUDA {
public:
  float gx1[1000];
  float gy1[1000];
  float gz1[1000];
  float gphi1[1000];
  float gr1[1000];
  int   layerID1[1000];
  float gx2[1000];
  float gy2[1000];
  float gz2[1000];
  float gphi2[1000];
  float gr2[1000];
  int   layerID2[1000];

  int nPairs; 

private:

};

#endif  // CUDADataFormats_Muon_interface_MuonSegmentPairs_h
