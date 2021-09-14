#ifndef CUDADataFormats_Muon_interface_MuonSegmentPairsCUDA_h
#define CUDADataFormats_Muon_interface_MuonSegmentPairsCUDA_h


class MuonSegmentPairs {
public:
//  MuonSegmentPairs() = default;
//  MuonSegmentPairs(cudaStream_t stream) {};
//  ~MuonSegmentPairs() = default;

//  MuonSegmentPairs(const MuonSegmentPairs &) = delete;
//  MuonSegmentPairs &operator=(const MuonSegmentPairs &) = delete;
//  MuonSegmentPairs(MuonSegmentPairs &&) = default;
//  MuonSegmentPairs &operator=(MuonSegmentPairs &&) = default;



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

#endif  // CUDADataFormats_Muon_interface_MuonSegmentPairsCUDA_h
