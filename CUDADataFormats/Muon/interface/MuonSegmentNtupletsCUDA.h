#ifndef CUDADataFormats_Muon_interface_MuonSegmentNtuplets_h
#define CUDADataFormats_Muon_interface_MuonSegmentNtuplets_h


class MuonSegmentNtupletsCUDA {
public:

  float gx1[100];
  float gy1[100];
  float gz1[100];
  float gphi1[100];
  float gr1[100];
  int   layerID1[100];
  float gx2[100];
  float gy2[100];
  float gz2[100];
  float gphi2[100];
  float gr2[100];
  int   layerID2[100];
  float gx3[100];
  float gy3[100];
  float gz3[100];
  float gphi3[100];
  float gr3[100];
  int   layerID3[100];
  float gx4[100];
  float gy4[100];
  float gz4[100];
  float gphi4[100];
  float gr4[100];
  int   layerID4[100];

  int nNtuplets; 
  int segmentsInNtuplet[100]; 

private:

};

#endif  // CUDADataFormats_Muon_interface_MuonSegmentNtuplets_h
