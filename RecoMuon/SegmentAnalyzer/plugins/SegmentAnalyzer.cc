// -*- C++ -*-
//
// Package:    RecoMuon/SegmentAnalyzer
// Class:      SegmentAnalyzer
//
/**\class SegmentAnalyzer SegmentAnalyzer.cc RecoMuon/SegmentAnalyzer/plugins/SegmentAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jan-Frederik Schulte
//         Created:  Mon, 30 Aug 2021 19:21:06 GMT
//
//

// system include files
#include <memory>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"

// user include files



#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHitBuilder.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>

#include "CUDADataFormats/Muon/interface/MuonSegmentPairsCUDA.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentPairsHeterogeneous.h"

#include <DataFormats/MuonDetId/interface/DTChamberId.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/MuonDetId/interface/MuonSubdetId.h>
//
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class SegmentAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SegmentAnalyzer(const edm::ParameterSet&);
  ~SegmentAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsGetToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsGetToken_;
  edm::EDGetTokenT<reco::TrackCollection> l2MuonGetToken_;
  edm::EDGetTokenT<MuonSegmentPairsHeterogeneous> segmentPairGetToken_;

//  std::string theTrackerRecHitBuilderName;

  TH1F* nSegments;
  
  TH1F* nSegmentsPair1;
  TH1F* nSegmentsPair2;
  TH1F* nSegmentsPair3;
  TH1F* nSegmentsPair4;
  TH1F* nSegmentsPair5;
  TH1F* nSegmentsPair6;
  TH1F* nSegmentsPair7;
  TH1F* nSegmentsPair8;
  TH1F* nSegmentsPair9;
  TH1F* nSegmentsPair10;
  TH1F* nSegmentsPair11;
  TH1F* nSegmentsPair12;
  TH1F* nSegmentsPair13;
  TH1F* nSegmentsPair14;
  TH1F* nSegmentsPair15;
  TH1F* nSegmentsPair16;
  TH1F* nSegmentsPair17;
  TH1F* nSegmentsPair18;
  TH1F* nSegmentsPair19;
  TH1F* nSegmentsPair20;
  TH1F* nSegmentsPair21;
  TH1F* nSegmentsPair22;
  TH1F* nSegmentsPair23;

  TH1F* nSegmentsDT;
  TH1F* nSegmentsCSC;
  TH1F* zMB1;
  TH1F* zMB2;
  TH1F* zMB3;
  TH1F* zMB4;
  TH1F* zME1p;
  TH1F* zME2p;
  TH1F* zME3p;
  TH1F* zME4p;
  TH1F* zME1n;
  TH1F* zME2n;
  TH1F* zME3n;
  TH1F* zME4n;

  TH1F* rMB1;
  TH1F* rMB2;
  TH1F* rMB3;
  TH1F* rMB4;
  TH1F* rME1p;
  TH1F* rME2p;
  TH1F* rME3p;
  TH1F* rME4p;
  TH1F* rME1n;
  TH1F* rME2n;
  TH1F* rME3n;
  TH1F* rME4n;

  TH2F* rzMB1;
  TH2F* rzMB2;
  TH2F* rzMB3;
  TH2F* rzMB4;
  TH2F* rzME1p;
  TH2F* rzME2p;
  TH2F* rzME3p;
  TH2F* rzME4p;
  TH2F* rzME1n;
  TH2F* rzME2n;
  TH2F* rzME3n;
  TH2F* rzME4n;

  TH1F* distPair1;
  TH1F* distPair2;
  TH1F* distPair3;
  TH1F* distPair4;
  TH1F* distPair5;
  TH1F* distPair6;
  TH1F* distPair7;
  TH1F* distPair8;
  TH1F* distPair9;
  TH1F* distPair10;
  TH1F* distPair11;
  TH1F* distPair12;
  TH1F* distPair13;
  TH1F* distPair14;
  TH1F* distPair15;
  TH1F* distPair16;
  TH1F* distPair17;
  TH1F* distPair18;
  TH1F* distPair19;
  TH1F* distPair20;
  TH1F* distPair21;
  TH1F* distPair22;
  TH1F* distPair23;

  TH1F* dDirPair1;
  TH1F* dDirPair2;
  TH1F* dDirPair3;
  TH1F* dDirPair4;
  TH1F* dDirPair5;
  TH1F* dDirPair6;
  TH1F* dDirPair7;
  TH1F* dDirPair8;
  TH1F* dDirPair9;
  TH1F* dDirPair10;
  TH1F* dDirPair11;
  TH1F* dDirPair12;
  TH1F* dDirPair13;
  TH1F* dDirPair14;
  TH1F* dDirPair15;
  TH1F* dDirPair16;
  TH1F* dDirPair17;
  TH1F* dDirPair18;
  TH1F* dDirPair19;
  TH1F* dDirPair20;
  TH1F* dDirPair21;
  TH1F* dDirPair22;
  TH1F* dDirPair23;

  TH1F* drPair1;
  TH1F* drPair2;
  TH1F* drPair3;
  TH1F* drPair4;
  TH1F* drPair5;
  TH1F* drPair6;
  TH1F* drPair7;
  TH1F* drPair8;
  TH1F* drPair9;
  TH1F* drPair10;
  TH1F* drPair11;
  TH1F* drPair12;
  TH1F* drPair13;
  TH1F* drPair14;
  TH1F* drPair15;
  TH1F* drPair16;
  TH1F* drPair17;
  TH1F* drPair18;
  TH1F* drPair19;
  TH1F* drPair20;
  TH1F* drPair21;
  TH1F* drPair22;
  TH1F* drPair23;

  TH1F* dzPair1;
  TH1F* dzPair2;
  TH1F* dzPair3;
  TH1F* dzPair4;
  TH1F* dzPair5;
  TH1F* dzPair6;
  TH1F* dzPair7;
  TH1F* dzPair8;
  TH1F* dzPair9;
  TH1F* dzPair10;
  TH1F* dzPair11;
  TH1F* dzPair12;
  TH1F* dzPair13;
  TH1F* dzPair14;
  TH1F* dzPair15;
  TH1F* dzPair16;
  TH1F* dzPair17;
  TH1F* dzPair18;
  TH1F* dzPair19;
  TH1F* dzPair20;
  TH1F* dzPair21;
  TH1F* dzPair22;
  TH1F* dzPair23;

  TH1F* dPhiPair1;
  TH1F* dPhiPair2;
  TH1F* dPhiPair3;
  TH1F* dPhiPair4;
  TH1F* dPhiPair5;
  TH1F* dPhiPair6;
  TH1F* dPhiPair7;
  TH1F* dPhiPair8;
  TH1F* dPhiPair9;
  TH1F* dPhiPair10;
  TH1F* dPhiPair11;
  TH1F* dPhiPair12;
  TH1F* dPhiPair13;
  TH1F* dPhiPair14;
  TH1F* dPhiPair15;
  TH1F* dPhiPair16;
  TH1F* dPhiPair17;
  TH1F* dPhiPair18;
  TH1F* dPhiPair19;
  TH1F* dPhiPair20;
  TH1F* dPhiPair21;
  TH1F* dPhiPair22;
  TH1F* dPhiPair23;

  TH1F* pTPair1;
  TH1F* pTPair2;
  TH1F* pTPair3;
  TH1F* pTPair4;
  TH1F* pTPair5;
  TH1F* pTPair6;
  TH1F* pTPair7;
  TH1F* pTPair8;
  TH1F* pTPair9;
  TH1F* pTPair10;
  TH1F* pTPair11;
  TH1F* pTPair12;
  TH1F* pTPair13;
  TH1F* pTPair14;
  TH1F* pTPair15;
  TH1F* pTPair16;
  TH1F* pTPair17;
  TH1F* pTPair18;
  TH1F* pTPair19;
  TH1F* pTPair20;
  TH1F* pTPair21;
  TH1F* pTPair22;
  TH1F* pTPair23;

  TH1F* zSMB1;
  TH1F* zSMB2;
  TH1F* zSMB3;
  TH1F* zSMB4;
  TH1F* zSME1p;
  TH1F* zSME2p;
  TH1F* zSME3p;
  TH1F* zSME4p;
  TH1F* zSME1n;
  TH1F* zSME2n;
  TH1F* zSME3n;
  TH1F* zSME4n;

  TH1F* rSMB1;
  TH1F* rSMB2;
  TH1F* rSMB3;
  TH1F* rSMB4;
  TH1F* rSME1p;
  TH1F* rSME2p;
  TH1F* rSME3p;
  TH1F* rSME4p;
  TH1F* rSME1n;
  TH1F* rSME2n;
  TH1F* rSME3n;
  TH1F* rSME4n;

  TH1F* nSegmentsPairS1;
  TH1F* nSegmentsPairS2;
  TH1F* nSegmentsPairS3;
  TH1F* nSegmentsPairS4;
  TH1F* nSegmentsPairS5;
  TH1F* nSegmentsPairS6;
  TH1F* nSegmentsPairS7;
  TH1F* nSegmentsPairS8;
  TH1F* nSegmentsPairS9;
  TH1F* nSegmentsPairS10;
  TH1F* nSegmentsPairS11;
  TH1F* nSegmentsPairS12;
  TH1F* nSegmentsPairS13;
  TH1F* nSegmentsPairS14;
  TH1F* nSegmentsPairS15;
  TH1F* nSegmentsPairS16;
  TH1F* nSegmentsPairS17;
  TH1F* nSegmentsPairS18;
  TH1F* nSegmentsPairS19;
  TH1F* nSegmentsPairS20;
  TH1F* nSegmentsPairS21;
  TH1F* nSegmentsPairS22;
  TH1F* nSegmentsPairS23;

  TH1F* nSegmentsPairG1;
  TH1F* nSegmentsPairG2;
  TH1F* nSegmentsPairG3;
  TH1F* nSegmentsPairG4;
  TH1F* nSegmentsPairG5;
  TH1F* nSegmentsPairG6;
  TH1F* nSegmentsPairG7;
  TH1F* nSegmentsPairG8;
  TH1F* nSegmentsPairG9;
  TH1F* nSegmentsPairG10;
  TH1F* nSegmentsPairG11;
  TH1F* nSegmentsPairG12;
  TH1F* nSegmentsPairG13;
  TH1F* nSegmentsPairG14;
  TH1F* nSegmentsPairG15;
  TH1F* nSegmentsPairG16;
  TH1F* nSegmentsPairG17;
  TH1F* nSegmentsPairG18;
  TH1F* nSegmentsPairG19;
  TH1F* nSegmentsPairG20;
  TH1F* nSegmentsPairG21;
  TH1F* nSegmentsPairG22;
  TH1F* nSegmentsPairG23;


  TH2F* rzSMB1;
  TH2F* rzSMB2;
  TH2F* rzSMB3;
  TH2F* rzSMB4;
  TH2F* rzSME1p;
  TH2F* rzSME2p;
  TH2F* rzSME3p;
  TH2F* rzSME4p;
  TH2F* rzSME1n;
  TH2F* rzSME2n;
  TH2F* rzSME3n;
  TH2F* rzSME4n;

  TH1F* distPairS1;
  TH1F* distPairS2;
  TH1F* distPairS3;
  TH1F* distPairS4;
  TH1F* distPairS5;
  TH1F* distPairS6;
  TH1F* distPairS7;
  TH1F* distPairS8;
  TH1F* distPairS9;
  TH1F* distPairS10;
  TH1F* distPairS11;
  TH1F* distPairS12;
  TH1F* distPairS13;
  TH1F* distPairS14;
  TH1F* distPairS15;
  TH1F* distPairS16;
  TH1F* distPairS17;
  TH1F* distPairS18;
  TH1F* distPairS19;
  TH1F* distPairS20;
  TH1F* distPairS21;
  TH1F* distPairS22;
  TH1F* distPairS23;

  TH1F* drPairS1;
  TH1F* drPairS2;
  TH1F* drPairS3;
  TH1F* drPairS4;
  TH1F* drPairS5;
  TH1F* drPairS6;
  TH1F* drPairS7;
  TH1F* drPairS8;
  TH1F* drPairS9;
  TH1F* drPairS10;
  TH1F* drPairS11;
  TH1F* drPairS12;
  TH1F* drPairS13;
  TH1F* drPairS14;
  TH1F* drPairS15;
  TH1F* drPairS16;
  TH1F* drPairS17;
  TH1F* drPairS18;
  TH1F* drPairS19;
  TH1F* drPairS20;
  TH1F* drPairS21;
  TH1F* drPairS22;
  TH1F* drPairS23;

  TH1F* dzPairS1;
  TH1F* dzPairS2;
  TH1F* dzPairS3;
  TH1F* dzPairS4;
  TH1F* dzPairS5;
  TH1F* dzPairS6;
  TH1F* dzPairS7;
  TH1F* dzPairS8;
  TH1F* dzPairS9;
  TH1F* dzPairS10;
  TH1F* dzPairS11;
  TH1F* dzPairS12;
  TH1F* dzPairS13;
  TH1F* dzPairS14;
  TH1F* dzPairS15;
  TH1F* dzPairS16;
  TH1F* dzPairS17;
  TH1F* dzPairS18;
  TH1F* dzPairS19;
  TH1F* dzPairS20;
  TH1F* dzPairS21;
  TH1F* dzPairS22;
  TH1F* dzPairS23;

  TH1F* dPhiPairS1;
  TH1F* dPhiPairS2;
  TH1F* dPhiPairS3;
  TH1F* dPhiPairS4;
  TH1F* dPhiPairS5;
  TH1F* dPhiPairS6;
  TH1F* dPhiPairS7;
  TH1F* dPhiPairS8;
  TH1F* dPhiPairS9;
  TH1F* dPhiPairS10;
  TH1F* dPhiPairS11;
  TH1F* dPhiPairS12;
  TH1F* dPhiPairS13;
  TH1F* dPhiPairS14;
  TH1F* dPhiPairS15;
  TH1F* dPhiPairS16;
  TH1F* dPhiPairS17;
  TH1F* dPhiPairS18;
  TH1F* dPhiPairS19;
  TH1F* dPhiPairS20;
  TH1F* dPhiPairS21;
  TH1F* dPhiPairS22;
  TH1F* dPhiPairS23;


  TH1F* pTPairS1;
  TH1F* pTPairS2;
  TH1F* pTPairS3;
  TH1F* pTPairS4;
  TH1F* pTPairS5;
  TH1F* pTPairS6;
  TH1F* pTPairS7;
  TH1F* pTPairS8;
  TH1F* pTPairS9;
  TH1F* pTPairS10;
  TH1F* pTPairS11;
  TH1F* pTPairS12;
  TH1F* pTPairS13;
  TH1F* pTPairS14;
  TH1F* pTPairS15;
  TH1F* pTPairS16;
  TH1F* pTPairS17;
  TH1F* pTPairS18;
  TH1F* pTPairS19;
  TH1F* pTPairS20;
  TH1F* pTPairS21;
  TH1F* pTPairS22;
  TH1F* pTPairS23;

  TH1F* dDirPairS1;
  TH1F* dDirPairS2;
  TH1F* dDirPairS3;
  TH1F* dDirPairS4;
  TH1F* dDirPairS5;
  TH1F* dDirPairS6;
  TH1F* dDirPairS7;
  TH1F* dDirPairS8;
  TH1F* dDirPairS9;
  TH1F* dDirPairS10;
  TH1F* dDirPairS11;
  TH1F* dDirPairS12;
  TH1F* dDirPairS13;
  TH1F* dDirPairS14;
  TH1F* dDirPairS15;
  TH1F* dDirPairS16;
  TH1F* dDirPairS17;
  TH1F* dDirPairS18;
  TH1F* dDirPairS19;
  TH1F* dDirPairS20;
  TH1F* dDirPairS21;
  TH1F* dDirPairS22;
  TH1F* dDirPairS23;

  TH1F* distPairG1;
  TH1F* distPairG2;
  TH1F* distPairG3;
  TH1F* distPairG4;
  TH1F* distPairG5;
  TH1F* distPairG6;
  TH1F* distPairG7;
  TH1F* distPairG8;
  TH1F* distPairG9;
  TH1F* distPairG10;
  TH1F* distPairG11;
  TH1F* distPairG12;
  TH1F* distPairG13;
  TH1F* distPairG14;
  TH1F* distPairG15;
  TH1F* distPairG16;
  TH1F* distPairG17;
  TH1F* distPairG18;
  TH1F* distPairG19;
  TH1F* distPairG20;
  TH1F* distPairG21;
  TH1F* distPairG22;
  TH1F* distPairG23;

  TH1F* drPairG1;
  TH1F* drPairG2;
  TH1F* drPairG3;
  TH1F* drPairG4;
  TH1F* drPairG5;
  TH1F* drPairG6;
  TH1F* drPairG7;
  TH1F* drPairG8;
  TH1F* drPairG9;
  TH1F* drPairG10;
  TH1F* drPairG11;
  TH1F* drPairG12;
  TH1F* drPairG13;
  TH1F* drPairG14;
  TH1F* drPairG15;
  TH1F* drPairG16;
  TH1F* drPairG17;
  TH1F* drPairG18;
  TH1F* drPairG19;
  TH1F* drPairG20;
  TH1F* drPairG21;
  TH1F* drPairG22;
  TH1F* drPairG23;

  TH1F* dzPairG1;
  TH1F* dzPairG2;
  TH1F* dzPairG3;
  TH1F* dzPairG4;
  TH1F* dzPairG5;
  TH1F* dzPairG6;
  TH1F* dzPairG7;
  TH1F* dzPairG8;
  TH1F* dzPairG9;
  TH1F* dzPairG10;
  TH1F* dzPairG11;
  TH1F* dzPairG12;
  TH1F* dzPairG13;
  TH1F* dzPairG14;
  TH1F* dzPairG15;
  TH1F* dzPairG16;
  TH1F* dzPairG17;
  TH1F* dzPairG18;
  TH1F* dzPairG19;
  TH1F* dzPairG20;
  TH1F* dzPairG21;
  TH1F* dzPairG22;
  TH1F* dzPairG23;

  TH1F* dPhiPairG1;
  TH1F* dPhiPairG2;
  TH1F* dPhiPairG3;
  TH1F* dPhiPairG4;
  TH1F* dPhiPairG5;
  TH1F* dPhiPairG6;
  TH1F* dPhiPairG7;
  TH1F* dPhiPairG8;
  TH1F* dPhiPairG9;
  TH1F* dPhiPairG10;
  TH1F* dPhiPairG11;
  TH1F* dPhiPairG12;
  TH1F* dPhiPairG13;
  TH1F* dPhiPairG14;
  TH1F* dPhiPairG15;
  TH1F* dPhiPairG16;
  TH1F* dPhiPairG17;
  TH1F* dPhiPairG18;
  TH1F* dPhiPairG19;
  TH1F* dPhiPairG20;
  TH1F* dPhiPairG21;
  TH1F* dPhiPairG22;
  TH1F* dPhiPairG23;

  TH1F* dDirPairG1;
  TH1F* dDirPairG2;
  TH1F* dDirPairG3;
  TH1F* dDirPairG4;
  TH1F* dDirPairG5;
  TH1F* dDirPairG6;
  TH1F* dDirPairG7;
  TH1F* dDirPairG8;
  TH1F* dDirPairG9;
  TH1F* dDirPairG10;
  TH1F* dDirPairG11;
  TH1F* dDirPairG12;
  TH1F* dDirPairG13;
  TH1F* dDirPairG14;
  TH1F* dDirPairG15;
  TH1F* dDirPairG16;
  TH1F* dDirPairG17;
  TH1F* dDirPairG18;
  TH1F* dDirPairG19;
  TH1F* dDirPairG20;
  TH1F* dDirPairG21;
  TH1F* dDirPairG22;
  TH1F* dDirPairG23;


  TH1F* purityNum;
  TH1F* purityDenom;

  TH1F* effNum;
  TH1F* effDenom;
 

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

SegmentAnalyzer::SegmentAnalyzer(const edm::ParameterSet& iConfig)
    : dtSegmentsGetToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("srcDT"))),
      cscSegmentsGetToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("srcCSC"))),
      l2MuonGetToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("srcL2"))),
      segmentPairGetToken_(consumes<MuonSegmentPairsHeterogeneous>(iConfig.getParameter<edm::InputTag>("srcSegmentPairs")))
//      theTrackerRecHitBuilderName(iConfig.getParameter<std::string>("TrackerRecHitBuilder"))
{

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);
  nSegments    = fs->make<TH1F>("nSegments", "segments / event",  300, 0,  300);
  nSegmentsDT  = fs->make<TH1F>("nSegmentsDT", "DT segments / event",  300, 0,  300);
  nSegmentsCSC = fs->make<TH1F>("nSegmentsCSC", "CSC segments / event",  300, 0,  300);

  zMB1 = fs->make<TH1F>("zMB1", "z of segments in MB1",  3000, -1500,  1500);
  zMB2 = fs->make<TH1F>("zMB2", "z of segments in MB2",  3000, -1500,  1500);
  zMB3 = fs->make<TH1F>("zMB3", "z of segments in MB3",  3000, -1500,  1500);
  zMB4 = fs->make<TH1F>("zMB4", "z of segments in MB4",  3000, -1500,  1500);
  zME1p = fs->make<TH1F>("zME1p", "z of segments in ME1p",  3000, -1500,  1500);
  zME2p = fs->make<TH1F>("zME2p", "z of segments in ME2p",  3000, -1500,  1500);
  zME3p = fs->make<TH1F>("zME3p", "z of segments in ME3p",  3000, -1500,  1500);
  zME4p = fs->make<TH1F>("zME4p", "z of segments in ME4p",  3000, -1500,  1500);
  zME1n = fs->make<TH1F>("zME1n", "z of segments in ME1n",  3000, -1500,  1500);
  zME2n = fs->make<TH1F>("zME2n", "z of segments in ME2n",  3000, -1500,  1500);
  zME3n = fs->make<TH1F>("zME3n", "z of segments in ME3n",  3000, -1500,  1500);
  zME4n = fs->make<TH1F>("zME4n", "z of segments in ME4n",  3000, -1500,  1500);

  rMB1 = fs->make<TH1F>("rMB1", "r of segments in MB1",  3000, -1500,  1500);
  rMB2 = fs->make<TH1F>("rMB2", "r of segments in MB2",  3000, -1500,  1500);
  rMB3 = fs->make<TH1F>("rMB3", "r of segments in MB3",  3000, -1500,  1500);
  rMB4 = fs->make<TH1F>("rMB4", "r of segments in MB4",  3000, -1500,  1500);
  rME1p = fs->make<TH1F>("rME1p", "r of segments in ME1p",  3000, -1500,  1500);
  rME2p = fs->make<TH1F>("rME2p", "r of segments in ME2p",  3000, -1500,  1500);
  rME3p = fs->make<TH1F>("rME3p", "r of segments in ME3p",  3000, -1500,  1500);
  rME4p = fs->make<TH1F>("rME4p", "r of segments in ME4p",  3000, -1500,  1500);
  rME1n = fs->make<TH1F>("rME1n", "r of segments in ME1n",  3000, -1500,  1500);
  rME2n = fs->make<TH1F>("rME2n", "r of segments in ME2n",  3000, -1500,  1500);
  rME3n = fs->make<TH1F>("rME3n", "r of segments in ME3n",  3000, -1500,  1500);
  rME4n = fs->make<TH1F>("rME4n", "r of segments in ME4n",  3000, -1500,  1500);

  rzMB1 = fs->make<TH2F>("rzMB1", "rz of segments in MB1",  600, -1500,  1500, 300, 0, 1500);
  rzMB2 = fs->make<TH2F>("rzMB2", "rz of segments in MB2",  600, -1500,  1500, 300, 0, 1500);
  rzMB3 = fs->make<TH2F>("rzMB3", "rz of segments in MB3",  600, -1500,  1500, 300, 0, 1500);
  rzMB4 = fs->make<TH2F>("rzMB4", "rz of segments in MB4",  600, -1500,  1500, 300, 0, 1500);
  rzME1p = fs->make<TH2F>("rzME1p", "rz of segments in ME1p",  600, -1500,  1500, 300, 0, 1500);
  rzME2p = fs->make<TH2F>("rzME2p", "rz of segments in ME2p",  600, -1500,  1500, 300, 0, 1500);
  rzME3p = fs->make<TH2F>("rzME3p", "rz of segments in ME3p",  600, -1500,  1500, 300, 0, 1500);
  rzME4p = fs->make<TH2F>("rzME4p", "rz of segments in ME4p",  600, -1500,  1500, 300, 0, 1500);
  rzME1n = fs->make<TH2F>("rzME1n", "rz of segments in ME1n",  600, -1500,  1500, 300, 0, 1500);
  rzME2n = fs->make<TH2F>("rzME2n", "rz of segments in ME2n",  600, -1500,  1500, 300, 0, 1500);
  rzME3n = fs->make<TH2F>("rzME3n", "rz of segments in ME3n",  600, -1500,  1500, 300, 0, 1500);
  rzME4n = fs->make<TH2F>("rzME4n", "rz of segments in ME4n",  600, -1500,  1500, 300, 0, 1500);




  nSegmentsPair1    = fs->make<TH1F>("nSegmentsPair1", "segments / event",  300, 0,  300);
  nSegmentsPair2    = fs->make<TH1F>("nSegmentsPair2", "segments / event",  300, 0,  300);
  nSegmentsPair3    = fs->make<TH1F>("nSegmentsPair3", "segments / event",  300, 0,  300);
  nSegmentsPair4    = fs->make<TH1F>("nSegmentsPair4", "segments / event",  300, 0,  300);
  nSegmentsPair5    = fs->make<TH1F>("nSegmentsPair5", "segments / event",  300, 0,  300);
  nSegmentsPair6    = fs->make<TH1F>("nSegmentsPair6", "segments / event",  300, 0,  300);
  nSegmentsPair7    = fs->make<TH1F>("nSegmentsPair7", "segments / event",  300, 0,  300);
  nSegmentsPair8    = fs->make<TH1F>("nSegmentsPair8", "segments / event",  300, 0,  300);
  nSegmentsPair9    = fs->make<TH1F>("nSegmentsPair9", "segments / event",  300, 0,  300);
  nSegmentsPair10    = fs->make<TH1F>("nSegmentsPair10", "segments / event",  300, 0,  300);
  nSegmentsPair11    = fs->make<TH1F>("nSegmentsPair11", "segments / event",  300, 0,  300);
  nSegmentsPair12    = fs->make<TH1F>("nSegmentsPair12", "segments / event",  300, 0,  300);
  nSegmentsPair13    = fs->make<TH1F>("nSegmentsPair13", "segments / event",  300, 0,  300);
  nSegmentsPair14    = fs->make<TH1F>("nSegmentsPair14", "segments / event",  300, 0,  300);
  nSegmentsPair15    = fs->make<TH1F>("nSegmentsPair15", "segments / event",  300, 0,  300);
  nSegmentsPair16    = fs->make<TH1F>("nSegmentsPair16", "segments / event",  300, 0,  300);
  nSegmentsPair17    = fs->make<TH1F>("nSegmentsPair17", "segments / event",  300, 0,  300);
  nSegmentsPair18    = fs->make<TH1F>("nSegmentsPair18", "segments / event",  300, 0,  300);
  nSegmentsPair19    = fs->make<TH1F>("nSegmentsPair19", "segments / event",  300, 0,  300);
  nSegmentsPair20    = fs->make<TH1F>("nSegmentsPair20", "segments / event",  300, 0,  300);
  nSegmentsPair21    = fs->make<TH1F>("nSegmentsPair21", "segments / event",  300, 0,  300);
  nSegmentsPair22    = fs->make<TH1F>("nSegmentsPair22", "segments / event",  300, 0,  300);
  nSegmentsPair23    = fs->make<TH1F>("nSegmentsPair23", "segments / event",  300, 0,  300);

  distPair1 = fs->make<TH1F>("distPair1", "distPair1",  3000, -1500,  1500);
  distPair2 = fs->make<TH1F>("distPair2", "distPair2",  3000, -1500,  1500);
  distPair3 = fs->make<TH1F>("distPair3", "distPair3",  3000, -1500,  1500);
  distPair4 = fs->make<TH1F>("distPair4", "distPair4",  3000, -1500,  1500);
  distPair5 = fs->make<TH1F>("distPair5", "distPair5",  3000, -1500,  1500);
  distPair6 = fs->make<TH1F>("distPair6", "distPair6",  3000, -1500,  1500);
  distPair7 = fs->make<TH1F>("distPair7", "distPair7",  3000, -1500,  1500);
  distPair8 = fs->make<TH1F>("distPair8", "distPair8",  3000, -1500,  1500);
  distPair9 = fs->make<TH1F>("distPair9", "distPair9",  3000, -1500,  1500);
  distPair10 = fs->make<TH1F>("distPair10", "distPair10",  3000, -1500,  1500);
  distPair11 = fs->make<TH1F>("distPair11", "distPair11",  3000, -1500,  1500);
  distPair12 = fs->make<TH1F>("distPair12", "distPair12",  3000, -1500,  1500);
  distPair13 = fs->make<TH1F>("distPair13", "distPair13",  3000, -1500,  1500);
  distPair14 = fs->make<TH1F>("distPair14", "distPair14",  3000, -1500,  1500);
  distPair15 = fs->make<TH1F>("distPair15", "distPair15",  3000, -1500,  1500);
  distPair16 = fs->make<TH1F>("distPair16", "distPair16",  3000, -1500,  1500);
  distPair17 = fs->make<TH1F>("distPair17", "distPair17",  3000, -1500,  1500);
  distPair18 = fs->make<TH1F>("distPair18", "distPair18",  3000, -1500,  1500);
  distPair19 = fs->make<TH1F>("distPair19", "distPair19",  3000, -1500,  1500);
  distPair20 = fs->make<TH1F>("distPair20", "distPair20",  3000, -1500,  1500);
  distPair21 = fs->make<TH1F>("distPair21", "distPair21",  3000, -1500,  1500);
  distPair22 = fs->make<TH1F>("distPair22", "distPair22",  3000, -1500,  1500);
  distPair23 = fs->make<TH1F>("distPair23", "distPair23",  3000, -1500,  1500);

  dDirPair1 = fs->make<TH1F>("dDirPair1", "dDirPair1",  1000, -1,  1);
  dDirPair2 = fs->make<TH1F>("dDirPair2", "dDirPair2",  1000, -1,  1);
  dDirPair3 = fs->make<TH1F>("dDirPair3", "dDirPair3",  1000, -1,  1);
  dDirPair4 = fs->make<TH1F>("dDirPair4", "dDirPair4",  1000, -1,  1);
  dDirPair5 = fs->make<TH1F>("dDirPair5", "dDirPair5",  1000, -1,  1);
  dDirPair6 = fs->make<TH1F>("dDirPair6", "dDirPair6",  1000, -1,  1);
  dDirPair7 = fs->make<TH1F>("dDirPair7", "dDirPair7",  1000, -1,  1);
  dDirPair8 = fs->make<TH1F>("dDirPair8", "dDirPair8",  1000, -1,  1);
  dDirPair9 = fs->make<TH1F>("dDirPair9", "dDirPair9",  1000, -1,  1);
  dDirPair10 = fs->make<TH1F>("dDirPair10", "dDirPair10",  1000, -1,  1);
  dDirPair11 = fs->make<TH1F>("dDirPair11", "dDirPair11",  1000, -1,  1);
  dDirPair12 = fs->make<TH1F>("dDirPair12", "dDirPair12",  1000, -1,  1);
  dDirPair13 = fs->make<TH1F>("dDirPair13", "dDirPair13",  1000, -1,  1);
  dDirPair14 = fs->make<TH1F>("dDirPair14", "dDirPair14",  1000, -1,  1);
  dDirPair15 = fs->make<TH1F>("dDirPair15", "dDirPair15",  1000, -1,  1);
  dDirPair16 = fs->make<TH1F>("dDirPair16", "dDirPair16",  1000, -1,  1);
  dDirPair17 = fs->make<TH1F>("dDirPair17", "dDirPair17",  1000, -1,  1);
  dDirPair18 = fs->make<TH1F>("dDirPair18", "dDirPair18",  1000, -1,  1);
  dDirPair19 = fs->make<TH1F>("dDirPair19", "dDirPair19",  1000, -1,  1);
  dDirPair20 = fs->make<TH1F>("dDirPair20", "dDirPair20",  1000, -1,  1);
  dDirPair21 = fs->make<TH1F>("dDirPair21", "dDirPair21",  1000, -1,  1);
  dDirPair22 = fs->make<TH1F>("dDirPair22", "dDirPair22",  1000, -1,  1);
  dDirPair23 = fs->make<TH1F>("dDirPair23", "dDirPair23",  1000, -1,  1);

  pTPair1 = fs->make<TH1F>("pTPair1", "pTPair1",  3000, 0,  1500);
  pTPair2 = fs->make<TH1F>("pTPair2", "pTPair2",  3000, 0,  1500);
  pTPair3 = fs->make<TH1F>("pTPair3", "pTPair3",  3000, 0,  1500);
  pTPair4 = fs->make<TH1F>("pTPair4", "pTPair4",  3000, 0,  1500);
  pTPair5 = fs->make<TH1F>("pTPair5", "pTPair5",  3000, 0,  1500);
  pTPair6 = fs->make<TH1F>("pTPair6", "pTPair6",  3000, 0,  1500);
  pTPair7 = fs->make<TH1F>("pTPair7", "pTPair7",  3000, 0,  1500);
  pTPair8 = fs->make<TH1F>("pTPair8", "pTPair8",  3000, 0,  1500);
  pTPair9 = fs->make<TH1F>("pTPair9", "pTPair9",  3000, 0,  1500);
  pTPair10 = fs->make<TH1F>("pTPair10", "pTPair10",  3000, 0,  1500);
  pTPair11 = fs->make<TH1F>("pTPair11", "pTPair11",  3000, 0,  1500);
  pTPair12 = fs->make<TH1F>("pTPair12", "pTPair12",  3000, 0,  1500);
  pTPair13 = fs->make<TH1F>("pTPair13", "pTPair13",  3000, 0,  1500);
  pTPair14 = fs->make<TH1F>("pTPair14", "pTPair14",  3000, 0,  1500);
  pTPair15 = fs->make<TH1F>("pTPair15", "pTPair15",  3000, 0,  1500);
  pTPair16 = fs->make<TH1F>("pTPair16", "pTPair16",  3000, 0,  1500);
  pTPair17 = fs->make<TH1F>("pTPair17", "pTPair17",  3000, 0,  1500);
  pTPair18 = fs->make<TH1F>("pTPair18", "pTPair18",  3000, 0,  1500);
  pTPair19 = fs->make<TH1F>("pTPair19", "pTPair19",  3000, 0,  1500);
  pTPair20 = fs->make<TH1F>("pTPair20", "pTPair20",  3000, 0,  1500);
  pTPair21 = fs->make<TH1F>("pTPair21", "pTPair21",  3000, 0,  1500);
  pTPair22 = fs->make<TH1F>("pTPair22", "pTPair22",  3000, 0,  1500);
  pTPair23 = fs->make<TH1F>("pTPair23", "pTPair23",  3000, 0,  1500);


  dPhiPair1 = fs->make<TH1F>("dPhiPair1", "dPhiPair1",  300, -3.2,  3.2);
  dPhiPair2 = fs->make<TH1F>("dPhiPair2", "dPhiPair2",  300, -3.2,  3.2);
  dPhiPair3 = fs->make<TH1F>("dPhiPair3", "dPhiPair3",  300, -3.2,  3.2);
  dPhiPair4 = fs->make<TH1F>("dPhiPair4", "dPhiPair4",  300, -3.2,  3.2);
  dPhiPair5 = fs->make<TH1F>("dPhiPair5", "dPhiPair5",  300, -3.2,  3.2);
  dPhiPair6 = fs->make<TH1F>("dPhiPair6", "dPhiPair6",  300, -3.2,  3.2);
  dPhiPair7 = fs->make<TH1F>("dPhiPair7", "dPhiPair7",  300, -3.2,  3.2);
  dPhiPair8 = fs->make<TH1F>("dPhiPair8", "dPhiPair8",  300, -3.2,  3.2);
  dPhiPair9 = fs->make<TH1F>("dPhiPair9", "dPhiPair9",  300, -3.2,  3.2);
  dPhiPair10 = fs->make<TH1F>("dPhiPair10", "dPhiPair10",  300, -3.2,  3.2);
  dPhiPair11 = fs->make<TH1F>("dPhiPair11", "dPhiPair11",  300, -3.2,  3.2);
  dPhiPair12 = fs->make<TH1F>("dPhiPair12", "dPhiPair12",  300, -3.2,  3.2);
  dPhiPair13 = fs->make<TH1F>("dPhiPair13", "dPhiPair13",  300, -3.2,  3.2);
  dPhiPair14 = fs->make<TH1F>("dPhiPair14", "dPhiPair14",  300, -3.2,  3.2);
  dPhiPair15 = fs->make<TH1F>("dPhiPair15", "dPhiPair15",  300, -3.2,  3.2);
  dPhiPair16 = fs->make<TH1F>("dPhiPair16", "dPhiPair16",  300, -3.2,  3.2);
  dPhiPair17 = fs->make<TH1F>("dPhiPair17", "dPhiPair17",  300, -3.2,  3.2);
  dPhiPair18 = fs->make<TH1F>("dPhiPair18", "dPhiPair18",  300, -3.2,  3.2);
  dPhiPair19 = fs->make<TH1F>("dPhiPair19", "dPhiPair19",  300, -3.2,  3.2);
  dPhiPair20 = fs->make<TH1F>("dPhiPair20", "dPhiPair20",  300, -3.2,  3.2);
  dPhiPair21 = fs->make<TH1F>("dPhiPair21", "dPhiPair21",  300, -3.2,  3.2);
  dPhiPair22 = fs->make<TH1F>("dPhiPair22", "dPhiPair22",  300, -3.2,  3.2);
  dPhiPair23 = fs->make<TH1F>("dPhiPair23", "dPhiPair23",  300, -3.2,  3.2);

  drPair1 = fs->make<TH1F>("drPair1", "drPair1",  3000, -1500,  1500);
  drPair2 = fs->make<TH1F>("drPair2", "drPair2",  3000, -1500,  1500);
  drPair3 = fs->make<TH1F>("drPair3", "drPair3",  3000, -1500,  1500);
  drPair4 = fs->make<TH1F>("drPair4", "drPair4",  3000, -1500,  1500);
  drPair5 = fs->make<TH1F>("drPair5", "drPair5",  3000, -1500,  1500);
  drPair6 = fs->make<TH1F>("drPair6", "drPair6",  3000, -1500,  1500);
  drPair7 = fs->make<TH1F>("drPair7", "drPair7",  3000, -1500,  1500);
  drPair8 = fs->make<TH1F>("drPair8", "drPair8",  3000, -1500,  1500);
  drPair9 = fs->make<TH1F>("drPair9", "drPair9",  3000, -1500,  1500);
  drPair10 = fs->make<TH1F>("drPair10", "drPair10",  3000, -1500,  1500);
  drPair11 = fs->make<TH1F>("drPair11", "drPair11",  3000, -1500,  1500);
  drPair12 = fs->make<TH1F>("drPair12", "drPair12",  3000, -1500,  1500);
  drPair13 = fs->make<TH1F>("drPair13", "drPair13",  3000, -1500,  1500);
  drPair14 = fs->make<TH1F>("drPair14", "drPair14",  3000, -1500,  1500);
  drPair15 = fs->make<TH1F>("drPair15", "drPair15",  3000, -1500,  1500);
  drPair16 = fs->make<TH1F>("drPair16", "drPair16",  3000, -1500,  1500);
  drPair17 = fs->make<TH1F>("drPair17", "drPair17",  3000, -1500,  1500);
  drPair18 = fs->make<TH1F>("drPair18", "drPair18",  3000, -1500,  1500);
  drPair19 = fs->make<TH1F>("drPair19", "drPair19",  3000, -1500,  1500);
  drPair20 = fs->make<TH1F>("drPair20", "drPair20",  3000, -1500,  1500);
  drPair21 = fs->make<TH1F>("drPair21", "drPair21",  3000, -1500,  1500);
  drPair22 = fs->make<TH1F>("drPair22", "drPair22",  3000, -1500,  1500);
  drPair23 = fs->make<TH1F>("drPair23", "drPair23",  3000, -1500,  1500);

  dzPair1 = fs->make<TH1F>("dzPair1", "dzPair1",  3000, -1500,  1500);
  dzPair2 = fs->make<TH1F>("dzPair2", "dzPair2",  3000, -1500,  1500);
  dzPair3 = fs->make<TH1F>("dzPair3", "dzPair3",  3000, -1500,  1500);
  dzPair4 = fs->make<TH1F>("dzPair4", "dzPair4",  3000, -1500,  1500);
  dzPair5 = fs->make<TH1F>("dzPair5", "dzPair5",  3000, -1500,  1500);
  dzPair6 = fs->make<TH1F>("dzPair6", "dzPair6",  3000, -1500,  1500);
  dzPair7 = fs->make<TH1F>("dzPair7", "dzPair7",  3000, -1500,  1500);
  dzPair8 = fs->make<TH1F>("dzPair8", "dzPair8",  3000, -1500,  1500);
  dzPair9 = fs->make<TH1F>("dzPair9", "dzPair9",  3000, -1500,  1500);
  dzPair10 = fs->make<TH1F>("dzPair10", "dzPair10",  3000, -1500,  1500);
  dzPair11 = fs->make<TH1F>("dzPair11", "dzPair11",  3000, -1500,  1500);
  dzPair12 = fs->make<TH1F>("dzPair12", "dzPair12",  3000, -1500,  1500);
  dzPair13 = fs->make<TH1F>("dzPair13", "dzPair13",  3000, -1500,  1500);
  dzPair14 = fs->make<TH1F>("dzPair14", "dzPair14",  3000, -1500,  1500);
  dzPair15 = fs->make<TH1F>("dzPair15", "dzPair15",  3000, -1500,  1500);
  dzPair16 = fs->make<TH1F>("dzPair16", "dzPair16",  3000, -1500,  1500);
  dzPair17 = fs->make<TH1F>("dzPair17", "dzPair17",  3000, -1500,  1500);
  dzPair18 = fs->make<TH1F>("dzPair18", "dzPair18",  3000, -1500,  1500);
  dzPair19 = fs->make<TH1F>("dzPair19", "dzPair19",  3000, -1500,  1500);
  dzPair20 = fs->make<TH1F>("dzPair20", "dzPair20",  3000, -1500,  1500);
  dzPair21 = fs->make<TH1F>("dzPair21", "dzPair21",  3000, -1500,  1500);
  dzPair22 = fs->make<TH1F>("dzPair22", "dzPair22",  3000, -1500,  1500);
  dzPair23 = fs->make<TH1F>("dzPair23", "dzPair23",  3000, -1500,  1500);




  zSMB1 = fs->make<TH1F>("zSMB1", "z of segments in MB1",  3000, -1500,  1500);
  zSMB2 = fs->make<TH1F>("zSMB2", "z of segments in MB2",  3000, -1500,  1500);
  zSMB3 = fs->make<TH1F>("zSMB3", "z of segments in MB3",  3000, -1500,  1500);
  zSMB4 = fs->make<TH1F>("zSMB4", "z of segments in MB4",  3000, -1500,  1500);
  zSME1p = fs->make<TH1F>("zSME1p", "z of segments in ME1p",  3000, -1500,  1500);
  zSME2p = fs->make<TH1F>("zSME2p", "z of segments in ME2p",  3000, -1500,  1500);
  zSME3p = fs->make<TH1F>("zSME3p", "z of segments in ME3p",  3000, -1500,  1500);
  zSME4p = fs->make<TH1F>("zSME4p", "z of segments in ME4p",  3000, -1500,  1500);
  zSME1n = fs->make<TH1F>("zSME1n", "z of segments in ME1n",  3000, -1500,  1500);
  zSME2n = fs->make<TH1F>("zSME2n", "z of segments in ME2n",  3000, -1500,  1500);
  zSME3n = fs->make<TH1F>("zSME3n", "z of segments in ME3n",  3000, -1500,  1500);
  zSME4n = fs->make<TH1F>("zSME4n", "z of segments in ME4n",  3000, -1500,  1500);

  rSMB1 = fs->make<TH1F>("rSMB1", "r of segments in MB1",  3000, -1500,  1500);
  rSMB2 = fs->make<TH1F>("rSMB2", "r of segments in MB2",  3000, -1500,  1500);
  rSMB3 = fs->make<TH1F>("rSMB3", "r of segments in MB3",  3000, -1500,  1500);
  rSMB4 = fs->make<TH1F>("rSMB4", "r of segments in MB4",  3000, -1500,  1500);
  rSME1p = fs->make<TH1F>("rSME1p", "r of segments in ME1p",  3000, -1500,  1500);
  rSME2p = fs->make<TH1F>("rSME2p", "r of segments in ME2p",  3000, -1500,  1500);
  rSME3p = fs->make<TH1F>("rSME3p", "r of segments in ME3p",  3000, -1500,  1500);
  rSME4p = fs->make<TH1F>("rSME4p", "r of segments in ME4p",  3000, -1500,  1500);
  rSME1n = fs->make<TH1F>("rSME1n", "r of segments in ME1n",  3000, -1500,  1500);
  rSME2n = fs->make<TH1F>("rSME2n", "r of segments in ME2n",  3000, -1500,  1500);
  rSME3n = fs->make<TH1F>("rSME3n", "r of segments in ME3n",  3000, -1500,  1500);
  rSME4n = fs->make<TH1F>("rSME4n", "r of segments in ME4n",  3000, -1500,  1500);

  rzSMB1 = fs->make<TH2F>("rzSMB1", "rz of segments in MB1",  3000, -1500,  1500, 1500, 0, 1500);
  rzSMB2 = fs->make<TH2F>("rzSMB2", "rz of segments in MB2",  3000, -1500,  1500, 1500, 0, 1500);
  rzSMB3 = fs->make<TH2F>("rzSMB3", "rz of segments in MB3",  3000, -1500,  1500, 1500, 0, 1500);
  rzSMB4 = fs->make<TH2F>("rzSMB4", "rz of segments in MB4",  3000, -1500,  1500, 1500, 0, 1500);
  rzSME1p = fs->make<TH2F>("rzSME1p", "rz of segments in ME1p",  3000, -1500,  1500, 1500, 0, 1500);
  rzSME2p = fs->make<TH2F>("rzSME2p", "rz of segments in ME2p",  3000, -1500,  1500, 1500, 0, 1500);
  rzSME3p = fs->make<TH2F>("rzSME3p", "rz of segments in ME3p",  3000, -1500,  1500, 1500, 0, 1500);
  rzSME4p = fs->make<TH2F>("rzSME4p", "rz of segments in ME4p",  3000, -1500,  1500, 1500, 0, 1500);
  rzSME1n = fs->make<TH2F>("rzSME1n", "rz of segments in ME1n",  3000, -1500,  1500, 1500, 0, 1500);
  rzSME2n = fs->make<TH2F>("rzSME2n", "rz of segments in ME2n",  3000, -1500,  1500, 1500, 0, 1500);
  rzSME3n = fs->make<TH2F>("rzSME3n", "rz of segments in ME3n",  3000, -1500,  1500, 1500, 0, 1500);
  rzSME4n = fs->make<TH2F>("rzSME4n", "rz of segments in ME4n",  3000, -1500,  1500, 1500, 0, 1500);


  nSegmentsPairS1    = fs->make<TH1F>("nSegmentsPairS1", "segments / event",  300, 0,  300);
  nSegmentsPairS2    = fs->make<TH1F>("nSegmentsPairS2", "segments / event",  300, 0,  300);
  nSegmentsPairS3    = fs->make<TH1F>("nSegmentsPairS3", "segments / event",  300, 0,  300);
  nSegmentsPairS4    = fs->make<TH1F>("nSegmentsPairS4", "segments / event",  300, 0,  300);
  nSegmentsPairS5    = fs->make<TH1F>("nSegmentsPairS5", "segments / event",  300, 0,  300);
  nSegmentsPairS6    = fs->make<TH1F>("nSegmentsPairS6", "segments / event",  300, 0,  300);
  nSegmentsPairS7    = fs->make<TH1F>("nSegmentsPairS7", "segments / event",  300, 0,  300);
  nSegmentsPairS8    = fs->make<TH1F>("nSegmentsPairS8", "segments / event",  300, 0,  300);
  nSegmentsPairS9    = fs->make<TH1F>("nSegmentsPairS9", "segments / event",  300, 0,  300);
  nSegmentsPairS10    = fs->make<TH1F>("nSegmentsPairS10", "segments / event",  300, 0,  300);
  nSegmentsPairS11    = fs->make<TH1F>("nSegmentsPairS11", "segments / event",  300, 0,  300);
  nSegmentsPairS12    = fs->make<TH1F>("nSegmentsPairS12", "segments / event",  300, 0,  300);
  nSegmentsPairS13    = fs->make<TH1F>("nSegmentsPairS13", "segments / event",  300, 0,  300);
  nSegmentsPairS14    = fs->make<TH1F>("nSegmentsPairS14", "segments / event",  300, 0,  300);
  nSegmentsPairS15    = fs->make<TH1F>("nSegmentsPairS15", "segments / event",  300, 0,  300);
  nSegmentsPairS16    = fs->make<TH1F>("nSegmentsPairS16", "segments / event",  300, 0,  300);
  nSegmentsPairS17    = fs->make<TH1F>("nSegmentsPairS17", "segments / event",  300, 0,  300);
  nSegmentsPairS18    = fs->make<TH1F>("nSegmentsPairS18", "segments / event",  300, 0,  300);
  nSegmentsPairS19    = fs->make<TH1F>("nSegmentsPairS19", "segments / event",  300, 0,  300);
  nSegmentsPairS20    = fs->make<TH1F>("nSegmentsPairS20", "segments / event",  300, 0,  300);
  nSegmentsPairS21    = fs->make<TH1F>("nSegmentsPairS21", "segments / event",  300, 0,  300);
  nSegmentsPairS22    = fs->make<TH1F>("nSegmentsPairS22", "segments / event",  300, 0,  300);
  nSegmentsPairS23    = fs->make<TH1F>("nSegmentsPairS23", "segments / event",  300, 0,  300);



  nSegmentsPairG1    = fs->make<TH1F>("nSegmentsPairG1", "segments / event",  300, 0,  300);
  nSegmentsPairG2    = fs->make<TH1F>("nSegmentsPairG2", "segments / event",  300, 0,  300);
  nSegmentsPairG3    = fs->make<TH1F>("nSegmentsPairG3", "segments / event",  300, 0,  300);
  nSegmentsPairG4    = fs->make<TH1F>("nSegmentsPairG4", "segments / event",  300, 0,  300);
  nSegmentsPairG5    = fs->make<TH1F>("nSegmentsPairG5", "segments / event",  300, 0,  300);
  nSegmentsPairG6    = fs->make<TH1F>("nSegmentsPairG6", "segments / event",  300, 0,  300);
  nSegmentsPairG7    = fs->make<TH1F>("nSegmentsPairG7", "segments / event",  300, 0,  300);
  nSegmentsPairG8    = fs->make<TH1F>("nSegmentsPairG8", "segments / event",  300, 0,  300);
  nSegmentsPairG9    = fs->make<TH1F>("nSegmentsPairG9", "segments / event",  300, 0,  300);
  nSegmentsPairG10    = fs->make<TH1F>("nSegmentsPairG10", "segments / event",  300, 0,  300);
  nSegmentsPairG11    = fs->make<TH1F>("nSegmentsPairG11", "segments / event",  300, 0,  300);
  nSegmentsPairG12    = fs->make<TH1F>("nSegmentsPairG12", "segments / event",  300, 0,  300);
  nSegmentsPairG13    = fs->make<TH1F>("nSegmentsPairG13", "segments / event",  300, 0,  300);
  nSegmentsPairG14    = fs->make<TH1F>("nSegmentsPairG14", "segments / event",  300, 0,  300);
  nSegmentsPairG15    = fs->make<TH1F>("nSegmentsPairG15", "segments / event",  300, 0,  300);
  nSegmentsPairG16    = fs->make<TH1F>("nSegmentsPairG16", "segments / event",  300, 0,  300);
  nSegmentsPairG17    = fs->make<TH1F>("nSegmentsPairG17", "segments / event",  300, 0,  300);
  nSegmentsPairG18    = fs->make<TH1F>("nSegmentsPairG18", "segments / event",  300, 0,  300);
  nSegmentsPairG19    = fs->make<TH1F>("nSegmentsPairG19", "segments / event",  300, 0,  300);
  nSegmentsPairG20    = fs->make<TH1F>("nSegmentsPairG20", "segments / event",  300, 0,  300);
  nSegmentsPairG21    = fs->make<TH1F>("nSegmentsPairG21", "segments / event",  300, 0,  300);
  nSegmentsPairG22    = fs->make<TH1F>("nSegmentsPairG22", "segments / event",  300, 0,  300);
  nSegmentsPairG23    = fs->make<TH1F>("nSegmentsPairG23", "segments / event",  300, 0,  300);

  distPairS1 = fs->make<TH1F>("distPairS1", "distPairS1",  3000, -1500,  1500);
  distPairS2 = fs->make<TH1F>("distPairS2", "distPairS2",  3000, -1500,  1500);
  distPairS3 = fs->make<TH1F>("distPairS3", "distPairS3",  3000, -1500,  1500);
  distPairS4 = fs->make<TH1F>("distPairS4", "distPairS4",  3000, -1500,  1500);
  distPairS5 = fs->make<TH1F>("distPairS5", "distPairS5",  3000, -1500,  1500);
  distPairS6 = fs->make<TH1F>("distPairS6", "distPairS6",  3000, -1500,  1500);
  distPairS7 = fs->make<TH1F>("distPairS7", "distPairS7",  3000, -1500,  1500);
  distPairS8 = fs->make<TH1F>("distPairS8", "distPairS8",  3000, -1500,  1500);
  distPairS9 = fs->make<TH1F>("distPairS9", "distPairS9",  3000, -1500,  1500);
  distPairS10 = fs->make<TH1F>("distPairS10", "distPairS10",  3000, -1500,  1500);
  distPairS11 = fs->make<TH1F>("distPairS11", "distPairS11",  3000, -1500,  1500);
  distPairS12 = fs->make<TH1F>("distPairS12", "distPairS12",  3000, -1500,  1500);
  distPairS13 = fs->make<TH1F>("distPairS13", "distPairS13",  3000, -1500,  1500);
  distPairS14 = fs->make<TH1F>("distPairS14", "distPairS14",  3000, -1500,  1500);
  distPairS15 = fs->make<TH1F>("distPairS15", "distPairS15",  3000, -1500,  1500);
  distPairS16 = fs->make<TH1F>("distPairS16", "distPairS16",  3000, -1500,  1500);
  distPairS17 = fs->make<TH1F>("distPairS17", "distPairS17",  3000, -1500,  1500);
  distPairS18 = fs->make<TH1F>("distPairS18", "distPairS18",  3000, -1500,  1500);
  distPairS19 = fs->make<TH1F>("distPairS19", "distPairS19",  3000, -1500,  1500);
  distPairS20 = fs->make<TH1F>("distPairS20", "distPairS20",  3000, -1500,  1500);
  distPairS21 = fs->make<TH1F>("distPairS21", "distPairS21",  3000, -1500,  1500);
  distPairS22 = fs->make<TH1F>("distPairS22", "distPairS22",  3000, -1500,  1500);
  distPairS23 = fs->make<TH1F>("distPairS23", "distPairS23",  3000, -1500,  1500);

  pTPairS1 = fs->make<TH1F>("pTPairS1", "pTPairS1",  3000, 0,  1500);
  pTPairS2 = fs->make<TH1F>("pTPairS2", "pTPairS2",  3000, 0,  1500);
  pTPairS3 = fs->make<TH1F>("pTPairS3", "pTPairS3",  3000, 0,  1500);
  pTPairS4 = fs->make<TH1F>("pTPairS4", "pTPairS4",  3000, 0,  1500);
  pTPairS5 = fs->make<TH1F>("pTPairS5", "pTPairS5",  3000, 0,  1500);
  pTPairS6 = fs->make<TH1F>("pTPairS6", "pTPairS6",  3000, 0,  1500);
  pTPairS7 = fs->make<TH1F>("pTPairS7", "pTPairS7",  3000, 0,  1500);
  pTPairS8 = fs->make<TH1F>("pTPairS8", "pTPairS8",  3000, 0,  1500);
  pTPairS9 = fs->make<TH1F>("pTPairS9", "pTPairS9",  3000, 0,  1500);
  pTPairS10 = fs->make<TH1F>("pTPairS10", "pTPairS10",  3000, 0,  1500);
  pTPairS11 = fs->make<TH1F>("pTPairS11", "pTPairS11",  3000, 0,  1500);
  pTPairS12 = fs->make<TH1F>("pTPairS12", "pTPairS12",  3000, 0,  1500);
  pTPairS13 = fs->make<TH1F>("pTPairS13", "pTPairS13",  3000, 0,  1500);
  pTPairS14 = fs->make<TH1F>("pTPairS14", "pTPairS14",  3000, 0,  1500);
  pTPairS15 = fs->make<TH1F>("pTPairS15", "pTPairS15",  3000, 0,  1500);
  pTPairS16 = fs->make<TH1F>("pTPairS16", "pTPairS16",  3000, 0,  1500);
  pTPairS17 = fs->make<TH1F>("pTPairS17", "pTPairS17",  3000, 0,  1500);
  pTPairS18 = fs->make<TH1F>("pTPairS18", "pTPairS18",  3000, 0,  1500);
  pTPairS19 = fs->make<TH1F>("pTPairS19", "pTPairS19",  3000, 0,  1500);
  pTPairS20 = fs->make<TH1F>("pTPairS20", "pTPairS20",  3000, 0,  1500);
  pTPairS21 = fs->make<TH1F>("pTPairS21", "pTPairS21",  3000, 0,  1500);
  pTPairS22 = fs->make<TH1F>("pTPairS22", "pTPairS22",  3000, 0,  1500);
  pTPairS23 = fs->make<TH1F>("pTPairS23", "pTPairS23",  3000, 0,  1500);


  dPhiPairS1 = fs->make<TH1F>("dPhiPairS1", "dPhiPairS1",  300, -3.2,  3.2);
  dPhiPairS2 = fs->make<TH1F>("dPhiPairS2", "dPhiPairS2",  300, -3.2,  3.2);
  dPhiPairS3 = fs->make<TH1F>("dPhiPairS3", "dPhiPairS3",  300, -3.2,  3.2);
  dPhiPairS4 = fs->make<TH1F>("dPhiPairS4", "dPhiPairS4",  300, -3.2,  3.2);
  dPhiPairS5 = fs->make<TH1F>("dPhiPairS5", "dPhiPairS5",  300, -3.2,  3.2);
  dPhiPairS6 = fs->make<TH1F>("dPhiPairS6", "dPhiPairS6",  300, -3.2,  3.2);
  dPhiPairS7 = fs->make<TH1F>("dPhiPairS7", "dPhiPairS7",  300, -3.2,  3.2);
  dPhiPairS8 = fs->make<TH1F>("dPhiPairS8", "dPhiPairS8",  300, -3.2,  3.2);
  dPhiPairS9 = fs->make<TH1F>("dPhiPairS9", "dPhiPairS9",  300, -3.2,  3.2);
  dPhiPairS10 = fs->make<TH1F>("dPhiPairS10", "dPhiPairS10",  300, -3.2,  3.2);
  dPhiPairS11 = fs->make<TH1F>("dPhiPairS11", "dPhiPairS11",  300, -3.2,  3.2);
  dPhiPairS12 = fs->make<TH1F>("dPhiPairS12", "dPhiPairS12",  300, -3.2,  3.2);
  dPhiPairS13 = fs->make<TH1F>("dPhiPairS13", "dPhiPairS13",  300, -3.2,  3.2);
  dPhiPairS14 = fs->make<TH1F>("dPhiPairS14", "dPhiPairS14",  300, -3.2,  3.2);
  dPhiPairS15 = fs->make<TH1F>("dPhiPairS15", "dPhiPairS15",  300, -3.2,  3.2);
  dPhiPairS16 = fs->make<TH1F>("dPhiPairS16", "dPhiPairS16",  300, -3.2,  3.2);
  dPhiPairS17 = fs->make<TH1F>("dPhiPairS17", "dPhiPairS17",  300, -3.2,  3.2);
  dPhiPairS18 = fs->make<TH1F>("dPhiPairS18", "dPhiPairS18",  300, -3.2,  3.2);
  dPhiPairS19 = fs->make<TH1F>("dPhiPairS19", "dPhiPairS19",  300, -3.2,  3.2);
  dPhiPairS20 = fs->make<TH1F>("dPhiPairS20", "dPhiPairS20",  300, -3.2,  3.2);
  dPhiPairS21 = fs->make<TH1F>("dPhiPairS21", "dPhiPairS21",  300, -3.2,  3.2);
  dPhiPairS22 = fs->make<TH1F>("dPhiPairS22", "dPhiPairS22",  300, -3.2,  3.2);
  dPhiPairS23 = fs->make<TH1F>("dPhiPairS23", "dPhiPairS23",  300, -3.2,  3.2);

  drPairS1 = fs->make<TH1F>("drPairS1", "drPairS1",  3000, -1500,  1500);
  drPairS2 = fs->make<TH1F>("drPairS2", "drPairS2",  3000, -1500,  1500);
  drPairS3 = fs->make<TH1F>("drPairS3", "drPairS3",  3000, -1500,  1500);
  drPairS4 = fs->make<TH1F>("drPairS4", "drPairS4",  3000, -1500,  1500);
  drPairS5 = fs->make<TH1F>("drPairS5", "drPairS5",  3000, -1500,  1500);
  drPairS6 = fs->make<TH1F>("drPairS6", "drPairS6",  3000, -1500,  1500);
  drPairS7 = fs->make<TH1F>("drPairS7", "drPairS7",  3000, -1500,  1500);
  drPairS8 = fs->make<TH1F>("drPairS8", "drPairS8",  3000, -1500,  1500);
  drPairS9 = fs->make<TH1F>("drPairS9", "drPairS9",  3000, -1500,  1500);
  drPairS10 = fs->make<TH1F>("drPairS10", "drPairS10",  3000, -1500,  1500);
  drPairS11 = fs->make<TH1F>("drPairS11", "drPairS11",  3000, -1500,  1500);
  drPairS12 = fs->make<TH1F>("drPairS12", "drPairS12",  3000, -1500,  1500);
  drPairS13 = fs->make<TH1F>("drPairS13", "drPairS13",  3000, -1500,  1500);
  drPairS14 = fs->make<TH1F>("drPairS14", "drPairS14",  3000, -1500,  1500);
  drPairS15 = fs->make<TH1F>("drPairS15", "drPairS15",  3000, -1500,  1500);
  drPairS16 = fs->make<TH1F>("drPairS16", "drPairS16",  3000, -1500,  1500);
  drPairS17 = fs->make<TH1F>("drPairS17", "drPairS17",  3000, -1500,  1500);
  drPairS18 = fs->make<TH1F>("drPairS18", "drPairS18",  3000, -1500,  1500);
  drPairS19 = fs->make<TH1F>("drPairS19", "drPairS19",  3000, -1500,  1500);
  drPairS20 = fs->make<TH1F>("drPairS20", "drPairS20",  3000, -1500,  1500);
  drPairS21 = fs->make<TH1F>("drPairS21", "drPairS21",  3000, -1500,  1500);
  drPairS22 = fs->make<TH1F>("drPairS22", "drPairS22",  3000, -1500,  1500);
  drPairS23 = fs->make<TH1F>("drPairS23", "drPairS23",  3000, -1500,  1500);

  dzPairS1 = fs->make<TH1F>("dzPairS1", "dzPairS1",  3000, -1500,  1500);
  dzPairS2 = fs->make<TH1F>("dzPairS2", "dzPairS2",  3000, -1500,  1500);
  dzPairS3 = fs->make<TH1F>("dzPairS3", "dzPairS3",  3000, -1500,  1500);
  dzPairS4 = fs->make<TH1F>("dzPairS4", "dzPairS4",  3000, -1500,  1500);
  dzPairS5 = fs->make<TH1F>("dzPairS5", "dzPairS5",  3000, -1500,  1500);
  dzPairS6 = fs->make<TH1F>("dzPairS6", "dzPairS6",  3000, -1500,  1500);
  dzPairS7 = fs->make<TH1F>("dzPairS7", "dzPairS7",  3000, -1500,  1500);
  dzPairS8 = fs->make<TH1F>("dzPairS8", "dzPairS8",  3000, -1500,  1500);
  dzPairS9 = fs->make<TH1F>("dzPairS9", "dzPairS9",  3000, -1500,  1500);
  dzPairS10 = fs->make<TH1F>("dzPairS10", "dzPairS10",  3000, -1500,  1500);
  dzPairS11 = fs->make<TH1F>("dzPairS11", "dzPairS11",  3000, -1500,  1500);
  dzPairS12 = fs->make<TH1F>("dzPairS12", "dzPairS12",  3000, -1500,  1500);
  dzPairS13 = fs->make<TH1F>("dzPairS13", "dzPairS13",  3000, -1500,  1500);
  dzPairS14 = fs->make<TH1F>("dzPairS14", "dzPairS14",  3000, -1500,  1500);
  dzPairS15 = fs->make<TH1F>("dzPairS15", "dzPairS15",  3000, -1500,  1500);
  dzPairS16 = fs->make<TH1F>("dzPairS16", "dzPairS16",  3000, -1500,  1500);
  dzPairS17 = fs->make<TH1F>("dzPairS17", "dzPairS17",  3000, -1500,  1500);
  dzPairS18 = fs->make<TH1F>("dzPairS18", "dzPairS18",  3000, -1500,  1500);
  dzPairS19 = fs->make<TH1F>("dzPairS19", "dzPairS19",  3000, -1500,  1500);
  dzPairS20 = fs->make<TH1F>("dzPairS20", "dzPairS20",  3000, -1500,  1500);
  dzPairS21 = fs->make<TH1F>("dzPairS21", "dzPairS21",  3000, -1500,  1500);
  dzPairS22 = fs->make<TH1F>("dzPairS22", "dzPairS22",  3000, -1500,  1500);
  dzPairS23 = fs->make<TH1F>("dzPairS23", "dzPairS23",  3000, -1500,  1500);

  dDirPairS1 = fs->make<TH1F>("dDirPairS1", "dDirPairS1",  1000, -1,  1);
  dDirPairS2 = fs->make<TH1F>("dDirPairS2", "dDirPairS2",  1000, -1,  1);
  dDirPairS3 = fs->make<TH1F>("dDirPairS3", "dDirPairS3",  1000, -1,  1);
  dDirPairS4 = fs->make<TH1F>("dDirPairS4", "dDirPairS4",  1000, -1,  1);
  dDirPairS5 = fs->make<TH1F>("dDirPairS5", "dDirPairS5",  1000, -1,  1);
  dDirPairS6 = fs->make<TH1F>("dDirPairS6", "dDirPairS6",  1000, -1,  1);
  dDirPairS7 = fs->make<TH1F>("dDirPairS7", "dDirPairS7",  1000, -1,  1);
  dDirPairS8 = fs->make<TH1F>("dDirPairS8", "dDirPairS8",  1000, -1,  1);
  dDirPairS9 = fs->make<TH1F>("dDirPairS9", "dDirPairS9",  1000, -1,  1);
  dDirPairS10 = fs->make<TH1F>("dDirPairS10", "dDirPairS10",  1000, -1,  1);
  dDirPairS11 = fs->make<TH1F>("dDirPairS11", "dDirPairS11",  1000, -1,  1);
  dDirPairS12 = fs->make<TH1F>("dDirPairS12", "dDirPairS12",  1000, -1,  1);
  dDirPairS13 = fs->make<TH1F>("dDirPairS13", "dDirPairS13",  1000, -1,  1);
  dDirPairS14 = fs->make<TH1F>("dDirPairS14", "dDirPairS14",  1000, -1,  1);
  dDirPairS15 = fs->make<TH1F>("dDirPairS15", "dDirPairS15",  1000, -1,  1);
  dDirPairS16 = fs->make<TH1F>("dDirPairS16", "dDirPairS16",  1000, -1,  1);
  dDirPairS17 = fs->make<TH1F>("dDirPairS17", "dDirPairS17",  1000, -1,  1);
  dDirPairS18 = fs->make<TH1F>("dDirPairS18", "dDirPairS18",  1000, -1,  1);
  dDirPairS19 = fs->make<TH1F>("dDirPairS19", "dDirPairS19",  1000, -1,  1);
  dDirPairS20 = fs->make<TH1F>("dDirPairS20", "dDirPairS20",  1000, -1,  1);
  dDirPairS21 = fs->make<TH1F>("dDirPairS21", "dDirPairS21",  1000, -1,  1);
  dDirPairS22 = fs->make<TH1F>("dDirPairS22", "dDirPairS22",  1000, -1,  1);
  dDirPairS23 = fs->make<TH1F>("dDirPairS23", "dDirPairS23",  1000, -1,  1);

  distPairG1 = fs->make<TH1F>("distPairG1", "distPairG1",  3000, -1500,  1500);
  distPairG2 = fs->make<TH1F>("distPairG2", "distPairG2",  3000, -1500,  1500);
  distPairG3 = fs->make<TH1F>("distPairG3", "distPairG3",  3000, -1500,  1500);
  distPairG4 = fs->make<TH1F>("distPairG4", "distPairG4",  3000, -1500,  1500);
  distPairG5 = fs->make<TH1F>("distPairG5", "distPairG5",  3000, -1500,  1500);
  distPairG6 = fs->make<TH1F>("distPairG6", "distPairG6",  3000, -1500,  1500);
  distPairG7 = fs->make<TH1F>("distPairG7", "distPairG7",  3000, -1500,  1500);
  distPairG8 = fs->make<TH1F>("distPairG8", "distPairG8",  3000, -1500,  1500);
  distPairG9 = fs->make<TH1F>("distPairG9", "distPairG9",  3000, -1500,  1500);
  distPairG10 = fs->make<TH1F>("distPairG10", "distPairG10",  3000, -1500,  1500);
  distPairG11 = fs->make<TH1F>("distPairG11", "distPairG11",  3000, -1500,  1500);
  distPairG12 = fs->make<TH1F>("distPairG12", "distPairG12",  3000, -1500,  1500);
  distPairG13 = fs->make<TH1F>("distPairG13", "distPairG13",  3000, -1500,  1500);
  distPairG14 = fs->make<TH1F>("distPairG14", "distPairG14",  3000, -1500,  1500);
  distPairG15 = fs->make<TH1F>("distPairG15", "distPairG15",  3000, -1500,  1500);
  distPairG16 = fs->make<TH1F>("distPairG16", "distPairG16",  3000, -1500,  1500);
  distPairG17 = fs->make<TH1F>("distPairG17", "distPairG17",  3000, -1500,  1500);
  distPairG18 = fs->make<TH1F>("distPairG18", "distPairG18",  3000, -1500,  1500);
  distPairG19 = fs->make<TH1F>("distPairG19", "distPairG19",  3000, -1500,  1500);
  distPairG20 = fs->make<TH1F>("distPairG20", "distPairG20",  3000, -1500,  1500);
  distPairG21 = fs->make<TH1F>("distPairG21", "distPairG21",  3000, -1500,  1500);
  distPairG22 = fs->make<TH1F>("distPairG22", "distPairG22",  3000, -1500,  1500);
  distPairG23 = fs->make<TH1F>("distPairG23", "distPairG23",  3000, -1500,  1500);

  dPhiPairG1 = fs->make<TH1F>("dPhiPairG1", "dPhiPairG1",  300, -3.2,  3.2);
  dPhiPairG2 = fs->make<TH1F>("dPhiPairG2", "dPhiPairG2",  300, -3.2,  3.2);
  dPhiPairG3 = fs->make<TH1F>("dPhiPairG3", "dPhiPairG3",  300, -3.2,  3.2);
  dPhiPairG4 = fs->make<TH1F>("dPhiPairG4", "dPhiPairG4",  300, -3.2,  3.2);
  dPhiPairG5 = fs->make<TH1F>("dPhiPairG5", "dPhiPairG5",  300, -3.2,  3.2);
  dPhiPairG6 = fs->make<TH1F>("dPhiPairG6", "dPhiPairG6",  300, -3.2,  3.2);
  dPhiPairG7 = fs->make<TH1F>("dPhiPairG7", "dPhiPairG7",  300, -3.2,  3.2);
  dPhiPairG8 = fs->make<TH1F>("dPhiPairG8", "dPhiPairG8",  300, -3.2,  3.2);
  dPhiPairG9 = fs->make<TH1F>("dPhiPairG9", "dPhiPairG9",  300, -3.2,  3.2);
  dPhiPairG10 = fs->make<TH1F>("dPhiPairG10", "dPhiPairG10",  300, -3.2,  3.2);
  dPhiPairG11 = fs->make<TH1F>("dPhiPairG11", "dPhiPairG11",  300, -3.2,  3.2);
  dPhiPairG12 = fs->make<TH1F>("dPhiPairG12", "dPhiPairG12",  300, -3.2,  3.2);
  dPhiPairG13 = fs->make<TH1F>("dPhiPairG13", "dPhiPairG13",  300, -3.2,  3.2);
  dPhiPairG14 = fs->make<TH1F>("dPhiPairG14", "dPhiPairG14",  300, -3.2,  3.2);
  dPhiPairG15 = fs->make<TH1F>("dPhiPairG15", "dPhiPairG15",  300, -3.2,  3.2);
  dPhiPairG16 = fs->make<TH1F>("dPhiPairG16", "dPhiPairG16",  300, -3.2,  3.2);
  dPhiPairG17 = fs->make<TH1F>("dPhiPairG17", "dPhiPairG17",  300, -3.2,  3.2);
  dPhiPairG18 = fs->make<TH1F>("dPhiPairG18", "dPhiPairG18",  300, -3.2,  3.2);
  dPhiPairG19 = fs->make<TH1F>("dPhiPairG19", "dPhiPairG19",  300, -3.2,  3.2);
  dPhiPairG20 = fs->make<TH1F>("dPhiPairG20", "dPhiPairG20",  300, -3.2,  3.2);
  dPhiPairG21 = fs->make<TH1F>("dPhiPairG21", "dPhiPairG21",  300, -3.2,  3.2);
  dPhiPairG22 = fs->make<TH1F>("dPhiPairG22", "dPhiPairG22",  300, -3.2,  3.2);
  dPhiPairG23 = fs->make<TH1F>("dPhiPairG23", "dPhiPairG23",  300, -3.2,  3.2);

  drPairG1 = fs->make<TH1F>("drPairG1", "drPairG1",  3000, -1500,  1500);
  drPairG2 = fs->make<TH1F>("drPairG2", "drPairG2",  3000, -1500,  1500);
  drPairG3 = fs->make<TH1F>("drPairG3", "drPairG3",  3000, -1500,  1500);
  drPairG4 = fs->make<TH1F>("drPairG4", "drPairG4",  3000, -1500,  1500);
  drPairG5 = fs->make<TH1F>("drPairG5", "drPairG5",  3000, -1500,  1500);
  drPairG6 = fs->make<TH1F>("drPairG6", "drPairG6",  3000, -1500,  1500);
  drPairG7 = fs->make<TH1F>("drPairG7", "drPairG7",  3000, -1500,  1500);
  drPairG8 = fs->make<TH1F>("drPairG8", "drPairG8",  3000, -1500,  1500);
  drPairG9 = fs->make<TH1F>("drPairG9", "drPairG9",  3000, -1500,  1500);
  drPairG10 = fs->make<TH1F>("drPairG10", "drPairG10",  3000, -1500,  1500);
  drPairG11 = fs->make<TH1F>("drPairG11", "drPairG11",  3000, -1500,  1500);
  drPairG12 = fs->make<TH1F>("drPairG12", "drPairG12",  3000, -1500,  1500);
  drPairG13 = fs->make<TH1F>("drPairG13", "drPairG13",  3000, -1500,  1500);
  drPairG14 = fs->make<TH1F>("drPairG14", "drPairG14",  3000, -1500,  1500);
  drPairG15 = fs->make<TH1F>("drPairG15", "drPairG15",  3000, -1500,  1500);
  drPairG16 = fs->make<TH1F>("drPairG16", "drPairG16",  3000, -1500,  1500);
  drPairG17 = fs->make<TH1F>("drPairG17", "drPairG17",  3000, -1500,  1500);
  drPairG18 = fs->make<TH1F>("drPairG18", "drPairG18",  3000, -1500,  1500);
  drPairG19 = fs->make<TH1F>("drPairG19", "drPairG19",  3000, -1500,  1500);
  drPairG20 = fs->make<TH1F>("drPairG20", "drPairG20",  3000, -1500,  1500);
  drPairG21 = fs->make<TH1F>("drPairG21", "drPairG21",  3000, -1500,  1500);
  drPairG22 = fs->make<TH1F>("drPairG22", "drPairG22",  3000, -1500,  1500);
  drPairG23 = fs->make<TH1F>("drPairG23", "drPairG23",  3000, -1500,  1500);

  dzPairG1 = fs->make<TH1F>("dzPairG1", "dzPairG1",  3000, -1500,  1500);
  dzPairG2 = fs->make<TH1F>("dzPairG2", "dzPairG2",  3000, -1500,  1500);
  dzPairG3 = fs->make<TH1F>("dzPairG3", "dzPairG3",  3000, -1500,  1500);
  dzPairG4 = fs->make<TH1F>("dzPairG4", "dzPairG4",  3000, -1500,  1500);
  dzPairG5 = fs->make<TH1F>("dzPairG5", "dzPairG5",  3000, -1500,  1500);
  dzPairG6 = fs->make<TH1F>("dzPairG6", "dzPairG6",  3000, -1500,  1500);
  dzPairG7 = fs->make<TH1F>("dzPairG7", "dzPairG7",  3000, -1500,  1500);
  dzPairG8 = fs->make<TH1F>("dzPairG8", "dzPairG8",  3000, -1500,  1500);
  dzPairG9 = fs->make<TH1F>("dzPairG9", "dzPairG9",  3000, -1500,  1500);
  dzPairG10 = fs->make<TH1F>("dzPairG10", "dzPairG10",  3000, -1500,  1500);
  dzPairG11 = fs->make<TH1F>("dzPairG11", "dzPairG11",  3000, -1500,  1500);
  dzPairG12 = fs->make<TH1F>("dzPairG12", "dzPairG12",  3000, -1500,  1500);
  dzPairG13 = fs->make<TH1F>("dzPairG13", "dzPairG13",  3000, -1500,  1500);
  dzPairG14 = fs->make<TH1F>("dzPairG14", "dzPairG14",  3000, -1500,  1500);
  dzPairG15 = fs->make<TH1F>("dzPairG15", "dzPairG15",  3000, -1500,  1500);
  dzPairG16 = fs->make<TH1F>("dzPairG16", "dzPairG16",  3000, -1500,  1500);
  dzPairG17 = fs->make<TH1F>("dzPairG17", "dzPairG17",  3000, -1500,  1500);
  dzPairG18 = fs->make<TH1F>("dzPairG18", "dzPairG18",  3000, -1500,  1500);
  dzPairG19 = fs->make<TH1F>("dzPairG19", "dzPairG19",  3000, -1500,  1500);
  dzPairG20 = fs->make<TH1F>("dzPairG20", "dzPairG20",  3000, -1500,  1500);
  dzPairG21 = fs->make<TH1F>("dzPairG21", "dzPairG21",  3000, -1500,  1500);
  dzPairG22 = fs->make<TH1F>("dzPairG22", "dzPairG22",  3000, -1500,  1500);
  dzPairG23 = fs->make<TH1F>("dzPairG23", "dzPairG23",  3000, -1500,  1500);

  dDirPairG1 = fs->make<TH1F>("dDirPairG1", "dDirPairG1",  1000, -1,  1);
  dDirPairG2 = fs->make<TH1F>("dDirPairG2", "dDirPairG2",  1000, -1,  1);
  dDirPairG3 = fs->make<TH1F>("dDirPairG3", "dDirPairG3",  1000, -1,  1);
  dDirPairG4 = fs->make<TH1F>("dDirPairG4", "dDirPairG4",  1000, -1,  1);
  dDirPairG5 = fs->make<TH1F>("dDirPairG5", "dDirPairG5",  1000, -1,  1);
  dDirPairG6 = fs->make<TH1F>("dDirPairG6", "dDirPairG6",  1000, -1,  1);
  dDirPairG7 = fs->make<TH1F>("dDirPairG7", "dDirPairG7",  1000, -1,  1);
  dDirPairG8 = fs->make<TH1F>("dDirPairG8", "dDirPairG8",  1000, -1,  1);
  dDirPairG9 = fs->make<TH1F>("dDirPairG9", "dDirPairG9",  1000, -1,  1);
  dDirPairG10 = fs->make<TH1F>("dDirPairG10", "dDirPairG10",  1000, -1,  1);
  dDirPairG11 = fs->make<TH1F>("dDirPairG11", "dDirPairG11",  1000, -1,  1);
  dDirPairG12 = fs->make<TH1F>("dDirPairG12", "dDirPairG12",  1000, -1,  1);
  dDirPairG13 = fs->make<TH1F>("dDirPairG13", "dDirPairG13",  1000, -1,  1);
  dDirPairG14 = fs->make<TH1F>("dDirPairG14", "dDirPairG14",  1000, -1,  1);
  dDirPairG15 = fs->make<TH1F>("dDirPairG15", "dDirPairG15",  1000, -1,  1);
  dDirPairG16 = fs->make<TH1F>("dDirPairG16", "dDirPairG16",  1000, -1,  1);
  dDirPairG17 = fs->make<TH1F>("dDirPairG17", "dDirPairG17",  1000, -1,  1);
  dDirPairG18 = fs->make<TH1F>("dDirPairG18", "dDirPairG18",  1000, -1,  1);
  dDirPairG19 = fs->make<TH1F>("dDirPairG19", "dDirPairG19",  1000, -1,  1);
  dDirPairG20 = fs->make<TH1F>("dDirPairG20", "dDirPairG20",  1000, -1,  1);
  dDirPairG21 = fs->make<TH1F>("dDirPairG21", "dDirPairG21",  1000, -1,  1);
  dDirPairG22 = fs->make<TH1F>("dDirPairG22", "dDirPairG22",  1000, -1,  1);
  dDirPairG23 = fs->make<TH1F>("dDirPairG23", "dDirPairG23",  1000, -1,  1);


  purityNum   = fs->make<TH1F>("purityNum", "purityNum",  23, 0.5,  23.5);
  purityDenom = fs->make<TH1F>("purityDenom", "purityDenom",  23, 0.5,  23.5);


  effNum   = fs->make<TH1F>("effNum", "effNum",  23, 0.5,  23.5);
  effDenom = fs->make<TH1F>("effDenom", "effDenom",  23, 0.5,  23.5);





}
SegmentAnalyzer::~SegmentAnalyzer() {
  // do anything here that needs to be done at desctruction time
  //   // (e.g. close files, deallocate resources etc.)
  //     //
  //       // please remove this method altogether if it would be left empty
}

void SegmentAnalyzer::beginJob() {
}
void SegmentAnalyzer::endJob() {}
void SegmentAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::ESHandle<DTGeometry> dtGeomHandle;
  iSetup.get<MuonGeometryRecord>().get(dtGeomHandle);
  const DTGeometry* dtGeom = &*dtGeomHandle;

  edm::ESHandle<CSCGeometry> cscGeomHandle;
  iSetup.get<MuonGeometryRecord>().get(cscGeomHandle);
  const CSCGeometry* cscGeom = &*cscGeomHandle;

  const DTRecSegment4DCollection& dtSegments = iEvent.get(dtSegmentsGetToken_);
  const CSCSegmentCollection& cscSegments = iEvent.get(cscSegmentsGetToken_);
  const reco::TrackCollection& l2Muons = iEvent.get(l2MuonGetToken_);
  const MuonSegmentPairsHeterogeneous& segmentPairs = iEvent.get(segmentPairGetToken_);

//  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
//  iSetup.get<TransientRecHitRecord>().get(theTrackerRecHitBuilderName,theTrackerRecHitBuilder);

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);

   MuonTransientTrackingRecHitBuilder muonTransBuilder;

  nSegments    -> Fill( dtSegments.size() + cscSegments.size() );
  nSegmentsDT  -> Fill( dtSegments.size()  );
  nSegmentsCSC -> Fill( cscSegments.size() );

  int nPair1 = 0;
  int nPair2 = 0;
  int nPair3 = 0;
  int nPair4 = 0;
  int nPair5 = 0;
  int nPair6 = 0;
  int nPair7 = 0;
  int nPair8 = 0;
  int nPair9 = 0;
  int nPair10 = 0;
  int nPair11 = 0;
  int nPair12 = 0;
  int nPair13 = 0;
  int nPair14 = 0;
  int nPair15 = 0;
  int nPair16 = 0;
  int nPair17 = 0;
  int nPair18 = 0;
  int nPair19 = 0;
  int nPair20 = 0;
  int nPair21 = 0;
  int nPair22 = 0;
  int nPair23 = 0;


  int nPairS1 = 0;
  int nPairS2 = 0;
  int nPairS3 = 0;
  int nPairS4 = 0;
  int nPairS5 = 0;
  int nPairS6 = 0;
  int nPairS7 = 0;
  int nPairS8 = 0;
  int nPairS9 = 0;
  int nPairS10 = 0;
  int nPairS11 = 0;
  int nPairS12 = 0;
  int nPairS13 = 0;
  int nPairS14 = 0;
  int nPairS15 = 0;
  int nPairS16 = 0;
  int nPairS17 = 0;
  int nPairS18 = 0;
  int nPairS19 = 0;
  int nPairS20 = 0;
  int nPairS21 = 0;
  int nPairS22 = 0;
  int nPairS23 = 0;
   for (DTRecSegment4DCollection::const_iterator it = dtSegments.begin(); it != dtSegments.end(); it++) {

        DTChamberId id = (DTChamberId)(*it).chamberId();
        GlobalPoint gp = dtGeom->chamber(id)->toGlobal((*it).localPosition());
        double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);
	int stationID = (*it).chamberId().station();
	if (stationID == 1) zMB1 -> Fill(gp.z());
	if (stationID == 2) zMB2 -> Fill(gp.z());
	if (stationID == 3) zMB3 -> Fill(gp.z());
	if (stationID == 4) zMB4 -> Fill(gp.z());
	if (stationID == 1) rMB1 -> Fill(r);
	if (stationID == 2) rMB2 -> Fill(r);
	if (stationID == 3) rMB3 -> Fill(r);
	if (stationID == 4) rMB4 -> Fill(r);

	if (stationID == 1) rzMB1 -> Fill(gp.z(),r);
	if (stationID == 2) rzMB2 -> Fill(gp.z(),r);
	if (stationID == 3) rzMB3 -> Fill(gp.z(),r);
	if (stationID == 4) rzMB4 -> Fill(gp.z(),r);
	bool signalSegment = false;
        for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {

	        for (auto recHit : (*itL2).recHits()){
			if (!recHit->isValid()) continue;
                        TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
                        //TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*recHit);
			if ((tthit->globalPosition() - gp).mag() < 1e-5) signalSegment = true;
		}
        }
	if (signalSegment){
		if (stationID == 1) zSMB1 -> Fill(gp.z());
		if (stationID == 2) zSMB2 -> Fill(gp.z());
		if (stationID == 3) zSMB3 -> Fill(gp.z());
		if (stationID == 4) zSMB4 -> Fill(gp.z());

		if (stationID == 1) rSMB1 -> Fill(r);
		if (stationID == 2) rSMB2 -> Fill(r);
		if (stationID == 3) rSMB3 -> Fill(r);
		if (stationID == 4) rSMB4 -> Fill(r);

		if (stationID == 1) rzSMB1 -> Fill(gp.z(),r);
		if (stationID == 2) rzSMB2 -> Fill(gp.z(),r);
		if (stationID == 3) rzSMB3 -> Fill(gp.z(),r);
		if (stationID == 4) rzSMB4 -> Fill(gp.z(),r);

	}	

  }
   for (CSCSegmentCollection::const_iterator it = cscSegments.begin(); it != cscSegments.end(); it++) {


        CSCDetId id = (CSCDetId)(*it).cscDetId();
        const CSCChamber* cscChamber = cscGeom->chamber(id);
        GlobalPoint gp = cscChamber->toGlobal((*it).localPosition());
        double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
	int layerID = -1;
	if (id.zendcap() > 0) layerID = id.station() + 4;
	else layerID = id.station() + 8;
	if (layerID == 5) zME1p -> Fill(gp.z());
	if (layerID == 6) zME2p -> Fill(gp.z());
	if (layerID == 7) zME3p -> Fill(gp.z());
	if (layerID == 8) zME4p -> Fill(gp.z());
	if (layerID == 9) zME1n -> Fill(gp.z());
	if (layerID == 10) zME2n -> Fill(gp.z());
	if (layerID == 11) zME3n -> Fill(gp.z());
	if (layerID == 12) zME4n -> Fill(gp.z());

	if (layerID == 5) rME1p -> Fill(r);
	if (layerID == 6) rME2p -> Fill(r);
	if (layerID == 7) rME3p -> Fill(r);
	if (layerID == 8) rME4p -> Fill(r);
	if (layerID == 9) rME1n -> Fill(r);
	if (layerID == 10) rME2n -> Fill(r);
	if (layerID == 11) rME3n -> Fill(r);
	if (layerID == 12) rME4n -> Fill(r);

	if (layerID == 5) rzME1p -> Fill(gp.z(),r);
	if (layerID == 6) rzME2p -> Fill(gp.z(),r);
	if (layerID == 7) rzME3p -> Fill(gp.z(),r);
	if (layerID == 8) rzME4p -> Fill(gp.z(),r);
	if (layerID == 9) rzME1n -> Fill(gp.z(),r);
	if (layerID == 10) rzME2n -> Fill(gp.z(),r);
	if (layerID == 11) rzME3n -> Fill(gp.z(),r);
	if (layerID == 12) rzME4n -> Fill(gp.z(),r);
	bool signalSegment = false;
        for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {

	        for (auto recHit : (*itL2).recHits()){
			if (!recHit->isValid()) continue;
                        TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
                        //TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*recHit);
			if ((tthit->globalPosition() - gp).mag() < 1e-5) signalSegment = true;
		}
        }
	if (signalSegment){
		if (layerID == 5) zSME1p -> Fill(gp.z());
		if (layerID == 6) zSME2p -> Fill(gp.z());
		if (layerID == 7) zSME3p -> Fill(gp.z());
		if (layerID == 8) zSME4p -> Fill(gp.z());
		if (layerID == 9) zSME1n -> Fill(gp.z());
		if (layerID == 10) zSME2n -> Fill(gp.z());
		if (layerID == 11) zSME3n -> Fill(gp.z());
		if (layerID == 12) zSME4n -> Fill(gp.z());

		if (layerID == 5) rSME1p -> Fill(r);
		if (layerID == 6) rSME2p -> Fill(r);
		if (layerID == 7) rSME3p -> Fill(r);
		if (layerID == 8) rSME4p -> Fill(r);
		if (layerID == 9) rSME1n -> Fill(r);
		if (layerID == 10) rSME2n -> Fill(r);
		if (layerID == 11) rSME3n -> Fill(r);
		if (layerID == 12) rSME4n -> Fill(r);

		if (layerID == 5) rzSME1p -> Fill(gp.z(),r);
		if (layerID == 6) rzSME2p -> Fill(gp.z(),r);
		if (layerID == 7) rzSME3p -> Fill(gp.z(),r);
		if (layerID == 8) rzSME4p -> Fill(gp.z(),r);
		if (layerID == 9) rzSME1n -> Fill(gp.z(),r);
		if (layerID == 10) rzSME2n -> Fill(gp.z(),r);
		if (layerID == 11) rzSME3n -> Fill(gp.z(),r);
		if (layerID == 12) rzSME4n -> Fill(gp.z(),r);
	}

  }

   for (DTRecSegment4DCollection::const_iterator it = dtSegments.begin(); it != dtSegments.end(); it++) {
   	for (DTRecSegment4DCollection::const_iterator it2 = dtSegments.begin(); it2 != dtSegments.end(); it2++) {

		DTChamberId id = (DTChamberId)(*it).chamberId();
		GlobalPoint gp = dtGeom->chamber(id)->toGlobal((*it).localPosition());
		GlobalVector gv = dtGeom->chamber(id)->toGlobal((*it).localDirection());

		DTChamberId id2 = (DTChamberId)(*it2).chamberId();
		GlobalPoint gp2 = dtGeom->chamber(id2)->toGlobal((*it2).localPosition());
		GlobalVector gv2 = dtGeom->chamber(id2)->toGlobal((*it2).localDirection());
	 


		double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
		double r2 = pow(gp2.x()*gp2.x() + gp2.y()*gp2.y(),0.5);	
		double dist = std::abs(gp.z()*r2 - gp2.z()*r)/(r2-r);
		int stationID = (*it).chamberId().station();
		int stationID2 = (*it2).chamberId().station();
		double dDir = gv.dot(gv2)/(gv.mag()*gv2.mag());
		double dPhi = std::min(std::abs(gp.phi() - gp2.phi()), std::abs(gp2.phi() - gp.phi()));
                double minRadius2T4 = 4. * (87.78*87.78);
		double pT = dPhi*dPhi * (minRadius2T4 - r*r2); 

		if (stationID == 1 && stationID2 == 2){
			distPair1 -> Fill(dist);
			drPair1   -> Fill(r2-r);
			dPhiPair1 -> Fill(dPhi);
			dDirPair1 -> Fill(dDir);
			pTPair1   -> Fill(pT);
			dzPair1   -> Fill(gp2.z() - gp.z());
			nPair1++;
		}
		if (stationID == 2 && stationID2 == 3){
			distPair4 -> Fill(dist);
			dDirPair4 -> Fill(dDir);
			drPair4   -> Fill(r2-r);
			dPhiPair4  -> Fill(dPhi);
			pTPair4   -> Fill(pT);
			dzPair4   -> Fill(gp2.z() - gp.z());
			nPair4++;

		}
		if (stationID == 3 && stationID2 == 4){
			distPair9 -> Fill(dist);
			dDirPair9 -> Fill(dDir);
			drPair9   -> Fill(r2-r);
			dPhiPair9 -> Fill(dPhi);
			pTPair9   -> Fill(pT);
			dzPair9   -> Fill(gp2.z() - gp.z());
			nPair9++;

		}
		if (stationID == 1 && stationID2 == 3){
			distPair16 -> Fill(dist);
			dDirPair16 -> Fill(dDir);
			drPair16   -> Fill(r2-r);
			dPhiPair16 -> Fill(dPhi);
			pTPair16   -> Fill(pT);
			dzPair16   -> Fill(gp2.z() - gp.z());
			nPair16++;

		}
		if (stationID == 2 && stationID2 == 4){
			distPair17 -> Fill(dist);
			dDirPair17 -> Fill(dDir);
			drPair17   -> Fill(r2-r);
			dPhiPair17 -> Fill(dPhi);
			pTPair17   -> Fill(pT);
			dzPair17   -> Fill(gp2.z() - gp.z());
			nPair17++;

		}
		bool signalSegment1 = false;
		bool signalSegment2 = false;
		for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {

			for (auto recHit : (*itL2).recHits()){
				if (!recHit->isValid()) continue;
				TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
				//TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*recHit);
				if ((tthit->globalPosition() - gp).mag() < 10e-5) signalSegment1 = true;
				if ((tthit->globalPosition() - gp2).mag() < 10e-5) signalSegment2 = true;
			}
		}
		if (signalSegment1 && signalSegment2 && dPhi < 0.4){
			if (stationID == 1 && stationID2 == 2){
				distPairS1 -> Fill(dist);
				dDirPairS1 -> Fill(dDir);
				drPairS1   -> Fill(r2-r);
				dPhiPairS1 -> Fill(dPhi);
				pTPairS1   -> Fill(pT);
				dzPairS1   -> Fill(gp2.z() - gp.z());
				nPairS1++;
			}
			if (stationID == 2 && stationID2 == 3){
				distPairS4 -> Fill(dist);
				dDirPairS4 -> Fill(dDir);
				drPairS4   -> Fill(r2-r);
				dPhiPairS4  -> Fill(dPhi);
				pTPairS4   -> Fill(pT);
				dzPairS4   -> Fill(gp2.z() - gp.z());
				nPairS4++;

			}
			if (stationID == 3 && stationID2 == 4){
				distPairS9 -> Fill(dist);
				dDirPairS9 -> Fill(dDir);
				drPairS9   -> Fill(r2-r);
				dPhiPairS9 -> Fill(dPhi);
				pTPairS9   -> Fill(pT);
				dzPairS9   -> Fill(gp2.z() - gp.z());
				nPairS9++;

			}
			if (stationID == 1 && stationID2 == 3){
				distPairS16 -> Fill(dist);
				dDirPairS16 -> Fill(dDir);
				drPairS16   -> Fill(r2-r);
				dPhiPairS16 -> Fill(dPhi);
				pTPairS16   -> Fill(pT);
				dzPairS16   -> Fill(gp2.z() - gp.z());
				nPairS16++;

			}
			if (stationID == 2 && stationID2 == 4){
				distPairS17 -> Fill(dist);
				dDirPairS17 -> Fill(dDir);
				drPairS17   -> Fill(r2-r);
				dPhiPairS17 -> Fill(dPhi);
				pTPairS17   -> Fill(pT);
				dzPairS17   -> Fill(gp2.z() - gp.z());
				nPairS17++;

			}
		}
	}
  }
   for (DTRecSegment4DCollection::const_iterator it = dtSegments.begin(); it != dtSegments.end(); it++) {
   	for (CSCSegmentCollection::const_iterator it2 = cscSegments.begin(); it2 != cscSegments.end(); it2++) {

		DTChamberId id = (DTChamberId)(*it).chamberId();
		GlobalPoint gp = dtGeom->chamber(id)->toGlobal((*it).localPosition());
		GlobalVector gv = dtGeom->chamber(id)->toGlobal((*it).localDirection());

		CSCDetId id2 = (CSCDetId)(*it2).cscDetId();
		const CSCChamber* cscChamber = cscGeom->chamber(id2);
		GlobalPoint gp2 = cscChamber->toGlobal((*it2).localPosition());
		GlobalVector gv2 = cscChamber->toGlobal((*it2).localDirection());
		int stationID2 = -1;
		if (id2.zendcap() > 0) stationID2 = id2.station() + 4;
		else stationID2 = id2.station() + 8;


		double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
		double r2 = pow(gp2.x()*gp2.x() + gp2.y()*gp2.y(),0.5);	
		double dist = std::abs(gp.z()*r2 - gp2.z()*r)/(r2-r);
		int stationID = (*it).chamberId().station();
		double dDir = gv.dot(gv2)/(gv.mag()*gv2.mag());
		double dPhi = std::min(std::abs(gp.phi() - gp2.phi()), std::abs(gp2.phi() - gp.phi()));
                double minRadius2T4 = 4. * (87.78*87.78);
		double pT = dPhi*dPhi * (minRadius2T4 - r*r2); 

		if (stationID == 1 && stationID2 == 5){
			distPair2 -> Fill(dist);
			dDirPair2 -> Fill(dDir);
			drPair2   -> Fill(r2-r);
			dPhiPair2 -> Fill(dPhi);
			pTPair2   -> Fill(pT);
			dzPair2   -> Fill(gp2.z() - gp.z());
			nPair2++;

		}
		if (stationID == 1 && stationID2 == 9){
			distPair3 -> Fill(dist);
			dDirPair3 -> Fill(dDir);
			drPair3   -> Fill(r2-r);
			dPhiPair3 -> Fill(dPhi);
			pTPair3   -> Fill(pT);
			dzPair3   -> Fill(gp2.z() - gp.z());
			nPair3++;

		}
		if (stationID == 2 && stationID2 == 5){
			distPair5 -> Fill(dist);
			dDirPair5 -> Fill(dDir);
			drPair5   -> Fill(r2-r);
			dPhiPair5 -> Fill(dPhi);
			pTPair5   -> Fill(pT);
			dzPair5   -> Fill(gp2.z() - gp.z());
			nPair5++;

		}
		if (stationID == 2 && stationID2 == 9){
			distPair6 -> Fill(dist);
			dDirPair6 -> Fill(dDir);
			drPair6   -> Fill(r2-r);
			dPhiPair6 -> Fill(dPhi);
			pTPair6   -> Fill(pT);
			dzPair6   -> Fill(gp2.z() - gp.z());
			nPair6++;

		}
		if (stationID == 3 && stationID2 == 5){
			distPair10 -> Fill(dist);
			dDirPair10 -> Fill(dDir);
			drPair10   -> Fill(r2-r);
			dPhiPair10 -> Fill(dPhi);
			pTPair10   -> Fill(pT);
			dzPair10   -> Fill(gp2.z() - gp.z());
			nPair10++;

		}
		if (stationID == 3 && stationID2 == 9){
			distPair11 -> Fill(dist);
			dDirPair11 -> Fill(dDir);
			drPair11   -> Fill(r2-r);
			dPhiPair11 -> Fill(dPhi);
			pTPair11   -> Fill(pT);
			dzPair11   -> Fill(gp2.z() - gp.z());
			nPair11++;

		}
		if (stationID == 1 && stationID2 == 6){
			distPair18 -> Fill(dist);
			dDirPair18 -> Fill(dDir);
			drPair18   -> Fill(r2-r);
			dPhiPair18 -> Fill(dPhi);
			pTPair18   -> Fill(pT);
			dzPair18   -> Fill(gp2.z() - gp.z());
			nPair18++;

		}
		if (stationID == 1 && stationID2 == 10){
			distPair19 -> Fill(dist);
			dDirPair19 -> Fill(dDir);
			drPair19   -> Fill(r2-r);			
			dPhiPair19 -> Fill(dPhi);
			pTPair19   -> Fill(pT);
			dzPair19   -> Fill(gp2.z() - gp.z());
			nPair19++;

		}
		bool signalSegment1 = false;
		bool signalSegment2 = false;
		for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {

			for (auto recHit : (*itL2).recHits()){
				if (!recHit->isValid()) continue;
				TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
				//TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*recHit);
				if ((tthit->globalPosition() - gp).mag() < 1e-5) signalSegment1 = true;
				if ((tthit->globalPosition() - gp2).mag() < 1e-5) signalSegment2 = true;
			}
		}
		if (signalSegment1 && signalSegment2 && dPhi < 0.4){

			if (stationID == 1 && stationID2 == 5){
				distPairS2 -> Fill(dist);
				dDirPairS2 -> Fill(dDir);
				drPairS2   -> Fill(r2-r);
				dPhiPairS2 -> Fill(dPhi);
				pTPairS2   -> Fill(pT);
				dzPairS2   -> Fill(gp2.z() - gp.z());
				nPairS2++;

			}
			if (stationID == 1 && stationID2 == 9){
				distPairS3 -> Fill(dist);
				dDirPairS3 -> Fill(dDir);
				drPairS3   -> Fill(r2-r);
				dPhiPairS3 -> Fill(dPhi);
				pTPairS3   -> Fill(pT);
				dzPairS3   -> Fill(gp2.z() - gp.z());
				nPairS3++;

			}
			if (stationID == 2 && stationID2 == 5){
				distPairS5 -> Fill(dist);
				dDirPairS5 -> Fill(dDir);
				drPairS5   -> Fill(r2-r);
				dPhiPairS5 -> Fill(dPhi);
				pTPairS5   -> Fill(pT);
				dzPairS5   -> Fill(gp2.z() - gp.z());
				nPairS5++;

			}
			if (stationID == 2 && stationID2 == 9){
				distPairS6 -> Fill(dist);
				dDirPairS6 -> Fill(dDir);
				drPairS6   -> Fill(r2-r);
				dPhiPairS6 -> Fill(dPhi);
				pTPairS6   -> Fill(pT);
				dzPairS6   -> Fill(gp2.z() - gp.z());
				nPairS6++;

			}
			if (stationID == 3 && stationID2 == 5){
				distPairS10 -> Fill(dist);
				dDirPairS10 -> Fill(dDir);
				drPairS10   -> Fill(r2-r);
				dPhiPairS10 -> Fill(dPhi);
				pTPairS10   -> Fill(pT);
				dzPairS10   -> Fill(gp2.z() - gp.z());
				nPairS10++;

			}
			if (stationID == 3 && stationID2 == 9){
				distPairS11 -> Fill(dist);
				dDirPairS11 -> Fill(dDir);
				drPairS11   -> Fill(r2-r);
				dPhiPairS11 -> Fill(dPhi);
				pTPairS11   -> Fill(pT);
				dzPairS11   -> Fill(gp2.z() - gp.z());
				nPairS11++;

			}
			if (stationID == 1 && stationID2 == 6){
				distPairS18 -> Fill(dist);
				dDirPairS18 -> Fill(dDir);
				drPairS18   -> Fill(r2-r);
				dPhiPairS18 -> Fill(dPhi);
				pTPairS18   -> Fill(pT);
				dzPairS18   -> Fill(gp2.z() - gp.z());
				nPairS18++;

			}
			if (stationID == 1 && stationID2 == 10){
				distPairS19 -> Fill(dist);
				dDirPairS19 -> Fill(dDir);
				drPairS19   -> Fill(r2-r);			
				dPhiPairS19 -> Fill(dPhi);
				pTPairS19   -> Fill(pT);
				dzPairS19   -> Fill(gp2.z() - gp.z());
				nPairS19++;

			}
		}

	}
  }
   for (CSCSegmentCollection::const_iterator it = cscSegments.begin(); it != cscSegments.end(); it++) {
   	for (CSCSegmentCollection::const_iterator it2 = cscSegments.begin(); it2 != cscSegments.end(); it2++) {

		CSCDetId id = (CSCDetId)(*it).cscDetId();
		const CSCChamber* cscChamber = cscGeom->chamber(id);
		GlobalPoint gp = cscChamber->toGlobal((*it).localPosition());
		GlobalVector gv = cscChamber->toGlobal((*it).localDirection());
		int stationID = -1;
		if (id.zendcap() > 0) stationID = id.station() + 4;
		else stationID = id.station() + 8;

		CSCDetId id2 = (CSCDetId)(*it2).cscDetId();
		const CSCChamber* cscChamber2 = cscGeom->chamber(id2);
		GlobalPoint gp2 = cscChamber2->toGlobal((*it2).localPosition());
		GlobalVector gv2 = cscChamber2->toGlobal((*it2).localDirection());
		int stationID2 = -1;
		if (id2.zendcap() > 0) stationID2 = id2.station() + 4;
		else stationID2 = id2.station() + 8;


		double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
		double r2 = pow(gp2.x()*gp2.x() + gp2.y()*gp2.y(),0.5);	
		double dist = std::abs(gp.z()*r2 - gp2.z()*r)/(r2-r);
		double dDir = gv.dot(gv2)/(gv.mag()*gv2.mag());
		double dPhi = std::min(std::abs(gp.phi() - gp2.phi()), std::abs(gp2.phi() - gp.phi()));
                double minRadius2T4 = 4. * (87.78*87.78);
		double pT = dPhi*dPhi * (minRadius2T4 - r*r2); 

		if (stationID == 5 && stationID2 == 6){
			distPair7 -> Fill(dist);
			dDirPair7 -> Fill(dDir);
			drPair7   -> Fill(r2-r);
			dPhiPair7 -> Fill(dPhi);
			pTPair7   -> Fill(pT);
			dzPair7   -> Fill(gp2.z() - gp.z());
			nPair7++;

		}
		if (stationID == 9 && stationID2 == 10){
			distPair8 -> Fill(dist);
			dDirPair8 -> Fill(dDir);
			drPair8   -> Fill(r2-r);
			dPhiPair8 -> Fill(dPhi);
			pTPair8   -> Fill(pT);
			dzPair8   -> Fill(gp2.z() - gp.z());
			nPair8++;

		}
		if (stationID == 6 && stationID2 == 7){
			distPair12 -> Fill(dist);
			dDirPair12 -> Fill(dDir);
			drPair12   -> Fill(r2-r);
			dPhiPair12 -> Fill(dPhi);
			pTPair12   -> Fill(pT);
			dzPair12   -> Fill(gp2.z() - gp.z());
			nPair12++;

		}
		if (stationID == 10 && stationID2 == 11){
			distPair13 -> Fill(dist);
			dDirPair13 -> Fill(dDir);
			drPair13   -> Fill(r2-r);
			dPhiPair13 -> Fill(dPhi);
			pTPair13   -> Fill(pT);
			dzPair13   -> Fill(gp2.z() - gp.z());
			nPair13++;

		}
		if (stationID == 7 && stationID2 == 8){
			distPair14 -> Fill(dist);
			dDirPair14 -> Fill(dDir);
			drPair14   -> Fill(r2-r);
			dPhiPair14 -> Fill(dPhi);
			pTPair14   -> Fill(pT);
			dzPair14   -> Fill(gp2.z() - gp.z());
			nPair14++;

		}
		if (stationID == 11 && stationID2 == 12){
			distPair15 -> Fill(dist);
			dDirPair15 -> Fill(dDir);
			drPair15   -> Fill(r2-r);
			dPhiPair15 -> Fill(dPhi);
			pTPair15   -> Fill(pT);
			dzPair15   -> Fill(gp2.z() - gp.z());
			nPair15++;

		}
		if (stationID == 5 && stationID2 == 7){
			distPair20 -> Fill(dist);
			dDirPair20 -> Fill(dDir);
			drPair20   -> Fill(r2-r);
			dPhiPair20 -> Fill(dPhi);
			pTPair20   -> Fill(pT);
			dzPair20   -> Fill(gp2.z() - gp.z());
			nPair20++;

		}
		if (stationID == 6 && stationID2 == 8){
			distPair21 -> Fill(dist);
			dDirPair21 -> Fill(dDir);
			drPair21   -> Fill(r2-r);
			dPhiPair21 -> Fill(dPhi);
			pTPair21   -> Fill(pT);
			dzPair21   -> Fill(gp2.z() - gp.z());
			nPair21++;


		}
		if (stationID == 9 && stationID2 == 11){
			distPair22 -> Fill(dist);
			dDirPair22 -> Fill(dDir);
			drPair22   -> Fill(r2-r);
			dPhiPair22 -> Fill(dPhi);
			pTPair22   -> Fill(pT);
			dzPair22   -> Fill(gp2.z() - gp.z());
			nPair22++;

		}
		if (stationID == 10 && stationID2 == 12){
			distPair23 -> Fill(dist);
			dDirPair23 -> Fill(dDir);
			drPair23   -> Fill(r2-r);
			dPhiPair23 -> Fill(dPhi);
			pTPair23   -> Fill(pT);
			dzPair23   -> Fill(gp2.z() - gp.z());
			nPair23++;

		}
		bool signalSegment1 = false;
		bool signalSegment2 = false;
		for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {

			for (auto recHit : (*itL2).recHits()){
				if (!recHit->isValid()) continue;
				TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
				//TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*recHit);
				if ((tthit->globalPosition() - gp).mag() < 1e-5) signalSegment1 = true;
				if ((tthit->globalPosition() - gp2).mag() < 1e-5) signalSegment2 = true;
			}
		}
		if (signalSegment1 && signalSegment2 && dPhi < 0.4){
			if (stationID == 5 && stationID2 == 6){
				distPairS7 -> Fill(dist);
				dDirPairS7 -> Fill(dDir);
				drPairS7   -> Fill(r2-r);
				dPhiPairS7 -> Fill(dPhi);
				pTPairS7   -> Fill(pT);
				dzPairS7   -> Fill(gp2.z() - gp.z());
				nPairS7++;

			}
			if (stationID == 9 && stationID2 == 10){
				distPairS8 -> Fill(dist);
				dDirPairS8 -> Fill(dDir);
				drPairS8   -> Fill(r2-r);
				dPhiPairS8 -> Fill(dPhi);
				pTPairS8   -> Fill(pT);
				dzPairS8   -> Fill(gp2.z() - gp.z());
				nPairS8++;

			}
			if (stationID == 6 && stationID2 == 7){
				distPairS12 -> Fill(dist);
				dDirPairS12 -> Fill(dDir);
				drPairS12   -> Fill(r2-r);
				dPhiPairS12 -> Fill(dPhi);
				pTPairS12   -> Fill(pT);
				dzPairS12   -> Fill(gp2.z() - gp.z());
				nPairS12++;

			}
			if (stationID == 10 && stationID2 == 11){
				distPairS13 -> Fill(dist);
				dDirPairS13 -> Fill(dDir);
				drPairS13   -> Fill(r2-r);
				dPhiPairS13 -> Fill(dPhi);
				pTPairS13   -> Fill(pT);
				dzPairS13   -> Fill(gp2.z() - gp.z());
				nPairS13++;

			}
			if (stationID == 7 && stationID2 == 8){
				distPairS14 -> Fill(dist);
				dDirPairS14 -> Fill(dDir);
				drPairS14   -> Fill(r2-r);
				dPhiPairS14 -> Fill(dPhi);
				pTPairS14   -> Fill(pT);
				dzPairS14   -> Fill(gp2.z() - gp.z());
				nPairS14++;

			}
			if (stationID == 11 && stationID2 == 12){
				distPairS15 -> Fill(dist);
				dDirPairS15 -> Fill(dDir);
				drPairS15   -> Fill(r2-r);
				dPhiPairS15 -> Fill(dPhi);
				pTPairS15   -> Fill(pT);
				dzPairS15   -> Fill(gp2.z() - gp.z());
				nPairS15++;

			}
			if (stationID == 5 && stationID2 == 7){
				distPairS20 -> Fill(dist);
				dDirPairS20 -> Fill(dDir);
				drPairS20   -> Fill(r2-r);
				dPhiPairS20 -> Fill(dPhi);
				pTPairS20   -> Fill(pT);
				dzPairS20   -> Fill(gp2.z() - gp.z());
				nPairS20++;

			}
			if (stationID == 6 && stationID2 == 8){
				distPairS21 -> Fill(dist);
				dDirPairS21 -> Fill(dDir);
				drPairS21   -> Fill(r2-r);
				dPhiPairS21 -> Fill(dPhi);
				pTPairS21   -> Fill(pT);
				dzPairS21   -> Fill(gp2.z() - gp.z());
				nPairS21++;


			}
			if (stationID == 9 && stationID2 == 11){
				distPairS22 -> Fill(dist);
				dDirPairS22 -> Fill(dDir);
				drPairS22   -> Fill(r2-r);
				dPhiPairS22 -> Fill(dPhi);
				pTPairS22   -> Fill(pT);
				dzPairS22   -> Fill(gp2.z() - gp.z());
				nPairS22++;

			}
			if (stationID == 10 && stationID2 == 12){
				distPairS23 -> Fill(dist);
				dDirPairS23 -> Fill(dDir);
				drPairS23   -> Fill(r2-r);
				dPhiPairS23 -> Fill(dPhi);
				pTPairS23   -> Fill(pT);
				dzPairS23   -> Fill(gp2.z() - gp.z());
				nPairS23++;

			}

		}

	}
  }

  int nPairMatchedG1 = 0;
  int nPairMatchedG2 = 0;
  int nPairMatchedG3 = 0;
  int nPairMatchedG4 = 0;
  int nPairMatchedG5 = 0;
  int nPairMatchedG6 = 0;
  int nPairMatchedG7 = 0;
  int nPairMatchedG8 = 0;
  int nPairMatchedG9 = 0;
  int nPairMatchedG10 = 0;
  int nPairMatchedG11 = 0;
  int nPairMatchedG12 = 0;
  int nPairMatchedG13 = 0;
  int nPairMatchedG14 = 0;
  int nPairMatchedG15 = 0;
  int nPairMatchedG16 = 0;
  int nPairMatchedG17 = 0;
  int nPairMatchedG18 = 0;
  int nPairMatchedG19 = 0;
  int nPairMatchedG20 = 0;
  int nPairMatchedG21 = 0;
  int nPairMatchedG22 = 0;
  int nPairMatchedG23 = 0;



  int nPairG1 = 0;
  int nPairG2 = 0;
  int nPairG3 = 0;
  int nPairG4 = 0;
  int nPairG5 = 0;
  int nPairG6 = 0;
  int nPairG7 = 0;
  int nPairG8 = 0;
  int nPairG9 = 0;
  int nPairG10 = 0;
  int nPairG11 = 0;
  int nPairG12 = 0;
  int nPairG13 = 0;
  int nPairG14 = 0;
  int nPairG15 = 0;
  int nPairG16 = 0;
  int nPairG17 = 0;
  int nPairG18 = 0;
  int nPairG19 = 0;
  int nPairG20 = 0;
  int nPairG21 = 0;
  int nPairG22 = 0;
  int nPairG23 = 0;

  auto pairs = segmentPairs.get();

  for (int i = 0; i < pairs->nPairs; i++){
	float dz = pairs->gz2[i] - pairs->gz1[i];	
	float dr = pairs->gr2[i] - pairs->gr1[i];	
	double dist = std::abs(pairs->gz1[i]*pairs->gr2[i] - pairs->gz2[i]*pairs->gr1[i])/(pairs->gr2[i]-pairs->gr1[i]);
	double dPhi = std::min(std::abs(pairs->gphi1[i] - pairs->gphi2[i]), std::abs(pairs->gphi2[i] - pairs->gphi1[i]));
	bool matchedSegment1 = false;
	bool matchedSegment2 = false;
	bool matchedPair = false;
        for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {

                for (auto recHit : (*itL2).recHits()){
                        if (!recHit->isValid()) continue;
                        TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
			if (tthit->globalPosition().x() == pairs->gx1[i] && tthit->globalPosition().y() == pairs->gy1[i]) matchedSegment1 = true;
			if (tthit->globalPosition().x() == pairs->gx2[i] && tthit->globalPosition().y() == pairs->gy2[i]) matchedSegment2 = true;
         	}
        }
        if (matchedSegment1 && matchedSegment2) matchedPair = true;                


	if (pairs->layerID1[i] == 0 && pairs->layerID2[i] == 1){
		 nPairG1++;
		 if (matchedPair) nPairMatchedG1++;
		 purityDenom->Fill(1);
		 if (matchedPair) purityNum->Fill(1);
                 distPairG1 -> Fill(dist);
                 drPairG1   -> Fill(dr);
                 dPhiPairG1 -> Fill(dPhi);
                 dzPairG1   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 1 && pairs->layerID2[i] == 2){
		 nPairG4++;
		 if (matchedPair) nPairMatchedG4++;
		 purityDenom->Fill(4);
		 if (matchedPair) purityNum->Fill(4);

                 distPairG4 -> Fill(dist);
                 drPairG4   -> Fill(dr);
                 dPhiPairG4 -> Fill(dPhi);
                 dzPairG4   -> Fill(dz);
	}

	if (pairs->layerID1[i] == 2 && pairs->layerID2[i] == 3){
		 nPairG9++;
		 if (matchedPair) nPairMatchedG9++;		 
		 purityDenom->Fill(9);
		 if (matchedPair) purityNum->Fill(9);

                 distPairG9 -> Fill(dist);
                 drPairG9   -> Fill(dr);
                 dPhiPairG9 -> Fill(dPhi);
                 dzPairG9   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 0 && pairs->layerID2[i] == 2){
		 nPairG16++;
		 if (matchedPair) nPairMatchedG16++;
		 purityDenom->Fill(16);
		 if (matchedPair) purityNum->Fill(16);

                 distPairG16 -> Fill(dist);
                 drPairG16   -> Fill(dr);
                 dPhiPairG16 -> Fill(dPhi);
                 dzPairG16   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 1 && pairs->layerID2[i] == 3){
		 nPairG17++;
		 if (matchedPair) nPairMatchedG17++;
		 purityDenom->Fill(17);
		 if (matchedPair) purityNum->Fill(17);

                 distPairG17 -> Fill(dist);
                 drPairG17   -> Fill(dr);
                 dPhiPairG17 -> Fill(dPhi);
                 dzPairG17   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 0 && pairs->layerID2[i] == 4){
		 nPairG2++;
		 if (matchedPair) nPairMatchedG2++;
		 purityDenom->Fill(2);
		 if (matchedPair) purityNum->Fill(2);

                 distPairG2 -> Fill(dist);
                 drPairG2   -> Fill(dr);
                 dPhiPairG2 -> Fill(dPhi);
                 dzPairG2   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 0 && pairs->layerID2[i] == 8){
		 nPairG3++;
		 if (matchedPair) nPairMatchedG3++;
		 purityDenom->Fill(3);
		 if (matchedPair) purityNum->Fill(3);

                 distPairG3 -> Fill(dist);
                 drPairG3   -> Fill(dr);
                 dPhiPairG3 -> Fill(dPhi);
                 dzPairG3   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 1 && pairs->layerID2[i] == 4){
		 nPairG5++;
		 if (matchedPair) nPairMatchedG5++;
		 purityDenom->Fill(5);
		 if (matchedPair) purityNum->Fill(5);

                 distPairG5 -> Fill(dist);
                 drPairG5   -> Fill(dr);
                 dPhiPairG5 -> Fill(dPhi);
                 dzPairG5   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 1 && pairs->layerID2[i] == 8){
		 nPairG6++;
		 if (matchedPair) nPairMatchedG6++;
		 purityDenom->Fill(6);
		 if (matchedPair) purityNum->Fill(6);

                 distPairG6 -> Fill(dist);
                 drPairG6   -> Fill(dr);
                 dPhiPairG6 -> Fill(dPhi);
                 dzPairG6   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 2 && pairs->layerID2[i] == 4){
		 nPairG10++;
		 if (matchedPair) nPairMatchedG10++;
		 purityDenom->Fill(10);
		 if (matchedPair) purityNum->Fill(10);

                 distPairG10 -> Fill(dist);
                 drPairG10   -> Fill(dr);
                 dPhiPairG10 -> Fill(dPhi);
                 dzPairG10   -> Fill(dz);
	}
	if (pairs->layerID1[i] == 2 && pairs->layerID2[i] == 8){
		 nPairG11++;
		 if (matchedPair) nPairMatchedG11++;
		 purityDenom->Fill(11);
		 if (matchedPair) purityNum->Fill(11);

                 distPairG11 -> Fill(dist);
                 drPairG11   -> Fill(dr);
                 dPhiPairG11 -> Fill(dPhi);
                 dzPairG11   -> Fill(dz);
		 nPairG11++;
	}
	if (pairs->layerID1[i] == 0 && pairs->layerID2[i] == 5){
		 nPairG18++;
		 if (matchedPair) nPairMatchedG18++;
		 purityDenom->Fill(18);
		 if (matchedPair) purityNum->Fill(18);

                 distPairG18 -> Fill(dist);
                 drPairG18   -> Fill(dr);
                 dPhiPairG18 -> Fill(dPhi);
                 dzPairG18   -> Fill(dz);
		 nPairG18++;
	}
	if (pairs->layerID1[i] == 0 && pairs->layerID2[i] == 9){
		 nPairG19++;
		 if (matchedPair) nPairMatchedG19++;
		 purityDenom->Fill(19);
		 if (matchedPair) purityNum->Fill(19);

                 distPairG19 -> Fill(dist);
                 drPairG19   -> Fill(dr);
                 dPhiPairG19 -> Fill(dPhi);
                 dzPairG19   -> Fill(dz);
		 nPairG19++;
	}
	if (pairs->layerID1[i] == 4 && pairs->layerID2[i] == 5){
		 nPairG7++;
		 if (matchedPair) nPairMatchedG7++;
		 purityDenom->Fill(7);
		 if (matchedPair) purityNum->Fill(7);

                 distPairG7 -> Fill(dist);
                 drPairG7   -> Fill(dr);
                 dPhiPairG7 -> Fill(dPhi);
                 dzPairG7   -> Fill(dz);
		 nPairG7++;
	}
	if (pairs->layerID1[i] == 8 && pairs->layerID2[i] == 9){
		 nPairG8++;
		 if (matchedPair) nPairMatchedG8++;
		 purityDenom->Fill(8);
		 if (matchedPair) purityNum->Fill(8);

                 distPairG8 -> Fill(dist);
                 drPairG8   -> Fill(dr);
                 dPhiPairG8 -> Fill(dPhi);
                 dzPairG8   -> Fill(dz);
		 nPairG8++;
	}
	if (pairs->layerID1[i] == 5 && pairs->layerID2[i] == 6){
		 nPairG12++;
		 if (matchedPair) nPairMatchedG12++;
		 purityDenom->Fill(12);
		 if (matchedPair) purityNum->Fill(12);

                 distPairG12 -> Fill(dist);
                 drPairG12   -> Fill(dr);
                 dPhiPairG12 -> Fill(dPhi);
                 dzPairG12   -> Fill(dz);
		 nPairG12++;
	}
	if (pairs->layerID1[i] == 9 && pairs->layerID2[i] == 10){
		 nPairG13++;
		 if (matchedPair) nPairMatchedG13++;
		 purityDenom->Fill(13);
		 if (matchedPair) purityNum->Fill(13);

                 distPairG13 -> Fill(dist);
                 drPairG13   -> Fill(dr);
                 dPhiPairG13 -> Fill(dPhi);
                 dzPairG13   -> Fill(dz);
		 nPairG13++;
	}
	if (pairs->layerID1[i] == 6 && pairs->layerID2[i] == 7){
		 nPairG14++;
		 if (matchedPair) nPairMatchedG14++;
		 purityDenom->Fill(14);
		 if (matchedPair) purityNum->Fill(14);

                 distPairG14 -> Fill(dist);
                 drPairG14   -> Fill(dr);
                 dPhiPairG14 -> Fill(dPhi);
                 dzPairG14   -> Fill(dz);
		 nPairG14++;
	}
	if (pairs->layerID1[i] == 10 && pairs->layerID2[i] == 11){
		 nPairG15++;
		 if (matchedPair) nPairMatchedG15++;
		 purityDenom->Fill(15);
		 if (matchedPair) purityNum->Fill(15);

                 distPairG15 -> Fill(dist);
                 drPairG15   -> Fill(dr);
                 dPhiPairG15 -> Fill(dPhi);
                 dzPairG15   -> Fill(dz);
		 nPairG15++;
	}
	if (pairs->layerID1[i] == 4 && pairs->layerID2[i] == 6){
		 nPairG20++;
		 if (matchedPair) nPairMatchedG20++;
		 purityDenom->Fill(20);
		 if (matchedPair) purityNum->Fill(20);

                 distPairG20 -> Fill(dist);
                 drPairG20   -> Fill(dr);
                 dPhiPairG20 -> Fill(dPhi);
                 dzPairG20   -> Fill(dz);
		 nPairG20++;
	}
	if (pairs->layerID1[i] == 5 && pairs->layerID2[i] == 7){
		 nPairG21++;
		 if (matchedPair) nPairMatchedG21++;
		 purityDenom->Fill(21);
		 if (matchedPair) purityNum->Fill(21);

                 distPairG21 -> Fill(dist);
                 drPairG21   -> Fill(dr);
                 dPhiPairG21 -> Fill(dPhi);
                 dzPairG21   -> Fill(dz);
		 nPairG21++;
	}
	if (pairs->layerID1[i] == 8 && pairs->layerID2[i] == 10){
		 nPairG22++;
		 if (matchedPair) nPairMatchedG22++;
		 purityDenom->Fill(22);
		 if (matchedPair) purityNum->Fill(22);

                 distPairG22 -> Fill(dist);
                 drPairG22   -> Fill(dr);
                 dPhiPairG22 -> Fill(dPhi);
                 dzPairG22   -> Fill(dz);
		 nPairG22++;
	}
	if (pairs->layerID1[i] == 9 && pairs->layerID2[i] == 11){
		 nPairG23++;
		 if (matchedPair) nPairMatchedG23++;
		 purityDenom->Fill(23);
		 if (matchedPair) purityNum->Fill(23);

                 distPairG23 -> Fill(dist);
                 drPairG23   -> Fill(dr);
                 dPhiPairG23 -> Fill(dPhi);
                 dzPairG23   -> Fill(dz);
		 nPairG23++;
	}


  }

  int nRecHits = 0;

  for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {
	int j = 0;
	for (auto recHit : (*itL2).recHits()){
		nRecHits++;
		if (!recHit->isValid()) continue;
		TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
		int k = -1;
		bool foundNeighbor = false;
		for (auto recHit2 : (*itL2).recHits()){
			k++;
			if (!recHit2->isValid()) continue;
			//if (foundNeighbor) continue;
			DetId i1 = recHit->geographicalId();
			DetId i2 = recHit2->geographicalId();
			if (!((i1.subdetId() == MuonSubdetId::DT || i1.subdetId() == MuonSubdetId::CSC) && (i2.subdetId() == MuonSubdetId::DT || i2.subdetId() == MuonSubdetId::CSC) )) continue;
			if (j == k) continue;
			//if (j != k) foundNeighbor = true;
			TrajectoryMeasurement::ConstRecHitPointer tthit2(muonTransBuilder.build(recHit2, globalGeometry));
			int binIndex = 0;
			int layerID = 0;
			int layerID2 = 0;
			if (i1.subdetId() == MuonSubdetId::DT){
				DTChamberId id =  (DTChamberId) recHit->geographicalId();
				layerID = id.station();
			} 
			else{
				CSCDetId id = (CSCDetId) recHit->geographicalId();
                		if (id.zendcap() > 0) layerID = id.station() + 4;
                		else layerID = id.station() + 8;
			}
			if (i2.subdetId() == MuonSubdetId::DT){
				DTChamberId id2 =  (DTChamberId) recHit2->geographicalId();
				layerID2 = id2.station();
			} 
			else{
				CSCDetId id2 = (CSCDetId) recHit2->geographicalId();
                		if (id2.zendcap() > 0) layerID2 = id2.station() + 4;
                		else layerID2 = id2.station() + 8;
			}
			if (layerID == 1 && layerID2 == 2) binIndex = 1;
			if (layerID == 1 && layerID2 == 5) binIndex = 2;
			if (layerID == 1 && layerID2 == 9) binIndex = 3;
			if (layerID == 2 && layerID2 == 3) binIndex = 4;
			if (layerID == 2 && layerID2 == 5) binIndex = 5;
			if (layerID == 2 && layerID2 == 9) binIndex = 6;
			if (layerID == 5 && layerID2 == 6) binIndex = 7;
			if (layerID == 9 && layerID2 == 10) binIndex = 8;
			if (layerID == 3 && layerID2 == 4) binIndex = 9;
			if (layerID == 3 && layerID2 == 5) binIndex = 10;
			if (layerID == 3 && layerID2 == 9) binIndex = 11;
			if (layerID == 6 && layerID2 == 7) binIndex = 12;
			if (layerID == 10 && layerID2 == 11) binIndex = 13;
			if (layerID == 7 && layerID2 == 8) binIndex = 14;
			if (layerID == 11 && layerID2 == 12) binIndex = 15;
			if (layerID == 1 && layerID2 == 3) binIndex = 16;
			if (layerID == 2 && layerID2 == 4) binIndex = 17;
			if (layerID == 1 && layerID2 == 6) binIndex = 18;
			if (layerID == 1 && layerID2 == 10) binIndex = 19;
			if (layerID == 5 && layerID2 == 7) binIndex = 20;
			if (layerID == 6 && layerID2 == 8) binIndex = 21;
			if (layerID == 9 && layerID2 == 11) binIndex = 22;
			if (layerID == 10 && layerID2 == 12) binIndex = 23;
                                               
			effDenom->Fill(binIndex);
			if (binIndex == 10) std::cout << "found pair on 2" << std::endl;
		        bool matched = false;
			GlobalPoint gp = tthit->globalPosition();
			GlobalPoint gp2 = tthit2->globalPosition();
			double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
			double r2 = pow(gp2.x()*gp2.x() + gp2.y()*gp2.y(),0.5);	
			double dist = std::abs(gp.z()*r2 - gp2.z()*r)/(r2-r);
			//double dDir = gv.dot(gv2)/(gv.mag()*gv2.mag());
			double dPhi = std::min(std::abs(gp.phi() - gp2.phi()), std::abs(gp2.phi() - gp.phi()));
		
                        for (int i = 0; i < pairs->nPairs; i++){
				if (tthit->globalPosition().x() == pairs->gx1[i] && tthit->globalPosition().y() == pairs->gy1[i] && tthit2->globalPosition().x() == pairs->gx2[i] && tthit2->globalPosition().y() == pairs->gy2[i]) matched = true;
				if (tthit2->globalPosition().x() == pairs->gx1[i] && tthit2->globalPosition().y() == pairs->gy1[i] && tthit->globalPosition().x() == pairs->gx2[i] && tthit->globalPosition().y() == pairs->gy2[i]) matched = true;
			}
			if (matched) effNum->Fill(binIndex);
			if (!matched && binIndex == 10) std::cout << "missing pair with dR " << r2-r << " dist " << dist << " dPhi " << dPhi <<  std::endl;
		}
		j++;
			
	}
  }  







  nSegmentsPair1->Fill(nPair1);
  nSegmentsPair2->Fill(nPair2);
  nSegmentsPair3->Fill(nPair3);
  nSegmentsPair4->Fill(nPair4);
  nSegmentsPair5->Fill(nPair5);
  nSegmentsPair6->Fill(nPair6);
  nSegmentsPair7->Fill(nPair7);
  nSegmentsPair8->Fill(nPair8);
  nSegmentsPair9->Fill(nPair9);
  nSegmentsPair10->Fill(nPair10);
  nSegmentsPair11->Fill(nPair11);
  nSegmentsPair12->Fill(nPair12);
  nSegmentsPair13->Fill(nPair13);
  nSegmentsPair14->Fill(nPair14);
  nSegmentsPair15->Fill(nPair15);
  nSegmentsPair16->Fill(nPair16);
  nSegmentsPair17->Fill(nPair17);
  nSegmentsPair18->Fill(nPair18);
  nSegmentsPair19->Fill(nPair19);
  nSegmentsPair20->Fill(nPair20);
  nSegmentsPair21->Fill(nPair21);
  nSegmentsPair22->Fill(nPair22);
  nSegmentsPair23->Fill(nPair23);

  nSegmentsPairS1->Fill(nPairS1);
  nSegmentsPairS2->Fill(nPairS2);
  nSegmentsPairS3->Fill(nPairS3);
  nSegmentsPairS4->Fill(nPairS4);
  nSegmentsPairS5->Fill(nPairS5);
  nSegmentsPairS6->Fill(nPairS6);
  nSegmentsPairS7->Fill(nPairS7);
  nSegmentsPairS8->Fill(nPairS8);
  nSegmentsPairS9->Fill(nPairS9);
  nSegmentsPairS10->Fill(nPairS10);
  nSegmentsPairS11->Fill(nPairS11);
  nSegmentsPairS12->Fill(nPairS12);
  nSegmentsPairS13->Fill(nPairS13);
  nSegmentsPairS14->Fill(nPairS14);
  nSegmentsPairS15->Fill(nPairS15);
  nSegmentsPairS16->Fill(nPairS16);
  nSegmentsPairS17->Fill(nPairS17);
  nSegmentsPairS18->Fill(nPairS18);
  nSegmentsPairS19->Fill(nPairS19);
  nSegmentsPairS20->Fill(nPairS20);
  nSegmentsPairS21->Fill(nPairS21);
  nSegmentsPairS22->Fill(nPairS22);
  nSegmentsPairS23->Fill(nPairS23);


  nSegmentsPairG1->Fill(nPairG1);
  nSegmentsPairG2->Fill(nPairG2);
  nSegmentsPairG3->Fill(nPairG3);
  nSegmentsPairG4->Fill(nPairG4);
  nSegmentsPairG5->Fill(nPairG5);
  nSegmentsPairG6->Fill(nPairG6);
  nSegmentsPairG7->Fill(nPairG7);
  nSegmentsPairG8->Fill(nPairG8);
  nSegmentsPairG9->Fill(nPairG9);
  nSegmentsPairG10->Fill(nPairG10);
  nSegmentsPairG11->Fill(nPairG11);
  nSegmentsPairG12->Fill(nPairG12);
  nSegmentsPairG13->Fill(nPairG13);
  nSegmentsPairG14->Fill(nPairG14);
  nSegmentsPairG15->Fill(nPairG15);
  nSegmentsPairG16->Fill(nPairG16);
  nSegmentsPairG17->Fill(nPairG17);
  nSegmentsPairG18->Fill(nPairG18);
  nSegmentsPairG19->Fill(nPairG19);
  nSegmentsPairG20->Fill(nPairG20);
  nSegmentsPairG21->Fill(nPairG21);
  nSegmentsPairG22->Fill(nPairG22);
  nSegmentsPairG23->Fill(nPairG23);



}

void SegmentAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcDT", edm::InputTag("hltDt4DSegments"));
  desc.add<edm::InputTag>("srcCSC", edm::InputTag("hltCscSegments"));
  desc.add<edm::InputTag>("srcL2", edm::InputTag("hltL2Muons"));
  desc.add<edm::InputTag>("srcSegmentPairs", edm::InputTag("hltSegmentPairsSoA"));
  desc.add<std::string>("TrackerRecHitBuilder", "WithTrackAngle");
  descriptions.add("SegmentAnalyzer",desc);

}

DEFINE_FWK_MODULE(SegmentAnalyzer);
