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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>


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

  TH1F* nSegments;
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


  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

SegmentAnalyzer::SegmentAnalyzer(const edm::ParameterSet& iConfig)
    : dtSegmentsGetToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("srcDT"))),
      cscSegmentsGetToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("srcCSC"))),
      l2MuonGetToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("srcL2"))){

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

  nSegments    -> Fill( dtSegments.size() + cscSegments.size() );
  nSegmentsDT  -> Fill( dtSegments.size()  );
  nSegmentsCSC -> Fill( cscSegments.size() );

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
	bool signalSegment = false;
        /*for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {

	        for (auto recHit : (*itL2).recHits()){
			if ((recHit->globalPosition() - gp).mag() < 10e-5) signalSegment = true;
		}
        }*/
	if (signalSegment){
		if (stationID == 1) zSMB1 -> Fill(gp.z());
		if (stationID == 2) zSMB2 -> Fill(gp.z());
		if (stationID == 3) zSMB3 -> Fill(gp.z());
		if (stationID == 4) zSMB4 -> Fill(gp.z());

		if (stationID == 1) rSMB1 -> Fill(r);
		if (stationID == 2) rSMB2 -> Fill(r);
		if (stationID == 3) rSMB3 -> Fill(r);
		if (stationID == 4) rSMB4 -> Fill(r);
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

  }

   for (DTRecSegment4DCollection::const_iterator it = dtSegments.begin(); it != dtSegments.end(); it++) {
   	for (DTRecSegment4DCollection::const_iterator it2 = dtSegments.begin(); it2 != dtSegments.end(); it2++) {

		DTChamberId id = (DTChamberId)(*it).chamberId();
		GlobalPoint gp = dtGeom->chamber(id)->toGlobal((*it).localPosition());

		DTChamberId id2 = (DTChamberId)(*it2).chamberId();
		GlobalPoint gp2 = dtGeom->chamber(id2)->toGlobal((*it2).localPosition());
	 


		double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
		double r2 = pow(gp2.x()*gp2.x() + gp2.y()*gp2.y(),0.5);	
		double dist = std::abs(gp.z()*r2 - gp2.z()*r)/(r2-r);
		int stationID = (*it).chamberId().station();
		int stationID2 = (*it2).chamberId().station();

		if (stationID == 1 && stationID2 == 2){
			distPair1 -> Fill(dist);
			drPair1   -> Fill(r2-r);
		}
		if (stationID == 2 && stationID2 == 3){
			distPair4 -> Fill(dist);
			drPair4   -> Fill(r2-r);
		}
		if (stationID == 3 && stationID2 == 4){
			distPair9 -> Fill(dist);
			drPair9   -> Fill(r2-r);
		}
		if (stationID == 1 && stationID2 == 3){
			distPair16 -> Fill(dist);
			drPair16   -> Fill(r2-r);
		}
		if (stationID == 2 && stationID2 == 4){
			distPair17 -> Fill(dist);
			drPair17   -> Fill(r2-r);
		}
	}
  }
   for (DTRecSegment4DCollection::const_iterator it = dtSegments.begin(); it != dtSegments.end(); it++) {
   	for (CSCSegmentCollection::const_iterator it2 = cscSegments.begin(); it2 != cscSegments.end(); it2++) {

		DTChamberId id = (DTChamberId)(*it).chamberId();
		GlobalPoint gp = dtGeom->chamber(id)->toGlobal((*it).localPosition());

		CSCDetId id2 = (CSCDetId)(*it2).cscDetId();
		const CSCChamber* cscChamber = cscGeom->chamber(id2);
		GlobalPoint gp2 = cscChamber->toGlobal((*it2).localPosition());
		int stationID2 = -1;
		if (id2.zendcap() > 0) stationID2 = id2.station() + 4;
		else stationID2 = id2.station() + 8;


		double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
		double r2 = pow(gp2.x()*gp2.x() + gp2.y()*gp2.y(),0.5);	
		//double dist = gp.z()*r2 - gp2.z()*r;
		double dist = std::abs(gp.z()*r2 - gp2.z()*r)/(r2-r);
		int stationID = (*it).chamberId().station();

		if (stationID == 1 && stationID2 == 5){
			distPair2 -> Fill(dist);
			drPair2   -> Fill(r2-r);
		}
		if (stationID == 1 && stationID2 == 9){
			distPair3 -> Fill(dist);
			drPair3   -> Fill(r2-r);
		}
		if (stationID == 2 && stationID2 == 5){
			distPair5 -> Fill(dist);
			drPair5   -> Fill(r2-r);
		}
		if (stationID == 2 && stationID2 == 9){
			distPair6 -> Fill(dist);
			drPair6   -> Fill(r2-r);
		}
		if (stationID == 3 && stationID2 == 5){
			distPair10 -> Fill(dist);
			drPair10   -> Fill(r2-r);
		}
		if (stationID == 3 && stationID2 == 9){
			distPair11 -> Fill(dist);
			drPair11   -> Fill(r2-r);
		}
		if (stationID == 1 && stationID2 == 6){
			distPair18 -> Fill(dist);
			drPair18   -> Fill(r2-r);
		}
		if (stationID == 1 && stationID2 == 10){
			distPair19 -> Fill(dist);
			drPair19   -> Fill(r2-r);
		}
	}
  }
   for (CSCSegmentCollection::const_iterator it = cscSegments.begin(); it != cscSegments.end(); it++) {
   	for (CSCSegmentCollection::const_iterator it2 = cscSegments.begin(); it2 != cscSegments.end(); it2++) {

		CSCDetId id = (CSCDetId)(*it).cscDetId();
		const CSCChamber* cscChamber = cscGeom->chamber(id);
		GlobalPoint gp = cscChamber->toGlobal((*it).localPosition());
		int stationID = -1;
		if (id.zendcap() > 0) stationID = id.station() + 4;
		else stationID = id.station() + 8;

		CSCDetId id2 = (CSCDetId)(*it2).cscDetId();
		const CSCChamber* cscChamber2 = cscGeom->chamber(id2);
		GlobalPoint gp2 = cscChamber2->toGlobal((*it2).localPosition());
		int stationID2 = -1;
		if (id2.zendcap() > 0) stationID2 = id2.station() + 4;
		else stationID2 = id2.station() + 8;


		double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
		double r2 = pow(gp2.x()*gp2.x() + gp2.y()*gp2.y(),0.5);	
		//double dist = gp.z()*r2 - gp2.z()*r;
		double dist = std::abs(gp.z()*r2 - gp2.z()*r)/(r2-r);

		if (stationID == 5 && stationID2 == 6){
			distPair7 -> Fill(dist);
			drPair7   -> Fill(r2-r);
		}
		if (stationID == 9 && stationID2 == 10){
			distPair8 -> Fill(dist);
			drPair8   -> Fill(r2-r);
		}
		if (stationID == 6 && stationID2 == 7){
			distPair12 -> Fill(dist);
			drPair12   -> Fill(r2-r);
		}
		if (stationID == 10 && stationID2 == 11){
			distPair13 -> Fill(dist);
			drPair13   -> Fill(r2-r);
		}
		if (stationID == 7 && stationID2 == 8){
			distPair14 -> Fill(dist);
			drPair14   -> Fill(r2-r);
		}
		if (stationID == 11 && stationID2 == 12){
			distPair15 -> Fill(dist);
			drPair15   -> Fill(r2-r);
		}
		if (stationID == 5 && stationID2 == 7){
			distPair20 -> Fill(dist);
			drPair20   -> Fill(r2-r);
		}
		if (stationID == 6 && stationID2 == 8){
			distPair21 -> Fill(dist);
			drPair21   -> Fill(r2-r);
		}
		if (stationID == 9 && stationID2 == 11){
			distPair22 -> Fill(dist);
			drPair22   -> Fill(r2-r);
		}
		if (stationID == 10 && stationID2 == 12){
			distPair23 -> Fill(dist);
			drPair23   -> Fill(r2-r);
		}
	}
  }
}

void SegmentAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcDT", edm::InputTag("hltDt4DSegments"));
  desc.add<edm::InputTag>("srcCSC", edm::InputTag("hltCscSegments"));
  desc.add<edm::InputTag>("srcL2", edm::InputTag("hltL2Muons"));
  descriptions.add("SegmentAnalyzer",desc);

}

DEFINE_FWK_MODULE(SegmentAnalyzer);
