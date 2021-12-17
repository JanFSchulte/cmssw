// -*- C++ -*-
//
// Package:    RecoMuon/SegmentAnalyzerForML
// Class:      SegmentAnalyzerForML
//
/**\class SegmentAnalyzerForML SegmentAnalyzerForML.cc RecoMuon/SegmentAnalyzer/plugins/SegmentAnalyzerForML.cc

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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
//
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class SegmentAnalyzerForML : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SegmentAnalyzerForML(const edm::ParameterSet&);
  ~SegmentAnalyzerForML();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesGetToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsGetToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsGetToken_;
  edm::EDGetTokenT<reco::TrackCollection> l2MuonGetToken_;


  struct tree_t {
	float gen_pt[100];
	float gen_eta[100];
	float gen_phi[100];
	int gen_charge[100];
	int nL2s;
	int nSegments[100];
	float l2_pt[100];
	float l2_eta[100];
	float l2_phi[100];
	int segment_L2ID[100];
	int segment_layerID[100];
	int segment_previousLayerID[100];
	float segment_globalX[100];
	float segment_globalY[100];
	float segment_globalZ[100];
	float segment_globalDX[100];
	float segment_globalDY[100];
	float segment_globalDZ[100];
	float segment_globalR[100];
	float segment_phi[100];
	float segment_deltaDir[100];
	float segment_deltaPhi[100];
	float segment_phiBend[100];
  };

  tree_t t;
  TTree* tree;

};

SegmentAnalyzerForML::SegmentAnalyzerForML(const edm::ParameterSet& iConfig)
    : genParticlesGetToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("srcGen"))),
      dtSegmentsGetToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("srcDT"))),
      cscSegmentsGetToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("srcCSC"))),
      l2MuonGetToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("srcL2")))
{
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");

  tree->Branch("gen_pt", &t.gen_pt, "gen_pt[100]/F");
  tree->Branch("gen_eta", &t.gen_eta, "gen_eta[100]/F");
  tree->Branch("gen_phi", &t.gen_phi, "gen_phi[100]/F");
  tree->Branch("gen_charge", &t.gen_charge, "gen_charge[100]/I");
  tree->Branch("nL2s", &t.nL2s, "nL2s/I");
  tree->Branch("nSegments", &t.nSegments, "nSegments[100]/I");
  tree->Branch("l2_pt", &t.l2_pt, "l2_pt[100]/F");
  tree->Branch("l2_eta", &t.l2_eta, "l2_eta[100]/F");
  tree->Branch("l2_phi", &t.l2_phi, "l2_phi[100]/F");
  tree->Branch("segment_L2ID", &t.segment_L2ID, "segment_L2ID[100]/I");
  tree->Branch("segment_layerID", &t.segment_layerID, "segment_layerID[100]/I");
  tree->Branch("segment_previousLayerID", &t.segment_previousLayerID, "segment_previousLayerID[100]/I");
  tree->Branch("segment_globalX", &t.segment_globalX, "segment_globalX[100]/F");
  tree->Branch("segment_globalY", &t.segment_globalY, "segment_globalY[100]/F");
  tree->Branch("segment_globalZ", &t.segment_globalZ, "segment_globalZ[100]/F");
  tree->Branch("segment_globalDX", &t.segment_globalDX, "segment_globalX[100]/F");
  tree->Branch("segment_globalDY", &t.segment_globalDY, "segment_globalY[100]/F");
  tree->Branch("segment_globalDZ", &t.segment_globalDZ, "segment_globalZ[100]/F");
  tree->Branch("segment_globalR", &t.segment_globalR, "segment_globalR[100]/F");
  tree->Branch("segment_phi", &t.segment_phi, "segment_phi[100]/F");
  tree->Branch("segment_deltaDir", &t.segment_deltaDir, "segment_deltaDir[100]/F");
  tree->Branch("segment_deltaPhi", &t.segment_deltaPhi, "segment_deltaPhi[100]/F");
  tree->Branch("segment_phiBend", &t.segment_phiBend, "segment_phiBend[100]/F");

}
SegmentAnalyzerForML::~SegmentAnalyzerForML() {
  // do anything here that needs to be done at desctruction time
  //   // (e.g. close files, deallocate resources etc.)
  //     //
  //       // please remove this method altogether if it would be left empty
}

void SegmentAnalyzerForML::beginJob() {
}
void SegmentAnalyzerForML::endJob() {}
void SegmentAnalyzerForML::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::ESHandle<DTGeometry> dtGeomHandle;
  iSetup.get<MuonGeometryRecord>().get(dtGeomHandle);
  const DTGeometry* dtGeom = &*dtGeomHandle;

  edm::ESHandle<CSCGeometry> cscGeomHandle;
  iSetup.get<MuonGeometryRecord>().get(cscGeomHandle);
  const CSCGeometry* cscGeom = &*cscGeomHandle;

  const DTRecSegment4DCollection& dtSegments = iEvent.get(dtSegmentsGetToken_);
  const CSCSegmentCollection& cscSegments = iEvent.get(cscSegmentsGetToken_);
  const reco::TrackCollection& l2Muons = iEvent.get(l2MuonGetToken_);
  const reco::GenParticleCollection& genParticles = iEvent.get(genParticlesGetToken_);

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);

   MuonTransientTrackingRecHitBuilder muonTransBuilder;

   for (int j = 0; j < 100; j++){
        t.gen_pt[j] = -999;
        t.gen_eta[j] = -999;
        t.gen_phi[j]= -999;
        t.gen_charge[j]= -999;
        t.nSegments[j]= -999;
        t.l2_pt[j]= -999;
        t.l2_eta[j]= -999;
        t.l2_phi[j]= -999;
        t.segment_L2ID[j]= -999;
        t.segment_layerID[j]= -999;
        t.segment_previousLayerID[j]= -999;
        t.segment_globalX[j]= -999;
        t.segment_globalY[j]= -999;
        t.segment_globalZ[j]= -999;
        t.segment_globalDX[j]= -999;
        t.segment_globalDY[j]= -999;
        t.segment_globalDZ[j]= -999;
        t.segment_globalR[j]= -999;
        t.segment_phi[j]= -999;
        t.segment_deltaDir[j]= -999;
        t.segment_deltaPhi[j]= -999;
        t.segment_phiBend[j]= -999;
   }


   int L2ID = 0;
   for (reco::TrackCollection::const_iterator itL2 = l2Muons.begin(); itL2 != l2Muons.end(); itL2++) {
	bool genMatched = false;
	for (reco::GenParticleCollection::const_iterator itGen = genParticles.begin(); itGen != genParticles.end(); itGen++) {

		if (deltaR(itL2->eta(),itL2->phi(),itGen->eta(),itGen->phi()) < 0.3 && std::abs(itGen->pdgId()) == 13 ){
			genMatched = true;
			t.gen_pt[L2ID] = itGen->pt();
			t.gen_eta[L2ID] = itGen->eta();
			t.gen_phi[L2ID] = itGen->phi();
			t.gen_charge[L2ID] = itGen->charge();
		}
		
	}
	if (!genMatched) continue;
	//if (!(abs(t.gen_eta[L2ID]) < 0.8)) continue; 
	t.l2_pt[L2ID] = itL2->pt();
	t.l2_eta[L2ID] = itL2->eta();
	t.l2_phi[L2ID] = itL2->phi();
	double previousPhi = -999;
	GlobalVector previousGv;
	int previousLayer = -1;
   	int nFound = 0;
        for (DTRecSegment4DCollection::const_iterator it = dtSegments.begin(); it != dtSegments.end(); it++) {
		DTChamberId id = (DTChamberId)(*it).chamberId();
		GlobalPoint gp = dtGeom->chamber(id)->toGlobal((*it).localPosition());
		GlobalVector gv = dtGeom->chamber(id)->toGlobal((*it).localDirection());
		double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);
		int layerID = (*it).chamberId().station();

		bool signalSegment = false;
		for (auto recHit : (*itL2).recHits()){
			if (!recHit->isValid()) continue;
			TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
			//TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*recHit);
			if ((tthit->globalPosition() - gp).mag() < 1e-5) signalSegment = true;
		}

		if (signalSegment){
			t.segment_L2ID[nFound] = L2ID;
			t.segment_layerID[nFound] = layerID-1;
			t.segment_globalX[nFound] = gp.x();
			t.segment_globalY[nFound] = gp.y();
			t.segment_globalZ[nFound] = gp.z();
			t.segment_globalDX[nFound] = gv.x();
			t.segment_globalDY[nFound] = gv.y();
			t.segment_globalDZ[nFound] = gv.z();
			t.segment_globalR[nFound] = r;
			t.segment_phi[nFound] = gp.phi().value();
			if (!( previousPhi == -999)){
				t.segment_deltaPhi[nFound-1] = deltaPhi(previousPhi,gp.phi().value());
				t.segment_deltaDir[nFound-1] = previousGv.dot(gv)/(gv.mag()*previousGv.mag());
				t.segment_previousLayerID[nFound] = previousLayer;
			}
			previousPhi = gp.phi().value();
			previousGv = gv;
			previousLayer = layerID-1;
			t.segment_phiBend[nFound] = -999;
			nFound++;

		}


	}	
	for (CSCSegmentCollection::const_iterator it = cscSegments.begin(); it != cscSegments.end(); it++) {

		CSCDetId id = (CSCDetId)(*it).cscDetId();
		const CSCChamber* cscChamber = cscGeom->chamber(id);
		GlobalPoint gp = cscChamber->toGlobal((*it).localPosition());
		GlobalVector gv = cscChamber->toGlobal((*it).localDirection());
		double r = pow(gp.x()*gp.x() + gp.y()*gp.y(),0.5);	
		int layerID = -1;
	        if (id.zendcap() > 0) layerID = id.station() + 4;
        	else layerID = id.station() + 8;

		bool signalSegment = false;
		for (auto recHit : (*itL2).recHits()){
			if (!recHit->isValid()) continue;
			TrajectoryMeasurement::ConstRecHitPointer tthit(muonTransBuilder.build(recHit, globalGeometry));
			//TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*recHit);
			if ((tthit->globalPosition() - gp).mag() < 1e-5) signalSegment = true;
		}
		if (signalSegment){
			t.segment_L2ID[nFound] = L2ID;
			t.segment_layerID[nFound] = layerID-1;
			t.segment_globalX[nFound] = gp.x();
			t.segment_globalY[nFound] = gp.y();
			t.segment_globalZ[nFound] = gp.z();
			t.segment_globalDX[nFound] = gv.x();
			t.segment_globalDY[nFound] = gv.y();
			t.segment_globalDZ[nFound] = gv.z();
			t.segment_globalR[nFound] = r;
			t.segment_phi[nFound] = gp.phi().value();
			if (! (previousPhi == -999)){
				t.segment_deltaDir[nFound-1] = previousGv.dot(gv)/(gv.mag()*previousGv.mag());
				t.segment_deltaPhi[nFound-1] = deltaPhi(previousPhi, gp.phi().value());
				t.segment_previousLayerID[nFound] = previousLayer;
			}
			previousGv = gv;
			previousPhi = gp.phi().value();
			previousLayer = layerID-1;
			//std::cout << "---------------------------------------------------------" << std::endl;
			const std::vector<CSCRecHit2D>&  recHits = (*it).specificRecHits();
  			
			GlobalPoint recHitGp1 = cscChamber->toGlobal(recHits.front().localPosition());
  			
			GlobalPoint recHitGp2 = cscChamber->toGlobal(recHits.back().localPosition());
			t.segment_phiBend[nFound] =  (pow(recHitGp1.x()*recHitGp1.x() + recHitGp1.y()*recHitGp1.y(),0.5) -  pow(recHitGp2.x()*recHitGp2.x() + recHitGp2.y()*recHitGp2.y(),0.5));
			//for (auto recHit : recHits){
			//	GlobalPoint recHitGp = cscChamber->toGlobal(recHit.localPosition());
				//std::cout << pow(recHitGp.x()*recHitGp.x() + recHitGp.y()*recHitGp.y(),0.5) << std::endl;		
			//}
			nFound++;

		}

	}
   	t.nSegments[L2ID]=nFound;
	L2ID++;
   }
   t.nL2s=L2ID;
   tree->Fill();

}  




void SegmentAnalyzerForML::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcGen", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("srcDT", edm::InputTag("hltDt4DSegments"));
  desc.add<edm::InputTag>("srcCSC", edm::InputTag("hltCscSegments"));
  desc.add<edm::InputTag>("srcL2", edm::InputTag("hltL2Muons"));
  desc.add<std::string>("TrackerRecHitBuilder", "WithTrackAngle");
  descriptions.add("SegmentAnalyzerForML",desc);

}

DEFINE_FWK_MODULE(SegmentAnalyzerForML);
