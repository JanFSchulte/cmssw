#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/UserData.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataMerger.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


  template <typename T>
  class IDResultEmbedder : public edm::stream::EDProducer<> {

    public:

      explicit IDResultEmbedder(const edm::ParameterSet & iConfig) :
            src_(consumes<std::vector<T>>(iConfig.getParameter<edm::InputTag>("src"))),
            vtx_src_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtx_src"))),
            rhoToken_ (consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))) {
            produces<std::vector<T>>();
        }

      ~IDResultEmbedder() override {}

      void produce(edm::Event & iEvent, const edm::EventSetup& iSetup) override;

      void setIDVariables(T &lep, reco::VertexCollection vertices, const double rho) const;

      static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
          edm::ParameterSetDescription desc;
          desc.add<edm::InputTag>("src");
          desc.add<edm::InputTag>("vtx_src");
          desc.add<edm::InputTag>("rho");
          if (typeid(T) == typeid(pat::Muon)) {
            descriptions.add("muonsWithIDResults", desc);
          }
          if (typeid(T) == typeid(pat::Electron)) {
            descriptions.add("electronsWithIDResults", desc);
          }

      }

    private:
      // configurables
      edm::EDGetToken src_;
      edm::EDGetToken vtx_src_;
      edm::EDGetToken rhoToken_;
  };


template <>
void IDResultEmbedder<pat::Electron>::setIDVariables(pat::Electron &anElectron, reco::VertexCollection vertices, const double rho) const {

      bool passEmHadIso2018 = false;
      bool passHOverE2018 = false;
      bool passEmHadIso = false;
      bool passHOverE = false;
      bool passShowerShape = false;
      bool passSieie = false;
      bool passEcalDriven = false;
      bool passDEtaIn = false;
      bool passDPhiIn = false;
      bool passTrackIso = false;
      bool passMissingHits = false;
      bool passDXY = false;


      double sc_e = anElectron.superCluster()->energy();
      if (anElectron.hasUserFloat("ecalEnergyPostCorr"))  sc_e = anElectron.userFloat("ecalEnergyPostCorr");
      double sc_et = anElectron.superCluster()->energy()*sin( anElectron.theta() ) ;
      if (anElectron.hasUserFloat("ecalEnergyPostCorr"))  sc_et = anElectron.userFloat("ecalEnergyPostCorr")*sin( anElectron.theta() ) ;

      if (fabs(anElectron.superCluster()->eta()) < 1.4442){
		if (anElectron.full5x5_e1x5()/anElectron.full5x5_e5x5() > 0.83 || anElectron.full5x5_e2x5Max()/anElectron.full5x5_e5x5() > 0.94) passShowerShape = true;
		passEmHadIso = (anElectron.dr03HcalDepth1TowerSumEt() + anElectron.dr03EcalRecHitSumEt()) < (2 + 0.03*sc_et + 0.28*rho);
		passEmHadIso2018 = (anElectron.dr03HcalDepth1TowerSumEt() + anElectron.dr03EcalRecHitSumEt()) < (2 + 0.03*sc_et + 0.28*rho);

		passHOverE = anElectron.hadronicOverEm() < (1./sc_e + 0.05);
		passHOverE2018 = anElectron.hadronicOverEm() < (1./sc_e + 0.05);

		passSieie = true;
                passEcalDriven = int(anElectron.ecalDrivenSeed());
		passDEtaIn = fabs(anElectron.deltaEtaSeedClusterTrackAtVtx()) < 0.004;
	        passDPhiIn = fabs(anElectron.deltaPhiSuperClusterTrackAtVtx()) < 0.06;

		passTrackIso = anElectron.dr03TkSumPtHEEP() < 5;
	        passMissingHits = anElectron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) < 2;
		passDXY = fabs(anElectron.gsfTrack()->dxy(vertices.at( 0 ).position())) < 0.02;
      }
      else if (fabs(anElectron.superCluster()->eta()) > 1.566 && fabs(anElectron.superCluster()->eta()) < 2.5 ){
		passShowerShape = true;
		passEmHadIso = (anElectron.dr03HcalDepth1TowerSumEt() + anElectron.dr03EcalRecHitSumEt()) < (2.5 + std::max(0.,0.03*(sc_et-50)) + 0.28*rho);
		passEmHadIso2018 = (anElectron.dr03HcalDepth1TowerSumEt() + anElectron.dr03EcalRecHitSumEt()) < (2.5 + std::max(0.,0.03*(sc_et-50)) + (0.15+0.07*fabs(anElectron.superCluster()->eta()))*rho);

		passHOverE = anElectron.hadronicOverEm() < (5./sc_e + 0.05);
		passHOverE2018 = anElectron.hadronicOverEm() < ((-0.4+0.4*fabs(anElectron.superCluster()->eta()))*rho/sc_e + 0.05);
		passSieie = anElectron.full5x5_sigmaIetaIeta() < 0.03;
                passEcalDriven = int(anElectron.ecalDrivenSeed());

		passDEtaIn = fabs(anElectron.deltaEtaSeedClusterTrackAtVtx()) < 0.006;
	        passDPhiIn = fabs(anElectron.deltaPhiSuperClusterTrackAtVtx()) < 0.06;
		passTrackIso = anElectron.dr03TkSumPtHEEP() < 5;
	        passMissingHits = anElectron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) < 2;
		passDXY = fabs(anElectron.gsfTrack()->dxy(vertices.at( 0 ).position())) < 0.05;
      }	

      const bool passID = (sc_et > 35.) && passEmHadIso && passHOverE && passShowerShape && passSieie && passEcalDriven && passDEtaIn && passDPhiIn && passTrackIso && passMissingHits && passDXY;
      const bool passID2018 = (sc_et > 35.) && passEmHadIso2018 && passHOverE2018 && passShowerShape && passSieie && passEcalDriven && passDEtaIn && passDPhiIn && passTrackIso && passMissingHits && passDXY;

      anElectron.addUserInt("passHEEPEmHadIso", passEmHadIso);
      anElectron.addUserInt("passHEEPEmHadIso2018", passEmHadIso2018);
      anElectron.addUserInt("passHEEPHOverE", passHOverE);
      anElectron.addUserInt("passHEEPHOverE2018", passHOverE2018);
      anElectron.addUserInt("passHEEPShowershape", passShowerShape);
      anElectron.addUserInt("passHEEPSieie", passSieie);
      anElectron.addUserInt("passHEEPEcalDriven", passEcalDriven);
      anElectron.addUserInt("passHEEPTrackIso", passTrackIso);
      anElectron.addUserInt("passHEEPMissingHits", passMissingHits);
      anElectron.addUserInt("passHEEPDXY", passDXY);
      anElectron.addUserInt("passHEEPDEta", passDEtaIn);
      anElectron.addUserInt("passHEEPDPhi", passDPhiIn);
      anElectron.addUserInt("passHEEPID", passID);
      anElectron.addUserInt("passHEEPID2018", passID2018);
}

template <>
void IDResultEmbedder<pat::Muon>::setIDVariables(pat::Muon &aMuon, reco::VertexCollection vertices, const double rho) const {

      bool passTrkIso = false;
      bool passTrackerLayers = false;
      bool passPixelHits = false;
      bool passValidMuonHits = false;
      bool passDPtOverPt = false;

      if (!(aMuon.innerTrack().isNull())) passTrkIso = (aMuon.isolationR03().sumPt / aMuon.innerTrack()->pt()) < 0.10;
      if (!(aMuon.globalTrack().isNull())){
           passTrackerLayers = aMuon.globalTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
           passPixelHits = aMuon.globalTrack()->hitPattern().numberOfValidPixelHits() >= 1;
      }
      if (!(aMuon.globalTrack().isNull()) && !(aMuon.tunePMuonBestTrack().isNull())) passValidMuonHits = ( (aMuon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0) || (aMuon.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits() > 0) );
      if (!(aMuon.tunePMuonBestTrack().isNull())) passDPtOverPt = aMuon.tunePMuonBestTrack()->ptError()/aMuon.tunePMuonBestTrack()->pt() < 0.3;

      bool passMatchedStations = (( aMuon.numberOfMatchedStations() > 1 ) || ( aMuon.numberOfMatchedStations() == 1 && ( aMuon.expectedNnumberOfMatchedStations() < 2 || !(aMuon.stationMask() ==1 || aMuon.stationMask() ==16) || aMuon.numberOfMatchedRPCLayers() > 2)));

      aMuon.addUserInt("passDXY", fabs(aMuon.dB()) < 0.2);
      aMuon.addUserInt("passTrkIso", passTrkIso);
      aMuon.addUserInt("passTrackerLayers", passTrackerLayers);
      aMuon.addUserInt("passPixelHits", passPixelHits);
      aMuon.addUserInt("passValidMuonHits", passValidMuonHits);
      aMuon.addUserInt("passMatchedStations", passMatchedStations);
      aMuon.addUserInt("passDPtOverPt", passDPtOverPt);

}

template <typename T>
void IDResultEmbedder<T>::produce(edm::Event & iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<std::vector<T>> src;
    iEvent.getByToken(src_, src);

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtx_src_, vertices);
    edm::Handle<double> rho_;
    iEvent.getByToken(rhoToken_,rho_);
    const double rho = *rho_;


    std::unique_ptr<std::vector<T>> out(new std::vector<T>(*src));

    for (unsigned int i = 0, n = src->size(); i < n; ++i) {
        T &lep = (*out)[i]; 
        setIDVariables(lep, *vertices.product(), rho);
        out->push_back(lep);
    }

    iEvent.put(std::move(out));
}

typedef IDResultEmbedder<pat::Electron> PATHEEPIDResultEmbedder;
typedef IDResultEmbedder<pat::Muon> PATHighPtMuonResultEmbedder;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATHEEPIDResultEmbedder);
DEFINE_FWK_MODULE(PATHighPtMuonResultEmbedder);
