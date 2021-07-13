//
//  Muon L1 prefiring study. It requires as input AOD or RECO events
//  that cannot have early trigger based on L1 trigger rules (no early
//  active bunches). In such unprefirable events we need to identify a
//  tag to avoid a selection bias. For SingleMuon PD such objects are
//  muons that triggered one of the triggers making the PD,
//  i.e. HLT_IsoMu24 or HLT_Mu50 for example. These muons are in time
//  by selection requirements. Any other muon in the event may have L1
//  muon object with wrong timing.
//
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/transform.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include <boost/regex.hpp>

#include "TTree.h"

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
// == reco::LorentzVector
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LorentzVector;

namespace {
  struct LumiInfo{
    Long64_t run;
    Long64_t lumi;
    int nEvents;
  };
  struct MuonInfo {
    Long64_t run;
    Long64_t lumi;
    Long64_t event;

    int n_vertex;

    LorentzVector probe_p4;
    int probe_q;
    
    float probe_dz, probe_dxy;

    bool probe_is_good_muon;
    bool probe_is_loose;
    bool probe_is_medium;
    bool probe_is_tight;
    bool probe_is_isolated;
    bool probe_is_isolated_tight;
    bool probe_is_pf_isolated_loose;
    bool probe_is_pf_isolated_tight;
    
    bool probe_isGlobal;
    bool probe_isPF;
    bool probe_isTracker;
    bool probe_normalizedChi2;
    bool probe_validMuonHits;
    bool probe_matchedStations;
    bool probe_validPixelHits;
    bool probe_trackerLayers;

    std::vector<LorentzVector> l1_p4;
    std::vector<int>           l1_quality;
    std::vector<int>           l1_bx = {-2,-1,0,1,2};

    std::vector<double>        l1_tf_pt;
    std::vector<double>        l1_tf_eta;
    std::vector<double>        l1_tf_phi;
    std::vector<int>           l1_tf_quality;
    std::vector<int>           l1_tf_bx = {-2,-1,0,1,2};


    void reset(){
      n_vertex = 0;
      run = lumi = event = 0;
      probe_is_tight = probe_is_isolated = false;
      probe_is_loose = probe_is_medium = false;
      probe_is_isolated_tight = probe_is_pf_isolated_tight = probe_is_pf_isolated_loose = false;
      probe_p4 = LorentzVector();
      probe_q = 0;
      probe_dxy = probe_dz = 0;
      l1_p4.assign(5,LorentzVector());
      l1_quality.assign(5,-1);
      l1_tf_pt.assign(5,-999.);
      l1_tf_eta.assign(5,-999.);
      l1_tf_phi.assign(5,-999.);
      l1_tf_quality.assign(5,-1);
 
    }
  };
  struct L1ObjInfo{
    int       bx;
    l1t::Muon obj;
    float     dR;
    float     dPhi;
    L1ObjInfo():bx(-999),dR(999),dPhi(999){
    }
  };
  struct L1TFObjInfo{
    int       bx;
    l1t::RegionalMuonCand obj;
    float     dR;
    float     dPhi;
    L1TFObjInfo():bx(-999),dR(999),dPhi(999){
    }

  };
}

class PrefiringMuAnaEGamma : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns,edm::one::WatchLuminosityBlocks>  {
public:
  explicit PrefiringMuAnaEGamma(const edm::ParameterSet&);
  ~PrefiringMuAnaEGamma();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
  void endLuminosityBlock(edm::LuminosityBlock const&, const edm::EventSetup&) override;
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void endRun(edm::Run const &, edm::EventSetup const &) override {}
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  bool isTagMuon(const reco::Muon& muon);
  bool isGoodMuon(const reco::Muon& muon);
  bool isLooseMuon(const reco::Muon& mu,const reco::Vertex& vtx);
  bool isMediumMuon(const reco::Muon& mu,const reco::Vertex& vtx);
  bool isTightMuon(const reco::Muon& mu,const reco::Vertex& vtx);
  bool isIsolatedMuon(const reco::Muon& muon);
  bool isIsolatedMuonTight(const reco::Muon& muon);
  bool isPFIsolatedMuonTight(const reco::Muon& muon);
  bool isPFIsolatedMuonLoose(const reco::Muon& muon);

  std::optional<GlobalPoint> 
  getMuonDirection(const reco::MuonChamberMatch& chamberMatch,
		   const edm::ESHandle<GlobalTrackingGeometry>& geometry,
		   const DetId& chamberId);
  std::vector<L1ObjInfo>
  findL1Objects(const reco::Muon& muon,
		const edm::ESHandle<GlobalTrackingGeometry>& geometry) ;
  std::vector<L1TFObjInfo>
  findL1ObjectsByTrackFinder(const reco::Muon& muon,
		const edm::ESHandle<GlobalTrackingGeometry>& geometry) ;


  void analyzeTrigger(const edm::Event& iEvent,
		      const edm::EventSetup& iSetup,
		      const std::string& triggerName);
  void trigger_init (const edm::TriggerNames & triggerNames);
  bool triggered(const reco::Muon& muon);

  // data members
  edm::EDGetTokenT<int> triggerRuleToken_;
  edm::EDGetTokenT<reco::MuonCollection>   muonToken_;
  edm::EDGetTokenT<BXVector<l1t::Muon>>    l1Token_;
  edm::EDGetTokenT<BXVector<l1t::RegionalMuonCand>>    l1BMTFToken_;
  edm::EDGetTokenT<BXVector<l1t::RegionalMuonCand>>    l1OMTFToken_;
  edm::EDGetTokenT<BXVector<l1t::RegionalMuonCand>>    l1EMTFToken_;
  edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;
  edm::EDGetTokenT<reco::VertexCollection>      vtx_token;

  edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTPrescaleProvider                hltPrescaleProvider_;
  edm::Handle<BXVector<l1t::Muon>>   l1Handle_;
  edm::Handle<BXVector<l1t::RegionalMuonCand>>   l1BMTFHandle_;
  edm::Handle<BXVector<l1t::RegionalMuonCand>>   l1OMTFHandle_;
  edm::Handle<BXVector<l1t::RegionalMuonCand>>   l1EMTFHandle_;

  edm::ParameterSetID                trigger_names_id_;
  std::vector<unsigned int>          trigger_ids_;
  const std::vector<std::string>     tag_hlt_trigger_patterns_ = {"HLT_IsoMu24","HLT_IsoMu27","HLT_HIMu15"};
  
  TTree* muonTree_;
  TTree* lumiTree_;
  MuonInfo muon_;
  LumiInfo lumi_;
};

void PrefiringMuAnaEGamma::trigger_init(const edm::TriggerNames & triggerNames){
  trigger_ids_.clear();
  for (unsigned int i = 0; i < triggerNames.size(); ++i) {
    std::string triggerName = triggerNames.triggerName(i);
    for (auto tag_trigger_name: tag_hlt_trigger_patterns_){
      if (boost::regex_match(triggerName,boost::regex("^"+tag_trigger_name+"_v\\d+$"))){
	trigger_ids_.push_back(i);
      }
    }
  }
}

bool PrefiringMuAnaEGamma::triggered(const reco::Muon& muon){
  // Loop over reference triggers and for each accepted trigger find
  // all objects passing the last filter and use them to match to
  // offline object
  HLTConfigProvider const& hltConfig = hltPrescaleProvider_.hltConfigProvider();

  const unsigned int n(hltConfig.size());

  for (auto triggerIndex : trigger_ids_){
    if (triggerIndex>=n){
      printf("hltConfig.size(): %u\n",n);
      throw cms::Exception("ConfigurationError") << "HLTPrescaleProvider is not properly initialized for HLTConfigProvider to work" << "\n";
    }

    if (not triggerResultsHandle_->accept(triggerIndex)) continue;

    // modules on this trigger path
    const unsigned int m(hltConfig.size(triggerIndex));
    const vector<string>& moduleLabels(hltConfig.moduleLabels(triggerIndex));

    const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
    assert(moduleIndex < m);

    // Results from TriggerEvent product - Attention: must look only for
    // modules actually run in this path for this event!
    trigger::TriggerObjectCollection firing_objects;
    for (unsigned int j = 0; j <= moduleIndex; ++j) {
      const string& moduleLabel(moduleLabels[j]);
      const string moduleType(hltConfig.moduleType(moduleLabel));
      // check whether the module is packed up in TriggerEvent product
      const unsigned int filterIndex(triggerEventHandle_->filterIndex(edm::InputTag(moduleLabel, "", "HLT")));
      if (filterIndex < triggerEventHandle_->sizeFilters()) {
	auto filterIds = triggerEventHandle_->filterIds(filterIndex);
	auto filterKeys = triggerEventHandle_->filterKeys(filterIndex);
	assert(filterIds.size() == filterKeys.size());
	unsigned int n = filterIds.size();
	auto triggerObjects = triggerEventHandle_->getObjects();
	firing_objects.clear();
	for (unsigned i = 0; i <= n; ++i) {
	  firing_objects.push_back(triggerObjects[filterKeys[i]]);
	}
      }
    }
    for (auto obj: firing_objects){
      if (deltaR2(muon,obj)<0.3) return true;
    }
  }
  return false;
}


PrefiringMuAnaEGamma::PrefiringMuAnaEGamma(const edm::ParameterSet& iConfig):
  triggerRuleToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("triggerRule"))),
  muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  l1Token_(consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("l1Src"))),
  l1BMTFToken_(consumes<BXVector<l1t::RegionalMuonCand>>(iConfig.getParameter<edm::InputTag>("l1BMTFSrc"))),
  l1OMTFToken_(consumes<BXVector<l1t::RegionalMuonCand>>(iConfig.getParameter<edm::InputTag>("l1OMTFSrc"))),
  l1EMTFToken_(consumes<BXVector<l1t::RegionalMuonCand>>(iConfig.getParameter<edm::InputTag>("l1EMTFSrc"))),
  l1GtToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1GtSrc"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  triggerEventToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEvents"))),
  vtx_token( consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices")) ),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  muonTree_ = fs->make<TTree>("muonTree","Muon Information");
  muonTree_->Branch("run",                   &muon_.run);
  muonTree_->Branch("lumi",                  &muon_.lumi);
  muonTree_->Branch("event",                 &muon_.event);
  muonTree_->Branch("n_vertex",              &muon_.n_vertex);
  muonTree_->Branch("probe_p4",              &muon_.probe_p4);
  muonTree_->Branch("probe_q",               &muon_.probe_q);
  muonTree_->Branch("probe_l1_p4",           &muon_.l1_p4);
  muonTree_->Branch("probe_l1_quality",      &muon_.l1_quality);
  muonTree_->Branch("probe_l1_bx",           &muon_.l1_bx);
  muonTree_->Branch("probe_l1_tf_pt",        &muon_.l1_tf_pt);
  muonTree_->Branch("probe_l1_tf_eta",       &muon_.l1_tf_eta);
  muonTree_->Branch("probe_l1_tf_phi",       &muon_.l1_tf_phi);
  muonTree_->Branch("probe_l1_tf_quality",   &muon_.l1_tf_quality);
  muonTree_->Branch("probe_l1_tf_bx",        &muon_.l1_tf_bx);
  muonTree_->Branch("probe_good_muon",       &muon_.probe_is_good_muon);
  muonTree_->Branch("probe_loose",           &muon_.probe_is_loose);
  muonTree_->Branch("probe_medium",          &muon_.probe_is_medium);
  muonTree_->Branch("probe_tight",           &muon_.probe_is_tight);
  muonTree_->Branch("probe_isolated",        &muon_.probe_is_isolated);
  muonTree_->Branch("probe_isolated_tight",  &muon_.probe_is_isolated_tight);
  muonTree_->Branch("probe_pf_isolated_tight",&muon_.probe_is_pf_isolated_tight);
  muonTree_->Branch("probe_pf_isolated_loose",&muon_.probe_is_pf_isolated_loose);
  muonTree_->Branch("probe_isGlobal",        &muon_.probe_isGlobal);
  muonTree_->Branch("probe_isPF",            &muon_.probe_isPF);
  muonTree_->Branch("probe_isTracker",       &muon_.probe_isTracker);
  muonTree_->Branch("probe_normalizedChi2",  &muon_.probe_normalizedChi2);
  muonTree_->Branch("probe_validMuonHits",   &muon_.probe_validMuonHits);
  muonTree_->Branch("probe_matchedStations", &muon_.probe_matchedStations);
  muonTree_->Branch("probe_validPixelHits",  &muon_.probe_validPixelHits);
  muonTree_->Branch("probe_trackerLayers",   &muon_.probe_trackerLayers);
  muonTree_->Branch("probe_dz",              &muon_.probe_dz);
  muonTree_->Branch("probe_dxy",             &muon_.probe_dxy);


  lumiTree_ = fs->make<TTree>("lumiTree","Lumi Information");
  lumiTree_->Branch("run",                   &lumi_.run);
  lumiTree_->Branch("lumi",                  &lumi_.lumi);
  lumiTree_->Branch("nevents",               &lumi_.nEvents);
}


void PrefiringMuAnaEGamma::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed(true);
  hltPrescaleProvider_.init(iRun, iSetup, "HLT", changed);
}

void PrefiringMuAnaEGamma::beginLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& setup)
{
  lumi_.run = lumiBlock.run();
  lumi_.lumi = lumiBlock.luminosityBlock();
  lumi_.nEvents = 0;
}

void PrefiringMuAnaEGamma::endLuminosityBlock(edm::LuminosityBlock const& lumiBlock, const edm::EventSetup& setup)
{
  lumiTree_->Fill();
}

PrefiringMuAnaEGamma::~PrefiringMuAnaEGamma()
{
}

bool PrefiringMuAnaEGamma::isTagMuon(const reco::Muon& muon){
  if (not isGoodMuon(muon)) return false;
  if (muon.pt()<30) return false;
  if (not triggered(muon)) return false;
  return true;
}

bool PrefiringMuAnaEGamma::isGoodMuon(const reco::Muon& muon){
  if ( not muon::isLooseMuon(muon) ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false;
  // if ( muon.pt() < 20 ) return false; 
  return true;
}

bool PrefiringMuAnaEGamma::isLooseMuon(const reco::Muon& muon,const reco::Vertex& vtx){
  if (not muon::isLooseMuon(muon)) return false;
  return true;
}
bool PrefiringMuAnaEGamma::isMediumMuon(const reco::Muon& muon,const reco::Vertex& vtx){
  if (not muon::isMediumMuon(muon)) return false;
  return true;
}
bool PrefiringMuAnaEGamma::isTightMuon(const reco::Muon& muon,const reco::Vertex& vtx){
  if (not muon::isTightMuon(muon,vtx)) return false;
  return true;
}
bool PrefiringMuAnaEGamma::isIsolatedMuon(const reco::Muon& muon){
  if ( muon.isolationR03().sumPt / muon.pt() > 0.1 ) return false;
  return true;
}
bool PrefiringMuAnaEGamma::isIsolatedMuonTight(const reco::Muon& muon){
  if ( muon.isolationR03().sumPt / muon.pt() > 0.05 ) return false;
  return true;
}
bool PrefiringMuAnaEGamma::isPFIsolatedMuonTight(const reco::Muon& muon){
  if ((muon.pfIsolationR04().sumChargedHadronPt + max(0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt))/muon.pt() > 0.15 ) return false;
  return true;
}
bool PrefiringMuAnaEGamma::isPFIsolatedMuonLoose(const reco::Muon& muon){
  if ((muon.pfIsolationR04().sumChargedHadronPt + max(0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt))/muon.pt() > 0.25 ) return false;
  return true;
}
std::optional<GlobalPoint> 
PrefiringMuAnaEGamma::getMuonDirection(const reco::MuonChamberMatch& chamberMatch,
				 const edm::ESHandle<GlobalTrackingGeometry>& geometry,
				 const DetId& chamberId) {
  const GeomDet* chamberGeometry = geometry->idToDet(chamberId);
  if (chamberGeometry) {
    LocalPoint localPosition(chamberMatch.x, chamberMatch.y, 0);
    return std::optional<GlobalPoint>(std::in_place, chamberGeometry->toGlobal(localPosition));
  }
  return std::optional<GlobalPoint>();
}

std::vector<L1ObjInfo>
PrefiringMuAnaEGamma::findL1Objects(const reco::Muon& muon,
			      const edm::ESHandle<GlobalTrackingGeometry>& geometry) 
{
  std::vector<L1ObjInfo> matches;

  // L1 trigger object parameters are defined at MB2/ME2. Use the muon
  // chamber matching information to get the local direction of the
  // muon trajectory and convert it to a global direction to match the
  // trigger objects

  std::optional<GlobalPoint> muonPosition;
  // Loop over chambers
  // initialize muonPosition with any available match, just in case
  // the second station is missing - it's better folling back to
  // dR matching at IP
  for (const auto& chamberMatch : muon.matches()) {
    if (chamberMatch.id.subdetId() == MuonSubdetId::DT) {
      DTChamberId detId(chamberMatch.id.rawId());
      if (abs(detId.station()) > 3)
	continue;
      muonPosition = getMuonDirection(chamberMatch, geometry, detId);
      if (abs(detId.station()) == 2)
	break;
    }
    if (chamberMatch.id.subdetId() == MuonSubdetId::CSC) {
      CSCDetId detId(chamberMatch.id.rawId());
      if (abs(detId.station()) > 3)
	continue;
      muonPosition = getMuonDirection(chamberMatch, geometry, detId);
      if (abs(detId.station()) == 2)
	break;
    }
  }
  if (not muonPosition) return matches;

  for (int ibx = l1Handle_->getFirstBX(); ibx <= l1Handle_->getLastBX(); ++ibx) {
    L1ObjInfo best_match;
    for (auto it = l1Handle_->begin(ibx); it != l1Handle_->end(ibx); it++){
      L1ObjInfo match;
      match.bx = ibx;
      match.obj = *it;
      // dPhi match is always possible, but dR is not
      match.dPhi = deltaPhi(it->phi(), muonPosition->phi());
      if (fabs(it->eta()) < 0.001) {
	// L1 is defined in X-Y plain
	if (match.dPhi > 0.1) continue;
      } else {
	// 3D L1
	match.dR = deltaR(it->p4(), *muonPosition);
	if (match.dR > 0.15) continue;
      }
      if (best_match.bx<-3){
	best_match = match;
      } else {
	if (best_match.dR<1 and match.dR<1){
	  if (best_match.dR > match.dR) best_match = match;
	} else {
	  if (best_match.dPhi > match.dPhi) best_match = match;
	}
      }
    }
    if (best_match.bx>=-3){
      matches.push_back(best_match);
    }
  }
  return matches;
}

std::vector<L1TFObjInfo>
PrefiringMuAnaEGamma::findL1ObjectsByTrackFinder(const reco::Muon& muon,
			      const edm::ESHandle<GlobalTrackingGeometry>& geometry) 
{
  std::vector<L1TFObjInfo> matches;

  // L1 trigger object parameters are defined at MB2/ME2. Use the muon
  // chamber matching information to get the local direction of the
  // muon trajectory and convert it to a global direction to match the
  // trigger objects

  std::optional<GlobalPoint> muonPosition;
  // Loop over chambers
  // initialize muonPosition with any available match, just in case
  // the second station is missing - it's better folling back to
  // dR matching at IP
  for (const auto& chamberMatch : muon.matches()) {
    if (chamberMatch.id.subdetId() == MuonSubdetId::DT) {
      DTChamberId detId(chamberMatch.id.rawId());
      if (abs(detId.station()) > 3)
	continue;
      muonPosition = getMuonDirection(chamberMatch, geometry, detId);
      if (abs(detId.station()) == 2)
	break;
    }
    if (chamberMatch.id.subdetId() == MuonSubdetId::CSC) {
      CSCDetId detId(chamberMatch.id.rawId());
      if (abs(detId.station()) > 3)
	continue;
      muonPosition = getMuonDirection(chamberMatch, geometry, detId);
      if (abs(detId.station()) == 2)
	break;
    }
  }
  if (not muonPosition) return matches;

  for (int ibx = l1BMTFHandle_->getFirstBX(); ibx <= l1BMTFHandle_->getLastBX(); ++ibx) {
    L1TFObjInfo best_match;
    for (auto it = l1BMTFHandle_->begin(ibx); it != l1BMTFHandle_->end(ibx); it++){
      L1TFObjInfo match;
      match.bx = ibx;
      match.obj = *it;
      // dPhi match is always possible, but dR is not
      //
      match.dPhi = deltaPhi( normalizedPhi(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor()) * (2 * M_PI / 576)) , muonPosition->phi());
      if (fabs(it->hwEta() * 0.010875) < 0.001) {
	// L1 is defined in X-Y plain
	if (match.dPhi > 0.1) continue;
      } else {
	// 3D L1
	//match.dR = deltaR(it->p4(), *muonPosition);
	match.dR = pow( pow(match.dPhi,2) + pow(it->hwEta() * 0.010875 - muonPosition->eta(),2),0.5);
	if (match.dR > 0.15) continue;
      }
      if (best_match.bx<-3){
	best_match = match;
      } else {
	if (best_match.dR<1 and match.dR<1){
	  if (best_match.dR > match.dR) best_match = match;
	} else {
	  if (best_match.dPhi > match.dPhi) best_match = match;
	}
      }
      matches.push_back(best_match);
    }
  }
  for (int ibx = l1OMTFHandle_->getFirstBX(); ibx <= l1OMTFHandle_->getLastBX(); ++ibx) {
    L1TFObjInfo best_match;
    for (auto it = l1OMTFHandle_->begin(ibx); it != l1OMTFHandle_->end(ibx); it++){
      L1TFObjInfo match;
      match.bx = ibx;
      match.obj = *it;
      // dPhi match is always possible, but dR is not
      match.dPhi = deltaPhi( normalizedPhi(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor()) * (2 * M_PI / 576)), muonPosition->phi());
      if (fabs(it->hwEta() * 0.010875) < 0.001) {
	// L1 is defined in X-Y plain
	if (match.dPhi > 0.1) continue;
      } else {
	// 3D L1
	//match.dR = deltaR(it->p4(), *muonPosition);
	match.dR = pow( pow(match.dPhi,2) + pow(it->hwEta() * 0.010875 - muonPosition->eta(),2),0.5);
	if (match.dR > 0.15) continue;
      }
      if (best_match.bx<-3){
	best_match = match;
      } else {
	if (best_match.dR<1 and match.dR<1){
	  if (best_match.dR > match.dR) best_match = match;
	} else {
	  if (best_match.dPhi > match.dPhi) best_match = match;
	}
      }
    }
    if (best_match.bx>=-3){
      matches.push_back(best_match);
    }
  }
  for (int ibx = l1EMTFHandle_->getFirstBX(); ibx <= l1EMTFHandle_->getLastBX(); ++ibx) {
    L1TFObjInfo best_match;
    for (auto it = l1EMTFHandle_->begin(ibx); it != l1EMTFHandle_->end(ibx); it++){
      L1TFObjInfo match;
      match.bx = ibx;
      match.obj = *it;
      // dPhi match is always possible, but dR is not
      match.dPhi = deltaPhi( normalizedPhi(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor()) * (2 * M_PI / 576)), muonPosition->phi());
      if (fabs(it->hwEta() * 0.010875) < 0.001) {
	// L1 is defined in X-Y plain
	if (match.dPhi > 0.1) continue;
      } else {
	// 3D L1
	//match.dR = deltaR(it->p4(), *muonPosition);
	match.dR = pow( pow(match.dPhi,2) + pow(it->hwEta() * 0.010875 - muonPosition->eta(),2),0.5);
	if (match.dR > 0.15) continue;
      }
      if (best_match.bx<-3){
	best_match = match;
      } else {
	if (best_match.dR<1 and match.dR<1){
	  if (best_match.dR > match.dR) best_match = match;
	} else {
	  if (best_match.dPhi > match.dPhi) best_match = match;
	}
      }
    }
    if (best_match.bx>=-3){
      matches.push_back(best_match);
    }
  }
  return matches;
}
void
PrefiringMuAnaEGamma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  lumi_.nEvents++;

  // get the tracking Geometry
  edm::ESHandle<GlobalTrackingGeometry> geometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(geometry);
  if (!geometry.isValid())
    throw cms::Exception("FatalError") << "Unable to find GlobalTrackingGeometryRecord in event!\n";

  edm::Handle<int> triggerRuleHandle;
  iEvent.getByToken(triggerRuleToken_, triggerRuleHandle);
  // event_.triggerRule = *triggerRuleHandle;

  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  iEvent.getByToken(l1Token_, l1Handle_);
  iEvent.getByToken(l1BMTFToken_, l1BMTFHandle_);
  iEvent.getByToken(l1OMTFToken_, l1OMTFHandle_);
  iEvent.getByToken(l1EMTFToken_, l1EMTFHandle_);

  edm::Handle<BXVector<GlobalAlgBlk>> l1GtHandle;
  iEvent.getByToken(l1GtToken_, l1GtHandle);

  iEvent.getByToken(triggerEventToken_, triggerEventHandle_);

  iEvent.getByToken(triggerResultsToken_, triggerResultsHandle_);

  const auto& triggerNames = iEvent.triggerNames(*triggerResultsHandle_);
  if (trigger_names_id_ != triggerNames.parameterSetID()) {
    trigger_names_id_ = triggerNames.parameterSetID();
    trigger_init(triggerNames);
  }

  // Primary vertex
  edm::Handle<reco::VertexCollection > vtx_handle;
  iEvent.getByToken(vtx_token,vtx_handle);
  const reco::Vertex& pv(*(vtx_handle->begin()));


  int vertex_count = 0;
  for (reco::VertexCollection::const_iterator it = vtx_handle->begin(), ite = vtx_handle->end(); it != ite; ++it) {
    //if (it->ndof() > 4 && fabs(it->z()) <= 24 && fabs(it->position().rho()) <= 2) {
      ++vertex_count;
    //}
  }


   int nL1Muons = 0;
   bool hasL1Above20 = false;
   for (int ibx = l1Handle_->getFirstBX(); ibx <= l1Handle_->getLastBX(); ++ibx) {
     for (auto it = l1Handle_->begin(ibx); it != l1Handle_->end(ibx); it++){
	if ((*it).pt() > 20) hasL1Above20 = true;
	nL1Muons++;
    }
   }



  // analyzeTrigger(iEvent,iSetup,triggerName);
  // consider all pairs
  if (!hasL1Above20 && nL1Muons ==1){
    for ( unsigned int i=0; i!=muonHandle->size(); ++i){
      const auto& probe_muon = muonHandle->at(i);
     
      muon_.reset();
      muon_.n_vertex = vertex_count;
      muon_.probe_p4 = probe_muon.p4();
      muon_.probe_q  = probe_muon.charge();
      muon_.probe_is_good_muon = isGoodMuon(probe_muon);
      muon_.probe_is_loose = isLooseMuon(probe_muon,pv);
      muon_.probe_is_medium = isMediumMuon(probe_muon,pv);
      muon_.probe_is_tight = isTightMuon(probe_muon,pv);
      muon_.probe_is_isolated = isIsolatedMuon(probe_muon);
      muon_.probe_is_isolated_tight = isIsolatedMuonTight(probe_muon);
      muon_.probe_is_pf_isolated_tight = isPFIsolatedMuonTight(probe_muon);
      muon_.probe_is_pf_isolated_loose = isPFIsolatedMuonLoose(probe_muon);
      muon_.probe_dxy = probe_muon.muonBestTrack()->dxy(pv.position());
      muon_.probe_dz = probe_muon.muonBestTrack()->dz(pv.position());
      muon_.probe_isGlobal = probe_muon.isGlobalMuon();
      muon_.probe_isPF = probe_muon.isPFMuon();
      muon_.probe_isTracker = probe_muon.isTrackerMuon();
      if (probe_muon.isGlobalMuon()){
      	muon_.probe_normalizedChi2 = probe_muon.globalTrack()->normalizedChi2();
	muon_.probe_validMuonHits = probe_muon.globalTrack()->hitPattern().numberOfValidMuonHits();
      }	
      else{
	muon_.probe_normalizedChi2 = -999.;
	muon_.probe_validMuonHits = -999.;
      }
      muon_.probe_matchedStations = probe_muon.numberOfMatchedStations();
      if (probe_muon.isTrackerMuon()){
	      muon_.probe_validPixelHits = probe_muon.innerTrack()->hitPattern().numberOfValidPixelHits();
	      muon_.probe_trackerLayers = probe_muon.innerTrack()->hitPattern().trackerLayersWithMeasurement();
      }
      else{
	muon_.probe_validPixelHits = 999.;
	muon_.probe_trackerLayers = 999.;

      }
      auto l1muons = findL1Objects(probe_muon,geometry);
      
      if (not l1muons.empty()){
	for ( const auto& l1: l1muons ){
		  unsigned int index = l1.bx+2;
	  muon_.l1_p4.at(index) = l1.obj.p4();
	  muon_.l1_quality.at(index) = l1.obj.hwQual();
	}
      }
      auto l1tfmuons = findL1ObjectsByTrackFinder(probe_muon,geometry);
      
      if (not l1tfmuons.empty()){
	for ( const auto& l1tf: l1tfmuons ){
		  unsigned int index = l1tf.bx+2;
	  muon_.l1_tf_pt.at(index) = (l1tf.obj.hwPt() - 1) * 0.5;
	  muon_.l1_tf_eta.at(index) = l1tf.obj.hwEta() * 0.010875;
	  muon_.l1_tf_phi.at(index) = normalizedPhi(l1t::MicroGMTConfiguration::calcGlobalPhi(l1tf.obj.hwPhi(), l1tf.obj.trackFinderType(), l1tf.obj.processor()) * (2 * M_PI / 576));
	  muon_.l1_tf_quality.at(index) = l1tf.obj.hwQual();
	}
      }

      muonTree_->Fill();

    }
  }
}

void PrefiringMuAnaEGamma::analyzeTrigger(const edm::Event& iEvent,
				    const edm::EventSetup& iSetup,
				    const std::string& triggerName) {
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  cout << endl;

  HLTConfigProvider const& hltConfig = hltPrescaleProvider_.hltConfigProvider();

  const unsigned int n(hltConfig.size());
  const unsigned int triggerIndex(hltConfig.triggerIndex(triggerName));

  if (triggerIndex != iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName)){
    printf("hltConfig.size(): %u\n",n);
    printf("triggerIndex (config): %u\n",triggerIndex);
    printf("triggerIndex (result): %u\n",iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));
    return;
  }

  // abort on invalid trigger name
  if (triggerIndex >= n) {
    cout << "HLTEventAnalyzerAOD::analyzeTrigger: path " << triggerName << " - not found!"
                                       << endl;
    return;
  }

  const std::pair<int, int> prescales(hltPrescaleProvider_.prescaleValues(iEvent, iSetup, triggerName));
  cout << "HLTEventAnalyzerAOD::analyzeTrigger: path " << triggerName << " ["
                                     << triggerIndex << "] "
                                     << "prescales L1T,HLT: " << prescales.first << "," << prescales.second << endl;
  const std::pair<std::vector<std::pair<std::string, int> >, int> prescalesInDetail(
      hltPrescaleProvider_.prescaleValuesInDetail(iEvent, iSetup, triggerName));
  std::ostringstream message;
  for (unsigned int i = 0; i < prescalesInDetail.first.size(); ++i) {
    message << " " << i << ":" << prescalesInDetail.first[i].first << "/" << prescalesInDetail.first[i].second;
  }
  cout << "HLTEventAnalyzerAOD::analyzeTrigger: path " << triggerName << " ["
                                     << triggerIndex << "] " << endl
                                     << "prescales L1T: " << prescalesInDetail.first.size() << message.str() << endl
                                     << " prescale HLT: " << prescalesInDetail.second << endl;

  // modules on this trigger path
  const unsigned int m(hltConfig.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig.moduleLabels(triggerIndex));

  // Results from TriggerResults product
  cout << " Trigger path status:"
                                     << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
                                     << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
                                     << " Error =" << triggerResultsHandle_->error(triggerIndex) << endl;
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  cout << " Last active module - label/type: " << moduleLabels[moduleIndex] << "/"
                                     << hltConfig.moduleType(moduleLabels[moduleIndex]) << " [" << moduleIndex
                                     << " out of 0-" << (m - 1) << " on this path]" << endl;
  assert(moduleIndex < m);

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  TriggerObjectCollection firing_objects;
  for (unsigned int j = 0; j <= moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string moduleType(hltConfig.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel, "", "HLT")));
    if (filterIndex < triggerEventHandle_->sizeFilters()) {
      cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/"
                                         << moduleType << endl;
      const Vids& VIDS(triggerEventHandle_->filterIds(filterIndex));
      const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI == nK);
      const size_type n(max(nI, nK));
      cout << "   " << n << " accepted 'L3' objects found: " << endl;
      const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
      firing_objects.clear();
      for (size_type i = 0; i != n; ++i) {
        const TriggerObject& TO(TOC[KEYS[i]]);
	firing_objects.push_back(TO);
        cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": " << TO.id() << " "
                                           << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass() << endl;
      }
    }
  }
  if (triggerResultsHandle_->accept(triggerIndex)){
    cout << "Objects that passed the last filter:" << endl;
    for (auto obj: firing_objects){
      cout << "\t" << obj.id() << " " << obj.pt() << " " << obj.eta() << " " << obj.phi() << " " << obj.mass() << endl;
    }
  }
  return;
}


void
PrefiringMuAnaEGamma::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(PrefiringMuAnaEGamma);
