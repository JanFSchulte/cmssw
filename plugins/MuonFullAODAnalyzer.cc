//======================================================//
//                                                      //
//                                                      //
//                MUON POG TnP NTUPLIZER                //
//                                                      //
//                                                      //
//======================================================//

// Package:    MuonAnalyzer for Run 3
//             version 2.0
//
/**\class

 Description: Ntuplizer class for full AOD files
*/
//
// Original Author:
//                george karathanasis
//         Created:  Thu, 20 feb 2020 17:40:23 GMT
//
// Modified:
//                Andre Frankenthal (Sept. 2020)
//
// Modified:
//                Minseok Oh (Feb. 2021)
//
//

// system include files
#include <iostream>
#include <memory>
#include <random>
#include <queue>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <vector>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "CondFormats/DataRecord/interface/JetResolutionRcd.h"
#include "CondFormats/DataRecord/interface/JetResolutionScaleFactorRcd.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"

#include "KlFitter.h"
#include "MuonBranches.h"
#include "MuonGenAnalyzer.h"
#include "NtupleContent.h"
#include "helper.h"
#include "MuonMiniIsolation.h"
#include "JetsBranches.h"

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class MuonFullAODAnalyzer : public edm::one::EDAnalyzer<> {
public:
  typedef std::vector<std::pair<reco::Muon, reco::TransientTrack>> RecoMuonAndTransientTrkCollection;
  typedef std::pair<reco::Muon, reco::TransientTrack> RecoMuonAndTransientTrk;
  typedef std::pair<reco::Track, reco::TransientTrack> RecoTrkAndTransientTrk;
  explicit MuonFullAODAnalyzer(const edm::ParameterSet&);
  ~MuonFullAODAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;

  bool HLTaccept(const edm::Event&, NtupleContent&, std::vector<std::string>&);
  void fillHLTmuon(const edm::Event&,
                   std::vector<TString>&,
                   std::vector<float>&,
                   std::vector<float>&,
                   std::vector<float>&,
                   std::vector<std::string>&,
                   const int&);
  void embedTriggerMatching(const reco::Muon&,
                            std::vector<TString>&,
                            std::vector<float>&,
                            std::vector<float>&,
                            std::vector<float>&,
                            std::vector<std::string>&,
                            bool,
                            const int&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupSummaryToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetTokenT<edm::View<reco::Muon>> muonsViewToken_;
  edm::EDGetToken tracksToken_;
  edm::EDGetToken dSAToken_;
  edm::EDGetToken dglToken_;
  edm::EDGetToken cosmicToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneMatch> l1MatchesToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> l1MatchesQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> l1MatchesDeltaRToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneMatch> l1MatchesByQToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> l1MatchesByQQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> l1MatchesByQDeltaRToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genToken_;
  edm::EDGetToken PATPFCands_;
  edm::EDGetTokenT<double> rhoJetsNC_;
  edm::EDGetToken jetsToken_;
  edm::EDGetToken jetCorrectorToken_;
  edm::EDGetToken genJetsToken_;
  edm::EDGetToken deepCSVProbbToken_;
  edm::EDGetToken deepCSVProbbbToken_;
  //  edm::EDGetToken deepFlavProbbToken_;
  //  edm::EDGetToken deepFlavProbbbToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonSimInfo>> simInfoToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
  // Jet resolution corrections
  JME::JetResolution::Token t_jet_resolution_token;
  JME::JetResolutionScaleFactor::Token t_jet_resolutionSF_token;

  std::vector<std::string> HLTPaths_;      // trigger fired
  std::vector<std::string> tagFilters_;    // tag-trigger matching
  std::vector<std::string> probeFilters_;  // probe-trigger matching
  std::vector<std::string> probeSelectorNames_;
  std::vector<unsigned> probeSelectorBits_;

  std::mt19937 m_random_generator = std::mt19937(37428479);
  const bool isMC_, includeJets_;
  const std::string era_;
  const double trgDRwindow_;
  const unsigned int tagQual_;
  const StringCutObjectSelector<reco::Muon> tagSelection_;  // kinematic cuts for tag
  const bool HighPurity_;
  const StringCutObjectSelector<reco::Track> probeSelection_;  // kinematic cuts for probe
  const bool muonOnly_;
  const StringCutObjectSelector<reco::Muon> probeMuonSelection_;
  const double pairMassMin_;
  const double pairMassMax_;
  const double pairDz_;
  const bool RequireVtxCreation_;  // if true skip pairs that do not create
                                   // that do not have a vertex
  const double minSVtxProb_;       // min probability of a vertex to be kept. If <0 inactive
  const double maxdz_trk_mu_;
  const double maxpt_relative_dif_trk_mu_;
  const double maxdr_trk_mu_;
  const double maxdr_trk_dsa_;
  const unsigned momPdgId_;
  const double genRecoDrMatch_;
  const int debug_;
  PropagateToMuon prop1_;
  PropagateToMuonSetup propSetup1_;

  edm::Service<TFileService> fs;
  TTree* t1;
  NtupleContent nt;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonFullAODAnalyzer::MuonFullAODAnalyzer(const edm::ParameterSet& iConfig)
    :  // inputs
      genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
      rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho"))),
      pileupSummaryToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfo"))),
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      muonsViewToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      dSAToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dSAmuons"))),
      dglToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dGlmuons"))),
      cosmicToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("staCosmic"))),
      trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      trigobjectsToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
      l1MatchesToken_(consumes<pat::TriggerObjectStandAloneMatch>(iConfig.getParameter<edm::InputTag>("l1Matches"))),
      l1MatchesQualityToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("l1MatchesQuality"))),
      l1MatchesDeltaRToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("l1MatchesDeltaR"))),
      l1MatchesByQToken_(
          consumes<pat::TriggerObjectStandAloneMatch>(iConfig.getParameter<edm::InputTag>("l1MatchesByQ"))),
      l1MatchesByQQualityToken_(
          consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("l1MatchesByQQuality"))),
      l1MatchesByQDeltaRToken_(
          consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("l1MatchesByQDeltaR"))),
      genToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen"))),
      PATPFCands_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("PATPFCands"))),
      rhoJetsNC_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoJetsNC"))),
      jetsToken_(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      jetCorrectorToken_(consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("jetCorrector"))),
      genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
      deepCSVProbbToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("deepCSVProbb"))),
      deepCSVProbbbToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("deepCSVProbbb"))),
      simInfoToken_(consumes<edm::ValueMap<reco::MuonSimInfo>>(iConfig.getParameter<edm::InputTag>("muonSimInfo"))),
      HLTPaths_(iConfig.getParameter<std::vector<std::string>>("triggerPaths")),
      tagFilters_(iConfig.getParameter<std::vector<std::string>>("tagFilters")),
      probeFilters_(iConfig.getParameter<std::vector<std::string>>("probeFilters")),
      probeSelectorNames_(iConfig.getParameter<std::vector<std::string>>("probeSelectorNames")),
      probeSelectorBits_(iConfig.getParameter<std::vector<unsigned>>("probeSelectorBits")),
      isMC_(iConfig.getParameter<bool>("isMC")),
      includeJets_(iConfig.getParameter<bool>("includeJets")),
      era_(iConfig.getParameter<std::string>("era")),
      trgDRwindow_(iConfig.getParameter<double>("trgDRwindow")),
      tagQual_(iConfig.getParameter<unsigned>("tagQuality")),
      tagSelection_(iConfig.getParameter<std::string>("tagSelection")),
      HighPurity_(iConfig.getParameter<bool>("probeHPurity")),
      probeSelection_(iConfig.getParameter<std::string>("probeSelection")),
      muonOnly_(iConfig.getParameter<bool>("muonOnly")),
      probeMuonSelection_(iConfig.getParameter<std::string>("probeMuonSelection")),
      pairMassMin_(iConfig.getParameter<double>("pairMassMin")),
      pairMassMax_(iConfig.getParameter<double>("pairMassMax")),
      pairDz_(iConfig.getParameter<double>("pairDz")),
      RequireVtxCreation_(iConfig.getParameter<bool>("RequireVtxCreation")),
      minSVtxProb_(iConfig.getParameter<double>("minSVtxProb")),
      maxdz_trk_mu_(iConfig.getParameter<double>("maxDzProbeTrkMuon")),
      maxpt_relative_dif_trk_mu_(iConfig.getParameter<double>("maxRelPtProbeTrkMuon")),
      maxdr_trk_mu_(iConfig.getParameter<double>("maxDRProbeTrkMuon")),
      maxdr_trk_dsa_(iConfig.getParameter<double>("maxDRProbeTrkDSA")),
      momPdgId_(iConfig.getParameter<unsigned>("momPdgId")),
      genRecoDrMatch_(iConfig.getParameter<double>("genRecoDrMatch")),
      debug_(iConfig.getParameter<int>("debug")),
      propSetup1_(iConfig, consumesCollector()) {
  edm::ConsumesCollector iC = consumesCollector();
  magfieldToken_ = iC.esConsumes<MagneticField, IdealMagneticFieldRecord>();
  t_jet_resolution_token = esConsumes(edm::ESInputTag("", "AK4PFPuppi_pt"));  // Load jet resolution files
  t_jet_resolutionSF_token = esConsumes(edm::ESInputTag("", "AK4PFPuppi"));

  if (probeSelectorNames_.size() != probeSelectorBits_.size()) {
    throw cms::Exception("ParameterError")
        << "length of probeSelectorNames and probeSelectorBits should be identical\n";
  }
}

MuonFullAODAnalyzer::~MuonFullAODAnalyzer() {}

// MEMBER FUNCTIONS

//
//   ---- HLTaccept ----
//
//   Function that takes as input the HLT paths
//   It check if they are fired and save the information in the ntuple
//

bool MuonFullAODAnalyzer::HLTaccept(const edm::Event& iEvent, NtupleContent& nt, std::vector<std::string>& HLTPaths) {
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  edm::TriggerNames trigName;
  trigName = iEvent.triggerNames(*trigResults);
  bool EvtFire = false;
  unsigned int ipath = 0;
  for (auto path : HLTPaths) {
    bool TrgFire = false;
    for (unsigned int itrg = 0; itrg < trigResults->size(); ++itrg) {
      TString TrigPath = trigName.triggerName(itrg);
      if (!trigResults->accept(itrg))
        continue;
      if (!TrigPath.Contains(path))
        continue;
      EvtFire = true;
      TrgFire = true;
    }
    nt.trigger[ipath] = TrgFire;
    ipath++;
  }
  return EvtFire;
}

//
// ---- fillHLTmuon ----
//
// Fill branches for HLT muon objects
//

void MuonFullAODAnalyzer::fillHLTmuon(const edm::Event& iEvent,
                                      std::vector<TString>& trg_filter,
                                      std::vector<float>& trg_pt,
                                      std::vector<float>& trg_eta,
                                      std::vector<float>& trg_phi,
                                      std::vector<std::string>& HLTFilters,
                                      const int& debug_) {
  edm::Handle<trigger::TriggerEvent> triggerObjects;
  iEvent.getByToken(trigobjectsToken_, triggerObjects);
  trigger::TriggerObjectCollection allTriggerObjects = triggerObjects->getObjects();
  for (auto ifilter : HLTFilters) {
    size_t filterIndex = (*triggerObjects).filterIndex(edm::InputTag(ifilter, "", "HLT"));
    if (filterIndex < (*triggerObjects).sizeFilters()) {
      const trigger::Keys& keys = (*triggerObjects).filterKeys(filterIndex);
      for (size_t j = 0; j < keys.size(); j++) {
        trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
        if (fabs(foundObject.id()) != 13)
          continue;
        trg_filter.push_back(TString(ifilter));
        trg_pt.push_back(foundObject.pt());
        trg_eta.push_back(foundObject.eta());
        trg_phi.push_back(foundObject.phi());
        if (debug_ > 0) {
          std::cout << "Trg muon " << foundObject.pt() << std::endl;
        }
      }
    }
  }
}

//
// ---- embedTriggerMatching ----
//
// Function that takes as input the muon collection and trigger objects
// It performs a matching between muon and trigger objects
//

void MuonFullAODAnalyzer::embedTriggerMatching(const reco::Muon& mu,
                                               std::vector<TString>& trg_filter,
                                               std::vector<float>& trg_pt,
                                               std::vector<float>& trg_eta,
                                               std::vector<float>& trg_phi,
                                               std::vector<std::string>& HLTFilters,
                                               bool isTag,
                                               const int& debug_ = 0) {
  for (const auto& thefilter : HLTFilters) {
    TString thefilter_tstr = TString(thefilter);
    // temporary method to tag L2 filters for dSA paths
    bool isL2DSA =
        thefilter_tstr.BeginsWith("hltL2") && (thefilter_tstr.Contains("NoVtx") || thefilter_tstr.Contains("NoVertex"));

    bool matched = false;
    float matched_pt = -99;
    float matched_eta = -99;
    float matched_phi = -99;
    float matched_dr = 99;
    for (unsigned itrg = 0; itrg < trg_filter.size(); ++itrg) {
      TString filter_tstr = TString(trg_filter.at(itrg));
      if (!filter_tstr.Contains(thefilter_tstr))
        continue;
      float dR_tmp = deltaR(mu.eta(), mu.phi(), trg_eta.at(itrg), trg_phi.at(itrg));
      if ((dR_tmp < matched_dr && dR_tmp < trgDRwindow_) ||
          (isL2DSA && dR_tmp < matched_dr && dR_tmp < maxdr_trk_dsa_)) {
        matched = true;
        matched_pt = trg_pt.at(itrg);
        matched_eta = trg_eta.at(itrg);
        matched_phi = trg_phi.at(itrg);
        matched_dr = dR_tmp;

        if (debug_ > 0) {
          std::cout << "embedTriggerMatching: isTag=" << isTag << "  filter=" << thefilter_tstr << "  dR=" << dR_tmp
                    << "  matched=" << matched << std::endl;
        }
      }
    }

    if (isTag) {
      nt.tag_trg[&thefilter - &HLTFilters[0]] = matched;
      nt.tag_trg_pt[&thefilter - &HLTFilters[0]] = matched_pt;
      nt.tag_trg_eta[&thefilter - &HLTFilters[0]] = matched_eta;
      nt.tag_trg_phi[&thefilter - &HLTFilters[0]] = matched_phi;
      nt.tag_trg_dr[&thefilter - &HLTFilters[0]] = matched_dr;
    } else {
      nt.probe_trg[&thefilter - &HLTFilters[0]] = matched;
      nt.probe_trg_pt[&thefilter - &HLTFilters[0]] = matched_pt;
      nt.probe_trg_eta[&thefilter - &HLTFilters[0]] = matched_eta;
      nt.probe_trg_phi[&thefilter - &HLTFilters[0]] = matched_phi;
      nt.probe_trg_dr[&thefilter - &HLTFilters[0]] = matched_dr;
    }
  }

  return;
}

// ------------------------------------
// ------------ TnP method ------------
// ------------------------------------

void MuonFullAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  prop1_ = propSetup1_.init(iSetup);

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_, theBeamSpot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  // Skip evts if there are no vertices
  if (vertices->empty())
    return;

  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  edm::Handle<edm::View<reco::Muon>> muonsView;
  iEvent.getByToken(muonsViewToken_, muonsView);
  edm::Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  edm::Handle<std::vector<reco::Track>> dSAmuons;
  iEvent.getByToken(dSAToken_, dSAmuons);
  edm::Handle<std::vector<reco::Track>> dGlmuons;
  iEvent.getByToken(dglToken_, dGlmuons);
  edm::Handle<std::vector<reco::Track>> staCosmic;
  iEvent.getByToken(cosmicToken_, staCosmic);

  edm::ESHandle<MagneticField> bField;
  bField = iSetup.getHandle(magfieldToken_);

  // mini isolation
  edm::Handle<pat::PackedCandidateCollection> patpfcands;
  iEvent.getByToken(PATPFCands_, patpfcands);
  edm::Handle<double> rhoJetsNC;
  iEvent.getByToken(rhoJetsNC_, rhoJetsNC);
  // jets
  edm::Handle<std::vector<reco::PFJet>> jets;
  edm::Handle<reco::JetCorrector> jetCorrector;
  edm::Handle<reco::JetTagCollection> deepCSVProbb;
  edm::Handle<reco::JetTagCollection> deepCSVProbbb;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor resolution_sf;
  if (includeJets_) {
    resolution = JME::JetResolution::get(iSetup, t_jet_resolution_token);
    resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, t_jet_resolutionSF_token);
  }

  edm::Handle<trigger::TriggerEvent> triggerObjects;
  iEvent.getByToken(trigobjectsToken_, triggerObjects);
  edm::Handle<pat::TriggerObjectStandAloneMatch> l1Matches;
  iEvent.getByToken(l1MatchesToken_, l1Matches);
  edm::Handle<edm::ValueMap<int>> l1Qualities;
  iEvent.getByToken(l1MatchesQualityToken_, l1Qualities);
  edm::Handle<edm::ValueMap<float>> l1Drs;
  iEvent.getByToken(l1MatchesDeltaRToken_, l1Drs);
  edm::Handle<pat::TriggerObjectStandAloneMatch> l1MatchesByQ;
  iEvent.getByToken(l1MatchesByQToken_, l1MatchesByQ);
  edm::Handle<edm::ValueMap<int>> l1QualitiesByQ;
  iEvent.getByToken(l1MatchesByQQualityToken_, l1QualitiesByQ);
  edm::Handle<edm::ValueMap<float>> l1DrsByQ;
  iEvent.getByToken(l1MatchesByQDeltaRToken_, l1DrsByQ);

  // Information about run -----------------------------------
  nt.ClearBranches();
  nt.branches["lumi"] = (int)iEvent.luminosityBlock();
  nt.branches["run"] = (int)iEvent.id().run();
  nt.branches["event"] = (int)iEvent.id().event();
  nt.branches["fromFullAOD"] = (bool)true;
  nt.branches["BSpot_x"] = (float)theBeamSpot->x0();
  nt.branches["BSpot_y"] = (float)theBeamSpot->y0();
  nt.branches["BSpot_z"] = (float)theBeamSpot->z0();
  nt.branches["nVertices"] = (int)vertices->size();

  // Gen weights, sim info -----------------------------------
  bool simInfoIsAvailalbe = false;
  edm::Handle<edm::ValueMap<reco::MuonSimInfo>> simInfo;
  if (!iEvent.isRealData()) {
    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);
    nt.branches["genWeight"] = (float)genEventInfoHandle->weight();

    simInfoIsAvailalbe = iEvent.getByToken(simInfoToken_, simInfo);
  } else {
    nt.branches["genWeight"] = (float)1.;
  }

  // Pileup information -----------------------------------
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  nt.branches["rho"] = (double)*rhoHandle;

  float trueNumInteractions = -1;
  int puNumInteractions = -1;
  if (isMC_) {
    edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);

    for (auto PVI : *PupInfo) {
      int BX = PVI.getBunchCrossing();
      if (BX == 0) {
        trueNumInteractions = PVI.getTrueNumInteractions();
        puNumInteractions = PVI.getPU_NumInteractions();
        continue;
      }
    }
  }

  nt.branches["nTrueInteractions"] = (float)trueNumInteractions;
  nt.branches["nPUInteractions"] = (float)puNumInteractions;

  if (debug_ > 0)
    std::cout << "New Evt " << std::get<int>(nt.branches["run"].value) << std::endl;

  // Vertex reconstruction -----------------------------------
  reco::TrackBase::Point vertex_point;
  bool goodVtx = false;
  reco::Vertex const* pv;
  for (const reco::Vertex& vtx : *vertices) {
    if (vtx.isFake() || !vtx.isValid())
      continue;
    nt.branches["pv_x"] = (float)vtx.x();
    nt.branches["pv_y"] = (float)vtx.y();
    nt.branches["pv_z"] = (float)vtx.z();
    goodVtx = true;
    pv = &vtx;
    break;
  }
  if (!goodVtx)
    return;  // skipping in absence of good vertex
  vertex_point.SetCoordinates(std::get<float>(nt.branches["pv_x"].value),
                              std::get<float>(nt.branches["pv_y"].value),
                              std::get<float>(nt.branches["pv_z"].value));

  // check if path fired, if so save hlt muons -----------------------------------
  if (!HLTaccept(iEvent, nt, HLTPaths_))
    return;

  fillHLTmuon(iEvent, nt.trg_filter, nt.trg_pt, nt.trg_eta, nt.trg_phi, tagFilters_, debug_);
  fillHLTmuon(iEvent, nt.prb_filter, nt.prb_pt, nt.prb_eta, nt.prb_phi, probeFilters_, debug_);

  // gen information and matching -----------------------------------
  MuonGenAnalyzer genmu;
  std::vector<unsigned> matched_muon_idx;
  std::vector<unsigned> matched_track_idx;
  if (!iEvent.isRealData()) {
    genmu.SetInputs(iEvent, genToken_, momPdgId_);
    genmu.FillNtuple(nt);

    auto reco_match_genmu1 = MatchReco<reco::Muon>(*muons,
                                                   std::get<float>(nt.branches["genmu1_eta"].value),
                                                   std::get<float>(nt.branches["genmu1_phi"].value),
                                                   std::get<int>(nt.branches["genmu1_charge"].value),
                                                   genRecoDrMatch_);
    auto reco_match_genmu2 = MatchReco<reco::Muon>(*muons,
                                                   std::get<float>(nt.branches["genmu2_eta"].value),
                                                   std::get<float>(nt.branches["genmu2_phi"].value),
                                                   std::get<int>(nt.branches["genmu2_charge"].value),
                                                   genRecoDrMatch_);
    if (reco_match_genmu1.first)
      matched_muon_idx.push_back(reco_match_genmu1.second);

    if (reco_match_genmu2.first)
      matched_muon_idx.push_back(reco_match_genmu2.second);

    reco_match_genmu1 = MatchReco<reco::Track>(*tracks,
                                               std::get<float>(nt.branches["genmu1_eta"].value),
                                               std::get<float>(nt.branches["genmu1_phi"].value),
                                               std::get<int>(nt.branches["genmu1_charge"].value),
                                               genRecoDrMatch_);
    reco_match_genmu2 = MatchReco<reco::Track>(*tracks,
                                               std::get<float>(nt.branches["genmu2_eta"].value),
                                               std::get<float>(nt.branches["genmu2_phi"].value),
                                               std::get<int>(nt.branches["genmu2_charge"].value),
                                               genRecoDrMatch_);

    if (reco_match_genmu1.first)
      matched_track_idx.push_back(reco_match_genmu1.second);
    if (reco_match_genmu2.first)
      matched_track_idx.push_back(reco_match_genmu2.second);
  }

  // Match between HLT and Reco::muon -----------------------------------
  std::vector<unsigned> trg_idx;
  for (unsigned itrg = 0; itrg < nt.trg_pt.size(); ++itrg) {
    float minDR = 1000;
    unsigned idx = 0;
    for (auto& mu : *muons) {
      if (minDR < deltaR(nt.trg_eta[itrg], nt.trg_phi[itrg], mu.eta(), mu.phi()))
        continue;
      minDR = deltaR(nt.trg_eta[itrg], nt.trg_phi[itrg], mu.eta(), mu.phi());
      idx = &mu - &muons->at(0);
    }
    if (debug_ > 0) {
      std::cout << tagFilters_[itrg] << std::endl;
      std::cout << "Trg " << itrg << ", min DR " << minDR << std::endl;
    }
    if (minDR < trgDRwindow_)
      trg_idx.push_back(idx);
    if (minDR < trgDRwindow_ && debug_ > 0)
      std::cout << "Matched!" << std::endl;
  }
  nt.branches["nmuons"] = (int)muons->size();

  /**
     
     Tag muons
     ---------
     
     Loop over Reco::Muons and find the tag muons
     They must be at least looseID and fire a trigger path
     This conditions are further tightened later
   
  **/

  std::vector<unsigned> tag_muon_map;  // idx of tag muon in muons
  RecoMuonAndTransientTrkCollection tag_trkttrk;
  std::vector<bool> genmatched_tag;
  for (const auto& mu : *muons) {
    if (mu.selectors() != 0) {  // Only 9_4_X and later have selector bits
      if (!mu.passed(static_cast<uint64_t>(pow(2, tagQual_))))
        continue;
    } else {  // For 2016, assume loose ID on the tag (can be tightened at spark level)
      if (!muon::isLooseMuon(mu))
        continue;
    }
    if (!tagSelection_(mu))
      continue;
    if (std::find(trg_idx.begin(), trg_idx.end(), &mu - &muons->at(0)) == trg_idx.end())
      continue;
    tag_trkttrk.emplace_back(std::make_pair(mu, reco::TransientTrack(*mu.bestTrack(), &(*bField))));
    tag_muon_map.push_back(&mu - &muons->at(0));
    if (std::find(matched_muon_idx.begin(), matched_muon_idx.end(), &mu - &muons->at(0)) != matched_muon_idx.end())
      genmatched_tag.push_back(true);
    else
      genmatched_tag.push_back(false);
  }
  nt.branches["ntag"] = (int)tag_trkttrk.size();

  if (debug_ > 0)
    std::cout << "Tag muons " << tag_trkttrk.size() << std::endl;

  /** 

      Probe tracks
      ------------  
      
      Loop over Reco::Muon and Reco::Track collections and perform a matching
      Tracks are going to be used as probes
      
  **/

  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_muon_map;
  for (const auto& mu : *muons) {
    if (muonOnly_ && !probeMuonSelection_(mu))
      continue;
    float minDR = 1000;
    unsigned int idx_trk;
    if (debug_ > 1)
      std::cout << "New trk-muon map entry pt " << mu.pt() << " eta " << mu.eta() << " phi " << mu.phi() << std::endl;
    for (const reco::Track& trk : *tracks) {
      if (mu.innerTrack().isNonnull() && mu.innerTrack().isAvailable()) {
        if (fabs(mu.innerTrack()->vz() - trk.vz()) > maxdz_trk_mu_ && maxdz_trk_mu_ > 0)
          continue;
        if (fabs(mu.innerTrack()->pt() - trk.pt()) / mu.innerTrack()->pt() > maxpt_relative_dif_trk_mu_ &&
            maxpt_relative_dif_trk_mu_ > 0)
          continue;
      }
      float DR = deltaR(mu.eta(), mu.phi(), trk.eta(), trk.phi());
      if (debug_ > 1)
        std::cout << "   DR " << DR << "  " << mu.eta() << "  " << mu.phi() << "  " << trk.eta() << "  " << trk.phi()
                  << std::endl;
      if (minDR < DR)
        continue;
      minDR = DR;
      idx_trk = &trk - &tracks->at(0);
    }
    if (minDR > maxdr_trk_mu_)
      continue;
    trk_muon_map.first.push_back(idx_trk);
    trk_muon_map.second.push_back(&mu - &muons->at(0));
  }
  if (debug_ > 0)
    std::cout << "Matched trk-mu " << trk_muon_map.first.size() << std::endl;

  using map_type = std::pair<std::vector<unsigned>, std::vector<unsigned>>;

  /**

     Further map
     -----------

     Generic functions to map probe track with a second collection of tracks (e.g. dSA, cosmics, dGl, etc)
     
     Simple mapping (not currently used) -- just check for the closest match between collections
     If closest match has dR < config cone value, include match in mapping
     This is not flexible for later use in spark (see Inclusive function below)

  **/

  auto mapTrackCollectionsSimple =
      [&](const std::vector<reco::Track>& coll_1, const std::vector<reco::Track>& coll_2, map_type& coll_map) {
        for (const reco::Track& trk_1 : coll_1) {
          unsigned idx_1 = &trk_1 - &coll_1.at(0);
          float min_dR = 1000;
          unsigned min_idx_2;
          for (const reco::Track& trk_2 : coll_2) {
            unsigned idx_2 = &trk_2 - &coll_2.at(0);
            float dR = deltaR(trk_1.eta(), trk_1.phi(), trk_2.eta(), trk_2.phi());
            if (min_dR < dR)
              continue;
            min_dR = dR;
            min_idx_2 = idx_2;
          }
          if (min_dR > maxdr_trk_dsa_)
            continue;
          coll_map.first.push_back(min_idx_2);
          coll_map.second.push_back(idx_1);
        }
      };
  (void)mapTrackCollectionsSimple;  // prevent "unused method" error

  /**

     Inclusive map
     -------------
  
     Inclusive mapping (currently used for all) -- check and count all possible matches within a config dR cone value
     Also pick the pair with smallest dR, *even* if it is greater than config dR cone -- this allows to customize match criteria in spark later on
     Also this should work even for collimated muons where nmatched may be greater than 1

     Start with probe at first loop level and dSA/cosmic/dGl at second level instead, to ensure each probe that has a match is actually recorded 
     i.e. a dSA/cosmic/dGl track can be matched to more than one generalTrack probem, but all generalTrack probes that match should be pushed to the map

  **/

  auto mapTrackCollectionsInclusive = [&]<typename T>(const vector<T>& coll_1,
                                                      const vector<reco::Track>& coll_2,
                                                      map_type& coll_map,
                                                      vector<float>& min_dRs,
                                                      vector<unsigned>& n_matched) {
    for (const T& trk_1 : coll_1) {
      unsigned idx_1 = &trk_1 - &coll_1.at(0);
      unsigned matched = 0;
      float min_dR = 1000;
      unsigned min_idx_2 = 1000;
      for (const reco::Track& trk_2 : coll_2) {
        unsigned idx_2 = &trk_2 - &coll_2.at(0);
        float dR = deltaR(trk_1.eta(), trk_1.phi(), trk_2.eta(), trk_2.phi());
        // first add to count of matched items in coll 2 within X dR cone of track in coll 1
        if (dR < maxdr_trk_dsa_)
          matched++;
        // now check if this is the closest one so far
        if (dR < min_dR) {
          min_dR = dR;
          min_idx_2 = idx_2;
        }
      }
      if (min_idx_2 == 1000)  // corner case: coll 2 is empty
        continue;
      coll_map.first.push_back(idx_1);
      coll_map.second.push_back(min_idx_2);
      min_dRs.push_back(min_dR);
      n_matched.push_back(matched);
    }
  };

  // Initialize different maps

  map_type probe_dGl_map;  // Map probe and dGl
  vector<float> probe_dGl_dRs;
  vector<unsigned> probe_dGl_segmentmatches;  // Not dR, segment matching instead
  map_type probe_cosmic_map;                  // Map probe and cosmics
  vector<float> probe_cosmic_dRs;
  vector<unsigned> probe_cosmic_nmatched;
  mapTrackCollectionsInclusive(*tracks, *staCosmic, probe_cosmic_map, probe_cosmic_dRs, probe_cosmic_nmatched);
  map_type probe_dSA_map;  // Map probe and dSA
  vector<float> probe_dSA_dRs;
  vector<unsigned> probe_dSA_segmentmatches;  // Not dR, segment matching instead
  map_type tag_dSA_map;                       // Map tag and dSA
  vector<float> tag_dSA_dRs;
  vector<unsigned> tag_dSA_segmentmatches;  // Not dR, segment matching instead

  /**
     
     Displaced StandAlone matching
     -----------------------------

     [Adapted from displaced dimuon analysis]
     Segment match: consider probe reco::Muons and dSA tracks
     matched if the segments used to build the dSA are the
     same as or a subset of segments of the reco::Muons.
  
  **/

  for (const auto& track : *tracks) {
    // can only do segment matching on probe tracks that are matched to muons
    auto it = std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), &track - &tracks->at(0));
    if (it == trk_muon_map.first.end())
      continue;
    unsigned idx_map = std::distance(trk_muon_map.first.begin(), it);
    unsigned idx_track = trk_muon_map.first[idx_map];
    unsigned idx_muon = trk_muon_map.second[idx_map];
    auto muon = muons->at(idx_muon);
    // also make sure probes are arbitrated tracker muons
    if (!(muon.isTrackerMuon() && muon::isGoodMuon(muon, muon::TrackerMuonArbitrated)))
      continue;

    int max_nmatches = -1;
    float min_dR = +99.f;
    unsigned matched_idx;
    for (const auto& dsa : *dSAmuons) {
      unsigned idx_dsa = &dsa - &dSAmuons->at(0);

      float dR = deltaR(dsa.eta(), dsa.phi(), muon.eta(), muon.phi());
      // Don't waste time with far away dSA tracks and muons
      // Aug. 2021: comment out for now to test dSA-dGl correspondence
      // if (dR > 0.7)
      // continue;

      int nmatches = 0;
      for (auto& hit : dsa.recHits()) {
        if (!hit->isValid())
          continue;
        DetId id = hit->geographicalId();
        if (id.det() != DetId::Muon)
          continue;
        if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
          for (auto& chamber : muon.matches()) {
            if (chamber.id.rawId() != id.rawId())
              continue;
            for (auto& segment : chamber.segmentMatches) {
              if (fabs(segment.x - hit->localPosition().x()) < 1e-6 &&
                  fabs(segment.y - hit->localPosition().y()) < 1e-6) {
                if (debug_)
                  std::cout << "matched segment found!!! subdet "
                            << ((chamber.id.subdetId() == MuonSubdetId::DT) ? "DT" : "CSC")
                            << " id = " << chamber.id.rawId() << " x = " << segment.x << " y = " << segment.y
                            << std::endl;
                nmatches++;
                break;
              }
            }
          }
        }
      }
      if (nmatches > max_nmatches) {
        max_nmatches = nmatches;
        min_dR = dR;
        matched_idx = idx_dsa;
      } else if (nmatches == max_nmatches && dR < min_dR) {
        min_dR = dR;
        matched_idx = idx_dsa;
      }
    }
    if (max_nmatches > -1) {
      probe_dSA_map.first.push_back(idx_track);
      probe_dSA_map.second.push_back(matched_idx);
      probe_dSA_dRs.push_back(min_dR);
      probe_dSA_segmentmatches.push_back(max_nmatches);
    }
  }

  /**
      
     Displaced Global matching 
     -------------------------
    
     [Adapted from displaced dimuon analysis]
     Segment match: consider probe reco::Muons and dGl tracks
     matched if the segments used to build the dGl are the
     same as or a subset of segments of the reco::Muons.
  
  **/

  for (const auto& track : *tracks) {
    // can only do segment matching on probe tracks that are matched to muons
    auto it = std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), &track - &tracks->at(0));
    if (it == trk_muon_map.first.end())
      continue;
    unsigned idx_map = std::distance(trk_muon_map.first.begin(), it);
    unsigned idx_track = trk_muon_map.first[idx_map];
    unsigned idx_muon = trk_muon_map.second[idx_map];
    auto muon = muons->at(idx_muon);
    // also make sure probes are arbitrated tracker muons
    if (!(muon.isTrackerMuon() && muon::isGoodMuon(muon, muon::TrackerMuonArbitrated)))
      continue;

    int max_nmatches = -1;
    float min_dR = +99.f;
    unsigned matched_idx;
    for (const auto& dgl : *dGlmuons) {
      unsigned idx_dgl = &dgl - &dGlmuons->at(0);

      float dR = deltaR(dgl.eta(), dgl.phi(), muon.eta(), muon.phi());
      // Don't waste time with far away dGl tracks and muons
      // Aug. 2021: comment out for now to test dSA-dGl correspondence
      // if (dR > 0.7)
      // continue;

      int nmatches = 0;
      for (auto& hit : dgl.recHits()) {
        if (!hit->isValid())
          continue;
        DetId id = hit->geographicalId();
        if (id.det() != DetId::Muon)
          continue;
        if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
          for (auto& chamber : muon.matches()) {
            if (chamber.id.rawId() != id.rawId())
              continue;
            for (auto& segment : chamber.segmentMatches) {
              if (fabs(segment.x - hit->localPosition().x()) < 1e-6 &&
                  fabs(segment.y - hit->localPosition().y()) < 1e-6) {
                if (debug_)
                  std::cout << "matched segment found!!! subdet "
                            << ((chamber.id.subdetId() == MuonSubdetId::DT) ? "DT" : "CSC")
                            << " id = " << chamber.id.rawId() << " x = " << segment.x << " y = " << segment.y
                            << std::endl;
                nmatches++;
                break;
              }
            }
          }
        }
      }
      if (nmatches > max_nmatches) {
        max_nmatches = nmatches;
        min_dR = dR;
        matched_idx = idx_dgl;
      } else if (nmatches == max_nmatches && dR < min_dR) {
        min_dR = dR;
        matched_idx = idx_dgl;
      }
    }
    if (max_nmatches > -1) {
      probe_dGl_map.first.push_back(idx_track);
      probe_dGl_map.second.push_back(matched_idx);
      probe_dGl_dRs.push_back(min_dR);
      probe_dGl_segmentmatches.push_back(max_nmatches);
    }
  }

  /**
     
     Displaced Standalone Tag matching 
     ---------------------------------
    
     [Adapted from displaced dimuon analysis]
     Segment match: consider tag reco::Muons and dSA tracks
     matched if the segments used to build the dSA are the
     same as or a subset of segments of the reco::Muons.

  **/

  for (const auto& muon : *muons) {
    unsigned idx_muon = &muon - &muons->at(0);

    int max_nmatches = -1;
    float min_dR = +99.f;
    unsigned matched_idx;
    for (const auto& dsa : *dSAmuons) {
      unsigned idx_dsa = &dsa - &dSAmuons->at(0);

      float dR = deltaR(dsa.eta(), dsa.phi(), muon.eta(), muon.phi());
      // Don't waste time with far away dSA tracks and muons
      // Aug. 2021: comment out for now to test dSA-dGl correspondence
      // if (dR > 0.7)
      // continue;

      int nmatches = 0;
      for (auto& hit : dsa.recHits()) {
        if (!hit->isValid())
          continue;
        DetId id = hit->geographicalId();
        if (id.det() != DetId::Muon)
          continue;
        if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
          for (auto& chamber : muon.matches()) {
            if (chamber.id.rawId() != id.rawId())
              continue;
            for (auto& segment : chamber.segmentMatches) {
              if (fabs(segment.x - hit->localPosition().x()) < 1e-6 &&
                  fabs(segment.y - hit->localPosition().y()) < 1e-6) {
                if (debug_)
                  std::cout << "matched segment found!!! subdet "
                            << ((chamber.id.subdetId() == MuonSubdetId::DT) ? "DT" : "CSC")
                            << " id = " << chamber.id.rawId() << " x = " << segment.x << " y = " << segment.y
                            << std::endl;
                nmatches++;
                break;
              }
            }
          }
        }
      }
      if (nmatches > max_nmatches) {
        max_nmatches = nmatches;
        min_dR = dR;
        matched_idx = idx_dsa;
      } else if (nmatches == max_nmatches && dR < min_dR) {
        min_dR = dR;
        matched_idx = idx_dsa;
      }
    }
    if (max_nmatches > -1) {
      tag_dSA_map.first.push_back(idx_muon);
      tag_dSA_map.second.push_back(matched_idx);
      tag_dSA_dRs.push_back(min_dR);
      tag_dSA_segmentmatches.push_back(max_nmatches);
    }
  }

  // Muon collection for jet cleaning
  std::vector<reco::Muon> muForJetCleaning;
  for (const auto& mu : *muons) {
    if (!muon::isLooseMuon(mu))
      continue;
    muForJetCleaning.push_back(mu);
  }
  sort(muForJetCleaning.begin(), muForJetCleaning.end(), [](const auto& l, const auto& r) { return l.pt() > r.pt(); });

  // Fill Jet branches   /  If IncludeJets == True ------------
  std::vector<reco::PFJet> corrJets;
  std::vector<float> genJets_pt;
  std::vector<float> genJets_eta;
  std::vector<float> genJets_phi;
  std::vector<float> genJets_mass;
  std::vector<float> jets_bTag_deepCSV;
  if (includeJets_) {
    edm::Handle<std::vector<reco::GenJet>> genJets;

    if (isMC_) {
      // Gen Jet Info
      iEvent.getByToken(genJetsToken_, genJets);
      for (const auto& genJet : *genJets) {
        genJets_pt.push_back((float)genJet.pt());
        genJets_eta.push_back((float)genJet.eta());
        genJets_phi.push_back((float)genJet.phi());
        genJets_mass.push_back((float)genJet.mass());
      }
      nt.branches["genJets_pt"] = (std::vector<float>)genJets_pt;
      nt.branches["genJets_eta"] = (std::vector<float>)genJets_eta;
      nt.branches["genJets_phi"] = (std::vector<float>)genJets_phi;
      nt.branches["genJets_mass"] = (std::vector<float>)genJets_mass;
    }

    // Get selected jets and fill branches
    for (size_t i = 0; i < jets->size(); ++i) {
      reco::PFJetRef jet(jets, i);
      if (CrossClean(*jet, muForJetCleaning))
        continue;
      std::unique_ptr<reco::PFJet> corrJet(jet->clone());
      double jec = jetCorrector->correction(*jet);
      corrJet->scaleEnergy(jec);
      double smearFactor = 1.0;
      if (isMC_) {
        double jet_resolution = resolution.getResolution({{JME::Binning::JetPt, corrJet->pt()},
                                                          {JME::Binning::JetEta, corrJet->eta()},
                                                          {JME::Binning::Rho, *rhoHandle}});
        double jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetPt, corrJet->pt()},
                                                      {JME::Binning::JetEta, corrJet->eta()},
                                                      {JME::Binning::Rho, *rhoHandle}},
                                                     Variation::NOMINAL);
        // gen matching
        double min_dR = std::numeric_limits<double>::infinity();
        const reco::GenJet* matched_genJet = nullptr;
        for (const auto& genJet : *genJets) {
          double dR = deltaR(genJet, *corrJet);
          if (dR > min_dR)
            continue;
          if (dR >= 0.2)
            continue;
          min_dR = dR;
          matched_genJet = &genJet;
        }
        if (matched_genJet) {
          double dPt = corrJet->pt() - matched_genJet->pt();
          smearFactor = 1 + (jer_sf - 1.) * dPt / corrJet->pt();
        } else if (jer_sf > 1) {
          double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
          std::normal_distribution<> d(0, sigma);
          smearFactor = 1. + d(m_random_generator);
        }
        if (corrJet->pt() * smearFactor < 0.01) {
          smearFactor = 0.01 / corrJet->energy();
        }
      }
      corrJet->scaleEnergy(smearFactor);
      corrJets.push_back(*corrJet);
      FillJetBranches(*jet, *corrJet, nt, era_);
      if (deepCSVProbb->size() > i && deepCSVProbbb->size() > i) {
        jets_bTag_deepCSV.push_back((*deepCSVProbb)[i].second + (*deepCSVProbbb)[i].second);
      } else
        jets_bTag_deepCSV.push_back(-9999.);
    }
    nt.branches["jets_bTag_deepCSV"] = (std::vector<float>)jets_bTag_deepCSV;
  }

  /**

     Find TnP pairs
     --------------
 
     run over tracks and probes once prior to filling tree to determine ordering of pairs
     this is necessary to use tag-probe pair with highest "quality" later on in spark_tnp.
     Also calculate the total number of pairs per event to include in ntuple

     1 - Loop over tag muons and probe tracks
     2 - Each Tag-Probe pair must pass a minimum dz and mass cuts and have opposite charge 
         Pairs are discarded otherwise
     3 - The Tag and Probe tracks are used to perform a fit to the primary vertex
     4 - All the probe tracks associated to a tag are sorted as a funtion of the compatibility to the PV

  **/

  using t_pair_prob = std::pair<float, std::pair<int, int>>;
  std::priority_queue<t_pair_prob> pair_dPhi_muons;
  std::priority_queue<t_pair_prob, vector<t_pair_prob>, std::greater<t_pair_prob>> pair_dz_PV_SV;    // inverse sort
  std::priority_queue<t_pair_prob, vector<t_pair_prob>, std::greater<t_pair_prob>> pair_dM_Z_Mmumu;  // inverse sort
  std::map<int, int> map_tagIdx_nprobes;

  using pair_prob = std::pair<int, std::vector<std::pair<int, float>>>;
  std::vector<pair_prob> pair_vtx_probs;
  map<std::pair<int, int>, float> pair_rank_vtx_prob;
  map<std::pair<int, int>, int> pair_rank_vtx_prob_idx;
  map<int, std::vector<float>> probe_vtxP;

  int nprobes;
  // loop over tags
  for (auto& tag : tag_trkttrk) {
    auto tag_idx = &tag - &tag_trkttrk[0];
    nprobes = 0;
    std::vector<std::pair<int, float>> tmp_vtx_probs;
    // loop over probes
    for (const reco::Track& probe : *tracks) {
      auto probe_idx = &probe - &tracks->at(0);

      // apply cuts on probe
      if (HighPurity_ && !probe.quality(Track::highPurity))
        continue;
      if (!probeSelection_(probe))
        continue;

      // apply cuts on pairs
      if (tag.first.charge() == probe.charge())
        continue;
      if (fabs(tag.first.vz() - probe.vz()) > pairDz_ && pairDz_ > 0)
        continue;
      float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(), probe.pt(), probe.eta(), probe.phi());
      if (mass < pairMassMin_ || mass > pairMassMax_)
        continue;

      // compute vtx
      std::vector<reco::TransientTrack> trk_pair = {tag.second, reco::TransientTrack(probe, &(*bField))};
      KlFitter vtx(trk_pair);
      if (RequireVtxCreation_ && !vtx.status())
        continue;
      if (minSVtxProb_ > 0) {
        if (vtx.prob() < minSVtxProb_) {
          continue;
        }
      }

      auto it = std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), &probe - &tracks->at(0));
      if (muonOnly_ && it == trk_muon_map.first.end())
        continue;

      nprobes++;

      float dPhi_muons = reco::deltaPhi(tag.first.phi(), probe.phi());
      math::PtEtaPhiMLorentzVector mu1(tag.first.pt(), tag.first.eta(), tag.first.phi(), MU_MASS);
      math::PtEtaPhiMLorentzVector mu2(probe.pt(), probe.eta(), probe.phi(), MU_MASS);
      float dM_Z_Mmumu = abs(91.2 - (mu1 + mu2).mass());

      // save quantities to ordered heap
      auto pair_idx = std::make_pair(tag_idx, probe_idx);
      if (RequireVtxCreation_) {
        pair_dz_PV_SV.push(std::make_pair(vtx.dz_PV_SV(std::get<float>(nt.branches["pv_z"].value)), pair_idx));
        if (vtx.prob() > 0 && vtx.status() == 1) {
          tmp_vtx_probs.emplace_back(std::make_pair(probe_idx, vtx.prob()));
        } else {
          tmp_vtx_probs.emplace_back(std::make_pair(probe_idx, -1));
        }
      } else {
        pair_dz_PV_SV.push(std::make_pair(0., std::make_pair(-1, -1)));
        tmp_vtx_probs.emplace_back(std::make_pair(probe_idx, 0.));
      }
      probe_vtxP[probe_idx].emplace_back(vtx.prob());
      pair_dPhi_muons.push(std::make_pair(dPhi_muons, pair_idx));
      pair_dM_Z_Mmumu.push(std::make_pair(dM_Z_Mmumu, pair_idx));
    }
    map_tagIdx_nprobes.insert(std::pair<int, int>(tag_idx, nprobes));

    auto compare_vtx = [=](std::pair<int, float>& a, std::pair<int, float>& b) { return a.second > b.second; };
    std::sort(tmp_vtx_probs.begin(), tmp_vtx_probs.end(), compare_vtx);

    pair_vtx_probs.emplace_back(pair_prob(tag_idx, tmp_vtx_probs));

    for (size_t j = 0; j < tmp_vtx_probs.size(); j++) {
      pair_rank_vtx_prob[std::make_pair(tag_idx, tmp_vtx_probs[j].first)] = tmp_vtx_probs[j].second;
      pair_rank_vtx_prob_idx[std::make_pair(tag_idx, tmp_vtx_probs[j].first)] = j;
    }
  }
  nt.branches["npairs"] = (int)pair_vtx_probs.size();

  map<std::pair<int, int>, int> pair_rank_dz_PV_SV;
  while (!pair_dz_PV_SV.empty()) {
    pair_rank_dz_PV_SV[pair_dz_PV_SV.top().second] = pair_rank_dz_PV_SV.size();  // careful: RHS evaluated first
    pair_dz_PV_SV.pop();
  }
  map<std::pair<int, int>, int> pair_rank_dPhi_muons;
  while (!pair_dPhi_muons.empty()) {
    pair_rank_dPhi_muons[pair_dPhi_muons.top().second] = pair_rank_dPhi_muons.size();  // careful: RHS evaluated first
    pair_dPhi_muons.pop();
  }
  map<std::pair<int, int>, int> pair_rank_dM_Z_Mmumu;
  while (!pair_dM_Z_Mmumu.empty()) {
    pair_rank_dM_Z_Mmumu[pair_dM_Z_Mmumu.top().second] = pair_rank_dM_Z_Mmumu.size();  // careful: RHS evaluated first
    pair_dM_Z_Mmumu.pop();
  }

  // Initialize some variables
  nt.branches["iprobe"] = (int)0;

  /**  
       
       Fill branches 
       -------------
       
       Once we have created the Tag-Probe pairs, loop again and save all the information. 
       The final output is a flat ntuple filled with pairs, not events.
       
       1 - Loop over Tag muons
       2 - Loop over Probe tracks (only over those already matched to the tag muon)  
       3 - Apply selection again. Redundant. 
       4 - Fill Tag branches   
       5 - Check if probe is a real muon in Reco::muon 
       6 - If so, fill probe muon branches      
       7 - If not, fill probe track branches, probe muon branches filled with -99 
       8 - Check vertex sorting, duplication of probes, and best pairs 
       9 - Write ntuple content and go for the next pair
       
  **/

  for (const auto& tag_prob_vtx : pair_vtx_probs) {
    auto& tag = tag_trkttrk[tag_prob_vtx.first];
    for (const auto& probe_vtx_vec : tag_prob_vtx.second) {
      const reco::Track& probe = tracks->at(probe_vtx_vec.first);
      // apply cuts on probe
      if (HighPurity_ && !probe.quality(Track::highPurity))
        continue;
      if (debug_ > 2)
        std::cout << "Passes high purity" << std::endl;
      if (!probeSelection_(probe))
        continue;
      if (debug_ > 2)
        std::cout << "Passes selection" << std::endl;

      // apply cuts on pairs; selected will be saved
      if (tag.first.charge() == probe.charge())
        continue;
      if (debug_ > 2)
        std::cout << "Passes charge requirement" << std::endl;
      if (fabs(tag.first.vz() - probe.vz()) > pairDz_ && pairDz_ > 0)
        continue;
      if (debug_ > 2)
        std::cout << "Passes dz cut" << std::endl;

      float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(), probe.pt(), probe.eta(), probe.phi());
      if (mass < pairMassMin_ || mass > pairMassMax_)
        continue;
      if (debug_ > 2)
        std::cout << "Passes mass cut" << std::endl;

      std::vector<reco::TransientTrack> trk_pair = {tag.second, reco::TransientTrack(probe, &(*bField))};
      KlFitter vtx(trk_pair);
      if (RequireVtxCreation_ && !vtx.status())
        continue;
      if (debug_ > 2)
        std::cout << "Passes vertex status" << std::endl;

      if (minSVtxProb_ > 0) {
        if (vtx.prob() < minSVtxProb_)
          continue;
      }
      if (debug_ > 2)
        std::cout << "Passes vertex probability cut" << std::endl;

      auto it = std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), &probe - &tracks->at(0));
      if (muonOnly_ && it == trk_muon_map.first.end())
        continue;
      if (debug_ > 2)
        std::cout << "Passes muon/track selection" << std::endl;

      math::PtEtaPhiMLorentzVector mu1(tag.first.pt(), tag.first.eta(), tag.first.phi(), MU_MASS);
      math::PtEtaPhiMLorentzVector mu2(probe.pt(), probe.eta(), probe.phi(), MU_MASS);

      nt.branches["tag_isMatchedGen"] = (bool)genmatched_tag[&tag - &tag_trkttrk[0]];

      FillTagBranches<reco::Muon, reco::Track>(tag.first, *tracks, nt, *pv);
      std::pair<bool, Measurement1D> sip3d = IPTools::signedImpactParameter3D(
          reco::TransientTrack(*tag.first.bestTrack(), &(*bField)),
          GlobalVector(tag.first.bestTrack()->px(), tag.first.bestTrack()->py(), tag.first.bestTrack()->pz()),
          *pv);
      nt.branches["tag_SIP3D"] = (float)sip3d.second.value();
      nt.branches["tag_SIP3D_err"] = (float)(pv->isValid() ? sip3d.second.error() : -1.0);

      FillMiniIso<reco::Muon, pat::PackedCandidate>(patpfcands.product(), tag.first, *rhoJetsNC, nt, true);

      // Tag-trigger matching
      auto tagRef = muonsView->refAt(tag_muon_map[&tag - &tag_trkttrk[0]]);
      pat::TriggerObjectStandAloneRef tagl1Match = (*l1Matches)[tagRef];
      if (tagl1Match.isNonnull()) {
        nt.branches["tag_l1pt"] = (float)tagl1Match->pt();
        nt.branches["tag_l1q"] = (float)(*l1Qualities)[tagRef];
        nt.branches["tag_l1dr"] = (float)(*l1Drs)[tagRef];
      } else {
        nt.branches["tag_l1pt"] = (float)-99.;
        nt.branches["tag_l1q"] = (float)-99;
        nt.branches["tag_l1dr"] = (float)99.;
      }

      pat::TriggerObjectStandAloneRef tagl1MatchByQ = (*l1MatchesByQ)[tagRef];
      if (tagl1MatchByQ.isNonnull()) {
        nt.branches["tag_l1ptByQ"] = (float)tagl1MatchByQ->pt();
        nt.branches["tag_l1qByQ"] = (float)(*l1QualitiesByQ)[tagRef];
        nt.branches["tag_l1drByQ"] = (float)(*l1DrsByQ)[tagRef];
      } else {
        nt.branches["tag_l1ptByQ"] = (float)-99.;
        nt.branches["tag_l1qByQ"] = (float)-99;
        nt.branches["tag_l1drByQ"] = (float)99.;
      }
      embedTriggerMatching(tag.first, nt.trg_filter, nt.trg_pt, nt.trg_eta, nt.trg_phi, tagFilters_, true, debug_);

      if (!simInfoIsAvailalbe) {
        FillSimMatchingBranchesDummy(nt, true);
      } else {
        const auto& msi = (*simInfo)[tagRef];
        FillSimMatchingBranchesAOD(msi, nt, true);
      }

      auto itdsa = std::find(probe_dSA_map.first.begin(), probe_dSA_map.first.end(), &probe - &tracks->at(0));
      auto itdsa_tag = std::find(tag_dSA_map.first.begin(), tag_dSA_map.first.end(), &tag - &tag_trkttrk[0]);
      auto itdgl = std::find(probe_dGl_map.first.begin(), probe_dGl_map.first.end(), &probe - &tracks->at(0));
      auto itcosmic = std::find(probe_cosmic_map.first.begin(), probe_cosmic_map.first.end(), &probe - &tracks->at(0));

      // If Probe track is also a Reco::Muon object -----------------------------------
      if (it == trk_muon_map.first.end()) {
        if (debug_ > 0)
          std::cout << "  Unsuccessful probe " << std::endl;
        reco::Muon fakeMuon;
        fakeMuon.setP4(mu2);
        fakeMuon.setCharge(probe.charge());
        FillProbeBranches<reco::Muon, reco::Track>(fakeMuon, *tracks, nt, false, *pv);
        FillProbeBranchesSelector<reco::Muon>(fakeMuon, nt, probeSelectorBits_, false);
        FillMiniIso<reco::Muon, pat::PackedCandidate>(patpfcands.product(), fakeMuon, *rhoJetsNC, nt, false);
        if (includeJets_)
          FindJetProbePair<reco::PFJet, reco::Muon>(corrJets, fakeMuon, nt);

        // store dummy trigger variables if offline muon is not found
        for (const auto& path : probeFilters_) {
          nt.probe_trg[&path - &probeFilters_[0]] = false;
          nt.probe_trg_pt[&path - &probeFilters_[0]] = -99;
          nt.probe_trg_eta[&path - &probeFilters_[0]] = -99;
          nt.probe_trg_phi[&path - &probeFilters_[0]] = -99;
          nt.probe_trg_dr[&path - &probeFilters_[0]] = 99;
        }
        nt.branches["l1pt"] = (float)-99.;
        nt.branches["l1q"] = (int)-99;
        nt.branches["l1dr"] = (float)99.;
        nt.branches["l1ptByQ"] = (float)-99.;
        nt.branches["l1qByQ"] = (int)-99;
        nt.branches["l1drByQ"] = (float)99.;

        FillSimMatchingBranchesDummy(nt, false);

        FillTunePPairBranchesDummy(nt);

      } else {  // Probe track is not a real Reco::Muon -----------------------------------

        unsigned idx = std::distance(trk_muon_map.first.begin(), it);
        if (debug_ > 0)
          std::cout << "  Successful probe pt " << muons->at(trk_muon_map.second[idx]).pt() << " eta "
                    << muons->at(trk_muon_map.second[idx]).eta() << " phi " << muons->at(trk_muon_map.second[idx]).phi()
                    << std::endl;
        FillProbeBranches<reco::Muon, reco::Track>(muons->at(trk_muon_map.second[idx]), *tracks, nt, true, *pv);
        FillProbeBranchesSelector<reco::Muon>(muons->at(trk_muon_map.second[idx]), nt, probeSelectorBits_, true);
        std::pair<bool, Measurement1D> sip3d = IPTools::signedImpactParameter3D(
            reco::TransientTrack(*muons->at(trk_muon_map.second[idx]).bestTrack(), &(*bField)),
            GlobalVector(muons->at(trk_muon_map.second[idx]).bestTrack()->px(),
                         muons->at(trk_muon_map.second[idx]).bestTrack()->py(),
                         muons->at(trk_muon_map.second[idx]).bestTrack()->pz()),
            *pv);
        nt.branches["probe_SIP3D"] = (float)sip3d.second.value();
        nt.branches["probe_SIP3D_err"] = (float)(pv->isValid() ? sip3d.second.error() : -1.0);
        FillMiniIso<reco::Muon, pat::PackedCandidate>(
            patpfcands.product(), muons->at(trk_muon_map.second[idx]), *rhoJetsNC, nt, false);
        if (includeJets_)
          FindJetProbePair<reco::PFJet, pat::Muon>(corrJets, muons->at(trk_muon_map.second[idx]), nt);

        // Probe-trigger matching
        auto muRef = muonsView->refAt(trk_muon_map.second[idx]);
        pat::TriggerObjectStandAloneRef l1Match = (*l1Matches)[muRef];
        if (l1Match.isNonnull()) {
          nt.branches["l1pt"] = (float)l1Match->pt();
          nt.branches["l1q"] = (float)(*l1Qualities)[muRef];
          nt.branches["l1dr"] = (float)(*l1Drs)[muRef];
        } else {
          nt.branches["l1pt"] = (float)-99.;
          nt.branches["l1q"] = (float)-99;
          nt.branches["l1dr"] = (float)99.;
        }

        pat::TriggerObjectStandAloneRef l1MatchByQ = (*l1MatchesByQ)[muRef];
        if (l1MatchByQ.isNonnull()) {
          nt.branches["l1ptByQ"] = (float)l1MatchByQ->pt();
          nt.branches["l1qByQ"] = (float)(*l1QualitiesByQ)[muRef];
          nt.branches["l1drByQ"] = (float)(*l1DrsByQ)[muRef];
        } else {
          nt.branches["l1ptByQ"] = (float)-99.;
          nt.branches["l1qByQ"] = (float)-99;
          nt.branches["l1drByQ"] = (float)99.;
        }

        embedTriggerMatching(muons->at(trk_muon_map.second[idx]),
                             nt.prb_filter,
                             nt.prb_pt,
                             nt.prb_eta,
                             nt.prb_phi,
                             probeFilters_,
                             false,
                             debug_);

        if (!simInfoIsAvailalbe) {
          FillSimMatchingBranchesDummy(nt, false);
        } else {
          const auto& msi = (*simInfo)[muRef];
          FillSimMatchingBranchesAOD(msi, nt, false);
        }

        // TuneP pair branches
        if (tag.first.tunePMuonBestTrack().isNonnull() &&
            muons->at(trk_muon_map.second[idx]).tunePMuonBestTrack().isNonnull()) {
          const reco::TrackRef tag_tuneP = tag.first.tunePMuonBestTrack();
          const reco::TrackRef probe_tuneP = muons->at(trk_muon_map.second[idx]).tunePMuonBestTrack();
          FillTunePPairBranches<reco::Track, reco::Track>(*tag_tuneP, *probe_tuneP, nt);
          std::vector<reco::TransientTrack> ttrk_pair_tuneP = {reco::TransientTrack(*tag_tuneP, &(*bField)),
                                                               reco::TransientTrack(*probe_tuneP, &(*bField))};
          KlFitter vtx_tuneP(ttrk_pair_tuneP);
          vtx_tuneP.fillNtuple(nt, true);
        } else {
          FillTunePPairBranchesDummy(nt);
        }
      }

      //
      // Fill Displaced variables and check matching
      //

      if (itdsa == probe_dSA_map.first.end()) {
        nt.branches["probe_dsa_segmentMatches"] = (int)-1;
        FillProbeBranchesdSA<reco::Track>(probe, nt, false);
      } else {
        unsigned idx = std::distance(probe_dSA_map.first.begin(), itdsa);
        nt.branches["probe_dsa_segmentMatches"] = (int)probe_dSA_segmentmatches[idx];
        nt.branches["probe_dsa_minDR"] = (float)probe_dSA_dRs[idx];
        if (debug_ > 0)
          std::cout << "Successful probe dSA " << dSAmuons->at(probe_dSA_map.second[idx]).pt() << " eta "
                    << dSAmuons->at(probe_dSA_map.second[idx]).eta() << " phi "
                    << dSAmuons->at(probe_dSA_map.second[idx]).phi() << std::endl;
        FillProbeBranchesdSA<reco::Track>(dSAmuons->at(probe_dSA_map.second[idx]), nt, true);
      }

      if (itdsa_tag == tag_dSA_map.first.end()) {
        nt.branches["tag_dsa_segmentMatches"] = (int)-1;
        FillTagBranchesdSA<reco::Track>(probe, nt, false);
      } else {
        unsigned idx = std::distance(tag_dSA_map.first.begin(), itdsa_tag);
        nt.branches["tag_dsa_segmentMatches"] = (int)tag_dSA_segmentmatches[idx];
        nt.branches["tag_dsa_minDR"] = (float)tag_dSA_dRs[idx];
        if (debug_ > 0)
          std::cout << "Successful tag dSA " << dSAmuons->at(tag_dSA_map.second[idx]).pt() << " eta "
                    << dSAmuons->at(tag_dSA_map.second[idx]).eta() << " phi "
                    << dSAmuons->at(tag_dSA_map.second[idx]).phi() << std::endl;
        FillTagBranchesdSA<reco::Track>(dSAmuons->at(tag_dSA_map.second[idx]), nt, true);
      }

      if (itdgl == probe_dGl_map.first.end()) {
        nt.branches["probe_dgl_segmentMatches"] = (float)-1;
        FillProbeBranchesdgl<reco::Track>(probe, nt, false);
      } else {
        unsigned idx = std::distance(probe_dGl_map.first.begin(), itdgl);
        nt.branches["probe_dgl_segmentMatches"] = (float)probe_dGl_segmentmatches[idx];
        nt.branches["probe_dgl_minDR"] = (float)probe_dGl_dRs[idx];
        if (debug_ > 0)
          std::cout << "Successful probe displaced global " << dGlmuons->at(probe_dGl_map.second[idx]).pt() << " eta "
                    << dGlmuons->at(probe_dGl_map.second[idx]).eta() << " phi "
                    << dGlmuons->at(probe_dGl_map.second[idx]).phi() << std::endl;
        FillProbeBranchesdgl<reco::Track>(dGlmuons->at(probe_dGl_map.second[idx]), nt, true);
      }

      if (itcosmic == probe_cosmic_map.first.end()) {
        nt.branches["probe_ncosmic"] = (int)0;
        FillProbeBranchesCosmic<reco::Track>(probe, nt, false);
      } else {
        unsigned idx = std::distance(probe_cosmic_map.first.begin(), itcosmic);
        nt.branches["probe_ncosmic"] = (int)probe_cosmic_nmatched[idx];
        nt.branches["probe_cosmic_minDR"] = (float)probe_cosmic_dRs[idx];
        if (debug_ > 0)
          std::cout << "Successful probe cosmic " << staCosmic->at(probe_cosmic_map.second[idx]).pt() << " eta "
                    << staCosmic->at(probe_cosmic_map.second[idx]).eta() << " phi "
                    << staCosmic->at(probe_cosmic_map.second[idx]).phi() << std::endl;
        FillProbeBranchesCosmic<reco::Track>(staCosmic->at(probe_cosmic_map.second[idx]), nt, true);
      }

      //Fill vertex variables
      RecoTrkAndTransientTrk probe_pair = std::make_pair(probe, reco::TransientTrack(probe, &(*bField)));
      FillPairBranches<RecoMuonAndTransientTrk, RecoTrkAndTransientTrk>(tag, probe_pair, nt, prop1_);

      vtx.fillNtuple(nt);

      // Fill Gen variables
      auto it_genmatch = std::find(matched_track_idx.begin(), matched_track_idx.end(), &probe - &tracks->at(0));
      nt.branches["probe_isMatchedGen"] = (bool)(it_genmatch != matched_track_idx.end());

      // Rank probes, check duplicates, and best fits
      nt.branches["iprobe"] = (int)(std::get<int>(nt.branches["iprobe"].value) + 1);
      nt.branches["pair_rank_dz_PV_SV"] = (float)pair_rank_dz_PV_SV[{&tag - &tag_trkttrk[0], &probe - &tracks->at(0)}];
      nt.branches["pair_rank_dPhi_muons"] =
          (float)pair_rank_dPhi_muons[{&tag - &tag_trkttrk[0], &probe - &tracks->at(0)}];
      nt.branches["pair_rank_dM_Z_Mmumu"] =
          (float)pair_rank_dM_Z_Mmumu[{&tag - &tag_trkttrk[0], &probe - &tracks->at(0)}];
      nt.branches["probe_isHighPurity"] = (bool)probe.quality(Track::highPurity);

      nt.branches["pair_rank_vtx_prob"] =
          (int)pair_rank_vtx_prob_idx[{&tag - &tag_trkttrk[0], &probe - &tracks->at(0)}];
      std::vector vtxProbs = probe_vtxP[probe_vtx_vec.first];
      float maximum_vtx = std::max_element(vtxProbs.begin(), vtxProbs.end())[0];
      if (probe_vtxP[probe_vtx_vec.first].size() == 1) {
        nt.branches["probe_isDuplicated"] = (bool)false;
        nt.branches["probe_isBestPair"] = (bool)true;
      } else {
        nt.branches["probe_isDuplicated"] = (bool)true;
        if (pair_rank_vtx_prob[std::make_pair(tag_prob_vtx.first, probe_vtx_vec.first)] >= maximum_vtx ||
            !RequireVtxCreation_) {
          nt.branches["probe_isBestPair"] = (bool)true;
        } else {
          nt.branches["probe_isBestPair"] = (bool)false;
        }
      }

      t1->Fill();
    }
  }
}

// ------------ method called once each job just before starting event loop
// ------------
void MuonFullAODAnalyzer::beginJob() {
  t1 = fs->make<TTree>("Events", "Events");
  nt.SetTree(t1);
  nt.CreateBranches(HLTPaths_, probeSelectorNames_);
  if (!tagFilters_.empty())
    nt.CreateExtraTrgBranches(tagFilters_, true);
  if (!probeFilters_.empty())
    nt.CreateExtraTrgBranches(probeFilters_, false);
}

// ------------ method called once each job just after ending the event loop
// ------------
void MuonFullAODAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void MuonFullAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

///////////////////////

// define this as a plug-in
DEFINE_FWK_MODULE(MuonFullAODAnalyzer);
