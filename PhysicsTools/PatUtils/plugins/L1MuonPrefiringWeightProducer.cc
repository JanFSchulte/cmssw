// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TFile.h"
#include "TF1.h"

#include <iostream>
enum fluctuations { central = 0, up, down };

class L1MuonPrefiringWeightProducer : public edm::global::EDProducer<> {
public:
  explicit L1MuonPrefiringWeightProducer(const edm::ParameterSet&);
  ~L1MuonPrefiringWeightProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  double getPrefiringRate(double eta, double phi, double pt, fluctuations fluctuation) const;

  edm::EDGetTokenT<std::vector<pat::Muon> > muon_token_;

  TFile* file_prefiringparams_;
  TF1* parametrization0p0To0p2_;
  TF1* parametrization0p2To0p3_;
  TF1* parametrization0p3To0p55_;
  TF1* parametrization0p55To0p83_;
  TF1* parametrization0p83To1p24_;
  TF1* parametrization1p24To1p4_;
  TF1* parametrization1p4To1p6_;
  TF1* parametrization1p6To1p8_;
  TF1* parametrization1p8To2p1_;
  TF1* parametrization2p1To2p25_;
  TF1* parametrization2p25To2p4_;

  std::string dataera_;
  double prefiringRateSystUnc_;
  bool skipwarnings_;
};

L1MuonPrefiringWeightProducer::L1MuonPrefiringWeightProducer(const edm::ParameterSet& iConfig) {
  muon_token_ = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("TheMuons"));

  dataera_ = iConfig.getParameter<std::string>("DataEra");
  prefiringRateSystUnc_ = iConfig.getParameter<double>("PrefiringRateSystematicUncty");
  skipwarnings_ = iConfig.getParameter<bool>("SkipWarnings");

  std::string fname = iConfig.getParameter<std::string>("L1MuonParametrizations");
  edm::FileInPath mapsfilepath("PhysicsTools/PatUtils/data/" + fname);
  file_prefiringparams_ = new TFile(mapsfilepath.fullPath().c_str(), "read");
  if (file_prefiringparams_ == nullptr && !skipwarnings_)
    std::cout << "File with maps not found. All prefiring weights set to 0." << std::endl;

  TString paramName = "L1prefiring_muonparam_0.0To0.2_" + dataera_;
  parametrization0p0To0p2_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_0.2To0.3_" + dataera_;
  parametrization0p2To0p3_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_0.3To0.55_" + dataera_;
  parametrization0p3To0p55_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_0.55To0.83_" + dataera_;
  parametrization0p55To0p83_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_0.83To1.24_" + dataera_;
  parametrization0p83To1p24_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_1.24To1.4_" + dataera_;
  parametrization1p24To1p4_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_1.4To1.6_" + dataera_;
  parametrization1p4To1p6_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_1.6To1.8_" + dataera_;
  parametrization1p6To1p8_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_1.8To2.1_" + dataera_;
  parametrization1p8To2p1_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_2.1To2.25_" + dataera_;
  parametrization2p1To2p25_ = (TF1*)file_prefiringparams_->Get(paramName);
  paramName = "L1prefiring_muonparam_2.25To2.4_" + dataera_;
  parametrization2p25To2p4_ = (TF1*)file_prefiringparams_->Get(paramName);


  produces<double>("nonPrefiringProbMuon").setBranchAlias("nonPrefiringProbMuon");
  produces<double>("nonPrefiringProbMuonUp").setBranchAlias("nonPrefiringProbMuonUp");
  produces<double>("nonPrefiringProbMuonDown").setBranchAlias("nonPrefiringProbMuonDown");
}

L1MuonPrefiringWeightProducer::~L1MuonPrefiringWeightProducer() {
  delete file_prefiringparams_;
  delete parametrization0p0To0p2_;
  delete parametrization0p2To0p3_;
  delete parametrization0p3To0p55_;
  delete parametrization0p55To0p83_;
  delete parametrization0p83To1p24_;
  delete parametrization1p24To1p4_;
  delete parametrization1p4To1p6_;
  delete parametrization1p6To1p8_;
  delete parametrization1p8To2p1_;
  delete parametrization2p1To2p25_;
  delete parametrization2p25To2p4_;
}

void L1MuonPrefiringWeightProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  using namespace edm;

  //Muons
  edm::Handle<std::vector<pat::Muon> > theMuons;
  iEvent.getByToken(muon_token_, theMuons);

  //Probability for the event NOT to prefire, computed with the prefiring maps per object.
  //Up and down values correspond to the resulting value when shifting up/down all prefiring rates in prefiring maps.
  double nonPrefiringProba[3] = {1., 1., 1.};  //0: central, 1: up, 2: down

  for (const auto fluct : {fluctuations::central, fluctuations::up, fluctuations::down}) {
    for (const auto& muon : *theMuons) {
      double pt = muon.pt();
      double phi = muon.eta();
      double eta = muon.eta();

      double prefiringprob_gam = getPrefiringRate(eta, phi, pt, fluct);
      nonPrefiringProba[fluct] *= (1. - prefiringprob_gam);
    }
  }

  auto nonPrefiringProb = std::make_unique<double>(nonPrefiringProba[0]);
  auto nonPrefiringProbUp = std::make_unique<double>(nonPrefiringProba[1]);
  auto nonPrefiringProbDown = std::make_unique<double>(nonPrefiringProba[2]);
  iEvent.put(std::move(nonPrefiringProb), "nonPrefiringProbMuon");
  iEvent.put(std::move(nonPrefiringProbUp), "nonPrefiringProbMuonUp");
  iEvent.put(std::move(nonPrefiringProbDown), "nonPrefiringProbMuonDown");
}

double L1MuonPrefiringWeightProducer::getPrefiringRate(double eta,
						       double phi,
                                                       double pt,
                                                       fluctuations fluctuation) const {


  double prefrate;
  double statuncty;
 
  if (std::abs(eta) < 0.2){ 
	if (parametrization0p0To0p2_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization0p0To0p2_ == nullptr)
    		return 0.;
  	prefrate = parametrization0p0To0p2_->Eval(pt);
  	statuncty = parametrization0p0To0p2_->GetParError(2);
  }
  else if (std::abs(eta) < 0.3){ 
	if (parametrization0p2To0p3_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization0p2To0p3_ == nullptr)
    		return 0.;

        prefrate = parametrization0p2To0p3_->Eval(pt);
        statuncty = parametrization0p2To0p3_->GetParError(2);
  }
  else if (std::abs(eta) < 0.55){ 
	if (parametrization0p3To0p55_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization0p3To0p55_ == nullptr)
    		return 0.;

        prefrate = parametrization0p3To0p55_->Eval(pt);
        statuncty = parametrization0p3To0p55_->GetParError(2);
  }
  else if (std::abs(eta) < 0.83){ 
	if (parametrization0p55To0p83_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization0p55To0p83_ == nullptr)
    		return 0.;

        prefrate = parametrization0p55To0p83_->Eval(pt);
        statuncty = parametrization0p55To0p83_->GetParError(2);
  }
  else if (std::abs(eta) < 1.24){ 
	if (parametrization0p83To1p24_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization0p83To1p24_ == nullptr)
    		return 0.;

        prefrate = parametrization0p83To1p24_->Eval(pt);
        statuncty = parametrization0p83To1p24_->GetParError(2);
  }
  else if (std::abs(eta) < 1.4){ 
	if (parametrization1p24To1p4_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization1p24To1p4_ == nullptr)
    		return 0.;

        prefrate = parametrization1p24To1p4_->Eval(pt);
        statuncty = parametrization1p24To1p4_->GetParError(2);
  }
  else if (std::abs(eta) < 1.6){ 
	if (parametrization1p4To1p6_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization1p4To1p6_ == nullptr)
    		return 0.;

        prefrate = parametrization1p4To1p6_->Eval(pt);
        statuncty = parametrization1p4To1p6_->GetParError(2);
  }
  else if (std::abs(eta) < 1.8){
	if (parametrization1p6To1p8_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization1p6To1p8_ == nullptr)
    		return 0.;
 
        prefrate = parametrization1p6To1p8_->Eval(pt);
        statuncty = parametrization1p6To1p8_->GetParError(2);
  }
  else if (std::abs(eta) < 2.1){
	if (parametrization1p8To2p1_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization1p8To2p1_ == nullptr)
    		return 0.;

        prefrate = parametrization1p8To2p1_->Eval(pt);
        statuncty = parametrization1p8To2p1_->GetParError(2);
  }
  else if (std::abs(eta) < 2.25){
	if (parametrization2p1To2p25_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization2p1To2p25_ == nullptr)
    		return 0.;

        prefrate = parametrization2p1To2p25_->Eval(pt);
        statuncty = parametrization2p1To2p25_->GetParError(2);
  }
  else if (std::abs(eta) < 2.4){
	if (parametrization2p25To2p4_ == nullptr && !skipwarnings_)
    		std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  	if (parametrization2p25To2p4_ == nullptr)
    		return 0.;

        prefrate = parametrization2p25To2p4_->Eval(pt);
        statuncty = parametrization2p25To2p4_->GetParError(2);
  }
  else
    return 0.;


  double systuncty = prefiringRateSystUnc_ * prefrate;

  if (fluctuation == up)
    prefrate = std::max(0., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
  if (fluctuation == down)
    prefrate = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
  return prefrate;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1MuonPrefiringWeightProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("TheMuons", edm::InputTag("slimmedMuons"));
  desc.add<std::string>("L1MuonParametrizations", "L1MuonPrefiringParametriations.root");
  desc.add<std::string>("DataEra", "2016");
  desc.add<double>("PrefiringRateSystematicUncty", 0.2);
  desc.add<bool>("SkipWarnings", true);
  descriptions.add("l1MuonPrefiringWeightProducer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1MuonPrefiringWeightProducer);
