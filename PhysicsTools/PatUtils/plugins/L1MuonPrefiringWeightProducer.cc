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

  produces<double>("nonPrefiringProbMuon").setBranchAlias("nonPrefiringProbMuon");
  produces<double>("nonPrefiringProbMuonUp").setBranchAlias("nonPrefiringProbMuonUp");
  produces<double>("nonPrefiringProbMuonDown").setBranchAlias("nonPrefiringProbMuonDown");
}

L1MuonPrefiringWeightProducer::~L1MuonPrefiringWeightProducer() {
  delete file_prefiringparams_;
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


  
  TF1* parametrization;
  if (std::abs(eta) < 0.2){ 
	TString paramName = "L1prefiring_muonparam_0.0To0.2_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 0.3){ 
	TString paramName = "L1prefiring_muonparam_0.2To0.3_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 0.55){ 
	TString paramName = "L1prefiring_muonparam_0.3To0.55_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 0.83){ 
	TString paramName = "L1prefiring_muonparam_0.55To0.83_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 1.24){ 
	TString paramName = "L1prefiring_muonparam_0.83To1.24_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 1.4){ 
	TString paramName = "L1prefiring_muonparam_1.24To1.4_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 1.6){ 
	TString paramName = "L1prefiring_muonparam_1.4To1.6_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 1.8){ 
	TString paramName = "L1prefiring_muonparam_1.6To1.8_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 2.1){ 
	TString paramName = "L1prefiring_muonparam_1.8To2.1_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 2.25){ 
	TString paramName = "L1prefiring_muonparam_2.1To2.25_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else if (std::abs(eta) < 2.4){ 
	TString paramName = "L1prefiring_muonparam_2.25To2.4_" + dataera_;
        parametrization = (TF1*)file_prefiringparams_->Get(paramName);
  }
  else
    return 0.;

  if (parametrization == nullptr && !skipwarnings_)
    std::cout << "Prefiring parametrization not found, setting prefiring rate to 0 " << eta << " " << std::endl;
  if (parametrization == nullptr)
    return 0.;


  double prefrate = parametrization->Eval(pt);
  double statuncty = parametrization->GetParError(2);
  double systuncty = prefiringRateSystUnc_ * prefrate;


  if (fluctuation == up)
    prefrate = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
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
