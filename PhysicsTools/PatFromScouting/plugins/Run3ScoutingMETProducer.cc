// -*- C++ -*-
//
// Package:    PhysicsTools/PatFromScouting
// Class:      Run3ScoutingMETProducer
//
/**\class Run3ScoutingMETProducer Run3ScoutingMETProducer.cc PhysicsTools/PatFromScouting/plugins/Run3ScoutingMETProducer.cc

 Description: Creates pat::MET from scouting MET (pt, phi) stored in event

 Implementation:
     Reads precomputed MET pt and phi from hltScoutingPFPacker and creates pat::MET
*/
//
// Original Author:  Dmytro Kovalskyi
//         Created:  Thu, 05 Dec 2024 15:27:09 GMT
//
//

#include <memory>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

class Run3ScoutingMETProducer : public edm::stream::EDProducer<> {
public:
  explicit Run3ScoutingMETProducer(const edm::ParameterSet&);
  ~Run3ScoutingMETProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<double> metPtToken_;
  const edm::EDGetTokenT<double> metPhiToken_;
  // Opt-in, default off (empty InputTag): embedding a GenMET is only
  // needed/possible on MC and only if the caller also wires up the
  // genMetTrue-equivalent production chain (see ScoutingPuppiCalibration's
  // _addMETTables, which builds it from packedGenParticles since this is a
  // MiniAOD-derived workflow, not the AOD-level genParticles the standard
  // RecoMET/Configuration/python/GenMETParticles_cff.py recipe expects).
  // Default-off so every other user of this shared plugin (ScoutingNanoProduction,
  // DeepNTuples, ...) is completely unaffected -- same pattern as
  // Run3ScoutingParticleToPackedCandidateProducer's useImprovedVertexAssociation.
  const bool hasGenMET_;
  edm::EDGetTokenT<reco::GenMETCollection> genMETToken_;
};

Run3ScoutingMETProducer::Run3ScoutingMETProducer(const edm::ParameterSet& iConfig)
    : metPtToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("metPt"))),
      metPhiToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("metPhi"))),
      hasGenMET_(!iConfig.getParameter<edm::InputTag>("genMET").label().empty()) {
  if (hasGenMET_) {
    genMETToken_ = consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMET"));
  }
  produces<pat::METCollection>();
}

void Run3ScoutingMETProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto patMETs = std::make_unique<pat::METCollection>();

  double metPt = iEvent.get(metPtToken_);
  double metPhi = iEvent.get(metPhiToken_);

  double metPx = metPt * std::cos(metPhi);
  double metPy = metPt * std::sin(metPhi);

  // sumEt is not available in scouting, use metPt as approximation
  double sumEt = metPt;

  reco::MET::LorentzVector p4(metPx, metPy, 0.0, metPt);
  reco::MET::Point vtx(0.0, 0.0, 0.0);

  reco::MET recoMET(sumEt, p4, vtx);
  pat::MET patMET(recoMET);

  // Initialize MET corrections to make NanoAOD happy
  // Using the same values since scouting MET is already the "raw" PF MET
  patMET.setCorShift(metPx, metPy, sumEt, pat::MET::None);
  patMET.setCorShift(metPx, metPy, sumEt, pat::MET::T1);
  patMET.setCorShift(metPx, metPy, sumEt, pat::MET::Calo);
  patMET.setCorShift(metPx, metPy, sumEt, pat::MET::Chs);
  patMET.setCorShift(metPx, metPy, sumEt, pat::MET::Trk);

  if (hasGenMET_) {
    const auto& genMETs = iEvent.get(genMETToken_);
    if (!genMETs.empty()) {
      patMET.setGenMET(genMETs.front());
    }
  }

  patMETs->push_back(patMET);

  iEvent.put(std::move(patMETs));
}

void Run3ScoutingMETProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("metPt", edm::InputTag("hltScoutingPFPacker", "pfMetPt"));
  desc.add<edm::InputTag>("metPhi", edm::InputTag("hltScoutingPFPacker", "pfMetPhi"));
  desc.add<edm::InputTag>("genMET", edm::InputTag(""))
      ->setComment("Optional reco::GenMETCollection to embed via pat::MET::setGenMET(), for "
                   "NanoAOD's standard GenMET table (PhysicsTools/NanoAOD/python/met_cff.py's "
                   "metMCTable, which reads pat::MET::genMET()). Empty (default) skips this "
                   "entirely -- MC-only, and no genMET production chain exists upstream for this "
                   "HLT-scouting-derived MET, unlike a standard offline MiniAOD's slimmedMETs.");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(Run3ScoutingMETProducer);
