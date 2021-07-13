// -*- C++ -*-
//
// Package:    L1Trigger/SingleMuL1Filter
// Class:      SingleMuL1Filter
// 
/**\class SingleMuL1Filter SingleMuL1Filter.cc L1Trigger/SingleMuL1Filter/plugins/SingleMuL1Filter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jan-Frederik Schulte
//         Created:  Thu, 01 Jul 2021 14:13:21 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
//
// class declaration
//

class SingleMuL1Filter : public edm::stream::EDFilter<> {
   public:
      explicit SingleMuL1Filter(const edm::ParameterSet&);
      ~SingleMuL1Filter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      edm::EDGetTokenT<BXVector<l1t::Muon>>    l1Token_;
      edm::EDGetTokenT<reco::MuonCollection>   muonToken_;
      edm::EDGetTokenT<reco::VertexCollection>      vtx_token;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
SingleMuL1Filter::SingleMuL1Filter(const edm::ParameterSet& iConfig):
	l1Token_(consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("l1Src"))),
        muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
        vtx_token( consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices")) )
{
   //now do what ever initialization is needed

}


SingleMuL1Filter::~SingleMuL1Filter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SingleMuL1Filter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  using namespace reco;

  // Primary vertex
  edm::Handle<reco::VertexCollection > vtx_handle;
  iEvent.getByToken(vtx_token,vtx_handle);
  const reco::Vertex& pv(*(vtx_handle->begin()));


   edm::Handle<reco::MuonCollection> muonHandle;
   iEvent.getByToken(muonToken_, muonHandle);
   bool hasLooseMuon = false;
   if (muonHandle->size() == 0) return false;
   for ( unsigned int i=0; i!=muonHandle->size(); ++i){
      const auto& muon = muonHandle->at(i);
      if (muon::isTightMuon(muon,pv)) hasLooseMuon = true;
   }
   if (!hasLooseMuon) return false;
   edm::Handle<BXVector<l1t::Muon>>   l1Handle_;
   iEvent.getByToken(l1Token_, l1Handle_);   
   int nL1Muons = 0;
   for (int ibx = l1Handle_->getFirstBX(); ibx <= l1Handle_->getLastBX(); ++ibx) {
     for (auto it = l1Handle_->begin(ibx); it != l1Handle_->end(ibx); it++){
	if ((*it).pt() > 20) return false;
	nL1Muons++;
    }
   }

   if (nL1Muons != 1) return false;
   return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
SingleMuL1Filter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SingleMuL1Filter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
SingleMuL1Filter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
SingleMuL1Filter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SingleMuL1Filter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SingleMuL1Filter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SingleMuL1Filter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1Src",edm::InputTag("gmtStage2Digis:Muon"))->setComment("L1 muon collection");;
  desc.add<edm::InputTag>("muonSrc",edm::InputTag("muons"))->setComment("reco muon collection");;
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SingleMuL1Filter);
