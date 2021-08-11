// C++ includes
#include <memory>
#include <string>
#include <vector>

// CMSSW includes
#include "CUDADataFormats/Common/interface/Product.h"
#include "CUDADataFormats/Muon/interface/MuonSegmentsCUDA.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCChamber.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HeterogeneousCore/CUDACore/interface/ScopedContext.h"
#include "HeterogeneousCore/CUDAServices/interface/CUDAService.h"
#include "HeterogeneousCore/CUDAUtilities/interface/copyAsync.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_noncached_unique_ptr.h"


class MuonSegmentsToCUDA : public edm::global::EDProducer<> {
public:
  explicit MuonSegmentsToCUDA(const edm::ParameterSet& iConfig);
  ~MuonSegmentsToCUDA() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;
private:

  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsGetToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsGetToken_;
  edm::EDPutTokenT<cms::cuda::Product<MuonSegmentsCUDA>> segmentsPutToken_;

};

MuonSegmentsToCUDA::MuonSegmentsToCUDA(const edm::ParameterSet& iConfig)
    : dtSegmentsGetToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("srcDT"))),
      cscSegmentsGetToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("srcCSC"))),
      segmentsPutToken_(produces<cms::cuda::Product<MuonSegmentsCUDA>>()) {

}

void MuonSegmentsToCUDA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcDT", edm::InputTag("hltDTSegments"));
  desc.add<edm::InputTag>("srcCSC", edm::InputTag("hltCSCSegments"));
  descriptions.add("MuonSegmentsToCUDA",desc);
}

void MuonSegmentsToCUDA::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {

  edm::ESHandle<DTGeometry> dtGeomHandle;
  iSetup.get<MuonGeometryRecord>().get(dtGeomHandle);
  const DTGeometry* dtGeom = &*dtGeomHandle;

  edm::ESHandle<CSCGeometry> cscGeomHandle;
  iSetup.get<MuonGeometryRecord>().get(cscGeomHandle);
  const CSCGeometry* cscGeom = &*cscGeomHandle;

  cms::cuda::ScopedContextProduce ctx{streamID};

  const DTRecSegment4DCollection& dtSegments = iEvent.get(dtSegmentsGetToken_);
  const CSCSegmentCollection& cscSegments = iEvent.get(cscSegmentsGetToken_);


  int nSegments = dtSegments.size() + cscSegments.size();

  auto segmentsCUDA = MuonSegmentsCUDA(nSegments,ctx.stream());  
  segmentsCUDA.setNSegents(nSegments);

  int offsets[12] = {0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  int index = 0;
  for (DTRecSegment4DCollection::const_iterator it = dtSegments.begin(); it != dtSegments.end(); it++) {

	segmentsCUDA.fillLocalX(index,(*it).parameters()[2]);
	segmentsCUDA.fillLocalY(index,(*it).parameters()[3]);
	segmentsCUDA.fillLocalDXDZ(index,(*it).parameters()[0]);
	segmentsCUDA.fillLocalDYDZ(index,(*it).parameters()[1]);

	segmentsCUDA.fillLocalSigmaX(index,(*it).parametersError()[2][2]);
	segmentsCUDA.fillLocalSigmaY(index,(*it).parametersError()[3][3]);
	segmentsCUDA.fillLocalSigmaDXDZ(index,(*it).parametersError()[0][0]);
	segmentsCUDA.fillLocalSigmaDYDZ(index,(*it).parametersError()[1][1]);

        DTChamberId id = (DTChamberId)(*it).chamberId();
        GlobalPoint gp = dtGeom->chamber(id)->toGlobal((*it).localPosition());

	segmentsCUDA.fillGlobalX(index,gp.x());
	segmentsCUDA.fillGlobalY(index,gp.y());
	segmentsCUDA.fillGlobalZ(index,gp.z());
	segmentsCUDA.fillGlobalR(index,pow(gp.x()*gp.x() + gp.y()*gp.y(),2));


	segmentsCUDA.fillLayerID(index,(*it).chamberId().station());
	if (offsets[(*it).chamberId().station()-1] == -1) offsets[(*it).chamberId().station()-1] = index;	
	index++;
  }

  for (CSCSegmentCollection::const_iterator it = cscSegments.begin(); it != cscSegments.end(); it++) {

	segmentsCUDA.fillLocalX(index,(*it).parameters()[2]);
	segmentsCUDA.fillLocalY(index,(*it).parameters()[3]);
	segmentsCUDA.fillLocalDXDZ(index,(*it).parameters()[0]);
	segmentsCUDA.fillLocalDYDZ(index,(*it).parameters()[1]);

	segmentsCUDA.fillLocalSigmaX(index,(*it).parametersError()[2][2]);
	segmentsCUDA.fillLocalSigmaY(index,(*it).parametersError()[3][3]);
	segmentsCUDA.fillLocalSigmaDXDZ(index,(*it).parametersError()[0][0]);
	segmentsCUDA.fillLocalSigmaDYDZ(index,(*it).parametersError()[1][1]);

        CSCDetId id = (CSCDetId)(*it).cscDetId();
        const CSCChamber* cscChamber = cscGeom->chamber(id);
        GlobalPoint gp = cscChamber->toGlobal((*it).localPosition());
	segmentsCUDA.fillGlobalX(index,gp.x());
	segmentsCUDA.fillGlobalY(index,gp.y());
	segmentsCUDA.fillGlobalZ(index,gp.z());
	segmentsCUDA.fillGlobalR(index,pow(gp.x()*gp.x() + gp.y()*gp.y(),2));

	int layerID = -1;
	if (id.zendcap() > 0) layerID = id.station() + 3;
	else layerID = id.station() + 7;

	segmentsCUDA.fillLayerID(index,layerID);
	if (offsets[layerID] == -1) offsets[layerID] = index;
	index++;
  }

  for (int i = 0; i < 12; i++){
     if (offsets[i] == -1) offsets[i] = offsets[i-1];
     segmentsCUDA.fillOffsets(i,offsets[i]);
  }

  segmentsCUDA.fillViewAndCopy(ctx.stream());
  ctx.emplace(iEvent, segmentsPutToken_, std::move(segmentsCUDA));
} 

DEFINE_FWK_MODULE(MuonSegmentsToCUDA);
