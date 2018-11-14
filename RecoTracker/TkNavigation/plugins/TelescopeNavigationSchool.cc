#include "SimpleNavigationSchool.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"

#include <vector>

/** Concrete navigation school for the Tracker, connecting disks only for traversing tracks : moslty beam halo muon 
 */
class dso_hidden TelescopeNavigationSchool final : public SimpleNavigationSchool {
public:
  
  TelescopeNavigationSchool(const GeometricSearchTracker* theTracker,
			 const MagneticField* field);
  ~TelescopeNavigationSchool() override{ cleanMemory();}

 protected:
  //addon to SimpleNavigationSchool
  void linkOtherEndLayers( SymmetricLayerFinder& symFinder);
  void addInward(const DetLayer * det, const FDLC& news);
  void addInward(const DetLayer * det, const ForwardDetLayer * newF);
  void establishInverseRelations() override;
  FDLC reachableFromHorizontal();
  std::vector<FDLC> splitForwardLayers();
};


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/DetLayers/src/DetBelowZ.h"
#include "TrackingTools/DetLayers/src/DetLessZ.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"

#include "DiskLessInnerRadius.h"
#include "SymmetricLayerFinder.h"
#include "SimpleBarrelNavigableLayer.h"
#include "SimpleForwardNavigableLayer.h"
#include "SimpleNavigableLayer.h"

using namespace std;

TelescopeNavigationSchool::TelescopeNavigationSchool(const GeometricSearchTracker* theInputTracker,
					       const MagneticField* field) 
{
  edm::LogInfo("TelescopeNavigationSchool")<<"*********Running TelescopeNavigationSchool *********";
  theBarrelLength = 0;theField = field; theTracker = theInputTracker;
  theAllDetLayersInSystem=&theInputTracker->allLayers();
  theAllNavigableLayer.resize(theInputTracker->allLayers().size(),nullptr);




  // Get barrel layers
  /*sideways does not need barrels*/
  /*  vector<BarrelDetLayer*> blc = theTracker->barrelLayers(); 
      for ( vector<BarrelDetLayer*>::iterator i = blc.begin(); i != blc.end(); i++) {
      theBarrelLayers.push_back( (*i) );
      }*/

  // get forward layers
  for( auto const& l : theTracker->forwardLayers()) {
    theForwardLayers.push_back(l);
  }
  
  FDLI middle = find_if( theForwardLayers.begin(), theForwardLayers.end(),
			 not1(DetBelowZ(0)));
  theLeftLayers  = FDLC( theForwardLayers.begin(), middle);
  theRightLayers = FDLC( middle, theForwardLayers.end());
  
  SymmetricLayerFinder symFinder( theForwardLayers);

  // only work on positive Z side; negative by mirror symmetry later
  /*sideways does not need barrels*/
  //  linkBarrelLayers( symFinder);

  linkForwardLayers( symFinder);
 
  setState(navigableLayers());

  LogDebug("TelescopeNavigationSchool")<<"inverse relation";
  establishInverseRelations();


  //add the necessary inward links to end caps
  LogDebug("TelescopeNavigationSchool")<<"linkOtherEndLayer";

  linkOtherEndLayers( symFinder);

  //set checkCrossing = false to all layers
  SimpleNavigationSchool::StateType allLayers=navigableLayers();
  SimpleNavigationSchool::StateType::iterator layerIt=allLayers.begin();
  SimpleNavigationSchool::StateType::iterator layerIt_end=allLayers.end();
  for (;layerIt!=layerIt_end;++layerIt)
    {
      //convert to SimpleNavigableLayer
      SimpleNavigableLayer* snl=dynamic_cast<SimpleNavigableLayer*>(*layerIt);
      if (!snl){
	edm::LogError("TelescopeNavigationSchool")<<"navigable layer not casting to simplenavigablelayer.";
	continue;}
      snl->setCheckCrossingSide(false);
    }

}

void TelescopeNavigationSchool::establishInverseRelations() {

  // find for each layer which are the barrel and forward
  // layers that point to it
  typedef map<const DetLayer*, vector<BarrelDetLayer const*>, less<const DetLayer*> > BarrelMapType;
  typedef map<const DetLayer*, vector<ForwardDetLayer const*>, less<const DetLayer*> > ForwardMapType;


  BarrelMapType reachedBarrelLayersMap;
  ForwardMapType reachedForwardLayersMap;

  for ( BDLI bli = theBarrelLayers.begin();
        bli!=theBarrelLayers.end(); bli++) {
    DLC reachedLC = nextLayers(**bli, insideOut);
    for ( DLI i = reachedLC.begin(); i != reachedLC.end(); i++) {
      reachedBarrelLayersMap[*i].push_back( *bli);
    }
  }

  for ( FDLI fli = theForwardLayers.begin();
        fli!=theForwardLayers.end(); fli++) {
    DLC reachedLC = nextLayers(**fli, insideOut);
    for ( DLI i = reachedLC.begin(); i != reachedLC.end(); i++) {
      reachedForwardLayersMap[*i].push_back( *fli);
    }
  }


  for ( auto const i : theTracker->allLayers()) {
    SimpleNavigableLayer* navigableLayer =
     dynamic_cast<SimpleNavigableLayer*>(theAllNavigableLayer[i->seqNum()]);
    if (!navigableLayer) {edm::LogInfo("TelescopeNavigationSchool")<<"a detlayer does not have a navigable layer, which is normal in beam halo navigation.";}
    if (navigableLayer){navigableLayer->setInwardLinks( reachedBarrelLayersMap[i],reachedForwardLayersMap[i], TkLayerLess(outsideIn, i) );}
  }

}


void TelescopeNavigationSchool::
linkOtherEndLayers(  SymmetricLayerFinder& symFinder){

  LogDebug("TelescopeNavigationSchool")<<"reachable from horizontal";
  //generally, on the right side, what are the forward layers reachable from the horizontal
  FDLC reachableFL= reachableFromHorizontal();

  //even simpler navigation from end to end.
  //for each of them
  for (FDLI fl=reachableFL.begin();fl!=reachableFL.end();fl++)
    {
      LogDebug("TelescopeNavigationSchool")<<"adding inward from right";
      //link it inward to the mirror reachable from horizontal
      addInward(static_cast<DetLayer const*>(*fl),symFinder.mirror(*fl));
      
      LogDebug("TelescopeNavigationSchool")<<"adding inward from mirror of right (left?)";
      addInward(static_cast<DetLayer const*>(symFinder.mirror(*fl)),*fl);
    }




}

void TelescopeNavigationSchool::
addInward(const DetLayer * det, const ForwardDetLayer * newF){
  //get the navigable layer for this DetLayer
  SimpleNavigableLayer* navigableLayer =
    dynamic_cast<SimpleNavigableLayer*>(theAllNavigableLayer[(det)->seqNum()]);

  LogDebug("TelescopeNavigationSchool")<<"retreive the nextlayer outsidein";
  //get the inward reachable layers.
  DLC inwardsLayers(navigableLayer->nextLayers(outsideIn));
  
  LogDebug("TelescopeNavigationSchool")<<"split them barrel/forward";
  // split barrel and forward layers
  BDLC inwardsBarrel;
  FDLC inwardsForward;
  for ( DLC::iterator dli=inwardsLayers.begin();dli!=inwardsLayers.end();dli++)
    {
      if ((**dli).location()==GeomDetEnumerators::barrel)
        inwardsBarrel.push_back(static_cast<const BarrelDetLayer*>(*dli));
      else
        inwardsForward.push_back(static_cast<const ForwardDetLayer*>(*dli));
    }
  LogDebug("TelescopeNavigationSchool")<<"add the new ones";
  //add the other forward layers provided
  inwardsForward.push_back(newF);

  LogDebug("TelescopeNavigationSchool")<<"no duplicate please";
  sort(inwardsForward.begin(),inwardsForward.end()); //if you don't sort, unique will not work
  //  FDLI read = inwardsForward.begin();
  //  std::stringstream showMe;
  //  for (;read !=inwardsForward.end();++read)  showMe<<" layer p: "<<*read<<"\n";
  //  LogDebug("TelescopeNavigationSchool")<<"list of layer pointers: \n"<<showMe.str();

  FDLI new_end =unique(inwardsForward.begin(),inwardsForward.end());
  //  if (new_end!=inwardsForward.end()) LogDebug("TelescopeNavigationSchool")<<"removing duplicates here";
  inwardsForward.erase(new_end,inwardsForward.end());

  LogDebug("TelescopeNavigationSchool")<<"set back the inward links (no duplicate)";
  //  set them back to the navigable layer
  navigableLayer->setInwardLinks( inwardsBarrel, inwardsForward, TkLayerLess(outsideIn, det));
}

void TelescopeNavigationSchool::
addInward(const DetLayer * det, const FDLC& news){
  //get the navigable layer for this DetLayer
  SimpleNavigableLayer* navigableLayer =
    dynamic_cast<SimpleNavigableLayer*>(theAllNavigableLayer[(det)->seqNum()]);

  LogDebug("TelescopeNavigationSchool")<<"retreive the nextlayer outsidein";
  //get the inward reachable layers.
  DLC inwardsLayers(navigableLayer->nextLayers(outsideIn));

  LogDebug("TelescopeNavigationSchool")<<"split them barrel/forward";
  // split barrel and forward layers
  BDLC inwardsBarrel;
  FDLC inwardsForward;
  for ( DLC::iterator dli=inwardsLayers.begin();dli!=inwardsLayers.end();dli++)
    {
      if ((**dli).location()==GeomDetEnumerators::barrel)
	inwardsBarrel.push_back(static_cast<const BarrelDetLayer*>(*dli));
      else
	inwardsForward.push_back(static_cast<const ForwardDetLayer*>(*dli));
    }
  
  LogDebug("TelescopeNavigationSchool")<<"add the new ones";
  //add the other forward layers provided
  inwardsForward.insert( inwardsForward.end(), news.begin(), news.end());

  LogDebug("TelescopeNavigationSchool")<<"no duplicate please";
  FDLI new_end =unique(inwardsForward.begin(),inwardsForward.end());
  inwardsForward.erase(new_end,inwardsForward.end());

  LogDebug("TelescopeNavigationSchool")<<"set back the inward links (no duplicate)";
  //  set them back to the navigable layer
  navigableLayer->setInwardLinks( inwardsBarrel, inwardsForward, TkLayerLess(outsideIn, det));
}

TelescopeNavigationSchool::FDLC
TelescopeNavigationSchool::reachableFromHorizontal()
{    
  //determine which is the list of forward layers that can be reached from inside-out
  //at horizontal direction

  FDLC myRightLayers( theRightLayers);
  FDLI begin = myRightLayers.begin();
  FDLI end   = myRightLayers.end();

  //sort along Z to be sure
  sort(begin, end, DetLessZ());

  FDLC reachableFL;

  begin = myRightLayers.begin();
  end   = myRightLayers.end();

  //the first one is always reachable
  reachableFL.push_back(*begin);
  FDLI current = begin;
  for (FDLI i = begin+1; i!= end; i++)
    {
      //is the previous layer NOT masking this one
      //inner radius smaller OR outer radius bigger
      if ((**i).specificSurface().innerRadius() < (**current).specificSurface().innerRadius() ||
	  (**i).specificSurface().outerRadius() > (**current).specificSurface().outerRadius())
	{	  //not masked
	  reachableFL.push_back(*i);
	  current=i;
	}
    }
  return reachableFL;
}
vector<TelescopeNavigationSchool::FDLC> 
TelescopeNavigationSchool::splitForwardLayers() 
{
  // only work on positive Z side; negative by mirror symmetry later
  FDLC myRightLayers( theRightLayers);
  FDLI begin = myRightLayers.begin();
  FDLI end   = myRightLayers.end();

  // sort according to inner radius, but keeping the ordering in z!
  stable_sort ( begin, end, DiskLessInnerRadius());
  // partition in cylinders
  vector<FDLC> result;
  FDLC current;
  current.push_back( *begin);
  for ( FDLI i = begin+1; i != end; i++) {

#ifdef EDM_ML_DEBUG
    LogDebug("TkNavigation") << "(**i).specificSurface().innerRadius()      = "
			     << (**i).specificSurface().innerRadius() << endl
			     << "(**(i-1)).specificSurface().outerRadius()) = "
			     << (**(i-1)).specificSurface().outerRadius() ;
    LogDebug("TkNavigation") << "(**i).specificSurface().position().z()      = "
                            << (**i).specificSurface().position().z() << endl
                            << "(**(i-1)).specificSurface().position().z() = "
                            << (**(i-1)).specificSurface().position().z() ;
#endif

    // if inner radius of i is larger than outer radius of i-1 then split!
    // FIXME: The solution found for phase2 is a bit dirty, we can do better.
    // For phase2 we compare the EXTENDED pixel with the TID to get the assignment right!
      if ( (**i).specificSurface().innerRadius() > (**(i-1)).specificSurface().outerRadius()){
      LogDebug("TkNavigation") << "found break between groups" ;

      // sort layers in group along Z
      stable_sort ( current.begin(), current.end(), DetLessZ());

      result.push_back(current);
      current.clear();
    }
    current.push_back(*i);
  }
  result.push_back(current); // save last one too 
  // now sort subsets in Z
  for ( vector<FDLC>::iterator ivec = result.begin();
	ivec != result.end(); ivec++) {
    stable_sort( ivec->begin(), ivec->end(), DetLessZ());
  }

  return result;
}
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

 
#include "NavigationSchoolFactory.h"
#include "TrackingTools/DetLayers/interface/NavigationSchool.h"
DEFINE_EDM_PLUGIN(NavigationSchoolFactory, TelescopeNavigationSchool, "TelescopeNavigationSchool");

