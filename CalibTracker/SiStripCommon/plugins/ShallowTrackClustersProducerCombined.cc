#include "CalibTracker/SiStripCommon/interface/ShallowTrackClustersProducerCombined.h"

#include "CalibTracker/SiStripCommon/interface/ShallowTools.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiStripLorentzAngleRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "boost/foreach.hpp"
#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/SiStripDigi/interface/SiStripProcessedRawDigi.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CalibTracker/Records/interface/SiStripDependentRecords.h"
#include <map>
#include <iostream>
#include <fstream>
#include "TMath.h"

using namespace std;

unsigned int ntracks2 =0;

  ofstream myfile2;

ShallowTrackClustersProducerCombined::ShallowTrackClustersProducerCombined(const edm::ParameterSet& iConfig)
  :  tracks_token_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("Tracks"))),
     association_token_(consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("Tracks"))),
     clusters_token_( consumes< edmNew::DetSetVector<SiStripCluster> >( iConfig.getParameter<edm::InputTag>("Clusters") ) ),
     theVertexToken_(consumes<std::vector<reco::Vertex> >          (iConfig.getParameter<edm::InputTag>("vertices"))),
     theDigisToken_    (consumes<edm::DetSetVector<SiStripProcessedRawDigi> > (edm::InputTag("siStripProcessedRawDigis", ""))),
     Suffix       ( iConfig.getParameter<std::string>("Suffix")    ),
     Prefix       ( iConfig.getParameter<std::string>("Prefix") ),
     isData       ( iConfig.getParameter<bool>("isData") )
{
  produces <std::vector<unsigned> >      ( Prefix + "StripIdx"        );
  produces <std::vector<float> >         ( Prefix + "Amplitudes"        );
  produces <std::vector<unsigned> >      ( Prefix + "Idx"        );
  produces <std::vector<float> >         ( Prefix + "Charge"        );
  produces <std::vector<unsigned> >      ( Prefix + "Width"        );
  produces <std::vector<float> >         ( Prefix + "LocalTrackPhi"        );
  produces <std::vector<float> >         ( Prefix + "LocalTrackTheta"        );
  produces <std::vector<int> >           ( Prefix + "Subdetid"        );
  produces <std::vector<int> >           ( Prefix + "Layerwheel"        );
  produces <std::vector<unsigned> >      ( Prefix + "Detid"        );
  produces <std::vector<float> >         ( Prefix + "LocalTrackX"        );
  produces <std::vector<float> >         ( Prefix + "LocalTrackY"        );
  produces <std::vector<float> >         ( Prefix + "LocalTrackZ"        );
  produces <std::vector<float> >         ( Prefix + "Localpitch"        );
  produces <std::vector<float> >         ( Prefix + "SensorThickness"        );
  produces <std::vector<bool> >          ( Prefix + "Saturation"        );
  produces <std::vector<bool> >          ( Prefix + "Overlapping"        );
  produces <std::vector<bool> >          ( Prefix + "Farfromedge"        );
  produces <std::vector<float> >         ( Prefix + "Path"        );
  produces <std::vector<float> >         ( Prefix + "BdotX"        );
  produces <std::vector<float> >         ( Prefix + "BdotY"        );
  produces <std::vector<float> >         ( Prefix + "BdotZ"        );
  produces <std::vector<float> >         ( Prefix + "BdotMag"        );
  produces  <std::vector<unsigned int>> (  Prefix +"Event" );
  produces  <std::vector<unsigned int>> (  Prefix +"StripEvent" );
}

void ShallowTrackClustersProducerCombined::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  shallow::CLUSTERMAP clustermap = shallow::make_cluster_map(iEvent, clusters_token_);
  edm::Handle<edm::View<reco::Track> > tracks;	             iEvent.getByToken(tracks_token_, tracks);	  


  int size = clustermap.size();
  //cout << "track size " << tracks->size() << "clustermap size " << size <<  endl;

  auto StripIdx  = std::make_unique <std::vector<unsigned> >();
  auto Amplitudes  = std::make_unique <std::vector<float> >();
  auto Idx  = std::make_unique <std::vector<unsigned> >();
  auto Charge  = std::make_unique<std::vector<float> >();
  auto Width  = std::make_unique<std::vector<unsigned> >();
  auto LocalTrackPhi  = std::make_unique<std::vector<float> >();
  auto LocalTrackTheta  = std::make_unique<std::vector<float> >();
  auto Subdetid  = std::make_unique<std::vector<int> >();
  auto Layerwheel  = std::make_unique<std::vector<int> >();
  auto Detid  = std::make_unique<std::vector<unsigned> >();
  auto LocalTrackX  = std::make_unique<std::vector<float> >();
  auto LocalTrackY  = std::make_unique<std::vector<float> >();
  auto LocalTrackZ  = std::make_unique<std::vector<float> >();
  auto Localpitch  = std::make_unique<std::vector<float> >();
  auto SensorThickness  = std::make_unique<std::vector<float> >();
  auto Saturation  = std::make_unique<std::vector<bool> >();
  auto Overlapping  = std::make_unique<std::vector<bool> >();
  auto Farfromedge  = std::make_unique<std::vector<bool> >();
  auto Path  = std::make_unique<std::vector<float> >();
  auto BdotX  = std::make_unique<std::vector<float> >();
  auto BdotY  = std::make_unique<std::vector<float> >();
  auto BdotZ  = std::make_unique<std::vector<float> >();
  auto BdotMag  = std::make_unique<std::vector<float> >();
  auto Event = std::make_unique <std::vector<unsigned> >();
  auto StripEvent = std::make_unique <std::vector<unsigned> >();

  edm::ESHandle<TrackerGeometry> theTrackerGeometry;         iSetup.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );  
  edm::ESHandle<MagneticField> magfield;		     iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  edm::ESHandle<SiStripLorentzAngle> SiStripLorentzAngle;    iSetup.get<SiStripLorentzAngleDepRcd>().get(SiStripLorentzAngle);

  edm::Handle<TrajTrackAssociationCollection> associations;  iEvent.getByToken(association_token_, associations);

    edm::Handle<edm::DetSetVector<SiStripProcessedRawDigi> > rawProcessedDigis;
    iEvent.getByToken(theDigisToken_,rawProcessedDigis);


  edm::Handle<std::vector<reco::Vertex> > vtx;
  iEvent.getByToken(theVertexToken_, vtx); 

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  const TrackerTopology* const tTopo = tTopoHandle.product();


  TrajectoryStateCombiner combiner;

	size_t ontrk_cluster_idx=0;
  std::map<size_t, std::vector<size_t> > mapping; //cluster idx --> on trk cluster idx (multiple)


  for( TrajTrackAssociationCollection::const_iterator association = associations->begin(); 
       association != associations->end(); association++) {
    const Trajectory*  traj  = association->key.get();
    const reco::Track* track = association->val.get();
		int trk_idx = shallow::findTrackIndex(tracks, track); 
		size_t trk_strt_idx = ontrk_cluster_idx;

    bool isGoodTrack = trackFilter(track);
    if(!isGoodTrack)
        continue;
   

    BOOST_FOREACH( const TrajectoryMeasurement measurement, traj->measurements() ) {
      const TrajectoryStateOnSurface tsos = measurement.updatedState();
      const TrajectoryStateOnSurface unbiased = combiner(measurement.forwardPredictedState(), measurement.backwardPredictedState());

      const TrackingRecHit*         hit        = measurement.recHit()->hit();
      const SiStripRecHit1D*        hit1D      = dynamic_cast<const SiStripRecHit1D*>(hit);
      const SiStripRecHit2D*        hit2D      = dynamic_cast<const SiStripRecHit2D*>(hit);
      const SiStripMatchedRecHit2D* matchedhit = dynamic_cast<const SiStripMatchedRecHit2D*>(hit);

      for(unsigned h=0; h<2; h++) { //loop over possible Hit options (1D, 2D)
				const SiStripCluster* cluster_ptr;
				if(!matchedhit && h==1) continue; 
				else if( matchedhit && h==0) cluster_ptr = &matchedhit->monoCluster(); 
				else if( matchedhit && h==1) cluster_ptr = &matchedhit->stereoCluster(); 
				else if(hit2D) cluster_ptr = (hit2D->cluster()).get(); 
				else if(hit1D) cluster_ptr = (hit1D->cluster()).get(); 
				else continue;

				shallow::CLUSTERMAP::const_iterator cluster = clustermap.find( std::make_pair( hit->geographicalId().rawId(), cluster_ptr->firstStrip() ));
				if(cluster == clustermap.end() ) throw cms::Exception("Logic Error") << "Cluster not found: this could be a configuration error" << std::endl;
	
				unsigned i = cluster->second;


                                uint32_t id = hit->geographicalId();
                                const moduleVars moduleV(id, tTopo);
                                const SiStripClusterInfo info(*cluster_ptr, iSetup, id);
                                //const NearDigis digis = rawProcessedDigis.isValid() ? NearDigis(info, *rawProcessedDigis) : NearDigis(info);
                                

				//find if cluster was already assigned to a previous track
				auto already_visited = mapping.find(i);
				int nassociations = 1;
				if(already_visited != mapping.end()) {
					nassociations += already_visited->second.size();
					/*for(size_t idx : already_visited->second) {
						trackmulti->at(idx)++;
					}*/
					already_visited->second.push_back(ontrk_cluster_idx);
				}
				else { //otherwise store this 
					std::vector<size_t> single = {ontrk_cluster_idx};
					mapping.insert( std::make_pair(i, single) );
				}

				const StripGeomDetUnit* theStripDet = dynamic_cast<const StripGeomDetUnit*>( theTrackerGeometry->idToDet( hit->geographicalId() ) );
	
 			//LocalVector drift = shallow::drift( theStripDet, *magfield, *SiStripLorentzAngle);

                                uint32_t eventNr = iEvent.id().event();				
				Idx->push_back( ontrk_cluster_idx  );  //link: on trk cluster --> general cluster info
				Event->push_back( eventNr  );  //link: on trk cluster --> general cluster info
                                cout<< "cluster idx " << ontrk_cluster_idx << " event " << eventNr << endl;             
                      //float langle = (SiStripLorentzAngle.isValid()) ? SiStripLorentzAngle->getLorentzAngle(id) : 0.;
                      //lorentzAngle->push_back(langle);
                      Width->push_back(        cluster_ptr->amplitudes().size()                              );
		      //barystrip->push_back(    cluster_ptr->barycenter()                                     );
		      //middlestrip->push_back(  info.firstStrip() + info.width()/2.0                    );
		      Charge->push_back(       info.charge()                                           );
		      //noise->push_back(        info.noiseRescaledByGain()                              );
		      //ston->push_back(         info.signalOverNoise()                                  );
		      //seedstrip->push_back(    info.maxStrip()                                         );
		      //seedindex->push_back(    info.maxIndex()                                         );
		      //seedcharge->push_back(   info.maxCharge()                                        );
		      //seednoise->push_back(    info.stripNoisesRescaledByGain().at(info.maxIndex())   );
		      //seednoisepure->push_back(     info.stripNoises().at(info.maxIndex())                  );
		      //seedgain->push_back(     info.stripGains().at(info.maxIndex())                  );
	 

	              LocalTrackTheta->push_back(  (theStripDet->toLocal(tsos.globalDirection())).theta() ); 
		      LocalTrackPhi->push_back(    (theStripDet->toLocal(tsos.globalDirection())).phi() );   
		      LocalTrackX->push_back(      (theStripDet->toLocal(tsos.globalPosition())).x() );    
		      LocalTrackY->push_back(      (theStripDet->toLocal(tsos.globalPosition())).y() );    
		      LocalTrackZ->push_back(      (theStripDet->toLocal(tsos.globalPosition())).z() );    
		      BdotX->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).x() );
		      BdotY->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).y() );
		      BdotZ->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).z() );
		      BdotMag->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).mag() );


		      Detid->push_back(            id                 );
		      Subdetid->push_back(         moduleV.subdetid      );
		      //side->push_back(             moduleV.side          );
		      //module->push_back(           moduleV.module        );
		      Layerwheel->push_back(       moduleV.layerwheel    );
		      //stringringrod->push_back(    moduleV.stringringrod );
		      //petal->push_back(            moduleV.petal         );
		      //stereo->push_back(           moduleV.stereo        );


                     //stripLength->push_back(           theStripDet->specificTopology().stripLength() );
                     SensorThickness->push_back(             theStripDet->specificSurface().bounds().thickness() );
                     Localpitch->push_back(  (theStripDet->specificTopology()).localPitch(theStripDet->toLocal(tsos.globalPosition())) ); 

                      bool sat = false;
		      for( auto itAmpl = cluster_ptr->amplitudes().begin(); itAmpl != cluster_ptr->amplitudes().end(); ++itAmpl){
			  Amplitudes->push_back(*itAmpl); 
			  StripIdx->push_back( ontrk_cluster_idx );
			  StripEvent->push_back( eventNr );
                                cout<< "cluster strip idx " << ontrk_cluster_idx << " event  " << eventNr << endl;             
                          if(*itAmpl >=254)
                              sat =true;
		      }
                      Saturation->push_back( sat);
                      
                      bool overlapping = false;
                      auto FirstStrip = info.firstStrip();

		       if(FirstStrip==0                                  )overlapping=true;
		       if(FirstStrip==128                                )overlapping=true;
		       if(FirstStrip==256                                )overlapping=true;
		       if(FirstStrip==384                                )overlapping=true;
		       if(FirstStrip==512                                )overlapping=true;
		       if(FirstStrip==640                                )overlapping=true;

		       if(FirstStrip<=127 && FirstStrip+cluster_ptr->amplitudes().size()>127)overlapping=true;
		       if(FirstStrip<=255 && FirstStrip+cluster_ptr->amplitudes().size()>255)overlapping=true;
		       if(FirstStrip<=383 && FirstStrip+cluster_ptr->amplitudes().size()>383)overlapping=true;
		       if(FirstStrip<=511 && FirstStrip+cluster_ptr->amplitudes().size()>511)overlapping=true;
		       if(FirstStrip<=639 && FirstStrip+cluster_ptr->amplitudes().size()>639)overlapping=true;

		       if(FirstStrip+cluster_ptr->amplitudes().size()==127                   )overlapping=true;
		       if(FirstStrip+cluster_ptr->amplitudes().size()==255                   )overlapping=true;
		       if(FirstStrip+cluster_ptr->amplitudes().size()==383                   )overlapping=true;
		       if(FirstStrip+cluster_ptr->amplitudes().size()==511                   )overlapping=true;
		       if(FirstStrip+cluster_ptr->amplitudes().size()==639                   )overlapping=true;
	               if(FirstStrip+cluster_ptr->amplitudes().size()==767                   )overlapping=true;

       
                      Overlapping->push_back(overlapping);
                      Farfromedge->push_back(  IsFarFromBorder(&tsos ,id , &iSetup)    );

                      LocalVector             trackDirection = tsos.localDirection();
                      double cosine = trackDirection.z()/trackDirection.mag();
                      double path = (10.0*thickness(id,  &iSetup))/fabs(cosine);
                      Path->push_back( path );

                      ontrk_cluster_idx++;
				//strip->push_back(       (theStripDet->specificTopology()).strip(theStripDet->toLocal(tsos.globalPosition())) );
				//globaltheta->push_back( tsos.globalDirection().theta() );                       
				//globalphi->push_back(   tsos.globalDirection().phi() );                         
				//globalx->push_back(     tsos.globalPosition().x() );                            
				//globaly->push_back(     tsos.globalPosition().y() );                            
				//globalz->push_back(     tsos.globalPosition().z() );                            
				//transverseCurvature->push_back(     tsos.transverseCurvature() );      
                                //trackPt->push_back(   track->pt() );                         
                                //trackEta->push_back(   track->eta() );                         
				//projwidth->push_back(   tan(localtheta->at(ontrk_cluster_idx))*cos(localphi->at(ontrk_cluster_idx)) );         
      } //for(unsigned h=0; h<2; h++) { //loop over possible Hit options (1D, 2D)
    } //BOOST_FOREACH( const TrajectoryMeasurement measurement, traj->measurements() )

  } //for(TrajTrackAssociationCollection::const_iterator association = associations->begin();

  iEvent.put(std::move(StripIdx),       Prefix +  "StripIdx"        );
  iEvent.put(std::move(Amplitudes),       Prefix +  "Amplitudes"        );
  iEvent.put(std::move(Idx),        Prefix + "Idx"        );
  iEvent.put(std::move(Charge),       Prefix +  "Charge"        );
  iEvent.put(std::move(Width),       Prefix +  "Width"        );
  iEvent.put(std::move(LocalTrackPhi),        Prefix + "LocalTrackPhi"        );
  iEvent.put(std::move(LocalTrackTheta),       Prefix +  "LocalTrackTheta"        );
  iEvent.put(std::move(Subdetid),    Prefix +     "Subdetid"        );
  iEvent.put(std::move(Layerwheel),      Prefix +   "Layerwheel"        );
  iEvent.put(std::move(Detid),      Prefix +   "Detid"        );
  iEvent.put(std::move(LocalTrackX),       Prefix +  "LocalTrackX"        );
  iEvent.put(std::move(LocalTrackY),       Prefix +  "LocalTrackY"        );
  iEvent.put(std::move(LocalTrackZ),        Prefix + "LocalTrackZ"        );
  iEvent.put(std::move(Localpitch),        Prefix + "Localpitch"        );
  iEvent.put(std::move(SensorThickness),    Prefix +     "SensorThickness"        );
  iEvent.put(std::move(Saturation),      Prefix +   "Saturation"        );
  iEvent.put(std::move(Overlapping),      Prefix +   "Overlapping"        );
  iEvent.put(std::move(Farfromedge),       Prefix +  "Farfromedge"        );
  iEvent.put(std::move(Path),       Prefix +  "Path"        );
  iEvent.put(std::move(BdotX),     Prefix +    "BdotX"        );
  iEvent.put(std::move(BdotY),       Prefix +  "BdotY"        );
  iEvent.put(std::move(BdotZ),       Prefix +  "BdotZ"        );
  iEvent.put(std::move(BdotMag),       Prefix +  "BdotMag"        );
  iEvent.put(std::move(Event),  Prefix + "Event" );
  iEvent.put(std::move(StripEvent),  Prefix + "StripEvent" );
}

bool ShallowTrackClustersProducerCombined::trackFilter(const reco::Track* trk)
{
  //if (trk->pt() < 0.8) return false;
  //if (trk->p()  < 2.0) return false;
  if (trk->hitPattern().numberOfValidTrackerHits()  <= 6) return false;
  if (trk->normalizedChi2() > 10.0) return false;
  //check PV compatibility ??
  return true;
}


ShallowTrackClustersProducerCombined::NearDigis::
NearDigis(const SiStripClusterInfo& info) {
  max =  info.maxCharge();
  left =           info.maxIndex()    > uint16_t(0)                ? info.stripCharges()[info.maxIndex()-1]      : 0 ;
  Lleft =          info.maxIndex()    > uint16_t(1)                ? info.stripCharges()[info.maxIndex()-2]      : 0 ;
  right=  unsigned(info.maxIndex()+1) < info.stripCharges().size() ? info.stripCharges()[info.maxIndex()+1]      : 0 ;
  Rright= unsigned(info.maxIndex()+2) < info.stripCharges().size() ? info.stripCharges()[info.maxIndex()+2]      : 0 ;
  first = info.stripCharges()[0];
  last =  info.stripCharges()[info.width()-1];
}

ShallowTrackClustersProducerCombined::NearDigis::
NearDigis(const SiStripClusterInfo& info, const edm::DetSetVector<SiStripProcessedRawDigi>& rawProcessedDigis) {
  edm::DetSetVector<SiStripProcessedRawDigi>::const_iterator digiframe = rawProcessedDigis.find(info.detId());
  if( digiframe != rawProcessedDigis.end()) {
    max =                                                            digiframe->data.at(info.maxStrip()).adc()       ;
    left =            info.maxStrip()    > uint16_t(0)             ? digiframe->data.at(info.maxStrip()-1).adc() : 0 ;
    Lleft =           info.maxStrip()    > uint16_t(1)             ? digiframe->data.at(info.maxStrip()-2).adc() : 0 ;
    right =  unsigned(info.maxStrip()+1) < digiframe->data.size()  ? digiframe->data.at(info.maxStrip()+1).adc() : 0 ;
    Rright = unsigned(info.maxStrip()+2) < digiframe->data.size()  ? digiframe->data.at(info.maxStrip()+2).adc() : 0 ;
    first = digiframe->data.at(info.firstStrip()).adc();
    last = digiframe->data.at(info.firstStrip()+info.width() - 1).adc();
  } else {
    *this = NearDigis(info);
  }
}

ShallowTrackClustersProducerCombined::moduleVars::
moduleVars(uint32_t detid, const TrackerTopology* tTopo) {
  SiStripDetId subdet(detid);
  subdetid = subdet.subDetector();
  if( SiStripDetId::TIB == subdetid ) {
    
    module        = tTopo->tibModule(detid); 
    side          = tTopo->tibIsZMinusSide(detid)?-1:1;  
    layerwheel    = tTopo->tibLayer(detid); 
    stringringrod = tTopo->tibString(detid); 
    stereo        = tTopo->tibIsStereo(detid) ? 1 : 0;
  } else
  if( SiStripDetId::TID == subdetid ) {
    
    module        = tTopo->tidModule(detid); 
    side          = tTopo->tidIsZMinusSide(detid)?-1:1;  
    layerwheel    = tTopo->tidWheel(detid); 
    stringringrod = tTopo->tidRing(detid); 
    stereo        = tTopo->tidIsStereo(detid) ? 1 : 0;
  } else
  if( SiStripDetId::TOB == subdetid ) {
    
    module        = tTopo->tobModule(detid); 
    side          = tTopo->tobIsZMinusSide(detid)?-1:1;  
    layerwheel    = tTopo->tobLayer(detid); 
    stringringrod = tTopo->tobRod(detid); 
    stereo        = tTopo->tobIsStereo(detid) ? 1 : 0;
  } else
  if( SiStripDetId::TEC == subdetid ) {
    
    module        = tTopo->tecModule(detid); 
    side          = tTopo->tecIsZMinusSide(detid)?-1:1;  
    layerwheel    = tTopo->tecWheel(detid); 
    stringringrod = tTopo->tecRing(detid); 
    petal         = tTopo->tecPetalNumber(detid); 
    stereo        = tTopo->tecIsStereo(detid) ? 1 : 0;
  } else {
    module = 0;
    side = 0;
    layerwheel=-1;
    stringringrod = -1;
    petal=-1;
  }
}

bool ShallowTrackClustersProducerCombined::IsFarFromBorder( const TrajectoryStateOnSurface* trajState, const uint32_t detid, const edm::EventSetup* iSetup)
{ 
  edm::ESHandle<TrackerGeometry> tkGeom; iSetup->get<TrackerDigiGeometryRecord>().get( tkGeom );

  LocalPoint  HitLocalPos   = trajState->localPosition();
  LocalError  HitLocalError = trajState->localError().positionError() ;

  const GeomDetUnit* it = tkGeom->idToDetUnit(DetId(detid));
  if (dynamic_cast<const StripGeomDetUnit*>(it)==nullptr ) {
     std::cout << "this detID doesn't seem to belong to the Tracker" << std::endl;
     return false;
  }

  const BoundPlane plane = it->surface();
  const TrapezoidalPlaneBounds* trapezoidalBounds( dynamic_cast<const TrapezoidalPlaneBounds*>(&(plane.bounds())));
  const RectangularPlaneBounds* rectangularBounds( dynamic_cast<const RectangularPlaneBounds*>(&(plane.bounds())));

  double DistFromBorder = 1.0;    
  double HalfLength     = it->surface().bounds().length() /2.0;

  if(trapezoidalBounds)
  {
      std::array<const float, 4> const & parameters = (*trapezoidalBounds).parameters();
     HalfLength     = parameters[3];
  }else if(rectangularBounds){
     HalfLength     = it->surface().bounds().length() /2.0;
  }else{return false;}

  if (fabs(HitLocalPos.y())+HitLocalError.yy() >= (HalfLength - DistFromBorder) ) return false;

  return true;
}

double ShallowTrackClustersProducerCombined::thickness(DetId id,  const edm::EventSetup* iSetup)
{
   edm::ESHandle<TrackerGeometry> tkGeom; iSetup->get<TrackerDigiGeometryRecord>().get( tkGeom );
   
 map<DetId,double>::iterator th=m_thicknessMap.find(id);
 if(th!=m_thicknessMap.end())
   return (*th).second;
 else {
   double detThickness=1.;
//compute thickness normalization
   const GeomDetUnit* it = tkGeom->idToDetUnit(DetId(id));
   //bool isPixel = dynamic_cast<const PixelGeomDetUnit*>(it)!=nullptr;
  // bool isStrip = dynamic_cast<const StripGeomDetUnit*>(it)!=nullptr;
/*if (!isPixel && ! isStrip) {
     edm::LogWarning("DeDxHitsProducer") << "\t\t this detID doesn't seem to belong to the Tracker";
      detThickness = 1.;*/
      detThickness = it->surface().bounds().thickness();

   m_thicknessMap[id]=detThickness;//computed value
   return detThickness;
}

}
