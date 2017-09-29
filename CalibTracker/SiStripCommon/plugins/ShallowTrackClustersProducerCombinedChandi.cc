#include "CalibTracker/SiStripCommon/interface/ShallowTrackClustersProducerCombinedChandi.h"

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


#include "CalibTracker/Records/interface/SiStripDependentRecords.h"
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

uint32_t ev3=0;
uint32_t ntracks3 =0;


ShallowTrackClustersProducerCombinedChandi::ShallowTrackClustersProducerCombinedChandi(const edm::ParameterSet& iConfig)
  :  tracks_token_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("Tracks"))),
     association_token_(consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("Tracks"))),
     clusters_token_( consumes< edmNew::DetSetVector<SiStripCluster> >( iConfig.getParameter<edm::InputTag>("Clusters") ) ),
     theVertexToken_(consumes<std::vector<reco::Vertex> >          (iConfig.getParameter<edm::InputTag>("vertices"))),
     theDigisToken_    (consumes<edm::DetSetVector<SiStripProcessedRawDigi> > (edm::InputTag("siStripProcessedRawDigis", ""))),
     Suffix       ( iConfig.getParameter<std::string>("Suffix")    ),
     Prefix       ( iConfig.getParameter<std::string>("Prefix") ),
     //lorentzAngleName(iConfig.getParameter<std::string>("LorentzAngle")),
     isData       ( iConfig.getParameter<bool>("isData") ), //1. @MJ@ I have exchanged isdata with edmstripdigisimlink_Token_ (I had problem with compilation)
     edmstripdigisimlink_Token_( consumes<edm::DetSetVector<StripDigiSimLink>>( iConfig.getParameter<edm::InputTag>("StripDigiSimLinkTag")))
     //lowBound       ( iConfig.getParameter<int32_t>("lowBound") ),
     //highBound       ( iConfig.getParameter<int32_t>("highBound") ),
     //filename       ( iConfig.getParameter<std::string>("filename") )
{
  produces<std::vector<int> > ( Suffix + "clusterIdx"    ); //link: on trk cluster --> general cluster info 
  produces<std::vector<int> > ( Suffix + "onTrkClusterIdx" ); //link: general cluster info --> on track cluster
  produces <std::vector<int> > ( Suffix + "onTrkClustersBegin"  ); //link: track --> onTrkInfo (range)
  produces <std::vector<int> > ( Suffix + "onTrkClustersEnd" ); //link: track --> onTrkInfo (range)
  produces <std::vector<int> > ( Suffix + "trackindex" ); //link: on trk cluster --> track index

  produces <std::vector<unsigned int> > ( Suffix + "trackmulti"   );
  produces <std::vector<float> >        ( Suffix + "localtheta"   );  
  produces <std::vector<float> >        ( Suffix + "localeta"   );  
  produces <std::vector<float> >        ( Suffix + "localphi"     );  
  produces <std::vector<float> >        ( Suffix + "localpitch");

  produces <std::vector<float> >        ( Suffix + "localx"       );  
  produces <std::vector<float> >        ( Suffix + "localy"       );  
  produces <std::vector<float> >        ( Suffix + "localz"       );  
  produces <std::vector<float> >        ( Suffix + "strip"        );  
  produces <std::vector<float> >        ( Suffix + "globaltheta"  );  
  produces <std::vector<float> >        ( Suffix + "globaleta"  );  
  produces <std::vector<float> >        ( Suffix + "globalphi"    );
  produces <std::vector<float> >        ( Suffix + "globalx"      );
  produces <std::vector<float> >        ( Suffix + "globaly"      );
  produces <std::vector<float> >        ( Suffix + "globalz"      );
  produces <std::vector<float> >        ( Suffix + "transverseCurvature"      );
  produces <std::vector<float> >        ( Suffix + "trackPt"      );
  produces <std::vector<float> >        ( Suffix + "trackEta"      );
  produces <std::vector<float> >        ( Suffix + "insidistance" );
  produces <std::vector<float> >        ( Suffix + "covered"      );
  produces <std::vector<float> >        ( Suffix + "projwidth"    );
  produces <std::vector<float> >        ( Suffix + "BdotX"        );
  produces <std::vector<float> >        ( Suffix + "BdotY"        );
  produces <std::vector<float> >        ( Suffix + "BdotZ"        );
  produces <std::vector<float> >        ( Suffix + "BdotMag"        );

  produces <std::vector<float> >        ( Suffix + "rhlocalx"      );   
  produces <std::vector<float> >        ( Suffix + "rhlocaly"      );   
  produces <std::vector<float> >        ( Suffix + "rhlocalxerr"   );   
  produces <std::vector<float> >        ( Suffix + "rhlocalyerr"   );   
  produces <std::vector<float> >        ( Suffix + "rhglobalx"     );   
  produces <std::vector<float> >        ( Suffix + "rhglobaly"     );   
  produces <std::vector<float> >        ( Suffix + "rhglobalz"     );
  produces <std::vector<float> >        ( Suffix + "rhstrip"       );   
  produces <std::vector<float> >        ( Suffix + "rhmerr"        );   
  produces <std::vector<float> >        ( Suffix + "ubstrip"       );   
  produces <std::vector<float> >        ( Suffix + "ubmerr"        );   

  produces <std::vector<float> >       ( Suffix + "driftx"         );
  produces <std::vector<float> >       ( Suffix + "drifty"         );
  produces <std::vector<float> >       ( Suffix + "driftz"         );
  produces <std::vector<float> >       ( Suffix + "globalZofunitlocalY"  );            


  produces <std::vector<float> >    ( Prefix + "lorentzAngle"       );
  produces <std::vector<unsigned> >    ( Prefix + "number"       );
  produces <std::vector<unsigned> >    ( Prefix + "width"        );
  produces <std::vector<float> >       ( Prefix + "variance"     );
  produces <std::vector<float> >       ( Prefix + "barystrip"    );
  produces <std::vector<float> >       ( Prefix + "middlestrip"  );
  produces <std::vector<unsigned> >    ( Prefix + "charge"       );
  produces <std::vector<float> >       ( Prefix + "noise"        );
  produces <std::vector<float> >       ( Prefix + "ston"         );
  produces <std::vector<unsigned> >    ( Prefix + "seedstrip"    );
  produces <std::vector<unsigned> >    ( Prefix + "seedindex"    );
  produces <std::vector<unsigned> >    ( Prefix + "seedcharge"   );
  produces <std::vector<float> >       ( Prefix + "seednoise"    );
  produces <std::vector<float> >       ( Prefix + "seedgain"     );
  produces <std::vector<unsigned> >    ( Prefix + "qualityisbad" );

  produces <std::vector<float> >       ( Prefix + "rawchargeC"   );
  produces <std::vector<float> >       ( Prefix + "rawchargeL"   );
  produces <std::vector<float> >       ( Prefix + "rawchargeR"   );
  produces <std::vector<float> >       ( Prefix + "rawchargeLL"   );
  produces <std::vector<float> >       ( Prefix + "rawchargeRR"   );
  produces <std::vector<float> >       ( Prefix + "eta"          );
  produces <std::vector<float> >       ( Prefix + "foldedeta"    );
  produces <std::vector<float> >       ( Prefix + "etaX"         );
  produces <std::vector<float> >       ( Prefix + "etaasymm"     );
  produces <std::vector<float> >       ( Prefix + "outsideasymm");
  produces <std::vector<float> >       ( Prefix + "neweta");
  produces <std::vector<float> >       ( Prefix + "newetaerr");

  produces <std::vector<unsigned> >    ( Prefix + "detid"         );
  produces <std::vector<int> >         ( Prefix + "subdetid"      );
  produces <std::vector<int> >         ( Prefix + "module"        );
  produces <std::vector<int> >         ( Prefix + "side"          );
  produces <std::vector<int> >         ( Prefix + "layerwheel"    );
  produces <std::vector<int> >         ( Prefix + "stringringrod" );
  produces <std::vector<int> >         ( Prefix + "petal"         );
  produces <std::vector<int> >         ( Prefix + "stereo"        );
  produces <std::vector<float> >         ( Prefix + "stripLength"        );
  produces <std::vector<float> >         ( Prefix + "sensorThickness"        );
  produces <std::vector<float> >         ( Prefix + "stripCharge"        );
  produces <std::vector<float> >         ( Prefix + "stripChargeTotCharge"        );
  produces <std::vector<unsigned> >         ( Prefix + "stripChargeTotWidth"        );
  produces <std::vector<int> >         ( Prefix + "stripChargeStripNr"        );
  produces <std::vector<float> >         ( Prefix + "stripChargeLocalTrackPhi"        );
  produces <std::vector<float> >         ( Prefix + "stripChargeGlobalTrackPhi"        );
  produces <std::vector<float> >         ( Prefix + "stripChargeLocalTrackTheta"        );
  produces <std::vector<float> >         ( Prefix + "stripChargeGlobalTrackTheta"        );
  produces <std::vector<int> >         ( Prefix + "stripChargeSubdetid"        );
  produces <std::vector<int> >         ( Prefix + "stripChargeLayerwheel"        );
  produces <std::vector<unsigned> >         ( Prefix + "stripChargeDetid"        );
  produces <std::vector<float> >         ( Prefix + "stripChargeLocalX"        );
  produces <std::vector<float> >         ( Prefix + "stripChargeLocalY"        );
  produces <std::vector<float> >         ( Prefix + "stripChargetrackPt"        );
  produces <std::vector<float> >         ( Prefix + "stripChargesensorThickness"        );
  produces <std::vector<float> >         ( Prefix + "stripChargeBdotY"        );
  produces <std::vector<float> >        ( "PU"       );
  produces <std::vector<unsigned int> > ( "bx" );

  produces <std::vector<unsigned int> > ( "channel"       );
  produces <std::vector<float> > ( "fraction"       );
  produces <std::vector<int> > ( "bunchcrossing"       );


  //
}

void ShallowTrackClustersProducerCombinedChandi::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  shallow::CLUSTERMAP clustermap = shallow::make_cluster_map(iEvent, clusters_token_);
  edm::Handle<edm::View<reco::Track> > tracks;             iEvent.getByToken(tracks_token_, tracks);  



  int size = clustermap.size();
  //cout << "track size " << tracks->size() << "clustermap size " << size <<  endl;

  //links
  auto clusterIdx        = std::make_unique<std::vector<int>>(); //link: on trk cluster --> general cluster info 
  auto onTrkClusterIdx   = std::make_unique<std::vector<int>>(size,   -1); //link: general cluster info --> on track cluster
  auto onTrkClustersBegin = std::make_unique<std::vector<int>>( tracks->size(), -1 ); //link: track --> on trk cluster
  auto onTrkClustersEnd   = std::make_unique<std::vector<int>>( tracks->size(), -1 ); //link: track --> on trk cluster
  auto trackindex    = std::make_unique<std::vector<int>>(); //link: on track cluster --> track
  clusterIdx->reserve(size);
  trackindex->reserve(size);
  
  auto  trackmulti   = std::make_unique<std::vector<unsigned int>>(); trackmulti  ->reserve(size);
  auto  localtheta   = std::make_unique<std::vector<float>>(); localtheta  ->reserve(size);
  auto  localeta   = std::make_unique<std::vector<float>>(); localeta  ->reserve(size);
  auto  localphi     = std::make_unique<std::vector<float>>(); localphi    ->reserve(size);
  auto  localpitch   = std::make_unique<std::vector<float>>(); localpitch  ->reserve(size);
  auto  localx       = std::make_unique<std::vector<float>>(); localx      ->reserve(size);
  auto  localy       = std::make_unique<std::vector<float>>(); localy      ->reserve(size);
  auto  localz       = std::make_unique<std::vector<float>>(); localz      ->reserve(size);
  auto  strip        = std::make_unique<std::vector<float>>(); strip       ->reserve(size);
  auto  globaltheta  = std::make_unique<std::vector<float>>(); globaltheta ->reserve(size);
  auto  globaleta  = std::make_unique<std::vector<float>>(); globaleta ->reserve(size);
  auto  globalphi    = std::make_unique<std::vector<float>>(); globalphi   ->reserve(size);
  auto  globalx      = std::make_unique<std::vector<float>>(); globalx     ->reserve(size);
  auto  globaly      = std::make_unique<std::vector<float>>(); globaly     ->reserve(size);
  auto  globalz      = std::make_unique<std::vector<float>>(); globalz     ->reserve(size);
  auto  transverseCurvature      = std::make_unique<std::vector<float>>(); transverseCurvature     ->reserve(size);
  auto  trackPt            = std::make_unique<std::vector<float>>(); trackPt     ->reserve(size);
  auto  trackEta            = std::make_unique<std::vector<float>>(); trackEta     ->reserve(size);
  auto  insidistance = std::make_unique<std::vector<float>>(); insidistance->reserve(size);
  auto  projwidth    = std::make_unique<std::vector<float>>(); projwidth   ->reserve(size);
  auto  BdotX        = std::make_unique<std::vector<float>>(); BdotX       ->reserve(size);
  auto  BdotY        = std::make_unique<std::vector<float>>(); BdotY       ->reserve(size);
  auto  BdotZ        = std::make_unique<std::vector<float>>(); BdotZ       ->reserve(size);
  auto  BdotMag        = std::make_unique<std::vector<float>>(); BdotMag       ->reserve(size);
  auto  covered      = std::make_unique<std::vector<float>>(); covered     ->reserve(size);
  auto  rhlocalx     = std::make_unique<std::vector<float>>(); rhlocalx   ->reserve(size);
  auto  rhlocaly     = std::make_unique<std::vector<float>>(); rhlocaly   ->reserve(size);  
  auto  rhlocalxerr  = std::make_unique<std::vector<float>>(); rhlocalxerr->reserve(size);  
  auto  rhlocalyerr  = std::make_unique<std::vector<float>>(); rhlocalyerr->reserve(size);    
  auto  rhglobalx    = std::make_unique<std::vector<float>>(); rhglobalx  ->reserve(size);  
  auto  rhglobaly    = std::make_unique<std::vector<float>>(); rhglobaly  ->reserve(size);  
  auto  rhglobalz    = std::make_unique<std::vector<float>>(); rhglobalz  ->reserve(size);  
  auto  rhstrip      = std::make_unique<std::vector<float>>(); rhstrip    ->reserve(size);  
  auto  rhmerr       = std::make_unique<std::vector<float>>(); rhmerr     ->reserve(size);  
  auto  ubstrip      = std::make_unique<std::vector<float>>(); ubstrip    ->reserve(size);  
  auto  ubmerr       = std::make_unique<std::vector<float>>(); ubmerr     ->reserve(size);  
  auto  driftx       = std::make_unique<std::vector<float>>(); driftx     ->reserve(size);
  auto  drifty       = std::make_unique<std::vector<float>>(); drifty     ->reserve(size);
  auto  driftz       = std::make_unique<std::vector<float>>(); driftz     ->reserve(size);
  auto  globalZofunitlocalY = std::make_unique<std::vector<float>>(); globalZofunitlocalY->reserve(size);
  auto       lorentzAngle         = std::make_unique<std::vector<float>>();
  auto       number        = std::make_unique<std::vector<unsigned>>(7,0);
  auto       width         = std::make_unique<std::vector<unsigned>>();
  auto       variance      = std::make_unique<std::vector<float>>();
  auto       barystrip     = std::make_unique<std::vector<float>>();
  auto       middlestrip   = std::make_unique<std::vector<float>>();
  auto       charge        = std::make_unique<std::vector<unsigned>>();
  auto       noise         = std::make_unique<std::vector<float>>();
  auto       ston          = std::make_unique<std::vector<float>>();
  auto       seedstrip     = std::make_unique<std::vector<unsigned>>();
  auto       seedindex     = std::make_unique<std::vector<unsigned>>();
  auto       seedcharge    = std::make_unique<std::vector<unsigned>>();
  auto       seednoise     = std::make_unique<std::vector<float>>();
  auto       seedgain      = std::make_unique<std::vector<float>>();
  auto       qualityisbad  = std::make_unique<std::vector<unsigned>>();

  auto       rawchargeC    = std::make_unique<std::vector<float>>();
  auto       rawchargeL    = std::make_unique<std::vector<float>>();
  auto       rawchargeR    = std::make_unique<std::vector<float>>();
  auto       rawchargeLL   = std::make_unique<std::vector<float>>();
  auto       rawchargeRR   = std::make_unique<std::vector<float>>();
  auto       etaX          = std::make_unique<std::vector<float>>();
  auto       eta           = std::make_unique<std::vector<float>>();
  auto       foldedeta     = std::make_unique<std::vector<float>>();
  auto       etaasymm      = std::make_unique<std::vector<float>>();
  auto       outsideasymm  = std::make_unique<std::vector<float>>();
  auto       neweta        = std::make_unique<std::vector<float>>();
  auto       newetaerr     = std::make_unique<std::vector<float>>();

  auto       detid         = std::make_unique<std::vector<unsigned>>();
  auto       subdetid      = std::make_unique<std::vector<int>>();
  auto       side          = std::make_unique<std::vector<int>>();
  auto       module        = std::make_unique<std::vector<int>>();
  auto       layerwheel    = std::make_unique<std::vector<int>>();
  auto       stringringrod = std::make_unique<std::vector<int>>();
  auto       petal         = std::make_unique<std::vector<int>>();
  auto       stereo        = std::make_unique<std::vector<int>>();
  auto       stripLength   = std::make_unique<std::vector<float>>();
  auto       sensorThickness     = std::make_unique<std::vector<float>>();
  auto       stripCharge     = std::make_unique<std::vector<float>>();
  auto       stripChargeTotCharge     = std::make_unique<std::vector<float>>();
  auto       stripChargeTotWidth     = std::make_unique<std::vector<unsigned>>();
  auto       stripChargeStripNr     = std::make_unique<std::vector<int>>();
  auto       stripChargeLocalTrackPhi     = std::make_unique<std::vector<float>>();
  auto       stripChargeGlobalTrackPhi     = std::make_unique<std::vector<float>>();
  auto       stripChargeLocalTrackTheta     = std::make_unique<std::vector<float>>();
  auto       stripChargeGlobalTrackTheta     = std::make_unique<std::vector<float>>();
  auto       stripChargeSubdetid     = std::make_unique<std::vector<int>>();
  auto       stripChargeLayerwheel     = std::make_unique<std::vector<int>>();
  auto       stripChargeDetid     = std::make_unique<std::vector<unsigned>>();
  auto       stripChargeLocalX     = std::make_unique<std::vector<float>>();
  auto       stripChargeLocalY     = std::make_unique<std::vector<float>>();
  auto       stripChargetrackPt     = std::make_unique<std::vector<float>>();
  auto       stripChargesensorThickness    = std::make_unique<std::vector<float>>();
  auto       stripChargeBdotY   = std::make_unique<std::vector<float>>();


  auto       PU      = std::make_unique<std::vector<float>>();
  auto       bx      = std::make_unique<std::vector<unsigned int>>();

   auto       channel      = std::make_unique<std::vector<unsigned int>>();
  auto       fraction     = std::make_unique<std::vector<float>>();
  auto       bunchcrossing      = std::make_unique<std::vector<int>>();

  edm::ESHandle<TrackerGeometry> theTrackerGeometry;         iSetup.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );  
  edm::ESHandle<MagneticField> magfield;     iSetup.get<IdealMagneticFieldRecord>().get(magfield);
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


  //myfile2.open(filename);
  float PU_=0;
  PU_=vtx->size();
  ev3++;
  PU->push_back(PU_);
  bx->push_back(iEvent.bunchCrossing());

  for( TrajTrackAssociationCollection::const_iterator association = associations->begin(); 
       association != associations->end(); association++) {
    const Trajectory*  traj  = association->key.get();
    const reco::Track* track = association->val.get();
    int trk_idx = shallow::findTrackIndex(tracks, track); 
    size_t trk_strt_idx = ontrk_cluster_idx;

    bool isGoodTrack = trackFilter(track);
    if(!isGoodTrack)
      continue;
    ntracks3++;

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
	//std::cout << "geografical id " << id << std::endl;
	const SiStripClusterInfo info(*cluster_ptr, iSetup, id);
	const NearDigis digis = rawProcessedDigis.isValid() ? NearDigis(info, *rawProcessedDigis) : NearDigis(info);
                                

	//find if cluster was already assigned to a previous track
	auto already_visited = mapping.find(i);
	int nassociations = 1;
	if(already_visited != mapping.end()) {
	  nassociations += already_visited->second.size();
	  for(size_t idx : already_visited->second) {
	    trackmulti->at(idx)++;
	  }
	  already_visited->second.push_back(ontrk_cluster_idx);
	}
	else { //otherwise store this 
	  std::vector<size_t> single = {ontrk_cluster_idx};
	  mapping.insert( std::make_pair(i, single) );
	}

	const StripGeomDetUnit* theStripDet = dynamic_cast<const StripGeomDetUnit*>( theTrackerGeometry->idToDet( hit->geographicalId() ) );
	LocalVector drift = shallow::drift( theStripDet, *magfield, *SiStripLorentzAngle);
	if(nassociations == 1) onTrkClusterIdx->at(i) = ontrk_cluster_idx; //link: general cluster info --> on track cluster
	clusterIdx->push_back(  i );  //link: on trk cluster --> general cluster info
	//if(PU_>lowBound && PU_<highBound)
	//{
	(number->at(0))++;
	(number->at(moduleV.subdetid))++;
                     
	float langle = (SiStripLorentzAngle.isValid()) ? SiStripLorentzAngle->getLorentzAngle(id) : 0.;;
	lorentzAngle->push_back(langle);
	width->push_back(        cluster_ptr->amplitudes().size()                              );
	barystrip->push_back(    cluster_ptr->barycenter()                                     );
	variance->push_back(     info.variance()                                         );
	middlestrip->push_back(  info.firstStrip() + info.width()/2.0                    );
	charge->push_back(       info.charge()                                           );
	noise->push_back(        info.noiseRescaledByGain()                              );
	ston->push_back(         info.signalOverNoise()                                  );
	seedstrip->push_back(    info.maxStrip()                                         );
	seedindex->push_back(    info.maxIndex()                                         );
	seedcharge->push_back(   info.maxCharge()                                        );
	seednoise->push_back(    info.stripNoisesRescaledByGain().at(info.maxIndex())   );
	seedgain->push_back(     info.stripGains().at(info.maxIndex())                  );
	qualityisbad->push_back( info.IsAnythingBad()                                    );
	 
	rawchargeC->push_back(   digis.max            );
	rawchargeL->push_back(   digis.left           );
	rawchargeR->push_back(   digis.right          );
	rawchargeLL->push_back(  digis.Lleft          );
	rawchargeRR->push_back(  digis.Rright         );
	etaX->push_back(         digis.etaX()         );
	eta->push_back(          digis.eta()          );
	etaasymm->push_back(     digis.etaasymm()     );
	outsideasymm->push_back( digis.outsideasymm() );
	neweta->push_back(       (digis.last-digis.first)/info.charge() );
	newetaerr->push_back(    (sqrt(digis.last+digis.first))/pow(info.charge(),1.5) );



	detid->push_back(            id                 );
	subdetid->push_back(         moduleV.subdetid      );
	side->push_back(             moduleV.side          );
	module->push_back(           moduleV.module        );
	layerwheel->push_back(       moduleV.layerwheel    );
	stringringrod->push_back(    moduleV.stringringrod );
	petal->push_back(            moduleV.petal         );
	stereo->push_back(           moduleV.stereo        );


	stripLength->push_back(           theStripDet->specificTopology().stripLength() );
	sensorThickness->push_back(             theStripDet->specificSurface().bounds().thickness() );


        //3. @MJ@ @PutTheCodeHere@
        //here I suggest to put the digisimlink code. The simlink works like this (probably, did not study that so far): You have your cluster and the link gives you fraction of the charge contribution on the given strip(channel) and also tells you from which event this contribution came. Probably you will have more contributions for one strip and of course if cluster width>1 even more strips. This means you will have many numbers per one cluster and you need to deal with it. There are few options:
        //a) Either you do a vector of vector (as cluster is already stored in vector), it is not that much intuitive, but I use it sometimes
        //b) Or you can somehow double information as I did below for the amplitudes in cluster
        //c) you can think about something more clever


	uint32_t strips=0;
	for( auto itAmpl = cluster_ptr->amplitudes().begin(); itAmpl != cluster_ptr->amplitudes().end(); ++itAmpl){
          stripCharge->push_back(*itAmpl); 
          stripChargeTotCharge->push_back(info.charge());
          stripChargeTotWidth->push_back(cluster_ptr->amplitudes().size());
          stripChargeStripNr->push_back(info.firstStrip() + strips);
          stripChargeLocalTrackPhi->push_back(  (theStripDet->toLocal(tsos.globalDirection())).phi());
          stripChargeGlobalTrackPhi->push_back(  tsos.globalDirection().phi());
          stripChargeLocalTrackTheta->push_back(  (theStripDet->toLocal(tsos.globalDirection())).theta());
          stripChargeGlobalTrackTheta->push_back(  tsos.globalDirection().theta());
          stripChargeSubdetid->push_back(moduleV.subdetid);
          stripChargeLayerwheel->push_back(moduleV.layerwheel);
          stripChargeDetid->push_back(id); 
          stripChargeLocalX->push_back((theStripDet->toLocal(tsos.globalPosition())).x());
          stripChargeLocalY->push_back((theStripDet->toLocal(tsos.globalPosition())).y());
          stripChargetrackPt->push_back(   track->pt() );                         
	  stripChargesensorThickness->push_back(             theStripDet->specificSurface().bounds().thickness() );
	  stripChargeBdotY->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).y() );
          ++strips; 
	}


	trackmulti->push_back(  nassociations );
	trackindex->push_back(  trk_idx );
	localtheta->push_back(  (theStripDet->toLocal(tsos.globalDirection())).theta() ); 
	localeta->push_back(  (theStripDet->toLocal(tsos.globalDirection())).eta() ); 
	localphi->push_back(    (theStripDet->toLocal(tsos.globalDirection())).phi() );   
	localpitch->push_back(  (theStripDet->specificTopology()).localPitch(theStripDet->toLocal(tsos.globalPosition())) ); 
	localx->push_back(      (theStripDet->toLocal(tsos.globalPosition())).x() );    
	localy->push_back(      (theStripDet->toLocal(tsos.globalPosition())).y() );    
	localz->push_back(      (theStripDet->toLocal(tsos.globalPosition())).z() );    
	strip->push_back(       (theStripDet->specificTopology()).strip(theStripDet->toLocal(tsos.globalPosition())) );
	globaltheta->push_back( tsos.globalDirection().theta() );                       
	globaleta->push_back( tsos.globalDirection().eta() );                       
	globalphi->push_back(   tsos.globalDirection().phi() );                         
	globalx->push_back(     tsos.globalPosition().x() );                            
	globaly->push_back(     tsos.globalPosition().y() );                            
	globalz->push_back(     tsos.globalPosition().z() );                            
	transverseCurvature->push_back(     tsos.transverseCurvature() );      

	//cout << "track pt" << track->pt() << endl;                      
	trackPt->push_back(   track->pt() );                         
	trackEta->push_back(   track->eta() );                         
	insidistance->push_back(1./fabs(cos(localtheta->at(ontrk_cluster_idx))) );                      
	projwidth->push_back(   tan(localtheta->at(ontrk_cluster_idx))*cos(localphi->at(ontrk_cluster_idx)) );         
	BdotX->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).x() );
	BdotY->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).y() );
	BdotZ->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).z() );
	BdotMag->push_back(       (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).mag() );
	covered->push_back(     drift.z()/localpitch->at(ontrk_cluster_idx) * fabs(projwidth->at(ontrk_cluster_idx) - drift.x()/drift.z()) );
	rhlocalx->push_back(    hit->localPosition().x() );
	rhlocaly->push_back(    hit->localPosition().y() );
	rhlocalxerr->push_back( sqrt(hit->localPositionError().xx()) );
	rhlocalyerr->push_back( sqrt(hit->localPositionError().yy()) );
	rhglobalx->push_back(   theStripDet->toGlobal(hit->localPosition()).x() );
	rhglobaly->push_back(   theStripDet->toGlobal(hit->localPosition()).y() );
	rhglobalz->push_back(   theStripDet->toGlobal(hit->localPosition()).z() );
	rhstrip->push_back(     theStripDet->specificTopology().strip(hit->localPosition()) );
	rhmerr->push_back(      sqrt(theStripDet->specificTopology().measurementError(hit->localPosition(), hit->localPositionError()).uu()) );
	ubstrip->push_back(     theStripDet->specificTopology().strip(unbiased.localPosition()) );
	ubmerr->push_back(      sqrt(theStripDet->specificTopology().measurementError(unbiased.localPosition(), unbiased.localError().positionError()).uu()) );
	driftx->push_back(      drift.x() );
	drifty->push_back(      drift.y() );
	driftz->push_back(      drift.z() );
	globalZofunitlocalY->push_back( (theStripDet->toGlobal(LocalVector(0,1,0))).z() );

	
	ontrk_cluster_idx++;
      } //for(unsigned h=0; h<2; h++) { //loop over possible Hit options (1D, 2D)
    } //BOOST_FOREACH( const TrajectoryMeasurement measurement, traj->measurements() )

    onTrkClustersBegin->at(trk_idx) = trk_strt_idx;
    onTrkClustersEnd->at(trk_idx)   = ontrk_cluster_idx;

    
  //inserted digisimlink
    edm::Handle<edmNew::DetSetVector<SiStripCluster>> clustercoll; //2. @MJ@ why you have one new cluster collection? You already have one above and you need to work with that one because you want clusters ON TRACK. The code above works like this: You have tracks, the you get hits on the track and then you link hits with clusters. Then your work starts: you need to take the clusters above and link them with simulated digis. I suggest to put the code somewhere around @PutTheCodeHere@ (this mark is in code).
    iEvent.getByToken(clusters_token_,clustercoll);

    edm::Handle< edm::DetSetVector<StripDigiSimLink> > stripdigilinkContainer;
    iEvent.getByToken( edmstripdigisimlink_Token_, stripdigilinkContainer);

    std::vector<int> map_subdets_cluster_charge_IT;
    std::vector<int> map_subdets_cluster_width_IT;
    std::vector<int> map_subdets_cluster_charge_OOT;
    std::vector<int> map_subdets_cluster_width_OOT;

    for(edmNew::DetSetVector<SiStripCluster>::const_iterator det = clustercoll->begin();det != clustercoll->end(); ++det) {

      DetId detid(det->detId());
      unsigned int subdet = detid.subdetId();

      edm::DetSetVector<StripDigiSimLink>::const_iterator idigilinksearch = stripdigilinkContainer->find(detid); 
      if(idigilinksearch != stripdigilinkContainer->end()) {  //if it is not empty
	edm::DetSet<StripDigiSimLink> stripdigilink_detset = (*idigilinksearch);
    
	for(edmNew::DetSet<SiStripCluster>::const_iterator clus=det->begin();clus!=det->end();++clus) {
	  int clus_size = clus->amplitudes().size();
	  int clus_firststrip  = clus->firstStrip();     
	  int clus_laststrip   = clus_firststrip + clus_size;
		  
	  for(edm::DetSet<StripDigiSimLink>::const_iterator simlinkiter = stripdigilink_detset.data.begin(), 
		simlinkerEnd = stripdigilink_detset.data.end(); simlinkiter != simlinkerEnd; ++simlinkiter){      
	    int channel_ = (int)(simlinkiter->channel()); 

	    //	    channel->push_back(channel_);
	    fraction->push_back(simlinkiter->fraction());
	    unsigned int clus_bunchcrossing;
	    if( channel_ >= clus_firststrip  && channel_ < clus_laststrip ){
	      clus_bunchcrossing=simlinkiter->eventId().bunchCrossing();

	    }else{
	      clus_bunchcrossing = -1;
	    }

	      bunchcrossing->push_back(clus_bunchcrossing);	 
	    if(clus_bunchcrossing==0){
	      map_subdets_cluster_charge_IT.push_back(clus->charge());
	      map_subdets_cluster_width_IT.push_back(clus->amplitudes().size());
	    }else{
	      map_subdets_cluster_charge_OOT.push_back(clus->charge()); 
	      map_subdets_cluster_width_OOT.push_back(clus->amplitudes().size()); 
	    }
	  }//Loop over clusters in each sub-setector
	}//Loop over edmstripsimdigilink detsetvector
      }//Loop over cluster collection over all SiStrip Clusters
    } //end digisimlink
   

  } //for(TrajTrackAssociationCollection::const_iterator association = associations->begin();

  iEvent.put(std::move(clusterIdx        ), Suffix + "clusterIdx"  );
  iEvent.put(std::move(onTrkClusterIdx   ), Suffix + "onTrkClusterIdx"  );
  iEvent.put(std::move(onTrkClustersBegin), Suffix + "onTrkClustersBegin"  );
  iEvent.put(std::move(onTrkClustersEnd  ), Suffix + "onTrkClustersEnd"  );
  iEvent.put(std::move(trackindex),  Suffix + "trackindex"   );

  iEvent.put(std::move(trackmulti),  Suffix + "trackmulti"   );
  iEvent.put(std::move(localtheta),  Suffix + "localtheta"   );
  iEvent.put(std::move(localeta),  Suffix + "localeta"   );
  iEvent.put(std::move(localphi),    Suffix + "localphi"     );
  iEvent.put(std::move(localpitch),  Suffix + "localpitch"   );
  iEvent.put(std::move(localx),      Suffix + "localx"       );
  iEvent.put(std::move(localy),      Suffix + "localy"       );
  iEvent.put(std::move(localz),      Suffix + "localz"       );
  iEvent.put(std::move(strip),       Suffix + "strip"        );
  iEvent.put(std::move(globaltheta), Suffix + "globaltheta"  );
  iEvent.put(std::move(globaleta), Suffix + "globaleta"  );
  iEvent.put(std::move(globalphi),   Suffix + "globalphi"    );
  iEvent.put(std::move(globalx),     Suffix + "globalx"      );
  iEvent.put(std::move(globaly),     Suffix + "globaly"      );
  iEvent.put(std::move(globalz),     Suffix + "globalz"      );  
  iEvent.put(std::move(transverseCurvature),     Suffix + "transverseCurvature"      );  
  iEvent.put(std::move(trackPt),     Suffix + "trackPt"      );  
  iEvent.put(std::move(trackEta),     Suffix + "trackEta"      );  
  iEvent.put(std::move(insidistance),Suffix + "insidistance" );
  iEvent.put(std::move(covered),     Suffix + "covered"      );
  iEvent.put(std::move(projwidth),   Suffix + "projwidth"    );
  iEvent.put(std::move(BdotX),       Suffix + "BdotX"        );
  iEvent.put(std::move(BdotY),       Suffix + "BdotY"        );
  iEvent.put(std::move(BdotZ),       Suffix + "BdotZ"        );
  iEvent.put(std::move(BdotMag),       Suffix + "BdotMag"        );
  iEvent.put(std::move(rhlocalx),    Suffix + "rhlocalx"     );   
  iEvent.put(std::move(rhlocaly),    Suffix + "rhlocaly"     );   
  iEvent.put(std::move(rhlocalxerr), Suffix + "rhlocalxerr"  );   
  iEvent.put(std::move(rhlocalyerr), Suffix + "rhlocalyerr"  );   
  iEvent.put(std::move(rhglobalx),   Suffix + "rhglobalx"    );   
  iEvent.put(std::move(rhglobaly),   Suffix + "rhglobaly"    );   
  iEvent.put(std::move(rhglobalz),   Suffix + "rhglobalz"    );   
  iEvent.put(std::move(rhstrip),     Suffix + "rhstrip"      );   
  iEvent.put(std::move(rhmerr),      Suffix + "rhmerr"       );   
  iEvent.put(std::move(ubstrip),     Suffix + "ubstrip"      );   
  iEvent.put(std::move(ubmerr),      Suffix + "ubmerr"       );   
  iEvent.put(std::move(driftx),      Suffix + "driftx"       );
  iEvent.put(std::move(drifty),      Suffix + "drifty"       );
  iEvent.put(std::move(driftz),      Suffix + "driftz"       );
  iEvent.put(std::move(globalZofunitlocalY), Suffix + "globalZofunitlocalY"  );


  iEvent.put(std::move(lorentzAngle),      Prefix + "lorentzAngle"       );
  iEvent.put(std::move(number),      Prefix + "number"       );
  iEvent.put(std::move(width),       Prefix + "width"        );
  iEvent.put(std::move(variance),    Prefix + "variance"     );
  iEvent.put(std::move(barystrip),   Prefix + "barystrip"    );
  iEvent.put(std::move(middlestrip), Prefix + "middlestrip"  );
  iEvent.put(std::move(charge),      Prefix + "charge"       );
  iEvent.put(std::move(noise),       Prefix + "noise"        );
  iEvent.put(std::move(ston),        Prefix + "ston"         );
  iEvent.put(std::move(seedstrip),   Prefix + "seedstrip"    );
  iEvent.put(std::move(seedindex),   Prefix + "seedindex"    );
  iEvent.put(std::move(seedcharge),  Prefix + "seedcharge"   );
  iEvent.put(std::move(seednoise),   Prefix + "seednoise"    );
  iEvent.put(std::move(seedgain),    Prefix + "seedgain"     );
  iEvent.put(std::move(qualityisbad),Prefix + "qualityisbad" );

  iEvent.put(std::move(rawchargeC),  Prefix + "rawchargeC"   );
  iEvent.put(std::move(rawchargeL),  Prefix + "rawchargeL"   );
  iEvent.put(std::move(rawchargeR),  Prefix + "rawchargeR"   );
  iEvent.put(std::move(rawchargeLL), Prefix + "rawchargeLL"  );
  iEvent.put(std::move(rawchargeRR), Prefix + "rawchargeRR"  );
  iEvent.put(std::move(etaX),        Prefix + "etaX"         );
  iEvent.put(std::move(eta),         Prefix + "eta"          );
  iEvent.put(std::move(foldedeta),   Prefix + "foldedeta"    );
  iEvent.put(std::move(etaasymm),    Prefix + "etaasymm"     );
  iEvent.put(std::move(outsideasymm),Prefix + "outsideasymm" );
  iEvent.put(std::move(neweta),      Prefix + "neweta"       );
  iEvent.put(std::move(newetaerr),   Prefix + "newetaerr"    );

  iEvent.put(std::move(detid),        Prefix + "detid"         );
  iEvent.put(std::move(subdetid),     Prefix + "subdetid"      );
  iEvent.put(std::move(module),       Prefix + "module"        );
  iEvent.put(std::move(side),         Prefix + "side"          );
  iEvent.put(std::move(layerwheel),   Prefix + "layerwheel"    );
  iEvent.put(std::move(stringringrod),Prefix + "stringringrod" );
  iEvent.put(std::move(petal),        Prefix + "petal"         );
  iEvent.put(std::move(stereo),       Prefix + "stereo"        );
  iEvent.put(std::move(stripLength),       Prefix + "stripLength"        );
  iEvent.put(std::move(sensorThickness),       Prefix + "sensorThickness"        );
  iEvent.put(std::move(stripCharge),       Prefix + "stripCharge"        );
  iEvent.put(std::move(stripChargeTotCharge),       Prefix + "stripChargeTotCharge"        );
  iEvent.put(std::move(stripChargeTotWidth),       Prefix + "stripChargeTotWidth"        );
  iEvent.put(std::move(stripChargeStripNr),       Prefix + "stripChargeStripNr"        );
  iEvent.put(std::move(stripChargeLocalTrackPhi),       Prefix + "stripChargeLocalTrackPhi"        );
  iEvent.put(std::move(stripChargeGlobalTrackPhi),       Prefix + "stripChargeGlobalTrackPhi"        );
  iEvent.put(std::move(stripChargeLocalTrackTheta),       Prefix + "stripChargeLocalTrackTheta"        );
  iEvent.put(std::move(stripChargeGlobalTrackTheta),       Prefix + "stripChargeGlobalTrackTheta"        );
  iEvent.put(std::move(stripChargeSubdetid),     Prefix + "stripChargeSubdetid"      );
  iEvent.put(std::move(stripChargeLayerwheel),   Prefix + "stripChargeLayerwheel"    );
  iEvent.put(std::move(stripChargeDetid),   Prefix + "stripChargeDetid"    );
  iEvent.put(std::move(stripChargeLocalX),   Prefix + "stripChargeLocalX"    );
  iEvent.put(std::move(stripChargeLocalY),   Prefix + "stripChargeLocalY"    );
  iEvent.put(std::move(stripChargetrackPt),   Prefix + "stripChargetrackPt"    );
  iEvent.put(std::move(stripChargesensorThickness),   Prefix + "stripChargesensorThickness"    );
  iEvent.put(std::move(stripChargeBdotY),   Prefix + "stripChargeBdotY"    );
  
  iEvent.put(std::move(PU),       "PU"        );
  iEvent.put(std::move(bx),       "bx"        );
  iEvent.put(std::move(channel),       "channel"        );
  iEvent.put(std::move(fraction),       "fraction"        );
  iEvent.put(std::move(bunchcrossing),       "bunchcrossing"        );

  
}

	     bool ShallowTrackClustersProducerCombinedChandi::trackFilter(const reco::Track* trk)
	     {
	       if (trk->pt() < 0.8) return false;
	       if (trk->p()  < 2.0) return false;
	       if (trk->hitPattern().numberOfValidTrackerHits()  <= 6) return false;
	       if (trk->normalizedChi2() > 10.0) return false;
	       return true;
	     }


	     ShallowTrackClustersProducerCombinedChandi::NearDigis::
	     NearDigis(const SiStripClusterInfo& info) {
	       max =  info.maxCharge();
	       left =           info.maxIndex()    > uint16_t(0)                ? info.stripCharges()[info.maxIndex()-1]      : 0 ;
	       Lleft =          info.maxIndex()    > uint16_t(1)                ? info.stripCharges()[info.maxIndex()-2]      : 0 ;
	       right=  unsigned(info.maxIndex()+1) < info.stripCharges().size() ? info.stripCharges()[info.maxIndex()+1]      : 0 ;
	       Rright= unsigned(info.maxIndex()+2) < info.stripCharges().size() ? info.stripCharges()[info.maxIndex()+2]      : 0 ;
	       first = info.stripCharges()[0];
	       last =  info.stripCharges()[info.width()-1];
	     }
	     ShallowTrackClustersProducerCombinedChandi::NearDigis::
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

	     ShallowTrackClustersProducerCombinedChandi::moduleVars::
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
	     
