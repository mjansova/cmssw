#ifndef SHALLOW_TRACKCLUSTERS_PRODUCER_COMBINED
#define SHALLOW_TRACKCLUSTERS_PRODUCER_COMBINED

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"

#include "CalibTracker/SiStripCommon/interface/ShallowTools.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

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


#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"



#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "DataFormats/TrackReco/interface/TrackDeDxHits.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
//#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"

#include <ext/hash_map>

class SiStripClusterInfo;
class SiStripProcessedRawDigi;
class TrackerTopology;
class SiStripLorentzAngle;
//class StripClusterParameterEstimator;

class ShallowTrackClustersProducerCombined : public edm::EDProducer {
public:
  explicit ShallowTrackClustersProducerCombined(const edm::ParameterSet&);
private:
  const edm::EDGetTokenT<edm::View<reco::Track> > tracks_token_;
  const edm::EDGetTokenT<TrajTrackAssociationCollection> association_token_;
  const edm::EDGetTokenT< edmNew::DetSetVector<SiStripCluster> > clusters_token_;
  edm::EDGetTokenT<std::vector<reco::Vertex> >          theVertexToken_;
  edm::EDGetTokenT<edm::DetSetVector<SiStripProcessedRawDigi> > theDigisToken_;
  const edm::EDGetTokenT<DetIdCollection> zsdigis_token_;
  edm::EDGetTokenT<edm::TriggerResults> theTriggerToken_;
  std::string Suffix;
  std::string Prefix;
  int32_t lowBound;
  int32_t highBound;
  std::string filename;
  edm::ESHandle<SiStripLorentzAngle> lorentzAngleHandle;
  const std::string lorentzAngleName;
  bool isData;

  void produce( edm::Event &, const edm::EventSetup & );
  bool trackFilter(const reco::Track* trk);
  double thickness(DetId id , const edm::EventSetup* iSetup);
  bool IsFarFromBorder(const TrajectoryStateOnSurface* trajState, const uint32_t detid, const edm::EventSetup* iSetup);

  std::map<DetId,double> m_thicknessMap;
  struct moduleVars {
    moduleVars(uint32_t, const TrackerTopology*);
    int subdetid, side, layerwheel, stringringrod, petal, stereo;
    uint32_t module;
  };

  struct NearDigis { 
    NearDigis(const SiStripClusterInfo&);
    NearDigis(const SiStripClusterInfo&, const edm::DetSetVector<SiStripProcessedRawDigi>&);
    float max, left, right, first, last, Lleft, Rright; 
    float etaX() const {return ((left+right)/max)/2.;}
    float eta()  const {return right>left ? max/(max+right) : left/(left+max);}
    float etaasymm() const {return right>left ? (right-max)/(right+max) : (max-left)/(max+left);}
    float outsideasymm() const {return (last-first)/(last+first);}
  };
};
#endif
