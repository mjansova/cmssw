import FWCore.ParameterSet.Config as cms

shallowTrackClusters = cms.EDProducer("ShallowTrackClustersProducerCombined",
                                      Tracks=cms.InputTag("generalTracks",""),
                                      Clusters=cms.InputTag("siStripClusters"),
                                      vertices=cms.InputTag("offlinePrimaryVertices"),
                                      Prefix=cms.string("cluster"),
                                      Suffix=cms.string("Tracking"))
