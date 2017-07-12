from CalibTracker.SiStripCommon.shallowTree_test_template import *
process.TFileService.fileName = 'test_shallowTrackClustersProducer.root'

#process.source.fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/CosmicFall16PhaseIDR/TKCosmics_38T/GEN-SIM-RECO/DECO_90X_upgrade2017cosmics_realistic_deco_v18-v1/00000/0A229457-9122-E711-8E68-0CC47A78A4A6.root') #3.8 cosmics MC 28300ev
#process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
)


from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter

process.load('CalibTracker.SiStripCommon.ShallowTrackClustersProducer_cfi')
process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
process.tracksRefit = TrackRefitter.clone()
process.shallowTrackClusters.Tracks = 'tracksRefit'

process.testTree = cms.EDAnalyzer(
   "ShallowTree",
   outputCommands = cms.untracked.vstring(
      'drop *',
      'keep *_shallowTrackClusters_*_*',
      )
   )
process.p = cms.Path(process.MeasurementTrackerEvent*process.tracksRefit*process.shallowTrackClusters*process.testTree)
