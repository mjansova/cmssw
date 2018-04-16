from CalibTracker.SiStripCommon.shallowTree_test_template import *
process.TFileService.fileName = 'EricsTree.root'
import sys

#process.source.fileNames = cms.untracked.vstring('file:output/step3_RAW2DIGI_L1Reco_RECO_19.root') #test_shallowTrackAndClusterFullInfo2016run284078.root.root
process.source.fileNames = cms.untracked.vstring(
'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/mjansova/simu/CRUZET_VR2018_RECO/ParkingCosmicsVirginRaw1/Commissioning2018_v1_RAW/180403_074415/0000/step3HICMN_RAW2DIGI_L1Reco_RECO_70.root')
#process.source.fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/RelValMinBias_13/GEN-SIM-RECO/90X_mcRun2_asymptotic_v5-v1/00000/9A94CC7A-4B0F-E711-8AC9-0CC47A4D7616.root') #relval MC 90, 16500 ev #test_shallowTrackClusterNoPUMCrun2_mcTAG.root




process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

#from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter
process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
process.shallowTrackClustersCombined = cms.EDProducer("ShallowTrackClustersProducerCombined",
                                      Tracks=cms.InputTag("ctfWithMaterialTracksP5",""),
                                      Clusters=cms.InputTag("siStripClusters"),
                                      vertices=cms.InputTag("offlinePrimaryVertices"),
                                      LorentzAngle = cms.string(''),
                                      Prefix=cms.string("cluster"),
                                      Suffix=cms.string("tsos"),
                                      lowBound=cms.int32(0),
                                      highBound=cms.int32(1000),
                                      filename=cms.string("lowPUlogMC.txt"),
                                      isData=cms.bool(False))

#process.GlobalTag = GlobalTag(process.GlobalTag, '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2017_realistic_v20')
#process.GlobalTag = GlobalTag(process.GlobalTag, '90X_dataRun2_Prompt_v2')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data_promptlike')
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_dataRun2_Prompt_v3')


#"global_tag": "90X_dataRun2_Prompt_v3" for run 293492 (from 902)
#"global_tag": "90X_dataRun2_Prompt_v2" for data (before 902)
#MC global tag: 90X_upgrade2017_realistic_v20
#process.load('CalibTracker.SiStripCommon.ShallowTrackClustersProducer_cfi')
#auto:run2_data_promptlike run 2 data global tag

process.TrackRefitter.src = "ctfWithMaterialTracksP5"
#process.TrackRefitter.src = "generalTracks"
process.TrackRefitter.TTRHBuilder = "WithTrackAngle"
process.TrackRefitter.NavigationSchool = ""
process.shallowTrackClustersCombined.Tracks = 'TrackRefitter'

#process.load('CalibTracker.SiStripCommon.ShallowTrackClustersProducer_cfi')
process.testTree = cms.EDAnalyzer(
   "ShallowTree",
   outputCommands = cms.untracked.vstring(
      'drop *',
      'keep *_shallowTrackClustersCombined_*_*',
      )
   )


process.p = cms.Path(process.TrackRefitter*process.shallowTrackClustersCombined*process.testTree)
