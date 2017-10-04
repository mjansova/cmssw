from CalibTracker.SiStripCommon.shallowTree_test_template import *
process.TFileService.fileName = 'bla.root'
import sys

process.source.fileNames = cms.untracked.vstring('file:output/step3_RAW2DIGI_L1Reco_RECO_19.root') #test_shallowTrackAndClusterFullInfo2016run284078.root.root
#process.source.fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/RelValMinBias_13/GEN-SIM-RECO/90X_mcRun2_asymptotic_v5-v1/00000/9A94CC7A-4B0F-E711-8AC9-0CC47A4D7616.root') #relval MC 90, 16500 ev #test_shallowTrackClusterNoPUMCrun2_mcTAG.root

inputStr = sys.argv[2]
update = sys.argv[3]
found = inputStr.find("root")

if found > 0:
    process.source.fileNames = cms.untracked.vstring("file:"+inputStr) 
    print(input)
    process.TFileService.fileName = 'test_shallowTrackAndClusterSeedInfo.root' + inputStr
else:
    process.TFileService.fileName = 'test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root'



if update == "update":
    print("updating database conditions")
    process.GlobalTag.toGet = cms.VPSet(
      cms.PSet(record = cms.string("SiStripApvGainSimRcd"),
               tag = cms.string("myTag"),
               connect = cms.string("sqlite_file:/opt/sbg/data/data6/cms/mjansova/CMSSW_9_2_0/src/CalibTracker/SiStripCommon/test/dbConditions/G1_run284078.db")
               ),
       cms.PSet(record = cms.string("SiStripApvGainRcd"),
               tag = cms.string("myTag"),
               connect = cms.string("sqlite_file:/opt/sbg/data/data6/cms/mjansova/CMSSW_9_2_0/src/CalibTracker/SiStripCommon/test/dbConditions/G1_run284078.db")
               ),
       cms.PSet(record = cms.string("SiStripApvGain2Rcd"),
               tag = cms.string("myTag"),
               connect = cms.string("sqlite_file:/opt/sbg/data/data6/cms/mjansova/CMSSW_9_2_0/src/CalibTracker/SiStripCommon/test/dbConditions/G2_run284078.db")
               ),
       cms.PSet(record = cms.string("SiStripNoisesRcd"),
               tag = cms.string("myTag"),
               connect = cms.string("sqlite_file:/opt/sbg/data/data6/cms/mjansova/CMSSW_9_2_0/src/CalibTracker/SiStripCommon/test/dbConditions/Noise_run284078.db")
               ),
       cms.PSet(record = cms.string("SiStripBadChannelRcd"),
               tag = cms.string("myTag"),
               connect = cms.string("sqlite_file:/opt/sbg/data/data6/cms/mjansova/CMSSW_9_2_0/src/CalibTracker/SiStripCommon/test/dbConditions/BadChannel_run284078.db")
               ),
       cms.PSet(record = cms.string("SiStripBadFiberRcd"),
               tag = cms.string("myTag"),
               connect = cms.string("sqlite_file:/opt/sbg/data/data6/cms/mjansova/CMSSW_9_2_0/src/CalibTracker/SiStripCommon/test/dbConditions/BadFiber_run284078.db")
               )
    )



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter

process.shallowTrackClustersCombined = cms.EDProducer("ShallowTrackClustersProducerCombined",
                                      Tracks=cms.InputTag("generalTracks",""),
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
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

#"global_tag": "90X_dataRun2_Prompt_v3" for run 293492 (from 902)
#"global_tag": "90X_dataRun2_Prompt_v2" for data (before 902)
#MC global tag: 90X_upgrade2017_realistic_v20
#process.load('CalibTracker.SiStripCommon.ShallowTrackClustersProducer_cfi')
#auto:run2_data_promptlike run 2 data global tag

process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
process.tracksRefit = TrackRefitter.clone()
process.shallowTrackClustersCombined.Tracks = 'tracksRefit'


#process.load('CalibTracker.SiStripCommon.ShallowTrackClustersProducer_cfi')
process.testTree = cms.EDAnalyzer(
   "ShallowTree",
   outputCommands = cms.untracked.vstring(
      'drop *',
      'keep *_shallowTrackClustersCombined_*_*',
      )
   )


process.p = cms.Path(process.MeasurementTrackerEvent*process.tracksRefit*process.shallowTrackClustersCombined*process.testTree)
