from CalibTracker.SiStripCommon.shallowTree_test_template import *
process.TFileService.fileName = 'bla.root'
import sys

#process.source.fileNames = cms.untracked.vstring('file:step3CRUZET.root') #0T gen
#process.source.fileNames = cms.untracked.vstring('file:CTcruzet_0.805step3.root') #0T gen

process.source.fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2017B/Cosmics/RECO/PromptReco-v2/000/298/647/00000/14EC1D25-5E66-E711-8C7C-02163E01A3A5.root') #2017 runB craft data
#process.source.fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2016H/Cosmics/RECO/PromptReco-v2/000/281/214/00000/2CBD8474-6782-E611-9360-02163E0145CC.root') #2016 runH craft data

inputStr = sys.argv[2]
found = inputStr.find("root")

if found > 0:
    process.source.fileNames = cms.untracked.vstring("file:"+inputStr)
    process.TFileService.fileName = 'test_shallowTrackAndClusterFullInfo.root' + inputStr
    print(input)
else:
    process.TFileService.fileName = 'test_shallowTrackAndClusterCRAFT_tests.root'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter

process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
process.ShallowTrackClustersCombined = cms.EDProducer("ShallowTrackClustersProducerCombined",
                                      #Tracks=cms.InputTag("generalTracks",""),
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

#process.load('CalibTracker.SiStripCommon.ShallowTrackClustersProducer_cfi')
#process.tracksRefit = TrackRefitter.clone()


process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_cosmics')
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v14')
#process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v4')

#"global_tag": "90X_dataRun2_Prompt_v3" for run 293492 (from 902)
#"global_tag": "90X_dataRun2_Prompt_v2" for data (before 902)
#MC global tag: 90X_upgrade2017_realistic_v20
#

process.TrackRefitter.src = "ctfWithMaterialTracksP5"
#process.TrackRefitter.src = "generalTracks"
process.TrackRefitter.TTRHBuilder = "WithTrackAngle"
process.TrackRefitter.NavigationSchool = ""
process.ShallowTrackClustersCombined.Tracks = 'TrackRefitter'


#process.load('CalibTracker.SiStripCommon.ShallowTrackClustersProducer_cfi')
process.testTree = cms.EDAnalyzer(
   "ShallowTree",
   outputCommands = cms.untracked.vstring(
      'drop *',
      'keep *_ShallowTrackClustersCombined_*_*',
      )
   )


process.p = cms.Path(#process.MeasurementTrackerEvent*
                     process.TrackRefitter*
                     process.ShallowTrackClustersCombined*process.testTree)



#print process.dumpPython() 
