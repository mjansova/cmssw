from CalibTracker.SiStripCommon.shallowTree_test_template import *
process.TFileService.fileName = 'bla.root'
import sys

#process.source.fileNames = cms.untracked.vstring('file:output/step3_RAW2DIGI_L1Reco_RECO_19.root') #test_shallowTrackAndClusterFullInfo2016run284078.root.root
process.source.fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/RelValMinBias_13/GEN-SIM-RECO/90X_mcRun2_asymptotic_v5-v1/00000/9A94CC7A-4B0F-E711-8AC9-0CC47A4D7616.root') #relval MC 90, 16500 ev #test_shallowTrackClusterNoPUMCrun2_mcTAG.root

inputStr = sys.argv[2]
found = inputStr.find("root")

if found > 0:
    process.source.fileNames = cms.untracked.vstring("file:"+inputStr) 
    process.TFileService.fileName = 'test_shallowTrackAndClusterNarrowInfo.root' + inputStr
    print(input)
else:
    process.TFileService.fileName = 'test_shallowTrackAndClusterNarrowInfoNoPU2016_data.root'

#process.source.fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/RelValMinBias_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v20_resub-v1/00000/7E9D0F0A-CC10-E711-BDBD-0CC47A7C354A.root') #test_shallowTrackClusterNoPUMCupgrade2017.root


#fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/ZeroBias/RECO/90X_dataRun2_PromptLike_v6_RelVal_zb2016H-v1/00000/1A2F6458-5F20-E711-884E-0025905A608C.root') #run H, 900, 1800 #test_shallowTrackClusterPUData.root 

#fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/RelValNuGun/GEN-SIM-RECO/PU25ns_90X_upgrade2017_realistic_v20_resub-v1/00000/16689437-7C14-E711-A30F-0025905A60CE.root') #test_shallowTrackClusterPUMCupgrade2017.root

#fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/RelValNuGun_UP15/GEN-SIM-RECO/PU25ns_90X_mcRun2_asymptotic_v5-v1/00000/8EA7CE2A-2810-E711-A617-0CC47A4C8E1E.root') #900 PU MCs 2000 ev test_shallowTrackClusterPUMCrun2_mc.root


#fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/RelValMinBias_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v20_resub-v1/00000/7E9D0F0A-CC10-E711-BDBD-0CC47A7C354A.root') #test_shallowTrackClusterNoPUMCupgrade2017.root

#fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_0_0/RelValMinBias_13/GEN-SIM-RECO/90X_mcRun2_asymptotic_v5-v1/00000/9A94CC7A-4B0F-E711-8AC9-0CC47A4D7616.root') #relval MC 90, 16500 ev #test_shallowTrackClusterNoPUMCrun2_mcTAG.root


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
