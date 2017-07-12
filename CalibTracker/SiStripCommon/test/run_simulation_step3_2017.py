from CalibTracker.SiStripCommon.step3_RAW2DIGI_L1Reco_RECO_2017 import * #good name
import sys

process.source.fileNames = cms.untracked.vstring("file:"+ sys.argv[2] + ".root") 

process.RECOSIMoutput.fileName = cms.untracked.string(sys.argv[2] + "step3.root")

print(sys.argv[2])

