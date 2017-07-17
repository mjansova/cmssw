from CalibTracker.SiStripCommon.step2_DIGI_L1_DIGI2RAW_HLT_2016 import * #good name

import sys
#process.source.fileNames = cms.untracked.vstring(sys.argv[3] + "step1.root")
process.FEVTDEBUGHLToutput.fileName = sys.argv[2] + ".root" 


