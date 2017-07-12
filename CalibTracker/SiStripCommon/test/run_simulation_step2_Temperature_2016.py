from SimGeneral.MixingModule.SiStripSimParameters_cfi import SiStripSimBlock 


import sys
value = map(float, sys.argv[3:])
SiStripSimBlock.Temperature = cms.double(value[0])

from CalibTracker.SiStripCommon.step2_DIGI_L1_DIGI2RAW_HLT_2016 import * #good name

process.FEVTDEBUGHLToutput.fileName = sys.argv[2] + ".root" 

print(sys.argv[3])

