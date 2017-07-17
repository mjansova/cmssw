from SimGeneral.MixingModule.SiStripSimParameters_cfi import SiStripSimBlock 


import sys
value = map(float, sys.argv[3:])

#TIB
SiStripSimBlock.CouplingConstantDecIB1  = cms.vdouble(value[0], value[1], value[2])                   
SiStripSimBlock.CouplingConstantDecIB2  = cms.vdouble(value[0], value[1], value[2])                    
SiStripSimBlock.CouplingConstantPeakIB1 = cms.vdouble(value[0], value[1], value[2])                          
SiStripSimBlock.CouplingConstantPeakIB2 = cms.vdouble(value[0], value[1], value[2])                           
#TOB
SiStripSimBlock.CouplingConstantDecOB1  = cms.vdouble(value[0], value[1], value[2])                
SiStripSimBlock.CouplingConstantDecOB2  = cms.vdouble(value[0], value[1], value[2])                   
SiStripSimBlock.CouplingConstantPeakOB1 = cms.vdouble(value[0], value[1], value[2])                           
SiStripSimBlock.CouplingConstantPeakOB2 = cms.vdouble(value[0], value[1], value[2])

from CalibTracker.SiStripCommon.step2_DIGI_L1_DIGI2RAW_HLT_2017CosCRUZET import * #good name

process.FEVTDEBUGoutput.fileName = sys.argv[2] + ".root" 

print(sys.argv[3])

