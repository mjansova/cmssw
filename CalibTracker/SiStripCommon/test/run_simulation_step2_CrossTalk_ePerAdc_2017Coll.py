# cmsRun name CRUZET2017 eperADC xtalk
from SimGeneral.MixingModule.SiStripSimParameters_cfi import SiStripSimBlock 


import sys
value = map(float, sys.argv[3:])

SiStripSimBlock.electronPerAdcDec = cms.double(value[0])
#TIB
SiStripSimBlock.CouplingConstantDecIB1  = cms.vdouble(value[1], value[2], value[3])                   
SiStripSimBlock.CouplingConstantDecIB2  = cms.vdouble(value[1], value[2], value[3])                    
#TOB
SiStripSimBlock.CouplingConstantDecOB1  = cms.vdouble(value[1], value[2], value[3])                
SiStripSimBlock.CouplingConstantDecOB2  = cms.vdouble(value[1], value[2], value[3])                   
#TID
SiStripSimBlock.CouplingConstantDecW1a  = cms.vdouble(value[1], value[2], value[3])                      
SiStripSimBlock.CouplingConstantDecW2a  = cms.vdouble(value[1], value[2], value[3])                   
SiStripSimBlock.CouplingConstantDecW3a  = cms.vdouble(value[1], value[2], value[3])
#TEC
SiStripSimBlock.CouplingConstantDecW1b  = cms.vdouble(value[1], value[2], value[3])                      
SiStripSimBlock.CouplingConstantDecW2b  = cms.vdouble(value[1], value[2], value[3])                      
SiStripSimBlock.CouplingConstantDecW3b  = cms.vdouble(value[1], value[2], value[3])                      
SiStripSimBlock.CouplingConstantDecW4   = cms.vdouble(value[1], value[2], value[3])                     
SiStripSimBlock.CouplingConstantDecW5   = cms.vdouble(value[1], value[2], value[3])                   
SiStripSimBlock.CouplingConstantDecW6   = cms.vdouble(value[1], value[2], value[3])                      
SiStripSimBlock.CouplingConstantDecW7   = cms.vdouble(value[1], value[2], value[3])

from CalibTracker.SiStripCommon.step2_DIGI_L1_DIGI2RAW_2017Coll import * #good name

process.source.fileNames = cms.untracked.vstring("file:/opt/sbg/scratch1/cms/mjansova/store/tmp/MCtuning/Coll2017step1.root")
process.FEVTDEBUGHLToutput.fileName=  sys.argv[2] + ".root"

print(sys.argv[2])

