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


SiStripSimBlock.electronPerAdcDec = cms.double(value[3])

from CalibTracker.SiStripCommon.step2_DIGI_L1_DIGI2RAW_HLT_2016 import * #good name

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


#process.source.fileNames = cms.untracked.vstring(sys.argv[2] + "step1.root")
process.FEVTDEBUGHLToutput.fileName = sys.argv[2] + ".root" 

print(sys.argv[3])

