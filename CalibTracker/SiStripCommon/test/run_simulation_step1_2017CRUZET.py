from CalibTracker.SiStripCommon.UndergroundCosmicMu_cfi_py_GEN_SIM_2017CosCRUZET import * #good name
import sys
process.RAWSIMoutput.fileName = cms.untracked.string(sys.argv[2] + "step1.root")

print(sys.argv[2])

