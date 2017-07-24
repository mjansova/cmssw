# Auto generated configuration file
# using:
#7.2 Cosmics_UP17+Cosmics_UP17+DIGICOS_UP17+RECOCOS_UP17+ALCACOS_UP17+HARVESTCOS_UP17 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: UndergroundCosmicMu_cfi.py --conditions auto:phase1_2017_cosmics -s GEN,SIM -n 500 --era Run2_2017 --eventcontent RAWSIM --relval 666000,100000 --scenario cosmics --datatier GEN-SIM --beamspot Realistic8TeVCollision --fileout file:step1.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('SimGeneral.MixingModule.mixCosmics_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimNOBEAM_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('UndergroundCosmicMu_cfi.py nevts:500'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:step1CRAFT.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_cosmics', '')

process.cosmicInTracker = cms.EDFilter("CosmicGenFilterHelix",
    charges = cms.vint32(1, -1),
    doMonitor = cms.untracked.bool(False),
    maxZ = cms.double(212.0),
    minP = cms.double(0.0),
    minPt = cms.double(0.0),
    minZ = cms.double(-212.0),
    pdgIds = cms.vint32(-13, 13),
    propagator = cms.string('SteppingHelixPropagatorAlong'),
    radius = cms.double(80.0),
    src = cms.InputTag("generator","unsmeared")
)


process.cosmicInPixelLoose = cms.EDFilter("CosmicGenFilterHelix",
    charges = cms.vint32(1, -1),
    doMonitor = cms.untracked.bool(False),
    maxZ = cms.double(100.0),
    minP = cms.double(0.0),
    minPt = cms.double(0.0),
    minZ = cms.double(-100.0),
    pdgIds = cms.vint32(-13, 13),
    propagator = cms.string('SteppingHelixPropagatorAlong'),
    radius = cms.double(20.0),
    src = cms.InputTag("generator","unsmeared")
)


process.generator = cms.EDProducer("CosMuoGenProducer",
    AcptAllMu = cms.bool(False),
    ClayWidth = cms.double(50000.0),
    ElossScaleFactor = cms.double(1.0),
    MTCCHalf = cms.bool(False),
    MaxEnu = cms.double(10000.0),
    MaxP = cms.double(3000.0),
    MaxPhi = cms.double(360.0),
    MaxT0 = cms.double(12.5),
    MaxTheta = cms.double(75.0),
    MinEnu = cms.double(10.0),
    MinP = cms.double(10.0),
    MinP_CMS = cms.double(-1.0),
    MinPhi = cms.double(0.0),
    MinT0 = cms.double(-12.5),
    MinTheta = cms.double(0.0),
    MultiMuon = cms.bool(False),
    MultiMuonFileFirstEvent = cms.int32(1),
    MultiMuonFileName = cms.string('CORSIKAmultiMuon.root'),
    MultiMuonNmin = cms.int32(2),
    NuProdAlt = cms.double(7500000.0),
    PlugVx = cms.double(0.0),
    PlugVz = cms.double(-14000.0),
    RadiusOfTarget = cms.double(8000.0),
    RhoAir = cms.double(0.001214),
    RhoClay = cms.double(2.3),
    RhoPlug = cms.double(2.5),
    RhoRock = cms.double(2.5),
    RhoWall = cms.double(2.5),
    TIFOnly_constant = cms.bool(False),
    TIFOnly_linear = cms.bool(False),
    TrackerOnly = cms.bool(False),
    Verbosity = cms.bool(False),
    ZCentrOfTarget = cms.double(0.0),
    ZDistOfTarget = cms.double(15000.0)
)


process.SteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('SteppingHelixPropagatorAlong'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.ProductionFilterSequence = cms.Sequence(process.generator+process.cosmicInTracker)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
