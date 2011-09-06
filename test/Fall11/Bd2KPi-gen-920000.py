# Auto generated configuration file
# using: 
# Revision: 1.303.2.7 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: cd-gen-sim -s GEN,SIM --filetype EDM --conditions START42_V13::All --pileup NoPileUp --datamix NODATAMIXER --datatier GEN-SIM --eventcontent RAWSIM -n 10 --no_exec --filein file:genFragment.py

# -- standard configurations
import FWCore.ParameterSet.Config as cms
process = cms.Process('SIM')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


# ----------------------------------------------------------------------
# -- Process specific setup
# ----------------------------------------------------------------------
process.bfilter = cms.EDFilter(
        "PythiaFilter",
        MaxEta = cms.untracked.double(9999.),
        MinEta = cms.untracked.double(-9999.),
        ParticleID = cms.untracked.int32(511)
        )

process.decayfilter = cms.EDFilter(
        "PythiaDauVFilter",
	verbose         = cms.untracked.int32(1), 
	NumberDaughters = cms.untracked.int32(2), 
	ParticleID      = cms.untracked.int32(511),  
        DaughterIDs     = cms.untracked.vint32(321, -211),
	MinPt           = cms.untracked.vdouble(3.5, 3.5), 
	MinEta          = cms.untracked.vdouble(-2.5, -2.5), 
	MaxEta          = cms.untracked.vdouble( 2.5,  2.5)
        )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('Bd2KPi'),
    name = cms.untracked.string('Bd2KPi')
)

# Output definition
process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('Bd2KPi-gen-920000.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# -- Generator definitions
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    ExternalDecays = cms.PSet(
        EvtGen = cms.untracked.PSet(
            use_default_decay = cms.untracked.bool(False),
            decay_table = cms.FileInPath('GeneratorInterface/ExternalDecays/data/DECAY_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt.pdl'),
            user_decay_file = cms.FileInPath('HeavyFlavorAnalysis/Bs2MuMu/data/Bd_Kpi.dec'),
            list_forced_decays = cms.vstring('MyB0', 
                'Myanti-B0'),
            operates_on_particles = cms.vint32(0)
        ),
        parameterSets = cms.vstring('EvtGen')
    ),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(0.003),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    crossSection = cms.untracked.double(48440000000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.832 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.275 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        bbbarSettings = cms.vstring('MSEL = 1'),  
        parameterSets = cms.vstring('pythiaUESettings', 
            'bbbarSettings')
    )
)


process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        g4SimHits = cms.untracked.uint32(920000),
        mix = cms.untracked.uint32(920000),
        VtxSmeared = cms.untracked.uint32(920000),
        caloRecHits = cms.untracked.uint32(920000),
        MuonSimHits = cms.untracked.uint32(920000),
        muonCSCDigis = cms.untracked.uint32(920000),
        muonDTDigis = cms.untracked.uint32(920000),
        famosSimHits = cms.untracked.uint32(920000),
        famosPileUp = cms.untracked.uint32(920000),
        l1ParamMuons = cms.untracked.uint32(920000),
        paramMuons = cms.untracked.uint32(920000),
        muonRPCDigis = cms.untracked.uint32(920000),
        siTrackerGaussianSmearingRecHits = cms.untracked.uint32(920000),
        simMuonCSCDigis = cms.untracked.uint32(920000),
        simMuonDTDigis = cms.untracked.uint32(920000),
        simMuonRPCDigis = cms.untracked.uint32(920000),
        simSiPixelDigis = cms.untracked.uint32(920000),
        simSiStripDigis = cms.untracked.uint32(920000),
        simEcalUnsuppressedDigis = cms.untracked.uint32(920000),
        simHcalUnsuppressedDigis = cms.untracked.uint32(920000),
        ecalRecHit = cms.untracked.uint32(920000),
        ecalPreshowerRecHit = cms.untracked.uint32(920000),
        hbhereco = cms.untracked.uint32(920000),
        horeco = cms.untracked.uint32(920000),
        hfreco = cms.untracked.uint32(920000),
        generator = cms.untracked.uint32(920000)
    ), 
    sourceSeed = cms.untracked.uint32(920000)
)
from IOMC.RandomEngine.RandomServiceHelper import  RandomNumberServiceHelper
randHelper =  RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.populate() 
process.RandomNumberGeneratorService.saveFileName =  cms.untracked.string("RandomEngineState-920000.log")


# ----------------------------------------------------------------------
# -- rest of setup
# ----------------------------------------------------------------------
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000000)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Other statements
process.GlobalTag.globaltag = 'START42_V13::All'

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)

# filter all path with the production filter sequence
process.ProductionFilterSequence = cms.Sequence(process.generator+process.bfilter+process.decayfilter)
for path in process.paths:
        getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 
