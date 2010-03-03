# Auto generated configuration file
# using: 
# Revision: 1.151 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/PYTHIA6_BsToMuMu_7TeV_noPtCut_cff.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT -n 10 --conditions DESIGN_3X_V8B::All --datatier GEN-SIM-RAW --eventcontent RAWSIM --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')

##########################################################################################################
#              services
##########################################################################################################

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        g4SimHits = cms.untracked.uint32(311422),
        mix = cms.untracked.uint32(123215),
        VtxSmeared = cms.untracked.uint32(9823432),
        caloRecHits = cms.untracked.uint32(6123321),
        MuonSimHits = cms.untracked.uint32(951331),
        muonCSCDigis = cms.untracked.uint32(2514232),
        muonDTDigis = cms.untracked.uint32(667376),
        famosSimHits = cms.untracked.uint32(357319),
        famosPileUp = cms.untracked.uint32(98273),
        l1ParamMuons = cms.untracked.uint32(870926),
        paramMuons = cms.untracked.uint32(541225),
        muonRPCDigis = cms.untracked.uint32(2146964),
        siTrackerGaussianSmearingRecHits = cms.untracked.uint32(321480),
        simMuonCSCDigis = cms.untracked.uint32(121245),
        simMuonDTDigis = cms.untracked.uint32(123115),
        simMuonRPCDigis = cms.untracked.uint32(1215235),
        simSiPixelDigis = cms.untracked.uint32(1215235),
        simSiStripDigis = cms.untracked.uint32(1215235),
        simEcalUnsuppressedDigis = cms.untracked.uint32(1215235),
        simHcalUnsuppressedDigis = cms.untracked.uint32(1215235),
        ecalRecHit = cms.untracked.uint32(18734),
        ecalPreshowerRecHit = cms.untracked.uint32(18734),
        hbhereco = cms.untracked.uint32(18734),
        horeco = cms.untracked.uint32(18734),
        hfreco = cms.untracked.uint32(18734),
        generator = cms.untracked.uint32(10000),
        LHCTransport = cms.untracked.uint32(12345678)    
    ), 
    sourceSeed = cms.untracked.uint32(10000)
)
from IOMC.RandomEngine.RandomServiceHelper import  RandomNumberServiceHelper
randHelper =  RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.populate()
process.RandomNumberGeneratorService.saveFileName =  cms.untracked.string("Bs-RandomEngineState-10000.log")
#############################################################################

process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
#process.load('HLTrigger/Configuration/HLT_1E31_cff')
process.load('HLTrigger/Configuration/HLT_8E29_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.load("Configuration.StandardSequences.Generator_cff")
process.genParticles.abortOnUnknownPDGCode = False

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    annotation = cms.untracked.string('January10: Pythia6+EvtGen generation of Bs->MuMu, 7TeV, D6T tune'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/HeavyFlavorAnalysis/Bs2MuMu/test/Winter10/BsToMuMu_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_START.py,v $')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('OtherCMS', 
        'StdException', 
        'Unknown', 
        'BadAlloc', 
        'BadExceptionType', 
        'ProductNotFound', 
        'DictionaryNotFound', 
        'InsertFailure', 
        'Configuration', 
        'LogicError', 
        'UnimplementedFeature', 
        'InvalidReference', 
        'NullPointerError', 
        'NoProductSpecified', 
        'EventTimeout', 
        'EventCorruption', 
        'ScheduleExecutionFailure', 
        'EventProcessorFailure', 
        'FileInPathError', 
        'FileOpenError', 
        'FileReadError', 
        'FatalRootError', 
        'MismatchedInputFiles', 
        'ProductDoesNotSupportViews', 
        'ProductDoesNotSupportPtr', 
        'NotFound'
#    'ProductNotFound'
    )
)
# Input source
process.source = cms.Source("EmptySource")

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('BsToMuMu_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT-10000.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START3X_V21B::All'
process.MuFromBs = cms.EDFilter("PythiaFilter",
    MaxEta = cms.untracked.double(2.5),
    Status = cms.untracked.int32(1),
    MinEta = cms.untracked.double(-2.5),
    MotherID = cms.untracked.int32(531),
    ParticleID = cms.untracked.int32(13)
)
process.mumugenfilter = cms.EDFilter("MCParticlePairFilter",
    MaxEta = cms.untracked.vdouble(2.5, 2.5),
    Status = cms.untracked.vint32(1, 1),
    MinEta = cms.untracked.vdouble(-2.5, -2.5),
    ParticleID1 = cms.untracked.vint32(-13, 13),
    ParticleID2 = cms.untracked.vint32(-13, 13)
)
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    ExternalDecays = cms.PSet(
        EvtGen = cms.untracked.PSet(
            use_default_decay = cms.untracked.bool(False),
            decay_table = cms.FileInPath('GeneratorInterface/ExternalDecays/data/DECAY_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt.pdl'),
#            user_decay_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/Bs_mumu.dec'),
            user_decay_file = cms.FileInPath('HeavyFlavorAnalysis/Bs2MuMu/data/Bs_mumu.dec'),
            list_forced_decays = cms.vstring('MyB_s0', 
#                'Myanti-B_s0'),
                'anti-MyB_s0'),
            operates_on_particles = cms.vint32(0)
        ),
        parameterSets = cms.vstring('EvtGen')
    ),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(4e-05),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    crossSection = cms.untracked.double(55000000000.0),
    maxEventsToPrint = cms.untracked.int32(2),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(51)=10042     ! CTEQ6L1 structure function chosen', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model', 
            'MSTU(21)=1     ! Check on possible errors during program execution', 
            'PARP(82)=1.8387   ! pt cutoff for multiparton interactions', 
            'PARP(89)=1960. ! sqrts for which PARP82 is set', 
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter', 
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter', 
            'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
            'PARP(67)=2.5    ! amount of initial-state radiation', 
            'PARP(85)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(86)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(62)=1.25   ! ', 
            'PARP(64)=0.2    ! ', 
            'MSTP(91)=1     !', 
            'PARP(91)=2.1   ! kt distribution', 
            'PARP(93)=15.0  ! '),
        bbbarSettings = cms.vstring('PMAS(5,1)=4.8          ! b quark mass', 
            'MSEL=1                 ! Min Bias'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'bbbarSettings')
    )
)
process.ProductionFilterSequence = cms.Sequence(process.generator*process.MuFromBs*process.mumugenfilter)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.out_step])
# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq
