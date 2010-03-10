import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


# ----------------------------------------------------------------------
# -- Database configuration
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START3X_V21B::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10000.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10001.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10002.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10003.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10004.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10005.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10006.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10007.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10008.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10009.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10010.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10011.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10012.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10013.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10014.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10015.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10016.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10018.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10019.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10020.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10021.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10022.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10023.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10024.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10025.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10026.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10027.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10028.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10029.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10030.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10031.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10032.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10033.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10034.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10035.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10036.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10037.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10038.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10039.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10040.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10041.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10042.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10043.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10044.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10045.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10046.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10047.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10048.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10049.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10050.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10051.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10052.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10053.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10054.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10055.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10056.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10057.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10058.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10059.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10060.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10061.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10062.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10063.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10064.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10065.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10066.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10067.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10068.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10069.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10070.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10071.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10072.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10073.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10074.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10075.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10076.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10077.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10078.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10079.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10080.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10081.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10082.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10083.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10084.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10085.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10086.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10087.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10088.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10089.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10090.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10091.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10092.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10093.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10094.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10095.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10096.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10097.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10098.root", 
    "/store/user/ursl/production/Winter10/BdToMuMu/BdToMuMu_7TeV_RAW2DIGI_RECO_START-10099.root"
    )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# ----------------------------------------------------------------------
process.HepPDTESSource = cms.ESSource(
    "HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/particle.tbl')
    )
process.genParticles = cms.EDProducer(
    "GenParticleProducer",
    saveBarCodes          = cms.untracked.bool(True),
    src                   = cms.InputTag("generator"),
    abortOnUnknownPDGCode = cms.untracked.bool(False)
    )


# ----------------------------------------------------------------------
process.genDump = cms.EDAnalyzer(
    "HFDumpGenerator",
    generatorCandidates = cms.untracked.string('genParticles'),
    generatorEvent = cms.untracked.string('generator')
    )


# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "test-BdToMuMu.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )


# ----------------------------------------------------------------------
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.trkDump = cms.EDAnalyzer(
    "HFDumpTracks",
    verbose                        = cms.untracked.int32(0),
    tracksLabel                    = cms.untracked.InputTag('generalTracks'),
    primaryVertexCollectionLabel   = cms.untracked.InputTag('offlinePrimaryVertices'),
    generatorEventLabel            = cms.untracked.InputTag('generator'),
    muonsLabel                     = cms.untracked.InputTag("muons"),
    calomuonsLabel                 = cms.untracked.InputTag("calomuons"),
    trackingParticlesLabel         = cms.untracked.InputTag('trackingParticles'),
    associatorLabel                = cms.untracked.InputTag('TrackAssociatorByHits'),
    doTruthMatching                = cms.untracked.int32(3),
    simTracksLabel                 = cms.untracked.InputTag('allLayer1TrackCands')
    )

# ----------------------------------------------------------------------
process.stuffDump = cms.EDAnalyzer(
    "HFDumpStuff",
    verbose                  = cms.untracked.int32(0),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    PrimaryVertexTracksLabel = cms.untracked.InputTag("generalTracks")
    )


# ----------------------------------------------------------------------
process.muonDump = cms.EDAnalyzer(
    "HFDumpMuons",
    verbose         = cms.untracked.int32(0),
    muonsLabel      = cms.untracked.InputTag("muons"),
    calomuonsLabel  = cms.untracked.InputTag("calomuons"),
    doTruthMatching = cms.untracked.int32(0),
    )


# ----------------------------------------------------------------------
process.triggerDump = cms.EDAnalyzer(
    "HFDumpTrigger",
    verbose                 = cms.untracked.int32(0),
    L1GTReadoutRecordLabel  = cms.untracked.InputTag("gtDigis"), 
    hltL1GtObjectMap        = cms.untracked.InputTag("hltL1GtObjectMap"), 
    L1MuonsLabel            = cms.untracked.InputTag("hltL1extraParticles"), 
    HLTResultsLabel         = cms.untracked.InputTag("TriggerResults::HLT"), 
    TriggerEventLabel       = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"), 
    hltLabel                = cms.untracked.InputTag("TriggerResults::HLT"), 
    )

# ----------------------------------------------------------------------
process.bmmDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    type         = cms.untracked.int32(1313), 
    vertexing    = cms.untracked.int32(0), 
    massLow      = cms.untracked.double(0.5), 
    massHigh     = cms.untracked.double(12.0)
    )


# ----------------------------------------------------------------------
process.truthBsToMuMuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(80),
    GenType      = cms.untracked.int32(-80),
    daughtersID  = cms.untracked.vint32(13, -13)
    )

# ----------------------------------------------------------------------
process.truthBsToMuMuGaDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(81),
    GenType      = cms.untracked.int32(-81),
    daughtersID  = cms.untracked.vint32(13, -13, 22)
    )

# ----------------------------------------------------------------------
process.truthBsToKKDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(82),
    GenType      = cms.untracked.int32(-82),
    daughtersID  = cms.untracked.vint32(321, -321)
    )

# ----------------------------------------------------------------------
process.truthBsToKPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(83),
    GenType      = cms.untracked.int32(-83),
    daughtersID  = cms.untracked.vint32(321, -211)
    )

# ----------------------------------------------------------------------
process.truthBsToPiPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(84),
    GenType      = cms.untracked.int32(-84),
    daughtersID  = cms.untracked.vint32(211, -211)
    )

# ----------------------------------------------------------------------
process.truthBuTo3MuNuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(70),
    GenType      = cms.untracked.int32(-70),
    daughtersID  = cms.untracked.vint32(13, 13, -13, 14)
    )

# ----------------------------------------------------------------------
process.truthBdToMuMuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(90),
    GenType      = cms.untracked.int32(-90),
    daughtersID  = cms.untracked.vint32(13, -13)
    )

# ----------------------------------------------------------------------
process.truthBdToPiPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(91),
    GenType      = cms.untracked.int32(-91),
    daughtersID  = cms.untracked.vint32(211, -211)
    )

# ----------------------------------------------------------------------
process.truthLambdaBToPPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(5122),
    type         = cms.untracked.int32(60),
    GenType      = cms.untracked.int32(-60),
    daughtersID  = cms.untracked.vint32(2212, -211)
    )

# ----------------------------------------------------------------------
process.truthLambdaBToPKDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(5122),
    type         = cms.untracked.int32(61),
    GenType      = cms.untracked.int32(-61),
    daughtersID  = cms.untracked.vint32(2212, -321)
    )

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genParticles* 
    process.genDump*
    process.stuffDump*
    process.trkDump*
    process.muonDump*
    process.triggerDump*
    process.bmmDump*
    process.truthBsToMuMuDump*
    process.truthBdToMuMuDump*
    process.truthBsToKKDump*
    process.truthBsToKPiDump*
    process.truthBdToPiPiDump*
    process.truthBsToPiPiDump*
    process.truthLambdaBToPPiDump*
    process.truthLambdaBToPKDump*
    process.truthBuTo3MuNuDump*
    process.truthBsToMuMuGaDump*
    process.tree
)




