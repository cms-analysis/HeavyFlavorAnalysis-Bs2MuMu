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
process.GlobalTag.globaltag = "START3X_V21A::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
#        "/store/user/starodumov/bsmm/reco/reco-10036.root",
        "/store/user/starodumov/kplus/reco/reco-kplus-10022.root"
#        "/store/user/starodumov/phi/reco/reco-phi-10014.root",
#        "/store/user/starodumov/kstar/reco/reco-kstar-10019.root"
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
    rootFileName = "test-signals.root"

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
process.bmtDump = cms.EDAnalyzer(
    "HFMuonAndTrack",
    verbose       = cms.untracked.int32(0), 
    muonsLabel    = cms.untracked.InputTag("muons"),
    tracksLabel   = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt        = cms.untracked.double(1.0),
    trackPt       = cms.untracked.double(0.5),
    type          = cms.untracked.int32(1300), 
    massLow       = cms.untracked.double(0.5), 
    massHigh      = cms.untracked.double(12.0)
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
process.bupsikpDump = cms.EDAnalyzer(
    "HFBu2JpsiKp",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    psiMuons     = cms.untracked.int32(2),
    trackPt      = cms.untracked.double(0.5),
    deltaR       = cms.untracked.double(1.5)
    )

# ----------------------------------------------------------------------
process.bdpsikstarDump = cms.EDAnalyzer(
    "HFBd2JpsiKstar",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    psiMuons     = cms.untracked.int32(2),
    trackPt       = cms.untracked.double(0.5)
    )

# ----------------------------------------------------------------------
process.bspsiphiDump = cms.EDAnalyzer(
    "HFBs2JpsiPhi",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    psiMuons     = cms.untracked.int32(2),
    trackPt      = cms.untracked.double(0.5),
    deltaR       = cms.untracked.double(1.5)
    )

# ----------------------------------------------------------------------
process.b2muD0Dump = cms.EDAnalyzer(
    "HFB2muD0",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    muonPt       = cms.untracked.double(1.0),
    trackPt      = cms.untracked.double(0.1),
    deltaR       = cms.untracked.double(1.5),
    deltaMD0     = cms.untracked.double(0.3),
    deltaMDs     = cms.untracked.double(0.3)
    )

# ----------------------------------------------------------------------
process.truthBsDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(67),
    GenType      = cms.untracked.int32(-67),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 333, 321, -321)
    )

# ----------------------------------------------------------------------
process.truthBuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(68),
    GenType      = cms.untracked.int32(-68),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 321)
    )

# ----------------------------------------------------------------------
process.truthBdDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(69),
    GenType      = cms.untracked.int32(-69),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 313, 321, 211)
    )

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genParticles* 
    process.genDump*
    process.stuffDump*
    process.trkDump*
    process.muonDump*
    process.triggerDump*
    process.bmtDump*
    process.bmmDump*
    process.bupsikpDump*
    process.bspsiphiDump*
    process.bdpsikstarDump*
    process.truthBuDump*
    process.truthBsDump*
    process.truthBdDump*
    process.tree
)




