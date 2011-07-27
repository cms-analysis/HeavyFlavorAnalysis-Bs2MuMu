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
process.GlobalTag.globaltag = "GR_P_V20::All"

# ----------------------------------------------------------------------
# -- Input files
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# ----------------------------------------------------------------------
rootFileName = "bmt-MuOnia-Run2011A-PromptReco-v4.root"


process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )

# ----------------------------------------------------------------------
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMuonAndTrackCandidates_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFLambdas_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFSkipEvents_cff")


# ----------------------------------------------------------------------
process.mtDump1 = cms.EDAnalyzer(
    "HFMuonAndTrack",
    verbose       = cms.untracked.int32(0), 
    muonsLabel    = cms.untracked.InputTag("muons"),
    tracksLabel   = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt        = cms.untracked.double(5.0),
    trackPt       = cms.untracked.double(2.8),
    type          = cms.untracked.int32(1301), 
    maxDoca       = cms.untracked.double(0.1), 
    massLow       = cms.untracked.double(2.8), 
    massHigh      = cms.untracked.double(3.3)
    )

# ----------------------------------------------------------------------
process.mtDump2 = cms.EDAnalyzer(
    "HFMuonAndTrack",
    verbose       = cms.untracked.int32(0), 
    muonsLabel    = cms.untracked.InputTag("muons"),
    tracksLabel   = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt        = cms.untracked.double(6.0),
    trackPt       = cms.untracked.double(2.8),
    type          = cms.untracked.int32(1302), 
    maxDoca       = cms.untracked.double(0.1), 
    massLow       = cms.untracked.double(9.0), 
    massHigh      = cms.untracked.double(10.6)
    )


# ----------------------------------------------------------------------
process.mmDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(2.8),
    type               = cms.untracked.int32(1313), 
    vertexing          = cms.untracked.int32(1), 
    maxDoca            = cms.untracked.double(0.1), 
    massLow            = cms.untracked.double(4.7), 
    massHigh           = cms.untracked.double(6.2)
    )

# ----------------------------------------------------------------------
process.bupsikpDump = cms.EDAnalyzer(
    "HFBu2JpsiKp",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(2.8),
    psiMuons           = cms.untracked.int32(2),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.1),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0)
    )


# ----------------------------------------------------------------------
process.skipEvents = cms.EDFilter(
    "HFSkipEvents",
    verbose                      = cms.untracked.int32(0),
    filterOnPrimaryVertex        =  cms.untracked.int32(1),
    PrimaryVertexCollectionLabel = cms.untracked.InputTag('offlinePrimaryVertices'),
    filterOnTrackMaximum         =  cms.untracked.int32(-1),
    filterOnMuonMinimum          =  cms.untracked.int32(1),
    TrackCollectionLabel         = cms.untracked.InputTag('generalTracks'),
    MuonCollectionLabel          = cms.untracked.InputTag('muons')
    )



# ----------------------------------------------------------------------
process.p = cms.Path(
    process.skipEvents*
    process.seqPhysDeclBitSelection*
    process.recoStuffSequence*
    process.mtDump1*
    process.mtDump2*
    process.mmDump*
    process.bupsikpDump*
    process.tree
)
