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
process.GlobalTag.globaltag = "START3X_V25::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
        "/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10041.root"
#        "/store/user/starodumov/kplus/reco/reco-kplus-10022.root"
#        "/store/user/starodumov/phi/reco/reco-phi-10014.root",
#        "/store/user/starodumov/kstar/reco/reco-kstar-10019.root"
    )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )



# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "validation-XXXX.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )

# ----------------------------------------------------------------------
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")


# ----------------------------------------------------------------------
trackList = "generalTracks"
process.truthBuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel          = cms.untracked.InputTag(trackList),
    verbose              = cms.untracked.int32(3),
    partialDecayMatching = cms.untracked.bool(True),
    motherID             = cms.untracked.int32(521),
    type                 = cms.untracked.int32(68),
    GenType              = cms.untracked.int32(-68),
#    daughtersID          = cms.untracked.vint32(443, 13, -13, 321)
    daughtersID          = cms.untracked.vint32(443, 321)
    )


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.MCTruthSequence*
    process.truthAllSequence*
    process.truthBuDump*
    process.tree
)




