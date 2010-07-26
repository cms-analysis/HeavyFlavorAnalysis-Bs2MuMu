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
# POOLSOURCE


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
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2JpsiCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2MuCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBmmCandidates_cff")

#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMuonAndTrackCandidates_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFDimuonsCandidates_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.MCTruthSequence*
    process.recoStuffSequence*
    process.bmmDump*
    process.truthAllSequence*
    process.tree
)




