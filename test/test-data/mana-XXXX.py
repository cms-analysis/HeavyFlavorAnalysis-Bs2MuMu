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
process.GlobalTag.globaltag = "GR10_P_V4::All"


# ----------------------------------------------------------------------
# POOLSOURCE


# ----------------------------------------------------------------------
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTree_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2JpsiCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFCharm_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMuonAndTrackCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFSkipEvents_cff")

# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "mana-XXXX.root"

process.tree.fileName = rootFileName


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.seqPhysDeclBitSelection*
    process.skipEvents*
    process.MCTruthSequence*
    process.recoStuffSequence*
    process.truthAllSequence*
    process.mtDump*
    process.B2JPsiSequence*
    process.HFCharmSequence*
    process.tree
)
