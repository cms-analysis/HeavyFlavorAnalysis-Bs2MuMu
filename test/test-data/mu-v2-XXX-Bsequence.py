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
process.GlobalTag.globaltag = "GR10_P_V7::All"

# ----------------------------------------------------------------------
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# ----------------------------------------------------------------------
rootFileName = "bsequence-mu-v2-run-XXX.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )

# ----------------------------------------------------------------------
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2JpsiCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFSkipEvents_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.seqPhysDeclBitSelection*
    process.skipEvents*
    process.recoStuffSequence*
    process.B2JPsiSequence*
    process.tree
    )
