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
#process.GlobalTag.globaltag = "GR_P_V22A::All"
process.GlobalTag.globaltag = "FT_R_44_V9::All"

# ----------------------------------------------------------------------
# -- Input files
POOLSOURCE

# ----------------------------------------------------------------------
rootFileName = "bmm-2011-rereco-XXXX.root"


process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(1),
    requireCand  =  cms.untracked.bool(True),
    printFrequency = cms.untracked.int32(1000),
    fileName     = cms.untracked.string(rootFileName)
    )

# ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBmm_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFLambdas_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFSkipEvents_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
#    process.skipEvents*
    process.recoStuffSequence*
    process.bmmSequence*
#    process.HFLambdasSequence*
    process.tree
)
