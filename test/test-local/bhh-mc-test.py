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
process.GlobalTag.globaltag = "START44_V10::All"

# ----------------------------------------------------------------------
# POOLSOURCE
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
         "/store//user/ursl/production/Fall11/reco2e33/Bd2PiPi0/Bd2PiPi-reco-910000.root",
#         "/store//user/ursl/production/Fall11/reco2e33/Bd2PiPi0/Bd2KK-reco-930000.root",
#         "/store//user/ursl/production/Fall11/reco2e33/Bd2PiPi0/Bd2KPi-reco-920000.root",
#         "/store//user/ursl/production/Fall11/reco2e33/Bs2PiPi0/Bd2PiPi-reco-840000.root",
#         "/store//user/ursl/production/Fall11/reco2e33/Bs2PiPi0/Bd2KK-reco-820000.root",
#         "/store//user/ursl/production/Fall11/reco2e33/Bs2PiPi0/Bd2KPi-reco-830000.root",
 )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "bhh-mc-test.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )


# ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
# process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBmm_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBhh_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")

process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")

#process.dstarDump.trackPt               = cms.untracked.double(1.0)
#process.dstarDump.slowPionPt            = cms.untracked.double(0.1)

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.MCTruthSequence*
    process.recoStuffSequence*
#    process.bmmSequence*
#    process.truthBmmSequence*
#    process.bmtSequence*
#    process.truthBmtSequence*
#    process.truthDstarToD0PiToKPiPi*
    process.hhDump*
    process.truthBdToPiPiDump*
    process.tree
)
