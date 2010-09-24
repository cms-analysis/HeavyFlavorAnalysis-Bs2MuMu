# ######################################################################
# t3ui02
# /shome/ursl/dana/CMSSW_3_7_0_patch2/src/HeavyFlavorAnalysis/Bs2MuMu/test/test-data
# file list contains 93768 events
# mkPyFiles -f ./bla.3 -t dimuons-XXXX.py -r -e 200000 -s bla
# ./dimuons-bla-0139980-0000.py with 200000 events, skipEvents = 0
# ######################################################################
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('HLTrigReport')
process.MessageLogger.categories.append('L1GtTrigReport')
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
process.GlobalTag.globaltag = "GR10_P_V10::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
  "/store/data/Run2010A/Mu/RECO/v4/000/140/401/C270D8D6-3093-DF11-B83B-003048F118C6.root"
 )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# ----------------------------------------------------------------------
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTree_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBmm_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFSkipEvents_cff")


# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "bmm-test.root"

process.tree.fileName = rootFileName


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.seqPhysDeclBitSelection*
    process.skipEvents*
#    process.MCTruthSequence*
    process.recoStuffSequence*
#    process.truthAllSequence*
#    process.mtDump*
#    process.B2JPsiSequence*
#    process.HFCharmSequence*
#    process.HFMuCharmDump*
    process.BmmSequence* 
    process.tree
)

