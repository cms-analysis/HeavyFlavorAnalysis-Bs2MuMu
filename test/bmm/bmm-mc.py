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
process.GlobalTag.globaltag = "START37_V6::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
   "/store/mc/Summer10/JPsiToMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START37_V5_S09-v1/0001/3E427161-D884-DF11-8658-90E6BA19A1F9.root"
#  "/store/user/ursl/production/Winter10/BsToMuMu/BsToMuMu_7TeV_RAW2DIGI_RECO_START-10000.root"
 )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# ----------------------------------------------------------------------
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTree_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBmm_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFSkipEvents_cff")

# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "bmm-mc.root"

process.tree.fileName = rootFileName

# -- adjustments
process.triggerDump.HLTResultsLabel = cms.untracked.InputTag("TriggerResults::REDIGI37X")
process.hltrep.HLTResultsLabel = cms.untracked.InputTag("TriggerResults::REDIGI37X")
process.mcrecoStuffSequence     = cms.Sequence(process.stuffDump*process.trkDump*process.muonDump)

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.seqPhysDeclBitSelection*
    process.skipEvents*
    process.MCTruthSequence*
    process.mcrecoStuffSequence*
    process.truthBmmSequence*
    process.BmmSequence*
    process.tree
)
