# ######################################################################
# t3ui02
# /shome/ursl/dana/CMSSW_3_5_6/src/HeavyFlavorAnalysis/Bs2MuMu/test/test-data/goodcoll
# mkPyFiles -f 100413.t1 -t ../dana-XXXX.py -e 200000 -s v05 -r
# ./dana-v05-132605-0000.py with 186277 events
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
process.GlobalTag.globaltag = "GR10_P_V9::All"


# ----------------------------------------------------------------------

# Dataset: /MuOnia/Run2010B-PromptReco-v2/RECO

process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
	'/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/148/068/DC4E6C99-5DDB-DF11-8F19-001617C3B6DC.root',
	'/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/148/064/6439CE7A-48DB-DF11-B35E-003048F118C2.root',
	'/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/148/062/5E34943B-36DB-DF11-A373-001617C3B6CC.root',
	'/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/148/058/DEE05B92-6DDB-DF11-892A-001617C3B6FE.root',
	'/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/148/058/C60D8491-6DDB-DF11-99A3-001617C3B69C.root',
	'/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/148/058/BAEC55ED-6CDB-DF11-8028-000423D98950.root',
	'/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/148/058/B03BA070-6EDB-DF11-BBEA-001D09F24600.root',
	'/store/data/Run2010B/MuOnia/RECO/PromptReco-v2/000/148/058/6A5A8AEB-6CDB-DF11-9590-000423D9997E.root'
	)
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "lambdaB-XXXX.root"

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
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2JpsiCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFLambdas_cff")
process.HFLambdasDump.verbose = cms.untracked.int32(0)
process.HFLambdasDump.maxDoca = cms.untracked.double(.1)
process.HFLambdasDump.minPocaJpsi = cms.untracked.double(.1)
process.HFLambdasDump.minPocaL0 = cms.untracked.double(1)
process.HFLambdasDump.psiMuons = cms.untracked.int32(2)
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMuonAndTrackCandidates_cff")

#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2MuCandidates_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFDimuonsCandidates_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
#    process.MCTruthSequence*
    process.recoStuffSequence*
#    process.bmmDump*
#    process.truthAllSequence*
#    process.mtDump*
    process.HFLambdasSequence*
#    process.B2JPsiSequence*
    process.tree
)




