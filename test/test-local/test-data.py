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
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
 "/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/0051B6BC-AA41-DF11-91CB-003048673EBA.root",
 "/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/08CF64D9-A941-DF11-9082-001A6478AB88.root",
 "/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/0EF92B4A-B041-DF11-BDCF-003048BAA5A8.root",
 "/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/1052ACD9-AE41-DF11-8711-0025B3E05DBE.root",
 "/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/146AE7C2-AA41-DF11-B9E1-0025B3E05C6E.root",
 "/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/1E49FAA8-AA41-DF11-A895-00E08179183D.root",
 "/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/1EFC3B2E-AE41-DF11-973C-0025B3E0653E.root",
 "/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/2EE1EDBA-AC41-DF11-9097-00E081791899.root"
    )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "test-data.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )

# ----------------------------------------------------------------------
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2JpsiCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFCharm_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMuonAndTrackCandidates_cff")

#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2MuCandidates_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFDimuonsCandidates_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
#    process.MCTruthSequence*
    process.recoStuffSequence*
#    process.bmmDump*
#    process.truthAllSequence*
    process.mtDump*
    process.HFCharmSequence*
    process.B2JPsiSequence*
    process.tree
)




