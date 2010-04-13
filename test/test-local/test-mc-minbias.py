# ######################################################################
# t3ui02
# /shome/ursl/dana/CMSSW_3_5_6/src/HeavyFlavorAnalysis/Bs2MuMu/test/test-data/goodcoll
# mkPyFiles -f 100413.mc.minbias -t ../mana-XXXX.py -e 200000 -s v05 -r
# ./mana-v05-4-0000.py with 199047 events
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
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FEE68846-F73B-DF11-98C8-002618943959.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FCECCB5E-F73B-DF11-8F01-002618943946.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FCAD1144-F73B-DF11-97E1-00261894391B.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FC59028E-F73B-DF11-B3F1-003048678B7E.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FC549636-F73B-DF11-B453-003048678BAA.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FC09D952-F73B-DF11-A53D-003048678B36.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FAD9E1D5-F73B-DF11-908A-003048678FE0.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FA856428-F93B-DF11-A84F-00261894391B.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FA70EC50-F43B-DF11-8827-0026189438C1.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FA35A244-F43B-DF11-BD0B-00261894396A.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FA1488B8-F73B-DF11-9E79-003048678B36.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F8FCB342-F73B-DF11-AADA-00261894396E.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F689B943-F73B-DF11-B702-003048678B36.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F67EED24-F73B-DF11-AEEC-003048678F74.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F454EDDB-F83B-DF11-ACCC-00304867BFF2.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F44C3A4D-F73B-DF11-B42D-003048678B94.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F42B6EA3-F73B-DF11-B1DC-00304867906C.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F4008891-F73B-DF11-921A-00304867906C.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F2F94152-F43B-DF11-A4EC-0026189438B3.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F2F0CF88-F73B-DF11-99AC-003048678B94.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F26C3E85-F73B-DF11-AC31-002618FDA277.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/F2159C33-F73B-DF11-ADDA-0030486790B8.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/EEC95E41-F43B-DF11-8ABE-002618943977.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/EE7A3E62-F73B-DF11-9DA0-003048678FC4.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/EC7E0D5F-F73B-DF11-A06C-00261894389F.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/EA53544D-F93B-DF11-9442-003048678B12.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E8B45E9F-F73B-DF11-AA57-003048678B36.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E891C380-F73B-DF11-88CA-0026189438CE.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E887B64C-F73B-DF11-A968-003048678ED4.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E431045E-F73B-DF11-968C-00261894389F.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E28F1523-F73B-DF11-BCC4-0030486790A6.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E25D5984-F73B-DF11-8E4D-002618943884.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E25D466F-F73B-DF11-B3AE-003048678BAA.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E0D6FC8C-F73B-DF11-8C49-003048678B8E.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/E037AE46-F43B-DF11-B726-002618943886.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/DC0BC822-F73B-DF11-AE96-0030486790B8.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/DAB75888-F73B-DF11-939A-00261894389C.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/D84FF25D-F73B-DF11-996E-00261894390B.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/D827CE43-F93B-DF11-B42E-003048679294.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/D684F444-F43B-DF11-9ABC-00261894396A.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/D66B0B9E-F73B-DF11-B747-003048678B94.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/D2276889-F73B-DF11-8365-0030486790B8.root",
 "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/CED6124D-F73B-DF11-8396-0030486790A6.root"
    )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "test-mc-minbias.root"

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
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFCharm_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMuonAndTrackCandidates_cff")

#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2MuCandidates_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFDimuonsCandidates_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.MCTruthSequence*
    process.recoStuffSequence*
#    process.bmmDump*
    process.truthAllSequence*
    process.mtDump*
    process.HFCharmSequence*
    process.B2JPsiSequence*
    process.tree
)
