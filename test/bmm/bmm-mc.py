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
process.GlobalTag.globaltag = "START38_V13::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
#   "/store/mc/Summer10/JPsiToMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START37_V5_S09-v1/0001/3E427161-D884-DF11-8658-90E6BA19A1F9.root"
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0023/22C16EF5-30E2-DF11-B4E2-00215E21DF18.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0023/22C16EF5-30E2-DF11-B4E2-00215E21DF18.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0023/02E8990D-73E2-DF11-A46F-00215E21DD44.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0008/FE574703-46D8-DF11-8BA5-003048C56D1A.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0008/802873D6-B2D8-DF11-87F5-001F29C9953C.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0008/562997FF-43D8-DF11-B8F9-003048976B48.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0008/4CCE5948-7FD8-DF11-93DA-0019BBEBD32C.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0008/269489BC-77D8-DF11-87EB-002618943C0F.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0008/10E30E18-86D8-DF11-AED2-00215E2211AC.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0008/00EC4D95-6CD8-DF11-8C20-001A6489C1AA.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/ECAA93BC-D4D5-DF11-9C2B-00E0812C7F31.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/D633FBBB-EBD5-DF11-A1F0-003048C56E80.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/D455B89B-98D5-DF11-816C-001A4B0A277E.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/B0F518BA-ACD5-DF11-AC09-001A4B0A271A.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/AEE02825-E8D6-DF11-85B0-001F29C945FE.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/ACB2C5E9-B7D5-DF11-8BD8-001A644E94E2.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/AA5C7E22-09D7-DF11-9844-001F29C76F6A.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/90E24668-9FD5-DF11-B092-003048976AD8.root",
   "/store/mc/Fall10/BsToMuMu_2MuPtFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START38_V12-v1/0004/78496870-F4D5-DF11-B1C6-00215E21DC7E.root"
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
    rootFileName = "/shome/ursl/root/bmm-mc.root"

process.tree.fileName = rootFileName

# -- adjustments
#process.triggerDump.HLTResultsLabel = cms.untracked.InputTag("TriggerResults::REDIGI37X")
#process.hltrep.HLTResultsLabel = cms.untracked.InputTag("TriggerResults::REDIGI37X")
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
