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
process.GlobalTag.globaltag = "GR_P_V27::All"

# ----------------------------------------------------------------------
# POOLSOURCE
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
 "/store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/F869A3C8-0305-E111-A199-BCAEC518FF76.root",
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/F2AE041D-0F05-E111-B5BF-BCAEC5364C42.root   14369
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/F29EE774-3405-E111-BF80-003048D373AE.root   13493
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/F07BAD25-ED04-E111-86F1-001D09F24448.root   14663
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/E6F53E38-0105-E111-8789-001D09F242EF.root   13375
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/E61469EC-EF04-E111-A94D-003048D2C01E.root   13752
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/DC6A28ED-2405-E111-B41E-003048F1182E.root   14274
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/D67236F9-EA04-E111-83A6-003048D2C16E.root   13347
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/C87E7EB6-0005-E111-8754-003048D2C01A.root   14329
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/C221F453-F604-E111-A113-003048D2C0F0.root   12893
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/AEC94A90-1005-E111-AAAB-003048F118C2.root   13621
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/AE73DE67-EC04-E111-A750-001D09F24600.root   13815
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/AA81A030-1105-E111-8E72-E0CB4E5536AE.root   12729
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/A828704D-0705-E111-8248-001D09F241F0.root   14515
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/A22A4631-E804-E111-88F6-003048D2C16E.root   13499
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/8E1C4F03-E404-E111-916F-BCAEC53296F3.root   14622
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/8872EE1D-1B05-E111-8FC6-BCAEC518FF7A.root   13296
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/86DE4D15-1905-E111-B5FF-0025901D5C86.root   14013
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/80526FAE-E604-E111-9118-001D09F242EF.root   13989
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/7C86B7C8-0305-E111-843E-BCAEC532970D.root   13895
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/6A6C7898-1205-E111-BEA0-BCAEC5329718.root   12611
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/621DEAF6-ED05-E111-965E-BCAEC5329700.root   13625
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/467A6ED9-0F05-E111-BDA0-003048F1C58C.root   14163
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/445295D8-E604-E111-90EB-003048CF99BA.root   14007
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/3CD6CC4B-1F05-E111-99BE-003048F1C836.root   12851
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/38911571-0905-E111-A14E-0025901D623C.root   13474
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/322F3766-E704-E111-A7ED-003048D3756A.root   13408
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/2235867E-F304-E111-9334-0030486730C6.root   13368
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/1EB01258-0205-E111-87EE-003048D2C1C4.root   13206
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/1C9E791E-1B05-E111-98B1-BCAEC518FF56.root   13537
## /store/data/Run2011B/MuOnia/AOD/PromptReco-v1/000/180/250/0C050942-0C05-E111-8ED7-E0CB4E4408C4.root   13022


 )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "bhh-test.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )


# ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBmm_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBhh_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFPhysicsDeclared_cff")

#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
#process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")

#process.dstarDump.trackPt               = cms.untracked.double(1.0)
#process.dstarDump.slowPionPt            = cms.untracked.double(0.1)

# ----------------------------------------------------------------------
process.p = cms.Path(
#    process.MCTruthSequence*
    process.recoStuffSequence*
#    process.bmmSequence*
#    process.truthBmmSequence*
#    process.bmtSequence*
#    process.truthBmtSequence*
#    process.truthDstarToD0PiToKPiPi*
    process.hhDump*
#    process.truthBdToPiPiDump*
    process.tree
)
