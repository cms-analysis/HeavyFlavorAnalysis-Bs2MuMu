import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
HFLambdasDump = cms.EDAnalyzer(
    "HFLambdas",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    muonPt             = cms.untracked.double(0.7),
    pionPt             = cms.untracked.double(0.1),
    protonPt           = cms.untracked.double(0.2),
    psiMuons	       = cms.untracked.int32(2),
    psiWindow	       = cms.untracked.double(0.3),
    L0Window           = cms.untracked.double(0.3),
    ksWindow	       = cms.untracked.double(0.3),
    LbWindow           = cms.untracked.double(0.8),
    deltaR             = cms.untracked.double(99.),
    maxDoca            = cms.untracked.double(0.5),
    pAngle	       = cms.untracked.double(0.1),
    useV0              = cms.untracked.bool(False)
)

HFLambdasDumpV0 = cms.EDAnalyzer(
    "HFLambdas",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    muonPt             = cms.untracked.double(0.7),
    pionPt             = cms.untracked.double(0.1),
    protonPt           = cms.untracked.double(0.2),
    psiMuons	       = cms.untracked.int32(2),
    psiWindow	       = cms.untracked.double(0.3),
    ksWindow	       = cms.untracked.double(0.3),
    LbWindow           = cms.untracked.double(0.8),
    deltaR             = cms.untracked.double(99.),
    pAngle	       = cms.untracked.double(0.1),
    useV0              = cms.untracked.bool(True)
)

# ######################################################################
# Sequences
# ######################################################################
HFLambdasSequence     = cms.Sequence(HFLambdasDump*HFLambdasDumpV0)

