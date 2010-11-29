import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
HFLambdasDump = cms.EDAnalyzer(
    "HFLambdas",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    muonPt             = cms.untracked.double(1.0),
    pionPt             = cms.untracked.double(0.7),
    protonPt           = cms.untracked.double(0.7),
    psiMuons	       = cms.untracked.int32(2),
    psiWindow	       = cms.untracked.double(0.3),
    L0Window           = cms.untracked.double(0.3),
    ksWindow	       = cms.untracked.double(0.3),
    LbWindow           = cms.untracked.double(0.8),
    deltaR             = cms.untracked.double(99.),
    maxDoca            = cms.untracked.double(0.1),
    pAngle	       = cms.untracked.double(0.1)
    )

# ######################################################################
# Sequences
# ######################################################################
HFLambdasSequence     = cms.Sequence(HFLambdasDump)

