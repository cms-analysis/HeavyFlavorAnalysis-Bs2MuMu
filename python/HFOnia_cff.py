import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
psiDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(4.0),
    type               = cms.untracked.int32(1313), 
    vertexing          = cms.untracked.int32(1), 
    maxDoca            = cms.untracked.double(0.1), 
    massLow            = cms.untracked.double(2.8), 
    massHigh           = cms.untracked.double(3.4)
    )

# ----------------------------------------------------------------------
upsDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(4.0),
    type               = cms.untracked.int32(1313), 
    vertexing          = cms.untracked.int32(1), 
    maxDoca            = cms.untracked.double(0.1), 
    massLow            = cms.untracked.double(8.9), 
    massHigh           = cms.untracked.double(9.9)
    )


# ######################################################################
# Sequences
# ######################################################################
oniaSequence     = cms.Sequence(psiDump*upsDump)
