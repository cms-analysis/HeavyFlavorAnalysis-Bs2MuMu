import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
dimuonsDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(3.5),
    type               = cms.untracked.int32(1313), 
    massLow            = cms.untracked.double(4.2), 
    massHigh           = cms.untracked.double(6.7),
    maxDoca            = cms.untracked.double(0.1),
    vertexing          = cms.untracked.int32(1), 
    pvWeight           = cms.untracked.double(0.6)
    )


# ######################################################################
# Sequences
# ######################################################################
dimuonsSequence     = cms.Sequence(dimuonsDump)
