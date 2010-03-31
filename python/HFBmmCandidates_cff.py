import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
bmmDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(1.0),
    type               = cms.untracked.int32(51313), 
    vertexing          = cms.untracked.int32(1), 
    massLow            = cms.untracked.double(4.0), 
    massHigh           = cms.untracked.double(7.0)
    )
