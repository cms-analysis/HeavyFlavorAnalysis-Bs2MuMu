import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
dimuonsDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(1.0),
    type               = cms.untracked.int32(1313), 
    vertexing          = cms.untracked.int32(1), 
    massLow            = cms.untracked.double(0.5), 
    massHigh           = cms.untracked.double(9999.0)
    )
