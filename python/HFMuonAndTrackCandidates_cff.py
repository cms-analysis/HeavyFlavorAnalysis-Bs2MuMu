import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
mtDump = cms.EDAnalyzer(
    "HFMuonAndTrack",
    verbose       = cms.untracked.int32(0), 
    muonsLabel    = cms.untracked.InputTag("muons"),
    tracksLabel   = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt        = cms.untracked.double(1.0),
    trackPt       = cms.untracked.double(0.5),
    type          = cms.untracked.int32(1300), 
    massLow       = cms.untracked.double(0.5), 
    massHigh      = cms.untracked.double(12.0)
    )
