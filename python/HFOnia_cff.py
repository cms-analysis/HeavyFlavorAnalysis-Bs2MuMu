import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
psiDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonPt             = cms.untracked.double(4.0),
    type               = cms.untracked.int32(1313), 
    maxDoca            = cms.untracked.double(0.1), 
    massLow            = cms.untracked.double(2.7), 
    massHigh           = cms.untracked.double(4.5),
    pvWeight           = cms.untracked.double(0.6),
    nbrMuons           = cms.untracked.int32(2),
    closeToMuons       = cms.untracked.bool(False)
    )

# ----------------------------------------------------------------------
upsDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonPt             = cms.untracked.double(4.0),
    type               = cms.untracked.int32(1313), 
    maxDoca            = cms.untracked.double(0.1), 
    massLow            = cms.untracked.double(8.0), 
    massHigh           = cms.untracked.double(12.0),
    pvWeight           = cms.untracked.double(0.6),
    nbrMuons           = cms.untracked.int32(2),
    closeToMuons       = cms.untracked.bool(False)
    )

# ######################################################################
# Sequences
# ######################################################################
oniaSequence     = cms.Sequence(psiDump*upsDump)
