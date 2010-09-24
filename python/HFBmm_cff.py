import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
bmmDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(1.0),
    type               = cms.untracked.int32(1313), 
    vertexing          = cms.untracked.int32(1), 
    massLow            = cms.untracked.double(4.5), 
    massHigh           = cms.untracked.double(6.5)
    )


# ----------------------------------------------------------------------
bupsikpDump = cms.EDAnalyzer(
    "HFBu2JpsiKp",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(1.0),
    psiMuons           = cms.untracked.int32(2),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.1),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0)
    )

# ----------------------------------------------------------------------
bspsiphiDump = cms.EDAnalyzer(
    "HFBs2JpsiPhi",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(1.0),
    psiMuons           = cms.untracked.int32(2),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.1),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0)
    )

# ######################################################################
# Sequences
# ######################################################################
BmmSequence     = cms.Sequence(bmmDump*bupsikpDump*bspsiphiDump)
