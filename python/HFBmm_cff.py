import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
bmmDump = cms.EDAnalyzer(
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
    vertexing          = cms.untracked.int32(1) 
    )


# ----------------------------------------------------------------------
bupsikpDump = cms.EDAnalyzer(
    "HFBu2JpsiKp",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(3.5),
    psiMuons           = cms.untracked.int32(2),
    psiWindow          = cms.untracked.double(0.2),
    BuWindow           = cms.untracked.double(1.0),
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
    muonPt             = cms.untracked.double(3.5),
    psiMuons           = cms.untracked.int32(2),
    psiWindow          = cms.untracked.double(0.2),
    phiWindow          = cms.untracked.double(0.2),
    BsWindow           = cms.untracked.double(1.0),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.1),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0)
    )

# ######################################################################
# Sequences
# ######################################################################
bmmSequence     = cms.Sequence(bmmDump*bupsikpDump*bspsiphiDump)
