import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
bupsikpDump = cms.EDAnalyzer(
    "HFBu2JpsiKp",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    psiMuons     = cms.untracked.int32(1),
    trackPt      = cms.untracked.double(0.5),
    deltaR       = cms.untracked.double(99.0),
    maxDoca      = cms.untracked.double(0.1),
    maxD0        = cms.untracked.double(5.0),
    maxDz        = cms.untracked.double(25.0)
    )

# ----------------------------------------------------------------------
bdpsikstarDump = cms.EDAnalyzer(
    "HFBd2JpsiKstar",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    psiMuons     = cms.untracked.int32(1),
    trackPt       = cms.untracked.double(0.5),
    deltaR       = cms.untracked.double(99.0),
    maxDoca      = cms.untracked.double(0.1),
    maxD0        = cms.untracked.double(5.0),
    maxDz        = cms.untracked.double(25.0)
    )

# ----------------------------------------------------------------------
bspsiphiDump = cms.EDAnalyzer(
    "HFBs2JpsiPhi",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    psiMuons     = cms.untracked.int32(1),
    trackPt      = cms.untracked.double(0.5),
    deltaR       = cms.untracked.double(99.0),
    maxDoca      = cms.untracked.double(0.1),
    maxD0        = cms.untracked.double(5.0),
    maxDz        = cms.untracked.double(25.0)
    )

bdpsiksDump = cms.EDAnalyzer(
    "HFBd2JpsiKs",
    verbose = cms.untracked.int32(0),
    tracksLabel = cms.untracked.InputTag("generalTracks"),
    muonsLabel = cms.untracked.InputTag("muons"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt = cms.untracked.double(1.0),
    pionPt = cms.untracked.double(0.4),
    psiMuons = cms.untracked.int32(1),
    psiWindow = cms.untracked.double(0.3),
    ksWindow = cms.untracked.double(0.3),
    bdWindow = cms.untracked.double(0.8),
    fDelta = cms.untracked.double(99.0),
    vertexing = cms.untracked.int32(1)
)

# ######################################################################
# Sequences
# ######################################################################
#B2JPsiSequence     = cms.Sequence(bupsikpDump*bdpsikstarDump*bspsiphiDump*bdpsiksDump)
B2JPsiSequence     = cms.Sequence(bdpsikstarDump*bspsiphiDump*bdpsiksDump)
