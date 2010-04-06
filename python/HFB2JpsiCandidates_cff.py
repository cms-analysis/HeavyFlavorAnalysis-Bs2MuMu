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
    deltaR       = cms.untracked.double(1.5)
    )

# ----------------------------------------------------------------------
bdpsikstarDump = cms.EDAnalyzer(
    "HFBd2JpsiKstar",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    psiMuons     = cms.untracked.int32(2),
    trackPt       = cms.untracked.double(0.5)
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
    deltaR       = cms.untracked.double(1.5)
    )


# ######################################################################
# Sequences
# ######################################################################
B2JPsiSequence     = cms.Sequence(bupsikpDump*bdpsikstarDump*bspsiphiDump)
