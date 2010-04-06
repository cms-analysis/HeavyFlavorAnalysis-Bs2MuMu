import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
HFCharmDump = cms.EDAnalyzer(
    "HFCharm",
    verbose            = cms.untracked.int32(0), 
    maxTracks          = cms.untracked.int32(100), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    muonPt             = cms.untracked.double(1.0),
    useMuon            = cms.untracked.int32(0),
    phiWindow          = cms.untracked.double(0.2),
    DWindow            = cms.untracked.double(0.2),
    trackPt            = cms.untracked.double(0.8),
    kaonPt             = cms.untracked.double(1.1),
    pionPt             = cms.untracked.double(0.8),
    deltaR             = cms.untracked.double(99.),
    type               = cms.untracked.int32(1) 
    )

HFMuCharmDump = cms.EDAnalyzer(
    "HFCharm",
    verbose            = cms.untracked.int32(0), 
    maxTracks          = cms.untracked.int32(100), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    muonPt             = cms.untracked.double(1.0),
    useMuon            = cms.untracked.int32(1),
    phiWindow          = cms.untracked.double(0.3),
    DWindow            = cms.untracked.double(0.3),
    trackPt            = cms.untracked.double(0.5),
    kaonPt             = cms.untracked.double(0.8),
    pionPt             = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(1.5),
    type               = cms.untracked.int32(2) 
    )


# ######################################################################
# Sequences
# ######################################################################
HFCharmSequence     = cms.Sequence(HFCharmDump*HFMuCharmDump)
