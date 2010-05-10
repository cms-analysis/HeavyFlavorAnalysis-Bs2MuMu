import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
HFCharmDump = cms.EDAnalyzer(
    "HFCharm",
    verbose            = cms.untracked.int32(0), 
    maxTracks          = cms.untracked.int32(1000), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    useMuon            = cms.untracked.int32(0),
    phiWindow          = cms.untracked.double(99.0),
    DWindow            = cms.untracked.double(0.2),
    LcWindow           = cms.untracked.double(0.4),
    muonPt             = cms.untracked.double(1.0),
    protonPt           = cms.untracked.double(0.7),
    kaonPt             = cms.untracked.double(0.7),
    pionPt             = cms.untracked.double(0.7),
    trackPt            = cms.untracked.double(0.7),
    deltaR             = cms.untracked.double(99.),
    maxDoca            = cms.untracked.double(0.2),
    maxVtxChi2         = cms.untracked.double(10.),
    minVtxSigXY        = cms.untracked.double(1.),
    minVtxSig3d        = cms.untracked.double(-1.),
    minCosAngle        = cms.untracked.double(0.),
    type               = cms.untracked.int32(1) 
    )

HFMuCharmDump = cms.EDAnalyzer(
    "HFCharm",
    verbose            = cms.untracked.int32(0), 
    maxTracks          = cms.untracked.int32(1000), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    useMuon            = cms.untracked.int32(1),
    phiWindow          = cms.untracked.double(99.0),
    DWindow            = cms.untracked.double(0.3),
    LcWindow           = cms.untracked.double(0.4),
    muonPt             = cms.untracked.double(1.0),
    protonPt           = cms.untracked.double(0.5),
    kaonPt             = cms.untracked.double(0.5),
    pionPt             = cms.untracked.double(0.5),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.5),
    maxVtxChi2         = cms.untracked.double(100.),
    minVtxSigXY        = cms.untracked.double(-1.),
    minVtxSig3d        = cms.untracked.double(-1.),
    minCosAngle        = cms.untracked.double(-9.),
    type               = cms.untracked.int32(2) 
    )


# ######################################################################
# Sequences
# ######################################################################
HFCharmSequence     = cms.Sequence(HFCharmDump*HFMuCharmDump)
