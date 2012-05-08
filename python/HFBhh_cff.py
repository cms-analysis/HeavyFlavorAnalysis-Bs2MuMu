import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
hhDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    trackPt            = cms.untracked.double(4.0),
    track1Mass         = cms.untracked.double(0.1057),
    track2Mass         = cms.untracked.double(0.1057),
    massLow            = cms.untracked.double(4.0),
    massHigh           = cms.untracked.double(7.0),
    maxDoca            = cms.untracked.double(0.025),
    pvWeight           = cms.untracked.double(0.70),
    type              = cms.untracked.int32(211211)
    )

# ######################################################################
# Sequences
# ######################################################################
bhhSequence     = cms.Sequence(hhDump)
