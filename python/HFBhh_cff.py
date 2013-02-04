import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
hhDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
	BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    trackPt            = cms.untracked.double(4.0),
    track1Mass         = cms.untracked.double(0.1057),
    track2Mass         = cms.untracked.double(0.1057),
    massLow            = cms.untracked.double(4.5),
    massHigh           = cms.untracked.double(6.5),
    maxDoca            = cms.untracked.double(0.025),
    pvWeight           = cms.untracked.double(0.70),
    type               = cms.untracked.int32(211211),
	nbrMuons           = cms.untracked.int32(0),
	closeToMuons       = cms.untracked.bool(False)
    )

# ######################################################################
# Sequences
# ######################################################################
bhhSequence     = cms.Sequence(hhDump)
