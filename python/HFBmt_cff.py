import FWCore.ParameterSet.Config as cms



# ----------------------------------------------------------------------
mtDump1 = cms.EDAnalyzer(
    "HFMuonAndTrack",
    verbose       = cms.untracked.int32(0), 
    muonsLabel    = cms.untracked.InputTag("muons"),
    tracksLabel   = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt        = cms.untracked.double(5.0),
    trackPt       = cms.untracked.double(3.5),
    type          = cms.untracked.int32(1301), 
    maxDoca       = cms.untracked.double(0.04), 
    massLow       = cms.untracked.double(2.8), 
    massHigh      = cms.untracked.double(3.3)
    )

# ----------------------------------------------------------------------
mtDump2 = cms.EDAnalyzer(
    "HFMuonAndTrack",
    verbose       = cms.untracked.int32(0), 
    muonsLabel    = cms.untracked.InputTag("muons"),
    tracksLabel   = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt        = cms.untracked.double(6.0),
    trackPt       = cms.untracked.double(3.5),
    type          = cms.untracked.int32(1302), 
    maxDoca       = cms.untracked.double(0.04), 
    massLow       = cms.untracked.double(9.0), 
    massHigh      = cms.untracked.double(10.6)
    )


# ----------------------------------------------------------------------
mmDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt             = cms.untracked.double(3.5),
    type               = cms.untracked.int32(1313), 
    vertexing          = cms.untracked.int32(1), 
    maxDoca            = cms.untracked.double(0.1), 
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
    muonPt             = cms.untracked.double(3.5),
    psiMuons           = cms.untracked.int32(2),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.1),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0)
    )

# ----------------------------------------------------------------------
dstarDump = cms.EDAnalyzer(
    "HFDstar",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    trackPt            = cms.untracked.double(3.5),
    slowPionPt         = cms.untracked.double(0.3),
    D0Window           = cms.untracked.double(0.1),
    deltaM             = cms.untracked.double(0.03),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.2),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0)
    )

# ######################################################################
# Sequences
# ######################################################################
bmtSequence     = cms.Sequence(mmDump*mtDump1*mtDump2*bupsikpDump*dstarDump)
