import FWCore.ParameterSet.Config as cms

trackList = "generalTracks"

# ######################################################################
# Bs modes
# ######################################################################

# ----------------------------------------------------------------------
truthBsToMuMuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(80),
    GenType      = cms.untracked.int32(-80),
    daughtersID  = cms.untracked.vint32(13, -13)
    )


# ----------------------------------------------------------------------
truthBsToMuMuGaDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(81),
    GenType      = cms.untracked.int32(-81),
    daughtersID  = cms.untracked.vint32(13, -13, 22)
    )

# ----------------------------------------------------------------------
truthBsToKKDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(82),
    GenType      = cms.untracked.int32(-82),
    daughtersID  = cms.untracked.vint32(321, -321)
    )

# ----------------------------------------------------------------------
truthBsToKPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(83),
    GenType      = cms.untracked.int32(-83),
    daughtersID  = cms.untracked.vint32(321, -211)
    )

# ----------------------------------------------------------------------
truthBsToPiPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(84),
    GenType      = cms.untracked.int32(-84),
    daughtersID  = cms.untracked.vint32(211, -211)
    )


# ######################################################################
# Bd modes
# ######################################################################

# ----------------------------------------------------------------------
truthBdToMuMuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(90),
    GenType      = cms.untracked.int32(-90),
    daughtersID  = cms.untracked.vint32(13, -13)
    )

# ----------------------------------------------------------------------
truthBdToPiPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(91),
    GenType      = cms.untracked.int32(-91),
    daughtersID  = cms.untracked.vint32(211, -211)
    )



# ######################################################################
# Bu modes
# ######################################################################

# ----------------------------------------------------------------------
truthBuTo3MuNuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(70),
    GenType      = cms.untracked.int32(-70),
    daughtersID  = cms.untracked.vint32(13, 13, -13, 14)
    )



# ######################################################################
# LambdaB modes
# ######################################################################

# ----------------------------------------------------------------------
truthLambdaBToPPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(5122),
    type         = cms.untracked.int32(60),
    GenType      = cms.untracked.int32(-60),
    daughtersID  = cms.untracked.vint32(2212, -211)
    )

# ----------------------------------------------------------------------
truthLambdaBToPKDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(5122),
    type         = cms.untracked.int32(61),
    GenType      = cms.untracked.int32(-61),
    daughtersID  = cms.untracked.vint32(2212, -321)
    )



# ######################################################################
# B -> J/psi X
# ######################################################################

# ----------------------------------------------------------------------
truthBsDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(67),
    GenType      = cms.untracked.int32(-67),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 333, 321, -321)
    )

# ----------------------------------------------------------------------
truthBuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(68),
    GenType      = cms.untracked.int32(-68),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 321)
    )

# ----------------------------------------------------------------------
truthBd2JpsiKstarDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(69),
    GenType      = cms.untracked.int32(-69),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 313, 321, 211)
    )

# ----------------------------------------------------------------------
truthBd2JpsiKsDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(66),
    GenType      = cms.untracked.int32(-66),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 310, 211, 211)
    )


# ######################################################################
# Sequences
# ######################################################################

truthSignalsSequence     = cms.Sequence(truthBsToMuMuDump*truthBdToMuMuDump)
truthRareBsSequence      = cms.Sequence(truthBsToMuMuGaDump*truthBsToKKDump*truthBsToKPiDump*truthBsToPiPiDump)
truthRareBdSequence      = cms.Sequence(truthBdToPiPiDump)
truthRareBuSequence      = cms.Sequence(truthBuTo3MuNuDump)
truthRareLambdaBSequence = cms.Sequence(truthLambdaBToPPiDump*truthLambdaBToPKDump)
truthB2JpsiSequence      = cms.Sequence(truthBsDump*truthBd2JpsiKsDump*truthBd2JpsiKstarDump*truthBuDump)
truthAllSequence         = cms.Sequence(truthSignalsSequence*truthRareBsSequence*truthRareBdSequence*truthRareBuSequence*truthRareLambdaBSequence*truthB2JpsiSequence)
