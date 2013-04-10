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


# ----------------------------------------------------------------------
truthBsToPiMuNuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(85),
    GenType      = cms.untracked.int32(-85),
    daughtersID  = cms.untracked.vint32(211, -13, 14)
    )

# ----------------------------------------------------------------------
truthBsToKMuNuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(86),
    GenType      = cms.untracked.int32(-86),
    daughtersID  = cms.untracked.vint32(321, -13, 14)
    )

# ----------------------------------------------------------------------
truthBsToPhiMuMuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(87),
    GenType      = cms.untracked.int32(-87),
    daughtersID  = cms.untracked.vint32(333, 13, -13),
    partialDecayMatching = cms.untracked.bool(True)
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

# ----------------------------------------------------------------------
truthBdToKPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(92),
    GenType      = cms.untracked.int32(-92),
    daughtersID  = cms.untracked.vint32(321, -211)
    )

# ----------------------------------------------------------------------
truthBdToKKDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(93),
    GenType      = cms.untracked.int32(-93),
    daughtersID  = cms.untracked.vint32(321, -321)
    )

# ----------------------------------------------------------------------
truthBdToMuMuPi0Dump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(94),
    GenType      = cms.untracked.int32(-94),
    daughtersID  = cms.untracked.vint32(13, -13, 111)
    )

# ----------------------------------------------------------------------
truthBdToPiMuNuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(95),
    GenType      = cms.untracked.int32(-95),
    daughtersID  = cms.untracked.vint32(211, -13, 14)
    )

# ----------------------------------------------------------------------
truthBdToMuMuGaDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(96),
    GenType      = cms.untracked.int32(-96),
    daughtersID  = cms.untracked.vint32(13, -13, 22)
    )

# ----------------------------------------------------------------------
truthBdToMuMuK0Dump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(97),
    GenType      = cms.untracked.int32(-97),
#    daughtersID  = cms.untracked.vint32(13, -13, 311),
    daughtersID  = cms.untracked.vint32(13, -13),
    partialDecayMatching = cms.untracked.bool(True)
    )
    
# ----------------------------------0000------------------------------------
truthBdToRhoPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(98),
    GenType      = cms.untracked.int32(-98),
#    daughtersID  = cms.untracked.vint32(213, 211),
    daughtersID  = cms.untracked.vint32(213, 211, 111, 211),
    partialDecayMatching = cms.untracked.bool(True)
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

# ----------------------------------------------------------------------
truthBuToMuMuKDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(77),
    GenType      = cms.untracked.int32(-77),
    daughtersID  = cms.untracked.vint32(13, -13, 321),
    partialDecayMatching = cms.untracked.bool(True)
    )

# ----------------------------------------------------------------------
truthBuToMuMuPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(78),
    GenType      = cms.untracked.int32(-78),
    daughtersID  = cms.untracked.vint32(13, -13),
    partialDecayMatching = cms.untracked.bool(True)
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

# ----------------------------------------------------------------------
truthLambdaBToKMuNuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(5122),
    type         = cms.untracked.int32(62),
    GenType      = cms.untracked.int32(-62),
    daughtersID  = cms.untracked.vint32(2212, -13, 14)
    )


# ######################################################################
# B -> J/psi X
# ######################################################################

# ----------------------------------------------------------------------
truthBu2JpsiPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(66),
    GenType      = cms.untracked.int32(-66),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 211)
    )

# ----------------------------------------------------------------------
truthBs2JpsiPhiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(67),
    GenType      = cms.untracked.int32(-67),
    daughtersID  = cms.untracked.vint32(443, 333, 13, -13, 321, -321)
    )

# ----------------------------------------------------------------------
truthBu2JpsiKpDump = cms.EDAnalyzer(
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
truthBd2JpsiKstarAsBpDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(10),
    GenType      = cms.untracked.int32(-10),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 313, 321, 211)
    )

# ----------------------------------------------------------------------
truthBd2JpsiKstarAsBsDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(11),
    GenType      = cms.untracked.int32(-11),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 313, 321, 211)
    )

# ----------------------------------------------------------------------
truthBd2JpsiKsDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(64),
    GenType      = cms.untracked.int32(-64),
    daughtersID  = cms.untracked.vint32(443, 13, -13, 310, 211, 211)
    )

# ----------------------------------------------------------------------
truthBd2DstarPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0),
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(30),
    GenType      = cms.untracked.int32(-30),
    daughtersID  = cms.untracked.vint32(-413, 211, -421, -211, 321, 211)
    )


# ######################################################################
# Charm
# ######################################################################

# ----------------------------------------------------------------------
truthD0ToKPi = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(421),
    type         = cms.untracked.int32(50),
    GenType      = cms.untracked.int32(-50),
    daughtersID  = cms.untracked.vint32(-321, 211)
    )

# ----------------------------------------------------------------------
truthDpToKPiPi = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(411),
    type         = cms.untracked.int32(51),
    GenType      = cms.untracked.int32(-51),
    daughtersID  = cms.untracked.vint32(-321, 211, 211)
    )

# ----------------------------------------------------------------------
truthDpToKstarPi = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(411),
    type         = cms.untracked.int32(52),
    GenType      = cms.untracked.int32(-52),
    daughtersID  = cms.untracked.vint32(313, 321, -211, 211)
    )

# ----------------------------------------------------------------------
truthDsToPhiPi = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(431),
    type         = cms.untracked.int32(53),
    GenType      = cms.untracked.int32(-53),
    daughtersID  = cms.untracked.vint32(333, 321, -321, 211)
    )

# ----------------------------------------------------------------------
truthDstarToD0PiToKPiPi = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(413),
    type         = cms.untracked.int32(54),
    GenType      = cms.untracked.int32(-54),
    daughtersID  = cms.untracked.vint32(421, 211, -321, -211)
    )

# ----------------------------------------------------------------------
truthDpToKKPi = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(411),
    type         = cms.untracked.int32(55),
    GenType      = cms.untracked.int32(-55),
    daughtersID  = cms.untracked.vint32(321, -321, 211)
    )

# ----------------------------------------------------------------------
truthLambdaCToPrKPi = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(4122),
    type         = cms.untracked.int32(56),
    GenType      = cms.untracked.int32(-56),
    daughtersID  = cms.untracked.vint32(2212, -321, 211)
    )


# ######################################################################
# Charmonium
# ######################################################################

# ----------------------------------------------------------------------
truthPsiToMuMu = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(443),
    type         = cms.untracked.int32(40),
    GenType      = cms.untracked.int32(-40),
    daughtersID  = cms.untracked.vint32(13, -13)
    )

# ----------------------------------------------------------------------
truthPsi2SToMuMu = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(100443),
    type         = cms.untracked.int32(41),
    GenType      = cms.untracked.int32(-41),
    daughtersID  = cms.untracked.vint32(13, -13)
    )



# ######################################################################
# Bottomium
# ######################################################################

# ----------------------------------------------------------------------
truthUps1SToMuMu = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(553),
    type         = cms.untracked.int32(45),
    GenType      = cms.untracked.int32(-45),
    daughtersID  = cms.untracked.vint32(13, -13)
    )

# ----------------------------------------------------------------------
truthUps2SToMuMu = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(100553),
    type         = cms.untracked.int32(46),
    GenType      = cms.untracked.int32(-46),
    daughtersID  = cms.untracked.vint32(13, -13)
    )

# ----------------------------------------------------------------------
truthUps3SToMuMu = cms.EDAnalyzer(
    "HFTruthCandidate",
    tracksLabel  = cms.untracked.InputTag(trackList),
    motherID     = cms.untracked.int32(200553),
    type         = cms.untracked.int32(47),
    GenType      = cms.untracked.int32(-47),
    daughtersID  = cms.untracked.vint32(13, -13)
    )



# ######################################################################
# Sequences
# ######################################################################

truthSignalsSequence     = cms.Sequence(truthBsToMuMuDump*truthBdToMuMuDump)

truthRareBsSequence      = cms.Sequence(truthBsToMuMuGaDump
                                        *truthBsToKKDump
                                        *truthBsToKPiDump
                                        *truthBsToPiPiDump
                                        *truthBsToPiMuNuDump
                                        *truthBsToKMuNuDump
                                        *truthBsToPhiMuMuDump
                                        )

truthRareBdSequence      = cms.Sequence(truthBdToPiPiDump
                                        *truthBdToKPiDump
                                        *truthBdToKKDump
                                        *truthBdToMuMuPi0Dump
                                        *truthBdToPiMuNuDump
                                        *truthBdToMuMuGaDump
                                        *truthBdToMuMuK0Dump
                                        *truthBdToRhoPiDump
                                        )

truthRareBuSequence      = cms.Sequence(truthBuTo3MuNuDump
                                        *truthBuToMuMuKDump
                                        *truthBuToMuMuPiDump
                                        )

truthRareLambdaBSequence = cms.Sequence(truthLambdaBToPPiDump
                                        *truthLambdaBToPKDump
                                        *truthLambdaBToKMuNuDump)

truthB2DstarSequence     = cms.Sequence(truthBd2DstarPiDump)

truthB2JpsiSequence      = cms.Sequence(truthBs2JpsiPhiDump
                                        *truthBd2JpsiKsDump
                                        *truthBd2JpsiKstarDump
                                        *truthBd2JpsiKstarAsBpDump
                                        *truthBd2JpsiKstarAsBsDump
                                        *truthBu2JpsiKpDump
                                        *truthBu2JpsiPiDump)

truthOniaSequence        = cms.Sequence(truthPsiToMuMu
                                        *truthPsi2SToMuMu
                                        *truthUps1SToMuMu
                                        *truthUps2SToMuMu
                                        *truthUps3SToMuMu)

truthCharmSequence       = cms.Sequence(truthD0ToKPi
                                        *truthDpToKPiPi
                                        *truthDpToKstarPi
                                        *truthDsToPhiPi
                                        *truthDstarToD0PiToKPiPi
                                        *truthDpToKKPi
                                        *truthLambdaCToPrKPi)

truthBmmSequence         = cms.Sequence(truthSignalsSequence
                                        *truthB2JpsiSequence
                                        *truthRareBsSequence
                                        *truthRareBdSequence
                                        *truthRareBuSequence
                                        *truthRareLambdaBSequence
                                        )

truthBmtSequence         = cms.Sequence(truthSignalsSequence
                                        *truthB2JpsiSequence
                                        *truthD0ToKPi
                                        *truthDstarToD0PiToKPiPi
                                        *truthOniaSequence
                                        )

truthAllSequence         = cms.Sequence(truthSignalsSequence
                                        *truthRareBsSequence
                                        *truthRareBdSequence
                                        *truthRareBuSequence
                                        *truthRareLambdaBSequence
                                        *truthB2JpsiSequence
                                        *truthOniaSequence
                                        *truthCharmSequence
                                        *truthB2DstarSequence)
