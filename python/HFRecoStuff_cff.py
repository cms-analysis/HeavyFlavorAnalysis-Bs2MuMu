import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
stuffDump = cms.EDAnalyzer(
    "HFDumpStuff",
    verbose                  = cms.untracked.int32(0),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    PrimaryVertexTracksLabel = cms.untracked.InputTag("generalTracks")
    )


# ----------------------------------------------------------------------
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *

trkDump = cms.EDAnalyzer(
    "HFDumpTracks",
    verbose                        = cms.untracked.int32(0),
    tracksLabel                    = cms.untracked.InputTag('generalTracks'),
    primaryVertexCollectionLabel   = cms.untracked.InputTag('offlinePrimaryVertices'),
    generatorEventLabel            = cms.untracked.InputTag('generator'),
    muonsLabel                     = cms.untracked.InputTag("muons"),
    calomuonsLabel                 = cms.untracked.InputTag("calomuons"),
    trackingParticlesLabel         = cms.untracked.InputTag('trackingParticles'),
    associatorLabel                = cms.untracked.InputTag('TrackAssociatorByHits'),
    doTruthMatching                = cms.untracked.int32(3),
    simTracksLabel                 = cms.untracked.InputTag('allLayer1TrackCands')
    )

# ----------------------------------------------------------------------
muonDump = cms.EDAnalyzer(
    "HFDumpMuons",
    verbose         = cms.untracked.int32(0),
    muonsLabel      = cms.untracked.InputTag("muons"),
    calomuonsLabel  = cms.untracked.InputTag("calomuons"),
    doTruthMatching = cms.untracked.int32(0),
    )


# ----------------------------------------------------------------------
triggerDump = cms.EDAnalyzer(
    "HFDumpTrigger",
    verbose                 = cms.untracked.int32(0),
    HLTProcessName          = cms.untracked.string('HLT'), 
    L1GTReadoutRecordLabel  = cms.untracked.InputTag("gtDigis"), 
    hltL1GtObjectMap        = cms.untracked.InputTag("hltL1GtObjectMap"), 
    L1MuonsLabel            = cms.untracked.InputTag("hltL1extraParticles"), 
    HLTResultsLabel         = cms.untracked.InputTag("TriggerResults::HLT"), 
    TriggerEventLabel       = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"), 
    )


# ----------------------------------------------------------------------
hltrep = cms.EDAnalyzer(
    "HLTrigReport",
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
    )


l1trep = cms.EDAnalyzer(
    "L1GtTrigReport",
    UseL1GlobalTriggerRecord = cms.bool( False ),
    L1GtRecordInputTag = cms.InputTag( "gtDigis" )
    )


# ######################################################################
# Sequences
# ######################################################################
recoStuffSequence     = cms.Sequence(stuffDump*trkDump*muonDump*triggerDump*hltrep*l1trep)

