import os
import FWCore.ParameterSet.Config as cms

skipEvents = cms.EDFilter(
    "HFSkipEvents",
    verbose                      = cms.untracked.int32(1),
    filterOnPrimaryVertex        =  cms.untracked.int32(1),
    PrimaryVertexCollectionLabel = cms.untracked.InputTag('offlinePrimaryVertices'),
    filterOnTrackMaximum         =  cms.untracked.int32(300),
    TrackCollectionLabel         = cms.untracked.InputTag('generalTracks')
    )
