import os
import FWCore.ParameterSet.Config as cms

skipEvents = cms.EDFilter(
    "HFSkipEvents",
    verbose                      = cms.untracked.int32(0),
    filterOnPrimaryVertex        =  cms.untracked.int32(1),
    PrimaryVertexCollectionLabel = cms.untracked.InputTag('offlinePrimaryVertices'),
    filterOnTrackMaximum         =  cms.untracked.int32(-1),
    filterOnMuonMinimum          =  cms.untracked.int32(-1),
    TrackCollectionLabel         = cms.untracked.InputTag('generalTracks')
    )
