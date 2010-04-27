import os
import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "test-signals.root"

tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(1),
    requireCand  = cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )

