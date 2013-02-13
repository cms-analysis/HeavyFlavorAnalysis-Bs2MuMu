#!/usr/bin/python

import subprocess
import sys
import os

CMSRUNFILE='''
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


# ----------------------------------------------------------------------
# -- Database configuration
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "%s::All"

# ----------------------------------------------------------------------
# POOLSOURCE


# ----------------------------------------------------------------------
rootFileName = "bmm-%s-%sXXXX.root"

process.tree = cms.EDAnalyzer(
"HFTree",
verbose      = cms.untracked.int32(0),
requireCand  =  cms.untracked.bool(True),
fileName     = cms.untracked.string(rootFileName)
)

# ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFBmm_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
process.genDump*
process.recoStuffSequence*
process.bmmDump*
process.truthSignalsSequence*
process.tree
)
'''


WORKDIR="submitting"
RUNDIR="run"
CWD=os.getcwd()

# load the configuration file
print "opening config file '%s'" % sys.argv[1]
f = open(sys.argv[1])
while True:
    line = f.readline()
    if line=="":
        break
    if line[0] == '#':
        continue

    # process this configuration
    line = line.split();

    folder=line[0]
    name=line[1]
    globaltag=line[2]
    dataset=line[3]

    os.chdir(CWD)
    
    # remove the subdirectory
    CMD = "rm -rf " + WORKDIR + "/" + folder + "/" + name
    print CMD
    os.system(CMD)
    
    # create the subdirectory
    CMD = "mkdir -p " + WORKDIR + "/" + folder + "/" + name
    print CMD
    os.system(CMD)
    os.chdir(CWD + "/" + WORKDIR + "/" + folder + "/" + name)

    # get the files for the dataset
    CMD = "dbs search --query='find file,file.numevents where dataset=%s' > %s" % (dataset,name)
    print CMD
    os.system(CMD)

    # create the run directory
    CMD = "mkdir -p " + RUNDIR
    print CMD
    os.system(CMD)

    template = "bmm-%s-%sXXXX.py" % (name,globaltag)
    cmsPy = open(template,"w")
    cmsPy.write( CMSRUNFILE % (globaltag,name,globaltag) )
    cmsPy.close()
    CMD = "mkPyFiles -t %s -e 60000 -f %s -d %s" % (template, name, RUNDIR)
    print CMD
    os.system(CMD)

    CMD = "ln -s %s" % (CWD + "/grid.tar.gz")
    print CMD
    os.system(CMD)

    # create the SE entry
    CMD = "srmmkdir srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/naegelic/ntuples/%s" % name
    print CMD
    os.system(CMD)

    os.chdir(RUNDIR)
    CMD = "run -q all.q -c %s -t ../grid.tar.gz -m batch -r 'STORAGE1 srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/naegelic/ntuples/%s' bmm-*.py" % (CWD + "/shell/prod.csh", name)
    print CMD
    os.system(CMD)

f.close()
