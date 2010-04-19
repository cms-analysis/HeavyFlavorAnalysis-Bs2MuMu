# ######################################################################
# t3ui02
# /shome/naegelic/CMSSW_3_5_6/src/HeavyFlavorAnalysis/Bs2MuMu/test/test-grid
# mkPyFiles -t validation-XXXX.py -n 99999 -f ../Winter10/files/BuToJPsiKplus.list -d .
# ./validation--0000.py with  events
# ######################################################################
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
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START3X_V25::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
    "/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10001.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10002.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10003.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10004.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10005.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10006.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10007.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10008.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10009.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10010.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10011.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10012.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10013.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10014.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10015.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10016.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10017.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10018.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10019.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10020.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10021.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10022.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10023.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10024.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10025.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10026.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10027.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10028.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10029.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10030.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10031.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10032.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10033.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10034.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10035.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10036.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10037.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10038.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10039.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10040.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10041.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10042.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10043.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10044.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10045.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10046.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10047.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10048.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10049.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10050.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10051.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10052.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10053.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10056.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10057.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10058.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10059.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10060.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10061.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10062.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10063.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10064.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10065.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10066.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10067.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10068.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10069.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10070.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10071.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10072.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10073.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10074.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10076.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10077.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10078.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10079.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10080.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10081.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10082.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10083.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10084.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10085.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10086.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10087.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10088.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10089.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10090.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10091.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10092.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10093.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10095.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10096.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10097.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10098.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10099.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10100.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10101.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10102.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10103.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10104.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10105.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10106.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10107.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10108.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10109.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10110.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10111.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10112.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10113.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10114.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10115.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10116.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10117.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10118.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10119.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10120.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10121.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10122.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10123.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10124.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10125.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10126.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10127.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10128.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10129.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10130.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10131.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10132.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10133.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10134.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10135.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10136.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10137.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10138.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10139.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10140.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10141.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10142.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10143.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10144.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10145.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10146.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10147.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10148.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10149.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10150.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10151.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10152.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10153.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10154.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10155.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10156.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10157.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10158.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10159.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10160.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10161.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10162.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10163.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10164.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10165.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10166.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10167.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10168.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10169.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10170.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10171.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10172.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10173.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10174.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10175.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10176.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10177.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10178.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10179.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10180.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10181.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10182.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10183.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10184.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10185.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10186.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10187.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10188.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10189.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10190.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10191.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10192.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10193.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10194.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10195.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10196.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10197.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10198.root",
"/store/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10199.root"
    )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + ".root"
except KeyError:
    rootFileName = "test-BuToJPsiKplus.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose      = cms.untracked.int32(0),
    requireCand  =  cms.untracked.bool(True),
    fileName     = cms.untracked.string(rootFileName)
    )

# ----------------------------------------------------------------------
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFMCTruth_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFRecoStuff_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFTruthCandidates_cff")
process.load("HeavyFlavorAnalysis.Bs2MuMu.HFB2JpsiCandidates_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.MCTruthSequence*
    process.recoStuffSequence*
    process.bupsikpDump*
    process.truthBuDump*
    process.tree
)
