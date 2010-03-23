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
process.GlobalTag.globaltag = "START3X_V21B::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10001.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10002.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10003.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10004.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10005.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10006.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10007.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10008.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10009.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10010.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10011.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10012.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10013.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10014.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10015.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10016.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10017.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10018.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10019.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10020.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10021.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10022.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10023.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10024.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10025.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10026.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10027.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10028.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10029.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10030.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10031.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10032.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10033.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10034.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10035.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10036.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10037.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10038.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10039.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10040.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10041.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10042.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10043.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10044.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10045.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10046.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10047.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10048.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10049.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10050.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10051.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10052.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10053.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10056.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10057.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10058.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10059.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10060.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10061.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10062.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10063.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10064.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10065.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10066.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10067.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10068.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10069.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10070.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10071.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10072.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10073.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10074.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10076.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10077.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10078.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10079.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10080.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10081.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10082.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10083.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10084.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10085.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10086.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10087.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10088.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10089.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10090.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10091.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10092.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10093.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10095.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10096.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10097.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10098.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10099.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10100.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10101.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10102.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10103.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10104.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10105.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10106.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10107.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10108.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10109.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10110.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10111.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10112.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10113.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10114.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10115.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10116.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10117.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10118.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10119.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10120.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10121.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10122.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10123.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10124.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10125.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10126.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10127.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10128.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10129.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10130.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10131.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10132.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10133.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10134.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10135.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10136.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10137.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10138.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10139.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10140.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10141.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10142.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10143.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10144.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10145.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10146.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10147.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10148.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10149.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10150.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10151.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10152.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10153.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10154.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10155.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10156.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10157.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10158.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10159.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10160.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10161.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10162.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10163.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10164.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10165.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10166.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10167.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10168.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10169.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10170.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10171.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10172.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10173.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10174.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10175.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10176.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10177.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10178.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10179.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10180.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10181.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10182.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10183.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10184.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10185.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10186.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10187.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10188.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10189.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10190.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10191.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10192.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10193.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10194.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10195.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10196.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10197.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10198.root",
    "/story/user/ursl/production/Winter10/BuToJPsiKplus/BuToJPsiKplus_7TeV_RAW2DIGI_RECO_START-10199.root"
    )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# ----------------------------------------------------------------------
process.HepPDTESSource = cms.ESSource(
    "HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/particle.tbl')
    )
process.genParticles = cms.EDProducer(
    "GenParticleProducer",
    saveBarCodes          = cms.untracked.bool(True),
    src                   = cms.InputTag("generator"),
    abortOnUnknownPDGCode = cms.untracked.bool(False)
    )


# ----------------------------------------------------------------------
process.genDump = cms.EDAnalyzer(
    "HFDumpGenerator",
    generatorCandidates = cms.untracked.string('genParticles'),
    generatorEvent = cms.untracked.string('generator')
    )


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
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.trkDump = cms.EDAnalyzer(
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
process.stuffDump = cms.EDAnalyzer(
    "HFDumpStuff",
    verbose                  = cms.untracked.int32(0),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    PrimaryVertexTracksLabel = cms.untracked.InputTag("generalTracks")
    )


# ----------------------------------------------------------------------
process.muonDump = cms.EDAnalyzer(
    "HFDumpMuons",
    verbose         = cms.untracked.int32(0),
    muonsLabel      = cms.untracked.InputTag("muons"),
    calomuonsLabel  = cms.untracked.InputTag("calomuons"),
    doTruthMatching = cms.untracked.int32(0),
    )


# ----------------------------------------------------------------------
process.triggerDump = cms.EDAnalyzer(
    "HFDumpTrigger",
    verbose                 = cms.untracked.int32(0),
    L1GTReadoutRecordLabel  = cms.untracked.InputTag("gtDigis"), 
    hltL1GtObjectMap        = cms.untracked.InputTag("hltL1GtObjectMap"), 
    L1MuonsLabel            = cms.untracked.InputTag("hltL1extraParticles"), 
    HLTResultsLabel         = cms.untracked.InputTag("TriggerResults::HLT"), 
    TriggerEventLabel       = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"), 
    hltLabel                = cms.untracked.InputTag("TriggerResults::HLT"), 
    )

# ----------------------------------------------------------------------
process.bmmDump = cms.EDAnalyzer(
    "HFDimuons",
    verbose      = cms.untracked.int32(0), 
    muonsLabel   = cms.untracked.InputTag("muons"),
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    muonPt       = cms.untracked.double(1.0),
    type         = cms.untracked.int32(1313), 
    vertexing    = cms.untracked.int32(0), 
    massLow      = cms.untracked.double(0.5), 
    massHigh     = cms.untracked.double(12.0)
    )


# ----------------------------------------------------------------------
process.truthBsToMuMuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(80),
    GenType      = cms.untracked.int32(-80),
    daughtersID  = cms.untracked.vint32(13, -13)
    )

# ----------------------------------------------------------------------
process.truthBsToMuMuGaDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(81),
    GenType      = cms.untracked.int32(-81),
    daughtersID  = cms.untracked.vint32(13, -13, 22)
    )

# ----------------------------------------------------------------------
process.truthBsToKKDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(82),
    GenType      = cms.untracked.int32(-82),
    daughtersID  = cms.untracked.vint32(321, -321)
    )

# ----------------------------------------------------------------------
process.truthBsToKPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(83),
    GenType      = cms.untracked.int32(-83),
    daughtersID  = cms.untracked.vint32(321, -211)
    )

# ----------------------------------------------------------------------
process.truthBsToPiPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(84),
    GenType      = cms.untracked.int32(-84),
    daughtersID  = cms.untracked.vint32(211, -211)
    )

# ----------------------------------------------------------------------
process.truthBuTo3MuNuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(70),
    GenType      = cms.untracked.int32(-70),
    daughtersID  = cms.untracked.vint32(13, 13, -13, 14)
    )

# ----------------------------------------------------------------------
process.truthBdToMuMuDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(90),
    GenType      = cms.untracked.int32(-90),
    daughtersID  = cms.untracked.vint32(13, -13)
    )

# ----------------------------------------------------------------------
process.truthBdToPiPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(91),
    GenType      = cms.untracked.int32(-91),
    daughtersID  = cms.untracked.vint32(211, -211)
    )

# ----------------------------------------------------------------------
process.truthLambdaBToPPiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(5122),
    type         = cms.untracked.int32(60),
    GenType      = cms.untracked.int32(-60),
    daughtersID  = cms.untracked.vint32(2212, -211)
    )

# ----------------------------------------------------------------------
process.truthLambdaBToPKDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(5122),
    type         = cms.untracked.int32(61),
    GenType      = cms.untracked.int32(-61),
    daughtersID  = cms.untracked.vint32(2212, -321)
    )

# ----------------------------------------------------------------------
process.truthBdToJPsiKsDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(92),
    GenType      = cms.untracked.int32(-92),
    daughtersID  = cms.untracked.vint32(443, 310)
    )

# ----------------------------------------------------------------------
process.truthBdToJPsiKstarDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(93),
    GenType      = cms.untracked.int32(-93),
    daughtersID  = cms.untracked.vint32(443, 313)
    )

# ----------------------------------------------------------------------
process.truthBdToMuMuPi0Dump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(511),
    type         = cms.untracked.int32(94),
    GenType      = cms.untracked.int32(-94),
    daughtersID  = cms.untracked.vint32(13, -13, 111)
    )

# ----------------------------------------------------------------------
process.truthBsToJPsiPhiDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(531),
    type         = cms.untracked.int32(95),
    GenType      = cms.untracked.int32(-95),
    daughtersID  = cms.untracked.vint32(443, 333)
    )

# ----------------------------------------------------------------------
process.truthBuToJPsiKplusDump = cms.EDAnalyzer(
    "HFTruthCandidate",
    verbose      = cms.untracked.int32(0), 
    tracksLabel  = cms.untracked.InputTag('generalTracks'),
    motherID     = cms.untracked.int32(521),
    type         = cms.untracked.int32(96),
    GenType      = cms.untracked.int32(-96),
    daughtersID  = cms.untracked.vint32(443, 321)
    )

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genParticles* 
    process.genDump*
    process.stuffDump*
    process.trkDump*
    process.muonDump*
    process.triggerDump*
    process.bmmDump*
    process.truthBsToMuMuDump*
    process.truthBdToMuMuDump*
    process.truthBsToKKDump*
    process.truthBsToKPiDump*
    process.truthBdToPiPiDump*
    process.truthBsToPiPiDump*
    process.truthLambdaBToPPiDump*
    process.truthLambdaBToPKDump*
    process.truthBuTo3MuNuDump*
    process.truthBsToMuMuGaDump*
    process.truthBdToJPsiKsDump*
    process.truthBdToJPsiKstarDump*
    process.truthBdToMuMuPi0Dump*
    process.truthBsToJPsiPhiDump*
    process.truthBuToJPsiKplusDump*
    process.tree
)




