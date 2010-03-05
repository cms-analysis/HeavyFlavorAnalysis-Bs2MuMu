#! /bin/csh -f

setenv CMSSW       CMSSW_3_5_0
setenv SCRAM_ARCH  slc5_ia32_gcc434
setenv JOB         

setenv SRMCP       "srmcp -2 -globus_tcp_port_range 20000,25000"
setenv SRMCP       "srmcp -2"
setenv FILE1       $JOB.root
#setenv STORAGE1    srm://storage01.lcg.cscs.ch:8443/srm/managerv2\?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/ursl/production/Winter10/BdToMuMu_7TeV
setenv STORAGE1

echo "========================"
echo "====> GRID wrapper <===="
echo "========================"

echo "--> Running grid job wrapper"

# ----------------------------------------------------------------------
# -- The Basics
# ----------------------------------------------------------------------
echo "--> Environment"
date
hostname
cat /proc/cpuinfo
uname -a

echo $VO_CMS_SW_DIR
ls -l $VO_CMS_SW_DIR
source $VO_CMS_SW_DIR/cmsset_default.csh
echo "-> which edg-gridftp-ls"
which edg-gridftp-ls
echo "-> which globus-url-copy"
which globus-url-copy
echo "-> which srmcp"
which srmcp

pwd
echo "--> End of env testing"

# BATCH

# ----------------------------------------------------------------------
# -- Setup CMSSW
# ----------------------------------------------------------------------
echo "--> Setup CMSSW"
pwd
date
cmsrel $CMSSW
cd $CMSSW/src
eval `scramv1 runtime -csh`
pwd

echo "--> Extract tar file"
date
tar zxf ../../$JOB.tar.gz
cd AnalysisDataFormats/HeavyFlavorObjects
make 
cd - 
scramv1 b -j4
mv ../../$JOB.py .


# ----------------------------------------------------------------------
# -- Run cmsRun
# ----------------------------------------------------------------------
echo "--> Run cmsRun"
pwd
date
which cmsRun
echo "cmsRun $JOB.py "
cmsRun $JOB.py >& $JOB.log 
date
cat $JOB.log 

# ----------------------------------------------------------------------
# -- Save Output to SE
# ----------------------------------------------------------------------
echo "--> Save output to SE: $STORAGE1/$FILE1"

echo srmrm     "$STORAGE1/$FILE1"
srmrm          "$STORAGE1/$FILE1"
echo $SRMCP    file:///`pwd`/$FILE1 "$STORAGE1/$FILE1"
$SRMCP         file:///`pwd`/$FILE1 "$STORAGE1/$FILE1"
echo srmls     "$STORAGE1/$FILE1"
srmls          "$STORAGE1/$FILE1"

echo srmrm  "$STORAGE1/$JOB.log"
srmrm       "$STORAGE1/$JOB.log"
echo $SRMCP file:///`pwd`/$JOB.log "$STORAGE1/$JOB.log"
$SRMCP      file:///`pwd`/$JOB.log "$STORAGE1/$JOB.log"
echo        srmls     "$STORAGE1/$JOB.log"
srmls       "$STORAGE1/$JOB.log"

echo "runGrid: This is the end, my friend"
