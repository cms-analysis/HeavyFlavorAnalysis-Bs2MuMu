#! /bin/csh -f

setenv CMSSW       
setenv SCRAM_ARCH  
setenv SRMCP       

setenv JOB
setenv FILE1    $JOB.root
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

# BATCH START

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
cmsRun $JOB.py |& tee $JOB.log
date

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

date

# BATCH END

echo "runGrid: This is the end, my friend"
