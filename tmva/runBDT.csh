#! /bin/csh -f

setenv CMSSW
setenv SCRAM_ARCH
setenv SRMCP

setenv JOB
setenv EXECUTABLE
setenv FILE1 $JOB.root
setenv BATCH
setenv CHANNEL
setenv STORAGE1 srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/ursl/tmva/2011barrel/$BATCH

setenv VERSION 053200

echo "========================"
echo "====> GRID wrapper <===="
echo "========================"

echo "--> Running grid job wrapper"

# ----------------------------------------------------------------------
# -- The Basics
# ----------------------------------------------------------------------
echo "--> Environment:"
date
hostname
uname -a
echo "--> df -kl:"
df -kl

echo $VO_CMS_SW_DIR
ls -l $VO_CMS_SW_DIR
source $VO_CMS_SW_DIR/cmsset_default.csh

printenv
pwd
echo "--> End of env testing"

# BATCH START

# ----------------------------------------------------------------------
# -- Setup TreeReader
# ----------------------------------------------------------------------
echo "--> Setup TreeReader"
pwd
date
cmsrel $CMSSW
cd $CMSSW/src
eval `scramv1 runtime -csh`
setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:${LD_LIBRARY_PATH}
# -- use special root version
setenv ROOTSYS /shome/naegelic/root
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
setenv PATH ${ROOTSYS}/bin:${PATH}
echo "--> Unsetting SCRAM_ARCH"
unsetenv SCRAM_ARCH
pwd

echo "--> Extract tar file"
date
tar zxf ../../$JOB.tar.gz
cd HeavyFlavorAnalysis/Bs2MuMu/tmva
pwd

date
echo "--> running $EXECUTABLE "
$EXECUTABLE
echo "--> done with trainBDT"
date
ls -rtl 
cp plots/TMVA-$JOB*  /shome/ursl/ana/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/tmva/tmp-$JOB/
cp TMVA-$JOB.log      /shome/ursl/ana/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/tmva/tmp-$JOB/
date

echo cp TMVA-$JOB.log "/shome/ursl/ana/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/tmva/tmp-$JOB/TMVA-$JOB.log"
cp TMVA-$JOB.log /shome/ursl/ana/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/tmva/tmp-$JOB/TMVA-$JOB.log
cp TMVA-$JOB.tex /shome/ursl/ana/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/tmva/tmp-$JOB/TMVA-$JOB.tex

# BATCH END


# -- cleanup

echo "run: This is the end, my friend"
