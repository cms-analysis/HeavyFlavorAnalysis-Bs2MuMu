#! /bin/csh -f

setenv CMSSW
setenv SCRAM_ARCH
setenv SRMCP

setenv JOB
setenv EXECUTABLE
setenv FILE1 $JOB.root
setenv BATCH
setenv CHANNEL
setenv STORAGE1 srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/ursl/tmva/$BATCH

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
setenv ROOTSYS /swshare/ROOT/root_v32.00.patches42593_slc5_amd64
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
setenv PATH ${ROOTSYS}/bin:${PATH}
pwd

echo "--> Extract tar file"
date
tar zxf ../../$JOB.tar.gz
cd HeavyFlavorAnalysis/Bs2MuMu/tmva
pwd
ls -rtl 

date
echo "--> run runBDT.C"
root -b -q "runBDT.C($JOB, $CHANNEL); > TMVA-$JOB.log"
echo "--> done with runBDT.C"
date
ls -rtl 
cp plots/TMVA-$JOB*  /shome/ursl/tmva/${BATCH}/plots
cp TMVA-$JOB.log  /shome/ursl/tmva/${BATCH}/log
date

echo cp TMVA-$JOB.log "/shome/ursl/prod/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Bs2MuMu/macros2/${BATCH}/tmp-$JOB/TMVA-$JOB.log"
cp TMVA-$JOB.log /shome/ursl/prod/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Bs2MuMu/macros2/${BATCH}/tmp-$JOB/TMVA-$JOB.log

srmrm       "$STORAGE1/root/TMVA-even-$JOB.root"
$SRMCP      file:///`pwd`/TMVA-even-$JOB.root "$STORAGE1/root/TMVA-even-$JOB.root"
srmls       "$STORAGE1/root/TMVA-even-$JOB.root"

srmrm       "$STORAGE1/root/TMVA-odd-$JOB.root"
$SRMCP      file:///`pwd`/TMVA-odd-$JOB.root "$STORAGE1/root/TMVA-odd-$JOB.root"
srmls       "$STORAGE1/root/TMVA-odd-$JOB.root"

srmrm       "$STORAGE1/xml/TMVA-even-${JOB}_BDT.weights.xml"
$SRMCP      file:///`pwd`/weights/TMVA-even-${JOB}_BDT.weights.xml "$STORAGE1/xml/TMVA-even-${JOB}_BDT.weights.xml"
srmls       "$STORAGE1/xml/TMVA-even-${JOB}_BDT.weights.xml"

srmrm       "$STORAGE1/xml/TMVA-odd-${JOB}_BDT.weights.xml"
$SRMCP      file:///`pwd`/weights/TMVA-odd-${JOB}_BDT.weights.xml "$STORAGE1/xml/TMVA-odd-${JOB}_BDT.weights.xml"
srmls       "$STORAGE1/xml/TMVA-odd-${JOB}_BDT.weights.xml"

# BATCH END


# -- cleanup

echo "run: This is the end, my friend"
