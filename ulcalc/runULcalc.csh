#! /bin/csh -f

setenv CMSSW
setenv SCRAM_ARCH
setenv SRMCP

setenv JOB
setenv EXECUTABLE
setenv FILE1 $JOB.root
setenv STORAGE1

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
cat /proc/cpuinfo
uname -a
echo "--> df -kl:"
df -kl

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
# -- Setup TreeReader
# ----------------------------------------------------------------------
echo "--> Setup TreeReader"
pwd
date
cmsrel $CMSSW
cd $CMSSW/src
eval `scramv1 runtime -csh`
setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:${LD_LIBRARY_PATH}
# -- use CN's root version
setenv ROOTSYS /shome/naegelic/root/
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
setenv PATH ${ROOTSYS}/bin:${PATH}
pwd

echo "--> Extract tar file"
date
tar zxf ../../$JOB.tar.gz
####echo "--> Building HeavyFlavorObjects"
####cd AnalysisDataFormats/HeavyFlavorObjects
####make 
####cd - 
####echo "--> Building TreeReader"
cd HeavyFlavorAnalysis/Bs2MuMu/ulcalc
####make
pwd

bin/ulcalc -a "clb" -l 0.95 -o clb.$JOB $JOB
cp clb.$JOB /shome/ursl/prod/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Bs2MuMu/ulcalc/output

bin/ulcalc -a "bayes" -l 0.95 --disable-errors -o bayes.$JOB $JOB
cp bayes.$JOB /shome/ursl/prod/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Bs2MuMu/ulcalc/output

bin/ulcalc -a "cls" -l 0.95 -o cls.$JOB $JOB 
cp cls.$JOB /shome/ursl/prod/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Bs2MuMu/ulcalc/output

#bin/ulcalc -a "fc" -l 0.68 -r 0.5,4 -n 101 -o fc68.$JOB $JOB
#cp fc68.$JOB /shome/ursl/prod/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Bs2MuMu/ulcalc/output

#bin/ulcalc -a "fc" -l 0.95 -r 0.5,5 -n 101 -o fc95.$JOB $JOB
#cp fc95.$JOB /shome/ursl/prod/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Bs2MuMu/ulcalc/output

# BATCH END


# -- cleanup

echo "run: This is the end, my friend"
