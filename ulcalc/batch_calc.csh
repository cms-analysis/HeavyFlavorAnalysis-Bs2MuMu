#! /bin/csh -f

setenv CMSSW
setenv SCRAM_ARCH
setenv SRMCP

setenv OUTPUT
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

# ----------------------------------------------------------------------
echo "--> Setup ulcalc"
# -- use CN's root version
setenv ROOTSYS /shome/naegelic/root
setenv LD_LIBRARY_PATH /shome/naegelic/root/lib:$LD_LIBRARY_PATH

echo "--> Extract tar file"
date
tar zxf $JOB.tar.gz
####echo "--> Building HeavyFlavorObjects"
####cd AnalysisDataFormats/HeavyFlavorObjects
####make 
####cd - 
####echo "--> Building TreeReader"
mv $JOB ulcalc
cd ulcalc
ls
echo "--> Building ulcalc"
make
pwd

echo "$EXECUTABLE -w $JOB.root $JOB |& tee $JOB.log"
$EXECUTABLE -w $JOB.root $JOB |& tee $JOB.log
echo "lcg-del -b -D srmv2 -l $STORAGE1/$JOB.root"
lcg-del -b -D srmv2 -l "$STORAGE1/$JOB.root"
echo "lcg-cp -b -D srmv2 $JOB.root $STORAGE1/$JOB.root"
lcg-cp -b -D srmv2 $JOB.root "$STORAGE1/$JOB.root"
echo "lcg-ls -b -D srmv2 $STORAGE1/$JOB.root"
lcg-ls -b -D srmv2 "$STORAGE1/$JOB.root"

echo "lcg-del -b -D srmv2 -l $STORAGE1/$JOB.log"
lcg-del -b -D srmv2 -l "$STORAGE1/$JOB.log"
echo "lcg-cp -b -D srmv2 $JOB.log $STORAGE1/$JOB.log"
lcg-cp -b -D srmv2 $JOB.log "$STORAGE1/$JOB.log"
echo "lcg-ls -b -D srmv2 $STORAGE1/$JOB.log"
lcg-ls -b -D srmv2 "$STORAGE1/$JOB.log"

# BATCH END

echo "run: This is the end, my friend"
