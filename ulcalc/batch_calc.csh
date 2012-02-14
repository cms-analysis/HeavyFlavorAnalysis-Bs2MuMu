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
setenv ROOTSYS /swshare/ROOT/root_v32.00.patches42593_slc5_amd64
setenv LD_LIBRARY_PATH /swshare/ROOT/root_v32.00.patches42593_slc5_amd64//lib/root:$LD_LIBRARY_PATH

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

echo "$EXECUTABLE -w $JOB.root -o $JOB.log $JOB"
$EXECUTABLE -w $JOB.root -o $JOB.log $JOB
cp $JOB.root $STORAGE1
cp $JOB.log $STORAGE1

# BATCH END

echo "run: This is the end, my friend"
