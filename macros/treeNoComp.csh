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
pwd

echo "--> Extract tar file"
date
tar zxf ../../$JOB.tar.gz
####echo "--> Building HeavyFlavorObjects"
####cd AnalysisDataFormats/HeavyFlavorObjects
####make 
####cd - 
####echo "--> Building TreeReader"
cd HeavyFlavorAnalysis/Bs2MuMu/macros2
####make
mkdir -p chains
mv ../../../../../$JOB chains/
cd -

# ----------------------------------------------------------------------
# -- Run Treereader
# ----------------------------------------------------------------------

echo "--> Run Treereader"

date
cd HeavyFlavorAnalysis/Bs2MuMu/macros2
pwd

echo "$EXECUTABLE -c chains/$JOB  -o $FILE1 |& tee $JOB.log"
$EXECUTABLE -c chains/$JOB -o $FILE1 |& tee $JOB.log
date

# ----------------------------------------------------------------------
# -- Save Output to NFS, not the SE
# ----------------------------------------------------------------------

# several files possible if tree grows to big. copy them all
echo "--> Save output to SE: $STORAGE1/$FILE1"
echo $SRMCP 

set FILES=`ls $JOB*.root`
echo "Found the following output root files: $FILES"
foreach f ($FILES)
  echo lcg-del "$STORAGE1/$f"
  lcg-del -b -D srmv2 -l "$STORAGE1/$f"
  echo lcg-cp    file:///`pwd`/$f "$STORAGE1/$f"
  lcg-cp -b -D srmv2  file:///`pwd`/$f "$STORAGE1/$f"
  echo lcg-ls     "$STORAGE1/$f"
  lcg-ls -b -D srmv2  "$STORAGE1/$f"
end

# copy the log file.
echo lcg-del  "$STORAGE1/$JOB.log"
lcg-del -b -D srmv2 -l "$STORAGE1/$JOB.log"
echo lcg-cp file:///`pwd`/$JOB.log "$STORAGE1/$JOB.log"
lcg-cp -b -D srmv2     file:///`pwd`/$JOB.log "$STORAGE1/$JOB.log"
echo lcg-ls     "$STORAGE1/$JOB.log"
lcg-ls -b -D srmv2  "$STORAGE1/$JOB.log"

date

# BATCH END


# -- cleanup
#/bin/rm -rf /scratch/ursl/dana*

echo "run: This is the end, my friend"
