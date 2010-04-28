#! /bin/csh -f

# tree.csh
#   Branched by Christoph Naegeli <naegelic@phys.ethz.ch> on 31.3.10.
#   Usage: Use in conjunction with 'run'.
#
#   TREEREADER should be the name of the readerclass to run as
#   can be found in runTreeReaders.cc and is supposed to be set
#   using the replacement list of the 'run' command.
#   JOB is the chain file on which the TreeReader should operate.

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
echo "--> Building HeavyFlavorObjects"
cd AnalysisDataFormats/HeavyFlavorObjects
make 
cd - 
echo "--> Building TreeReader"
cd HeavyFlavorAnalysis/Bs2MuMu/macros
make
mkdir -p chains
mv ../../../../../$JOB chains/
cd -

# ----------------------------------------------------------------------
# -- Run Treereader
# ----------------------------------------------------------------------

echo "--> Run Treereader"

date
cd HeavyFlavorAnalysis/Bs2MuMu/macros
pwd

echo "$EXECUTABLE -c chains/$JOB  -o $FILE1 |& tee $JOB.log"
$EXECUTABLE -c chains/$JOB -o $FILE1 |& tee $JOB.log
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

echo "run: This is the end, my friend"
