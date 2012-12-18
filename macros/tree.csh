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

# several files possible if tree grows to big. copy them all
echo "--> Save output to SE: $STORAGE1/$FILE1"

set FILES=`ls $JOB*.root`
echo "Found the following output root files: $FILES"
foreach f ($FILES)
    echo lcg-del -b -D srmv2 -l "$STORAGE1/$f"
    lcg-del -b -D srmv2 -l "$STORAGE1/$f"
    echo lcg-cp -b -D srmv2 $f "$STORAGE1/$f"
    lcg-cp -b -D srmv2 $f "$STORAGE1/$f"
    echo lcg-ls -b -D srmv2 "$STORAGE1/$f"
    lcg-ls -b -D srmv2 "$STORAGE1/$f"
end

# copy the log file.
echo lcg-del -b -D srmv2 -l "$STORAGE1/$JOB.log"
lcg-del -b -D srmv2 -l "$STORAGE1/$JOB.log"
echo lcg-cp -b -D srmv2 $JOB.log "$STORAGE1/$JOB.log"
lcg-cp -b -D srmv2 $JOB.log "$STORAGE1/$JOB.log"
echo lcg-ls -b -D srmv2 "$STORAGE1/$JOB.log"
lcg-ls -b -D srmv2 "$STORAGE1/$JOB.log"

date


# BATCH END

echo "run: This is the end, my friend"
