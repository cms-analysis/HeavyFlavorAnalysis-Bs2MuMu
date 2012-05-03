#! /bin/tcsh -f

setenv CMSSW
setenv SCRAM_ARCH
setenv SRMCP

setenv OUTPUT
setenv JOB
setenv EXECUTABLE
setenv FILE1 $JOB.root
setenv STORAGE1
setenv SIGNALSTORE
setenv BKGSTORE

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

echo "--> Setup mva training"
pwd
date
cmsrel $CMSSW
cd $CMSSW/src
eval `scramv1 runtime -csh`
setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:${LD_LIBRARY_PATH}

echo "--> Extract tar file"
date
tar zxf ../../$JOB.tar.gz
cd HeavyFlavorAnalysis/Bs2MuMu/macros
mv ../../../../../$JOB .

# -- use CN's root version
setenv ROOTSYS /shome/naegelic/root
set path = ($ROOTSYS/bin $path)
if ($?LD_LIBRARY_PATH) then
   setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH      # Linux, ELF HP-UX                                                                                  
else
   setenv LD_LIBRARY_PATH $ROOTSYS/lib
endif

if ($?DYLD_LIBRARY_PATH) then
   setenv DYLD_LIBRARY_PATH $ROOTSYS/lib:$DYLD_LIBRARY_PATH  # Mac OS X                                                                                          
else
   setenv DYLD_LIBRARY_PATH $ROOTSYS/lib
endif

if ($?SHLIB_PATH) then
   setenv SHLIB_PATH $ROOTSYS/lib:$SHLIB_PATH                # legacy HP-UX                                                                                      
else
   setenv SHLIB_PATH $ROOTSYS/lib
endif

if ($?LIBPATH) then
   setenv LIBPATH $ROOTSYS/lib:$LIBPATH                      # AIX                                                                                               
else
   setenv LIBPATH $ROOTSYS/lib
endif

if ($?PYTHONPATH) then
   setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH
else
   setenv PYTHONPATH $ROOTSYS/lib
endif

if ($?MANPATH) then
   setenv MANPATH `dirname $ROOTSYS/man/man1`:$MANPATH
else
   setenv MANPATH `dirname $ROOTSYS/man/man1`:$default_manpath
endif

echo "--> Building ncBatch"
make ncBatch
pwd

SIGNALFILE=`basename $SIGNALSTORE`
BKGFILE=`basename $BKGSTORE`

echo "--> Copying training files to scratch"
echo "$SRMCP $SIGNALSTORE file:///`pwd`/$SIGNALFILE"
$SRMCP "$SIGNALSTORE" file:///`pwd`/$SIGNALFILE
echo "$SRMCP $BKGSTORE file:///`pwd`/$BKGFILE"
$SRMCP "$BKGSTORE" file:///`pwd`/$BKGFILE

MLPOPTIONS=`cat $JOB`
echo "$EXECUTABLE $MLPOPTIONS $SIGNALFILE $BKGFILE"
$EXECUTABLE $MLPOPTIONS $SIGNALFILE $BKGFILE

# copy to storage element
FILE_BARREL=${JOB}_B.root
FILE_ENDCAP=${JOB}_E.root

# Barrel
echo "srmrm $STORAGE1/$FILE_BARREL"
srmrm "$STORAGE1/$FILE_BARREL"
echo "$SRMCP file:///`pwd`/TMVA_0.root $STORAGE1/$FILE_BARREL"
$SRMCP file:///`pwd`/TMVA_0.root "$STORAGE1/$FILE_BARREL"
echo "srmls $STORAGE1/$FILE_BARREL"
srmls "$STORAGE1/$FILE_BARREL"

# Endcap
echo "srmrm $STORAGE1/$FILE_ENDCAP"
srmrm "$STORAGE1/$FILE_ENDCAP"
echo "$SRMCP file:///`pwd`/TMVA_1.root $STORAGE1/$FILE_ENDCAP"
$SRMCP file:///`pwd`/TMVA_1.root "$STORAGE1/$FILE_ENDCAP"
echo "srmls $STORAGE1/$FILE_ENDCAP"
srmls "$STORAGE1/$FILE_ENDCAP"

# BATCH END

echo "run: This is the end, my friend"
