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
cd HeavyFlavorAnalysis/Bs2MuMu/ana-naegeli
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

echo "--> Building rootutils"
cd ../rootutils; make; cd -

echo "--> Building ncBatch"
make
pwd

echo "SIGNALFILE=/bin/basename $SIGNALSTORE"
set SIGNALFILE = `/bin/basename "$SIGNALSTORE"`
echo "BKGFILE=/bin/basename $BKGSTORE"
set BKGFILE = `/bin/basename "$BKGSTORE"`

echo "--> Copying training files to scratch"
echo "$SRMCP $SIGNALSTORE file:///`pwd`/$SIGNALFILE"
$SRMCP "$SIGNALSTORE" file:///`pwd`/$SIGNALFILE
echo "$SRMCP $BKGSTORE file:///`pwd`/$BKGFILE"
$SRMCP "$BKGSTORE" file:///`pwd`/$BKGFILE

set MLPOPTIONS = `cat $JOB`
echo "$EXECUTABLE $MLPOPTIONS $SIGNALFILE $BKGFILE"
$EXECUTABLE $MLPOPTIONS $SIGNALFILE $BKGFILE

# copy to storage element
set FILE_EVEN_BARREL = ${JOB}_s0_c0.root
set FILE_ODD_BARREL = ${JOB}_s1_c0.root
set FILE_EVEN_ENDCAP = ${JOB}_s0_c1.root
set FILE_ODD_ENDCAP = ${JOB}_s1_c1.root

# Barrel Even
echo "srmrm $STORAGE1/$FILE_EVEN_BARREL"
srmrm "$STORAGE1/$FILE_EVEN_BARREL"
echo "$SRMCP file:///`pwd`/TMVA_s0_c0.root $STORAGE1/$FILE_EVEN_BARREL"
$SRMCP file:///`pwd`/TMVA_s0_c0.root "$STORAGE1/$FILE_EVEN_BARREL"
echo "srmls $STORAGE1/$FILE_EVEN_BARREL"
srmls "$STORAGE1/$FILE_EVEN_BARREL"

# Barrel Odd
echo "srmrm $STORAGE1/$FILE_ODD_BARREL"
srmrm "$STORAGE1/$FILE_ODD_BARREL"
echo "$SRMCP file:///`pwd`/TMVA_s1_c0.root $STORAGE1/$FILE_ODD_BARREL"
$SRMCP file:///`pwd`/TMVA_s1_c0.root "$STORAGE1/$FILE_ODD_BARREL"
echo "srmls $STORAGE1/$FILE_ODD_BARREL"
srmls "$STORAGE1/$FILE_ODD_BARREL"

# Endcap Even
echo "srmrm $STORAGE1/$FILE_EVEN_ENDCAP"
srmrm "$STORAGE1/$FILE_EVEN_ENDCAP"
echo "$SRMCP file:///`pwd`/TMVA_s0_c1.root $STORAGE1/$FILE_EVEN_ENDCAP"
$SRMCP file:///`pwd`/TMVA_s0_c1.root "$STORAGE1/$FILE_EVEN_ENDCAP"
echo "srmls $STORAGE1/$FILE_EVEN_ENDCAP"
srmls "$STORAGE1/$FILE_EVEN_ENDCAP"

# Endcap Odd
echo "srmrm $STORAGE1/$FILE_ODD_ENDCAP"
srmrm "$STORAGE1/$FILE_ODD_ENDCAP"
echo "$SRMCP file:///`pwd`/TMVA_s1_c1.root $STORAGE1/$FILE_ODD_ENDCAP"
$SRMCP file:///`pwd`/TMVA_s1_c1.root "$STORAGE1/$FILE_ODD_ENDCAP"
echo "srmls $STORAGE1/$FILE_ODD_ENDCAP"
srmls "$STORAGE1/$FILE_ODD_ENDCAP"


# BATCH END

echo "run: This is the end, my friend"
