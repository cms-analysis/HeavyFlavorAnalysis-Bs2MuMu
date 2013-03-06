#! /bin/csh -f


# ----------------------------------------------------------------------
# example submission: 
# -------------------
#  run -t ../../../../130131.tar.gz -m grid -D unl.edu -c ../osg.csh -r 'PFNS srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat%STORAGE1 /store/user/ursl/test%SITE T3_CH_PSI' osgtest.py
# ----------------------------------------------------------------------

setenv CMSSW       
setenv SCRAM_ARCH  
setenv SRMCP       

setenv JOB      
setenv FILE1    $JOB.root
setenv STORAGE1 
setenv PFNS 
setenv SITE

echo "========================"
echo "====> osg wrapper <===="
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
df -kl 

source $OSG_GRID/setup.csh
echo $OSG_APP
ls -l $OSG_APP/cmssoft/cms
echo source $OSG_APP/cmssoft/cms/cmsset_default.csh
source $OSG_APP/cmssoft/cms/cmsset_default.csh
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
# -- Setup CMSSW
# ----------------------------------------------------------------------
echo "--> Setup CMSSW"
pwd
date
cmsrel $CMSSW
cd $CMSSW
eval `scramv1 runtime -csh`
pwd

echo "--> Extract tar file"
date
tar zxf ../$JOB.tar.gz
cd src
mv ../../$JOB.py .
mv ../../data_replica.py .
chmod 755 data_replica.py

# ----------------------------------------------------------------------
# -- Run cmsRun
# ----------------------------------------------------------------------
echo "--> Run cmsRun"
pwd
date
which cmsRun
echo "cmsRun $JOB.py "
cmsRun $JOB.py |& tee $JOB.log
date
pwd
ls -rtl 

setenv LD_LIBRARY_PATH /lib64:${LD_LIBRARY_PATH}

setenv ROOTFILE `ls *.root`

# ----------------------------------------------------------------------
# -- Save Output to SE
# ----------------------------------------------------------------------
echo "--> Save output to SE: $PFNS/$STORAGE1/$FILE1"
echo " local rootfile: $ROOTFILE"
echo " job   rootfile: $FILE1"


echo lcg-del -b -D srmv2 -l  "$PFNS/$STORAGE1/$FILE1"
lcg-del -b -D srmv2 -l "$PFNS/$STORAGE1/$FILE1"

# -- switch to data_replica.py
ls `pwd`/$FILE1 > dr.list
echo "--> cat dr.list: " 
cat dr.list
echo "--> AM running data_replica.py: " 
./data_replica.py --from-site LOCAL --to-site $SITE dr.list  "$STORAGE1"

echo "--> lcg-ls : $PFNS/$STORAGE1/$FILE1" 
echo lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE1/$FILE1"
lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE1/$FILE1"

date

# BATCH END

echo "run: This is the end, my friend"
