#!/bin/bash
#################################
# PSI Tier-3 example batch Job  #
#################################

##### CONFIGURATION ##############################################
# Output files to be copied back to the User Interface
# (the file path must be given relative to the working directory)
#OUTFILES="UpperlimitCalc$1.txt"
if [ "$1" == "" ] || [ "$2" == "" ] 
then
	echo "Option 1 is the channel name"
	echo "Option 2 is the random seed"
	exit 0
fi

if [ "$1" == "comb11" ]
then
	OUTFILES="CMS_LHCb_S11_Comb_Bs_Results$2.txt CMS_LHCb_S11_Comb_Bs_Results$2.root CMS_LHCb_S11_Comb_Bs_Results$2.pdf"
elif [ "$1" == "comb12" ]
then
	OUTFILES="CMS_LHCb_W12_Comb_Bs_Results$2.txt CMS_LHCb_W12_Comb_Bs_Results$2.root CMS_LHCb_W12_Comb_Bs_Results$2.pdf"
elif [ "$1" == "comb12bd" ]
then
	OUTFILES="CMS_LHCb_W12_Comb_Bd_Results$2.txt CMS_LHCb_W12_Comb_Bd_Results$2.root CMS_LHCb_W12_Comb_Bd_Results$2.pdf"
elif [ "$1" == "cms11bs" ]
then
	OUTFILES="CMS_S11_Bs_Results$2.txt CMS_S11_Bs_Results$2.root CMS_S11_Bs_Results$2.pdf"
elif [ "$1" == "cms11bd" ]
then
	OUTFILES="CMS_S11_Bd_Results$2.txt CMS_S11_Bd_Results$2.root CMS_S11_Bd_Results$2.pdf"
elif [ "$1" == "cms12bs" ]
then
	OUTFILES="CMS_W12_Bs_Results$2.txt CMS_W12_Bs_Results$2.root CMS_W12_Bs_Results$2.pdf"
elif [ "$1" == "cms12bd" ]
then
	OUTFILES="CMS_W12_Bd_Results$2.txt CMS_W12_Bd_Results$2.root CMS_W12_Bd_Results$2.pdf"
elif [ "$1" == "lhcb" ]
then 
	OUTFILES="LHCb_10_11_Comb_Bs_Results$2.txt LHCb_10_11_Comb_Bs_Results$2.root LHCb_10_11_Comb_Bs_Results$2.pdf"
elif [ "$1" == "lhcb10" ]
then 
	OUTFILES="LHCb_10_Bs_Results$2.txt LHCb_10_Bs_Results$2.root LHCb_10_Bs_Results$2.pdf"
elif [ "$1" == "lhcb11" ]
then
	OUTFILES="LHCb_11_Bs_Results$2.txt LHCb_11_Bs_Results$2.root LHCb_11_Bs_Results$2.pdf"
elif [ "$1" == "lhcb12" ]
then
	OUTFILES="LHCb_12_Bs_Results$2.txt LHCb_12_Bs_Results$2.root LHCb_12_Bs_Results$2.pdf"
elif [ "$1" == "lhcb12bd" ]
then
	OUTFILES="LHCb_12_Bd_Results$2.txt LHCb_12_Bd_Results$2.root LHCb_12_Bd_Results$2.pdf"
fi

# Output files to be copied to the SE
# (as above the file path must be given relative to the working directory)
SEOUTFILES=""
#
# By default, the files will be copied to $USER_SRM_HOME/$username/$JOBDIR,
# but here you can define the subdirectory under your SE storage path
# to which files will be copied (uncomment option)
#SEUSERSUBDIR="mytestsubdir/somedir"
#
# User's CMS hypernews name (needed for user's SE storage home path
# USER_SRM_HOME below)
HN_NAME=laflores

# set DBG=1 for additional debug output in the job report files
# DBG=2 will also give detailed output on SRM operations
DBG=0

#### The following configurations you should not need to change
# The SE's user home area (SRMv2 URL)
USER_SRM_HOME="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/"

# Top working directory on worker node's local disk. The batch
# job working directory will be created below this
TOPWORKDIR=/scratch/`whoami`

# Basename of job sandbox (job workdir will be $TOPWORKDIR/$JOBDIR)
JOBDIR=sgejob-$JOB_ID
##################################################################


############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N upperlim_calc

# Run time soft and hard limits hh:mm:ss
# soft=CPU time; hard=Wallclock time

### Specify the queue on which to run
#$ -q all.q

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
#  #$ -o /shome/username/mydir/
#  #$ -e /shome/username/mydir/

##################################################################



##### MONITORING/DEBUG INFORMATION ###############################
DATE_START=`date +%s`
echo "Job started at " `date`
cat <<EOF
################################################################
## QUEUEING SYSTEM SETTINGS:
HOME=$HOME
USER=$USER
JOB_ID=$JOB_ID
JOB_NAME=$JOB_NAME
HOSTNAME=$HOSTNAME
TASK_ID=$TASK_ID
QUEUE=$QUEUE

EOF

echo "################################################################"

if test 0"$DBG" -gt 0; then
   echo "######## Environment Variables ##########"
   env
   echo "################################################################"
fi


##### SET UP WORKDIR AND ENVIRONMENT ######################################
STARTDIR=`pwd`
WORKDIR=$TOPWORKDIR/$JOBDIR
RESULTDIR=$STARTDIR/$JOBDIR
if test -e "$WORKDIR"; then
   echo "ERROR: WORKDIR ($WORKDIR) already exists! Aborting..." >&2
   exit 1
fi
mkdir -p $WORKDIR
if test ! -d "$WORKDIR"; then
   echo "ERROR: Failed to create workdir ($WORKDIR)! Aborting..." >&2
   exit 1
fi

cd $WORKDIR
cat <<EOF
################################################################
## JOB SETTINGS:
STARTDIR=$STARTDIR
WORKDIR=$WORKDIR
RESULTDIR=$RESULTDIR
EOF


###########################################################################
## YOUR FUNCTIONALITY CODE GOES HERE
CMSSW_DIR=/shome/lazo-flores/CMSSW_4_2_8
echo "source $VO_CMS_SW_DIR/cmsset_default.sh"
source $VO_CMS_SW_DIR/cmsset_default.sh
cd $CMSSW_DIR/src
#cmsenv
eval `scramv1 runtime -sh`

cd $WORKDIR

echo " about to run the job "
cd /shome/lazo-flores/CMSSW_4_2_8/src/HeavyFlavorAnalysis/Bs2MuMu/mclimit
bin/upperlimitCalc $1 $2 

for n in $OUTFILES; do
	mv $n $WORKDIR
done

#### RETRIEVAL OF OUTPUT FILES AND CLEANING UP ############################
cd $WORKDIR
if test 0"$DBG" -gt 0; then
    echo "########################################################"
    echo "############# Working directory contents ###############"
    echo "pwd: " `pwd`
    ls -Rl
    echo "########################################################"
    echo "YOUR OUTPUT WILL BE MOVED TO $RESULTDIR"
    echo "########################################################"
fi

if test x"$OUTFILES" != x; then
   mkdir -p $RESULTDIR
   if test ! -e "$RESULTDIR"; then
          echo "ERROR: Failed to create $RESULTDIR ...Aborting..." >&2
	  exit 1
   fi
   for n in $OUTFILES; do
       if test ! -e $WORKDIR/$n; then
          echo "WARNING: Cannot find output file $WORKDIR/$n. Ignoring it" >&2
       else
          cp -a $WORKDIR/$n $RESULTDIR/$n
          if test $? -ne 0; then
             echo "ERROR: Failed to copy $WORKDIR/$n to $RESULTDIR/$n" >&2
          fi
   fi
   done
fi

if test x"$SEOUTFILES" != x; then
   if test x"$SEUSERSUBDIR" = x; then
      SERESULTDIR=$USER_SRM_HOME/$HN_NAME/$JOBDIR
   else
      SERESULTDIR=$USER_SRM_HOME/$HN_NAME/$SEUSERSUBDIR
   fi
   if test 0"$DBG" -ge 2; then
      srmdebug="-debug"
   fi
   for n in $SEOUTFILES; do
       if test ! -e $WORKDIR/$n; then
          echo "WARNING: Cannot find output file $WORKDIR/$n. Ignoring it" >&2
       else
          srmcp $srmdebug -2 file:///$WORKDIR/$n $SERESULTDIR/$n
          if test $? -ne 0; then
             echo "ERROR: Failed to copy $WORKDIR/$n to $SERESULTDIR/$n" >&2
          fi
   fi
   done
fi

echo "Cleaning up $WORKDIR"
rm -rf $WORKDIR

###########################################################################
DATE_END=`date +%s`
RUNTIME=$((DATE_END-DATE_START))
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0
