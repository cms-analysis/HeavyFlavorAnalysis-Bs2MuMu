#!/bin/bash

# This is expected to be set by the site
#export VO_CMS_SW_DIR=/experiment-software/cms 


# Defaults:
CMSSW_VERSION=CMSSW_1_3_1_HLT6
BASEDIR=`pwd`  # BASEDIR must be the dir containing the input files
# tar -chzf libHeavyFlavorAnalysisBs2MuMu.tar.gz lib module src/HeavyFlavorAnalysis/Bs2MuMu/data/*.c*
TARBALL=libHeavyFlavorAnalysisBs2MuMu.tar.gz
CMSSW_CFG=mytest.cfg
OUTPUTFILE=output.root

usage() {
    cat <<EOF
synopsis: runAnalysis.sh -o srm-URL [options]

    -o srm-URL     :  SRM filename for where to copy the result file to

    -c filename    :  Name of the CMSSW configuration file [$CMSSW_CFG]
    -v string      :  Version string for the CMSSW version [$CMSSW_VERSION]
    -t filename    :  Name of a tarball to unpack into the working directory

The script expects to find the defined filenames in the directory from which it
gets started.
The script relies on the output file from the cmsRun job being called $OUTPUTFILE

EOF
}


# This makes error output easier
die() {
    echo "runAnalysis error: $1" >&2
    exit 1
}

getFileFromSE() {
    local URL=$1

    local protocol=`expr $URL : "\([^:]*\)://"`
    local base=`basename $URL`
    local localfile=$PWD/$base
    case x"$protocol" in
	xsrm)
          cmd="srmcp -retry_num=1 $URL file:///$localfile"
	  ;;
        xdcap)
          cmd="dccp $URL $localfile"
	  ;;
        xgsiftp)
          cmd="globus-url-copy $URL file:$localfile"
	  ;;
	x)
	  cmd="cp $URL $localfile"
	  exit
	  ;;
        *)
          die "Unknown protocol ($protocol) for URL $URL"
	  ;;
    esac

    echo "runAnalysis: $cmd"
    eval $cmd
    status=$?
    if test $status -ne 0; then
	die "Failed downloading $URL"
	return $status
    fi

    return 0
}


# Argument parsing

echo " "
echo "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
echo "                         ---> Start of grid job wrapper:                           "
echo "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
date

TEMP=`getopt -o c:ho:t:v: --long help -n 'runAnalysis.sh' -- "$@"`
if [ $? != 0 ] ; then usage ; echo "Terminating..." >&2 ; exit 1 ; fi
#echo "TEMP: $TEMP"
eval set -- "$TEMP"

while true; do
    case "$1" in
        --help|-h)
            usage
            exit
            ;;
	-o)
	    PFN=$2
	    shift 2
	    ;;
        -c)
	    CMSSW_CFG=$2
            shift 2
            ;;
        -t)
	    TARBALL=$2
            shift 2
            ;;
        -v)
	    CMSSW_VERSION=$2
            shift 2
            ;;
	--)
            shift;
            break;
            ;;
	*)
            echo "Internal error!"
            exit 1
            ;;
    esac
done

# Sanity checks
echo " "
echo "ooooooooooooooooooooooooooooooooo Sanity checks ooooooooooooooooooooooooooooooooo"
if test ! -r $CMSSW_CFG; then
    die "CMSSW configuration file ($CMSSW_CFG) cannot be read"
fi
if test x"$PFN" = x -o ${PFN#srm://} = $PFN; then
    usage
    echo  "$PFN  ${PFN#srm://}"
    die 'Need an SRM address (srm://....) as argument to the -o option!!!'
fi


if test ! -r $VO_CMS_SW_DIR/cmsset_default.sh; then
    die "No environment file $VO_CMS_SW_DIR/cmsset_default.sh found"
fi
source $VO_CMS_SW_DIR/cmsset_default.sh
if test $? -ne 0; then
    die "Error in sourcing $VO_CMS_SW_DIR/cmsset_default.sh found"
fi


WORKDIR=`pwd`/$CMSSW_VERSION
rm -rf $WORKDIR

# Create project release area
echo " "
echo "ooooooooooooooooooooooooooo Create project release area ooooooooooooooooooooooooooo"
export SCRAM_ARCH=slc3_ia32_gcc323
scramv1 project CMSSW $CMSSW_VERSION


if test ! -d $WORKDIR; then
    die "Work area $WORKDIR was not created"
fi

# unpack any supplied tarball into the working directory
echo " "
echo "ooooooooooooooooooooooooooooooooo Unpack tarball ooooooooooooooooooooooooooooooooo"
if test x"$TARBALL" != x; then
    cd $WORKDIR
    getFileFromSE $TARBALL 
    tarfile=$WORKDIR/$(basename $TARBALL)
    echo "pwd: "
    pwd
    echo "ls -l: "
    ls -l
    echo "Extracting library tarball $tarfile"
    tar tzf $tarfile
    tar xzf $tarfile
    if test $? -ne 0; then
	die "Failed unpacking tarball \($tarfile\)"
    fi
fi
    
cd $WORKDIR/src/
eval `scramv1 runtime -sh`
if test $? -ne 0; then
    die "scram failed to source the environment"
fi


# run CMSSW
echo " "
echo "oooooooooooooooooooooooooooooooooooo Run CMSSW oooooooooooooooooooooooooooooooooooo"
cd $WORKDIR
cmsRun $BASEDIR/$CMSSW_CFG


echo "oooooooooooooooooooooooooooooooo End of CMSSW run ooooooooooooooooooooooooooooooooo"
echo "pwd: "
pwd
echo "ls -l: "
ls -l

# Save the output to the SE
echo " "
echo "oooooooooooooooooooooooooooooooo Save output to SE oooooooooooooooooooooooooooooooo"
if test ! -e $OUTPUTFILE; then
    die "No outputfile $OUTPUTFILE was written"
fi

echo "Copy file $OUTPUTFILE to $PFN"
srmcp file:///`pwd`/$OUTPUTFILE $PFN
if test $? -ne 0; then
    die "srmcp failed to copy the file to $PFN"
    # Do some smart thing to save the file somewher
fi

echo "runAnalysis: saved job to $PFN"

echo " "
echo "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
echo "                         ---> End of grid job wrapper:                             "
echo "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
date

exit 0

