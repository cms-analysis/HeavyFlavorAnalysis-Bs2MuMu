#!/bin/csh -f

#################################################################################
# Usage: ./build_crab <XXX> <crab.cfg> <cmsRunConfig.py> <jsonfile> [<prefix>] ##
#   Builds a directory to run crab. All 'XXX' is replaced by the command text  ##
#   given in the <XXX> argument.                                               ##
#################################################################################

# command line parsing...
if ($# < 4) then
    echo "usage: $0 <XXX> <crab.cfg> <cmsRunConfig.py> <jsonfile> [<prefix>]"
    exit 1
endif

# print the current configuration
set RUNNUMBER = $1
set CRAB_INPUT = $2
set CMSRUN_INPUT = $3
set JSONFILE = $4
if ($# >= 5) then
    set WORKINGDIR = $5-$RUNNUMBER
else
    set WORKINGDIR = $RUNNUMBER
endif

echo "Working directory: $WORKINGDIR"
echo "Run number: $RUNNUMBER"
echo "Crab config file: $CRAB_INPUT"
echo "cmsRun config file: $CMSRUN_INPUT"
echo "json file: $JSONFILE"

# create the directory
mkdir $WORKINGDIR
if ($? != 0) then
    echo "Error with working directory. Aborting..."
    exit 1
endif

# copy the json file
cp $JSONFILE $WORKINGDIR/jsonfile

# create the cmsRun config file
set CMSRUN_OUTPUT = `ls $CMSRUN_INPUT | awk '{gsub("XXX","'$RUNNUMBER'");print $0}'`
awk '{gsub("XXX","'$RUNNUMBER'");print $0}' < $CMSRUN_INPUT > $WORKINGDIR/$CMSRUN_OUTPUT

# create the crab config file
awk '{gsub("XXX","'$RUNNUMBER'");print $0}' < $CRAB_INPUT > $WORKINGDIR/crab.cfg

echo "Successfully built the directory $WORKINGDIR"

