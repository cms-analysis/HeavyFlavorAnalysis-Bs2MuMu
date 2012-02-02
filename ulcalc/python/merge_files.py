#!/usr/bin/python

import sys
import os

if len(sys.argv) < 3:
    print "usage: " + sys.argv[0] + " <targetfile> <name> <sourcefiles>"
    sys.exit(0)


rootfilename="merge.C"
rootfilecontent='''
#include <TFile.h>
#include <RooStats/HypoTestInverterResult.h>

using namespace RooStats;

void merge(const char *filename, const char *resultname, const char *fromfile)
{
	TFile theFile(filename,"update");
	TFile fromF(fromfile);
	
	RooStats::HypoTestInverterResult *result = dynamic_cast<RooStats::HypoTestInverterResult*>(theFile.Get(resultname));
	if(!result) // not yet loaded
	  result = new RooStats::HypoTestInverterResult(resultname);

	RooWorkspace *wspace = dynamic_cast<RooWorkspace*> (fromF.Get("wspace"));
	RooStats::HypoTestInverterResult *toadd = dynamic_cast<RooStats::HypoTestInverterResult*> (wspace->obj(resultname));
	
	result->Add(*toadd);
	theFile.cd();
	result->Write(resultname,TObject::kOverwrite);
	theFile.Close();
} // merge_result()
'''

targetfile=sys.argv[1]
name=sys.argv[2]
sources=sys.argv[3:len(sys.argv)]

print "Targetfile: " + targetfile
print "Name: " + name
print "Sources: " + str(sources)

# create the temporary root file
rfile = open(rootfilename,"w")
rfile.write(rootfilecontent)
rfile.close()

for src in sources:
    cmd = 'root -l -b -q \'' + rootfilename + '("' + targetfile + '","' + name + '","' + src + '")\''
    print cmd
    os.system(cmd)

os.remove(rootfilename)
