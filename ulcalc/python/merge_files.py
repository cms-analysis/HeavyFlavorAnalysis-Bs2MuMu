#!/usr/bin/python

import sys
import os

if len(sys.argv) < 3:
    print "usage: " + sys.argv[0] + " <targetfile> <name> <sourcefiles>"
    sys.exit(0)


rootfilename="merge.C"
rootfilecontent='''
#include <TROOT.h>
#include <TFile.h>
#include <RooWorkspace.h>
#include <RooStats/HypoTestResult.h>
#include <RooStats/HypoTestInverterResult.h>

const char* result_names[] = {
	"Hybrid_CLb",
	"Hybrid_Ul",
	"Hybrid_Int"
};

void mergeSpaces(RooWorkspace *targetSpace, RooWorkspace *addSpace, std::string resultName)
{
	TObject *obj;
	RooStats::HypoTestResult *targetResult = NULL;
	RooStats::HypoTestResult *addResult = NULL;
	RooStats::HypoTestInverterResult *targetInverter = NULL;
	RooStats::HypoTestInverterResult *addInverter = NULL;
	
	// load the target Result
	obj = targetSpace->obj(resultName.c_str());
	if (obj) {
		if (obj->InheritsFrom("RooStats::HypoTestResult"))
			targetResult = (RooStats::HypoTestResult*)obj;
		if (obj->InheritsFrom("RooStats::HypoTestInverterResult"))
			targetInverter = (RooStats::HypoTestInverterResult*)obj;
	}
	
	obj = addSpace->obj(resultName.c_str());
	if (obj) {
		if (obj->InheritsFrom("RooStats::HypoTestResult"))
			addResult = (RooStats::HypoTestResult*)obj;
		if (obj->InheritsFrom("RooStats::HypoTestInverterResult"))
			addInverter = (RooStats::HypoTestInverterResult*)obj;
	}
	
	if (addResult) {
		if (targetResult)	targetResult->Append(addResult);
		else				targetResult = addResult;
		
		cout << resultName << " size = " << targetResult->GetNullDistribution()->GetSize() << endl;
		targetSpace->import(*targetResult, kTRUE);
	}
	
	if (addInverter) {
		if (targetInverter)	targetInverter->Add(*addInverter);
		else				targetInverter = addInverter;
		
		targetSpace->import(*targetInverter, kTRUE);
	}
} // mergeSpaces()

// void merge(const char *filename, const char *resultname, const char *fromfile)
void merge(const char *targetName, const char *addName)
{
	TFile *targetFile = new TFile(targetName, "update");
	TFile *addFile = TFile::Open(addName);
	const char *wname = "tmp-wspace-merge";
	unsigned j;
	
	// load the workspace
	RooWorkspace *targetSpace = (RooWorkspace*)targetFile->Get("wspace");
	RooWorkspace *addSpace = (RooWorkspace*)addFile->Get("wspace");
	
	if (addSpace) {
		if (targetSpace) {
			for (j = 0; j < sizeof(result_names)/sizeof(const char *); j++) {
				mergeSpaces(targetSpace, addSpace, std::string(Form("%s_mu_s",result_names[j])));
				mergeSpaces(targetSpace, addSpace, std::string(Form("%s_mu_s_SM",result_names[j])));
				mergeSpaces(targetSpace, addSpace, std::string(Form("%s_mu_d",result_names[j])));
				mergeSpaces(targetSpace, addSpace, std::string(Form("%s_mu_d_SM",result_names[j])));
			}
		} else {
			cout << Form("No target workspace found, copying from '%s'",addName) << endl;
			targetSpace = addSpace;
		}
	}
	
	if (targetSpace)
		targetSpace->writeToFile(wname, kTRUE);
	
	// close the files
	delete targetFile;
	delete addFile;
	
	// move the workspace in place
	if (targetSpace)
		gROOT->ProcessLine(Form(".!mv %s %s", wname, targetName));
} // merge()
'''

targetfile=sys.argv[1]
sources=sys.argv[2:len(sys.argv)]

print "Targetfile: " + targetfile
print "Sources: " + str(sources)

# create the temporary root file
rfile = open(rootfilename,"w")
rfile.write(rootfilecontent)
rfile.close()

for src in sources:
    cmd = 'root -l -b -q \'' + rootfilename + '("' + targetfile + '","' + src + '")\''
    print cmd
    os.system(cmd)

os.remove(rootfilename)
