#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>

#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TString.h"
#include "TRandom.h"
#include "TUnixSystem.h"

#include "bmmReader.hh"
#include "dpReader.hh"
#include "d0Reader.hh"
#include "massReader.hh"
#include "ksReader.hh"
#include "copyReader.hh"
#include "dumpReader.hh"
#include "muCharmReader.hh"
#include "triggerValidation.hh"
#include "genLevel.hh"

using namespace std;


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: ./runTreeReaders -f test.root
// %%        ./runTreeReaders -c chains/bg-test -D root
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "Running under process ID  " << processID << endl;

  string progName  = argv[0]; 
  string writeName, fileName;
  int file(0);
  int dirspec(0);
  int nevents(-1), start(-1);
  int randomSeed(processID);

  // -- Some defaults
  string dirBase("./");               // this could point to "/home/ursl/data/root/."
  string dirName("."); dirspec = 0;   // and this to, e.g. "bmm", "bee", "bem", ...
  string cutFile("tree.defaults.cuts");

  string treeName("T1");
  string evtClassName("TAna01Event");

  string readerName("bmmReader"); 
  TString histfile("");

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i],"-c"))  {fileName   = string(argv[++i]); file = 0; }     // file with chain definition
    if (!strcmp(argv[i],"-C"))  {cutFile    = string(argv[++i]);           }     // file with cuts
    if (!strcmp(argv[i],"-D"))  {dirName    = string(argv[++i]);  dirspec = 1; } // where to put the output
    if (!strcmp(argv[i],"-f"))  {fileName   = string(argv[++i]); file = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-n"))  {nevents    = atoi(argv[++i]); }                 // number of events to run 
    if (!strcmp(argv[i],"-r"))  {readerName = string(argv[++i]); }               // which tree reader class to run
    if (!strcmp(argv[i],"-s"))  {randomSeed = atoi(argv[++i]); }                 // set seed for random gen.
    if (!strcmp(argv[i],"-o"))  {histfile   = TString(argv[++i]); }              // set output file
  }


  // -- Prepare histfilename variation with (part of) cut file name
  TString fn(cutFile);
  fn.ReplaceAll("cuts/", "");
  fn.ReplaceAll(".cuts", "");
  fn.ReplaceAll("tree", "");

  // -- Determine filename for output histograms and 'final' small/reduced tree
  TString meta = fileName;
  if(histfile == "") {
    TString  barefile(fileName), chainFile, meta;
    if (file == 0) {
      // -- input from chain
      if (barefile.Contains("chains/")) {
	barefile.ReplaceAll("chains/", "");
	histfile = barefile + "." + fn + ".root";
	if (dirspec) {
	  if (dirName[0] == '/') {
	    histfile = dirName + "/" + histfile;
	  } else {
	    histfile = dirBase + "/" + dirName + "/" + histfile;
	  }
	}
      } else {
	histfile =  barefile + "." + fn + ".root";
	if (dirspec) {
	  if (dirName[0] == '/') {
	    histfile = dirName + "/" + histfile;
	  } else {
	    histfile = dirBase + "/" + dirName + "/" + histfile;
	  }
	}
      }
      // -- The following lines strip everything from the string up to and including the last '/'
      int fl = barefile.Last('/');
      TString bla(barefile);
      bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
      histfile =  bla + "." + fn + ".root";
      if (dirspec) {
	histfile = dirBase + "/" + dirName + "/" + histfile;
      }
    }  else if (file == 1) {
      // -- single file input
      // -- The following lines strip everything from the string up to and including the last '/'
      int fl = barefile.Last('/');
      TString bla(barefile);
      bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
      histfile =  bla;
      histfile.ReplaceAll(".root", "");
      histfile +=  "." + fn + ".root";
      if (dirspec) {
	if (dirName[0] == '/') {
	  histfile = dirName + "/" + histfile;
	} else {
	  histfile = dirBase + "/" + dirName + "/" + histfile;
	}
      }
    }
  }

  cout << "Opening " << histfile.Data() << " for output histograms" << endl;
  cout << "Opening " << fileName.c_str() << " for input" << endl;


  // -- Set up chain
  TChain *chain = new TChain(TString(treeName));
  cout << "Chaining ... " << treeName << endl;
  char pName[2000]; 
  int nentries; 
  if (file == 0) {
    // -- non-trivial chain input
    ifstream is(meta);  
    while(meta.ReadLine(is) && (!meta.IsNull())){ 
      nentries = -1;
      if (meta.Data()[0] == '#') continue; 
      sscanf(meta.Data(), "%s %d", pName, &nentries); 
      if (nentries > -1) {
        cout << pName << " -> " << nentries << " entries" << endl; 
        chain->Add(pName, nentries); 
      } else {
        cout << meta << endl;
        chain->Add(meta); 
      }
    }
    is.close();
  }
  else if (file == 1) {
    // -- single file input
    cout << fileName << endl;
    chain->Add(TString(fileName));
  }

  // -- Now instantiate the tree-analysis class object, initialize, and run it ...
  treeReader01 *a; 
  if (readerName == "bmmReader") a = new bmmReader(chain, TString(evtClassName));
  else if (readerName == "massReader") a = new massReader(chain,TString(evtClassName));
  else if (readerName == "ksReader") a = new ksReader(chain,TString(evtClassName));
  else if (readerName == "dumpReader") a = new dumpReader(chain,TString(evtClassName));
  else if (readerName == "dpReader") a = new dpReader(chain,TString(evtClassName));
  else if (readerName == "d0Reader") a = new d0Reader(chain,TString(evtClassName));
  else if (readerName == "copyReader") a = new copyReader(chain,TString(evtClassName));
  else if (readerName == "muCharmReader") a = new muCharmReader(chain,TString(evtClassName));
  else if (readerName == "triggerValidation") a = new triggerValidation(chain,TString(evtClassName));
  else if (readerName == "genLevel") a = new genLevel(chain,TString(evtClassName));
  else {
    cout << "please provide a class name to instantiate" << endl;
  }


  if (a) {
    a->openHistFile(histfile); 
    a->bookHist(); 
    a->readCuts(cutFile, 1);
    
    a->startAnalysis(); 
    a->loop(nevents, start);
    a->endAnalysis();
    a->closeHistFile(); 
  } else
    cerr << "Readerclass '" << readerName << "' not found" << endl;

  delete a; // so we can dump some information in the destructor
  
  return 0;
}

