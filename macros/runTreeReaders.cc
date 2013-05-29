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
#include "TRandom.h"
#include "TUnixSystem.h"

#include "massReader.hh"
#include "copyReader.hh"
#include "genLevel.hh"
#include "impactReader.hh"
#include "pixelReader.hh"

using namespace std;


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: ./runTreeReaders -f test.root
// %%        ./runTreeReaders -c chains/bg-test -D root
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "Running under process ID  " << processID << endl;

  string progName  = argv[0]; 
  string writeName, fileName, jsonName;
  int file(0), json(0);
  int dirspec(0);
  int nevents(-1), start(-1);
  int randomSeed(processID);
  int verbose(-99); 
  int blind(1); 
  int isMC(0);

  // Change the MaxTreeSize to 100 GB (default since root v5.26)
  TTree::SetMaxTreeSize(100000000000ll); // 100 GB

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
    if (!strcmp(argv[i],"-h")) {
	std::cout << "List of arguments:" << std::endl;
	cout << "-b {0,1}      run blind?" << endl; 
	std::cout << "-c filename   chain definition file" << std::endl;
	std::cout << "-C filename   file with cuts" << std::endl;
	std::cout << "-D path       where to put the output" << std::endl;
	std::cout << "-f filename   single file instead of chain" << std::endl;
	cout << "-m            use MC information" << endl;
	std::cout << "-n integer    number of events to run on" << std::endl;
	std::cout << "-r class name which tree reader class to run" << std::endl;
	std::cout << "-s number     seed for random number generator" << std::endl;
	cout << "-S start     starting event number" << endl;
	std::cout << "-o filename   set output file" << std::endl;
	cout << "-v level      set verbosity level" << endl;
	std::cout << "-h            prints this message and exits" << std::endl;
	return 0;
    }
    if (!strcmp(argv[i],"-b"))  {blind      = atoi(argv[++i]); }                 // run blind?
    if (!strcmp(argv[i],"-c"))  {fileName   = string(argv[++i]); file = 0; }     // file with chain definition
    if (!strcmp(argv[i],"-C"))  {cutFile    = string(argv[++i]);           }     // file with cuts
    if (!strcmp(argv[i],"-D"))  {dirName    = string(argv[++i]);  dirspec = 1; } // where to put the output
    if (!strcmp(argv[i],"-f"))  {fileName   = string(argv[++i]); file = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-j"))  {jsonName   = string(argv[++i]); json = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-m"))  {isMC       = 1; }                               // use MC information?
    if (!strcmp(argv[i],"-n"))  {nevents    = atoi(argv[++i]); }                 // number of events to run 
    if (!strcmp(argv[i],"-r"))  {readerName = string(argv[++i]); }               // which tree reader class to run
    if (!strcmp(argv[i],"-s"))  {randomSeed = atoi(argv[++i]); }                 // set seed for random gen.
    if (!strcmp(argv[i],"-S"))  {start = atoi(argv[++i]); }                      // set start event number
    if (!strcmp(argv[i],"-o"))  {histfile   = TString(argv[++i]); }              // set output file
    if (!strcmp(argv[i],"-v"))  {verbose    = atoi(argv[++i]); }                 // set verbosity level
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
  treeReader01 *a = NULL;
  if (readerName == "massReader") a = new massReader(chain,TString(evtClassName));
  else if (readerName == "copyReader") a = new copyReader(chain,TString(evtClassName));
  else if (readerName == "genLevel") a = new genLevel(chain,TString(evtClassName));
  else if (readerName == "impactReader") a = new impactReader(chain,TString(evtClassName));
  else if (readerName == "pixelReader") a = new pixelReader(chain,TString(evtClassName));
  else {
    cout << "please provide a class name to instantiate" << endl;
  }


  if (a) {
    if (verbose > -99) {
      cout << "Setting verbosity to " << verbose << endl;
      a->setVerbosity(verbose); 
    }
    a->openHistFile(histfile); 
    a->readCuts(cutFile, 1);
    a->bookHist();
    if (json) {
      a->setJSONFile(jsonName.c_str()); 
      a->forceJSON();
    }

    if (isMC) {
      blind = 0; 
      a->setMC(1);
    } else {
      a->setMC(0); 
    }
    if (blind) a->runBlind();
    
    a->startAnalysis(); 
    a->loop(nevents, start);
    a->endAnalysis();
    a->closeHistFile(); 
  } else
    cerr << "Readerclass '" << readerName << "' not found" << endl;

  delete a; // so we can dump some information in the destructor
  
  return 0;
}
