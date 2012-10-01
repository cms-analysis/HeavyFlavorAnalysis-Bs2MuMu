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

#include "plotBDT.hh"
#include "plotResults.hh"
#include "plotReducedOverlays.hh"

using namespace std;

// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  string progName  = argv[0]; 

  string sample("Cs"), 
    dir("default"), 
    cuts("default"), 
    files("anaBmm.overlays.files");
  bool restricted(false);
  int mode(0); 

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-s"))  {restricted = true; sample = string(argv[++i]);}     
    if (!strcmp(argv[i], "-d"))  {dir   = string(argv[++i]);}     
    if (!strcmp(argv[i], "-c"))  {cuts  = string(argv[++i]);}     
    if (!strcmp(argv[i], "-f"))  {files = string(argv[++i]);}     
    if (!strcmp(argv[i], "-m"))  {mode  = atoi(argv[++i]);}     
  }

  // -- BDT plots
  if (mode & 1) {
    plotBDT a(files.c_str(), dir.c_str(), cuts.c_str());
    a.makeAll(1);
  } 

  // -- results
  if (mode & 2) {
    plotResults a(files.c_str(), dir.c_str(), cuts.c_str());
    a.makeAll(1);
  } 

  // -- overlays
  if (mode & 4) {
    plotReducedOverlays *a;
    
    vector<string> chanlist; 
    chanlist.push_back("B"); 
    chanlist.push_back("E"); 
    
    vector<string> dolist; 
    
    if (restricted) {
      dolist.push_back(sample);
    } else {
      dolist.push_back("Cs"); 
      dolist.push_back("Sg"); 
      dolist.push_back("No"); 
    }
    
    string mode1, mode2; 
    for (int j = 0; j < chanlist.size(); ++j) {
      for (int i = 0; i < dolist.size(); ++i) {
	mode1 = dolist[i] + string("Data");
	mode2 = dolist[i] + string("Mc");
	a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str());
	a->makeSampleOverlay(mode1, mode2, chanlist[j], "Presel");
	//      a->makeOverlay(mode1, mode2, chanlist[j], "Presel");
	delete a;
      }
    }
    
    
    a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str());
    a->allSystematics(); 
    delete a;
  }

  return 0;
}

