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
#include "plotEfficiencies.hh"

using namespace std;

// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  string progName  = argv[0]; 

  string sample("nada"), 
    dir("nada"), 
    cuts("nada"), 
    files("nada");
  bool restricted(false), remove(false), doUseBDT(true);
  int year(2012), mode(255), suffixMode(0);   

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-cnc")){doUseBDT = false;}     
    if (!strcmp(argv[i], "-s"))  {restricted = true; sample = string(argv[++i]);}     
    if (!strcmp(argv[i], "-d"))  {dir   = string(argv[++i]);}     
    if (!strcmp(argv[i], "-c"))  {cuts  = string(argv[++i]);}     
    if (!strcmp(argv[i], "-f"))  {files = string(argv[++i]);}     
    if (!strcmp(argv[i], "-m"))  {mode  = atoi(argv[++i]);}     
    if (!strcmp(argv[i], "-x"))  {remove= true;}     
    if (!strcmp(argv[i], "-y"))  {year  = atoi(argv[++i]);}     
  }

  if (2012 == year) {
    if ("nada" == files) files = "ul-anaBmm.2012.files";
    if ("nada" == cuts)  cuts  = "2012";    
    if ("nada" == dir)   dir   = "2012";
  }

  if (2011 == year) {
    if ("nada" == files) files = "ul-anaBmm.2011.files";
    if ("nada" == cuts)  cuts  = "2011";    
    if ("nada" == dir)   dir   = "2011";
  }

  cout << "mode:       " << mode << endl;
  cout << "  1:  BDT " << endl;
  cout << "  2:  results " << endl;
  cout << "  4:  overlays filling " << endl;
  cout << "  8:  overlays sbs/overlays " << endl;
  cout << " 16:  overlays systematics " << endl;
  cout << " 32:  muon ID/trigger " << endl;
  cout << "dir:        " << dir << endl;
  cout << "cuts:       " << cuts << endl;
  cout << "files:      " << files << endl;
  cout << "mode:       " << mode << endl;
  cout << "restricted: " << restricted  << (restricted?(string(" to ") + sample):"") << endl;
  cout << "remove:     " << remove << endl;

  // -- cleanup 
  if (remove) {
    string rootfile(Form("anaBmm.plotBDT.%s.root", cuts.c_str())); 
    if (mode & 1) {
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
    rootfile = Form("anaBmm.plotResults.%s.root", cuts.c_str()); 
    if (mode & 2) {
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
    rootfile = Form("anaBmm.plotReducedOverlays.%s.root", cuts.c_str()); 
    if (mode & 4) {
      rootfile = Form("anaBmm.plotReducedOverlays.%s.root", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
    if (mode & 8) {
      rootfile = Form("anaBmm.plotReducedOverlaysSystematics.%s.root", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
  }

  // -- BDT plots
  if (mode & 1 && doUseBDT) {
    plotBDT a(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    a.makeAll(1);
  } 

  // -- results
  if (mode & 2) {
    plotResults a(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    if (!doUseBDT) a.fDoUseBDT = false; 
    a.makeAll(1);
  }
  if (mode & 2) {
    plotResults b(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    if (!doUseBDT) b.fDoUseBDT = false; 
    b.makeAll(2);
  } 

  // -- overlays histogram filling
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
    for (unsigned int j = 0; j < chanlist.size(); ++j) {
      for (unsigned int i = 0; i < dolist.size(); ++i) {
	mode1 = dolist[i] + string("Data");
	mode2 = dolist[i] + string("Mc");
	a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
	if (!doUseBDT) a->fDoUseBDT = false; 
    
  	a->makeSample(mode1, "Presel", chanlist[j]);
  	a->makeSample(mode2, "Presel", chanlist[j]);

	delete a;
      }
    }
  }

  // -- overlays histogram: sbs and overlay
  if (mode & 8) {
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
    for (unsigned int j = 0; j < chanlist.size(); ++j) {
      for (unsigned int i = 0; i < dolist.size(); ++i) {
	mode1 = dolist[i] + string("Data");
	mode2 = dolist[i] + string("Mc");
	a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
	if (!doUseBDT) a->fDoUseBDT = false; 

	a->makeOverlay(mode1, mode2, chanlist[j], "Presel"); 

	delete a;
      }
    }
    
  }


  // -- systematics
  if (mode & 16) {
    plotReducedOverlays *a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    if (!doUseBDT) a->fDoUseBDT = false; 
    a->allSystematics(); 
    delete a;
  }

  // -- muon id/trigger efficiency numbers
  if (mode & 32) {
    plotEfficiencies *a = new plotEfficiencies(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    if (!doUseBDT) a->fDoUseBDT = false; 
    a->makeAll(1); 
    delete a;
  }

  return 0;
}

