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
#include "plotMisc.hh"

using namespace std;

// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  string progName  = argv[0]; 

  string sample("nada"), 
    dir("nada"), 
    cuts("nada"), 
    files("nada");
  bool restricted(false), remove(false), doUseBDT(true), smallTree(true);
  int year(2012), mode(255), Mode(-1), suffixMode(0), nevents(-1);   

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-cnc")){doUseBDT = false;}     
    if (!strcmp(argv[i], "-s"))  {restricted = true; sample = string(argv[++i]);}     
    if (!strcmp(argv[i], "-d"))  {dir   = string(argv[++i]);}     
    if (!strcmp(argv[i], "-c"))  {cuts  = string(argv[++i]);}     
    if (!strcmp(argv[i], "-f"))  {files = string(argv[++i]);}     
    if (!strcmp(argv[i], "-m"))  {mode  = atoi(argv[++i]);}     
    if (!strcmp(argv[i], "-M"))  {Mode  = atoi(argv[++i]);}     
    if (!strcmp(argv[i], "-n"))  {nevents = atoi(argv[++i]);}     
    if (!strcmp(argv[i], "-nst")) {smallTree = false;}     
    if (!strcmp(argv[i], "-x"))  {remove= true;}     
    if (!strcmp(argv[i], "-y"))  {year  = atoi(argv[++i]);}     
  }

  if (2012 == year) {
    if ("nada" == files) files = "anaBmm.v16-2012.files";
    if ("nada" == cuts)  cuts  = "2012";    
    if ("nada" == dir)   dir   = "2012";
  }

  if (2011 == year) {
    if ("nada" == files) files = "anaBmm.v16-2011.files";
    if ("nada" == cuts)  cuts  = "2011";    
    if ("nada" == dir)   dir   = "2011";
  }

  cout << "mode:       " << mode << endl;
  cout << "Mode:       " << Mode << endl;
  cout << " -m 1:       BDT " << endl;
  cout << " -m 2 -M 1:  results tree loops" << endl;
  cout << " -m 2 -M 2:  results finalization" << endl;
  cout << " -m 4 -M 1:  overlays filling " << endl;
  cout << " -m 4 -M 2:  overlays sbs/overlays " << endl;
  cout << " -m 4 -M 3:  overlays systematics " << endl;
  cout << " -m 8:       PU overlays " << endl;
  cout << " 32:         muon ID/trigger " << endl;
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
      rootfile = Form("anaBmm.plotBDT.%s.tex", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
      rootfile = Form("anaBmm.plotBDT.%s.root", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
    if (mode & 2) {
      rootfile = Form("anaBmm.plotResults.%s.tex", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
      rootfile = Form("anaBmm.plotResults.%s.root", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
    if ((mode & 4) && (1 == Mode)) {
      rootfile = Form("anaBmm.plotReducedOverlays.%s.tex", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
      rootfile = Form("anaBmm.plotReducedOverlays.%s.root", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
    if ((mode & 4) && (2 == Mode)) {
      rootfile = Form("anaBmm.plotReducedOverlaysSystematics.%s.tex", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
      rootfile = Form("anaBmm.plotReducedOverlaysSystematics.%s.root", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
    if (mode & 32) {
      rootfile = Form("anaBmm.plotEfficiencies.%s.tex", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
      rootfile = Form("anaBmm.plotEfficiencies.%s.root", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
    if (mode & 64) {
      rootfile = Form("anaBmm.plotMisc.%s.tex", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
      rootfile = Form("anaBmm.plotMisc.%s.root", cuts.c_str()); 
      cout << "Removing " << dir.c_str() << "/" <<  rootfile.c_str() << endl;
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));
    }
  }

  // -- BDT plots
  if (mode & 1 && doUseBDT) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotBDT a(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    if (Mode > -1) {
      a.makeAll(Mode);
    } else {
      a.makeAll(1); 
    } 
  }

  // -- results
  if (mode & 2) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotResults a(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    if (!doUseBDT) a.fDoUseBDT = false; 
    a.fSaveSmallTree = smallTree; 
    if (Mode > -1) {
      a.makeAll(Mode, nevents);
    } else {
      a.makeAll(1, nevents);
      
      gROOT->Clear();  gROOT->DeleteAll();
      plotResults b(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
      if (!doUseBDT) b.fDoUseBDT = false; 
      b.makeAll(2);
    } 
  }

  // -- overlays histogram filling
  if (mode & 4) {
    if (1 == Mode) {
      gROOT->Clear();  gROOT->DeleteAll();
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
	  gROOT->Clear();  gROOT->DeleteAll();
	  a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
	  if (!doUseBDT) a->fDoUseBDT = false; 
	  
	  a->makeSample(mode1, "Presel", chanlist[j], nevents);
	  a->makeSample(mode2, "Presel", chanlist[j], nevents);
	  
	  delete a;
	}
      }
    }

    // -- overlays histogram: sbs and overlay
    if (2 == Mode) {
      gROOT->Clear();  gROOT->DeleteAll();
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
	  gROOT->Clear();  gROOT->DeleteAll();
	  a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
	  if (!doUseBDT) a->fDoUseBDT = false; 
	  
	  a->makeOverlay(mode1, mode2, chanlist[j], "Presel"); 
	  
	  delete a;
	}
      }
      
    }
    

    // -- systematics
    if (3 == Mode) {
      gROOT->Clear();  gROOT->DeleteAll();
      plotReducedOverlays *a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
      if (!doUseBDT) a->fDoUseBDT = false; 
      a->allSystematics(); 
      delete a;
    }
  }


  // -- PU overlays
  if (mode & 8) {

    vector<string> chanlist; 
    chanlist.push_back("BLoPU"); 
    chanlist.push_back("BHiPU"); 
    chanlist.push_back("ELoPU"); 
    chanlist.push_back("EHiPU"); 
    
    chanlist.push_back("BClosePV"); 
    chanlist.push_back("BFarPV"); 
    chanlist.push_back("EClosePV"); 
    chanlist.push_back("EFarPV"); 

    vector<string> dolist; 
    if (restricted) {
      dolist.push_back(sample);
    } else {
      dolist.push_back("SgMc"); 
      dolist.push_back("NoMc"); 
    }

    plotReducedOverlays *a;
    if (1 == Mode) {
      gROOT->Clear();  gROOT->DeleteAll();
      
      for (unsigned int j = 0; j < chanlist.size(); ++j) {
	for (unsigned int i = 0; i < dolist.size(); ++i) {
	  gROOT->Clear();  gROOT->DeleteAll();
	  a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
	  if (!doUseBDT) a->fDoUseBDT = false; 
	  a->makeSample(dolist[i], "Presel", chanlist[j], nevents);
	  delete a;
	}
      }
    }
    
    if (2 == Mode) {
      string rootfile = Form("%s/anaBmm.plotReducedOverlaysSbs.%s.root", dir.c_str(), cuts.c_str()); 
      system(Form("/bin/rm -f %s/%s", dir.c_str(), rootfile.c_str()));

      rootfile = Form("%s/anaBmm.plotReducedOverlays.%s.root", dir.c_str(), cuts.c_str()); 
      string mode1, mode2; 
      // -- sbs
      gROOT->Clear();  gROOT->DeleteAll();
      a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
      for (unsigned int j = 0; j < chanlist.size(); ++j) {
	for (unsigned int i = 0; i < dolist.size(); ++i) {
	  if (!doUseBDT) a->fDoUseBDT = false; 
	  a->sbsSingleFile(rootfile, dolist[i], chanlist[j], "Presel"); 
	}
      }
      delete a; 
      
      gROOT->Clear();  gROOT->DeleteAll();
      a = new plotReducedOverlays(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
      for (unsigned int i = 0; i < dolist.size(); ++i) {
	if (!doUseBDT) a->fDoUseBDT = false; 
	rootfile = Form("%s/anaBmm.plotReducedOverlaysSbs.%s.root", dir.c_str(), cuts.c_str()); 
	a->overlay2Files(rootfile, dolist[i], rootfile, dolist[i], "BLoPU", "BHiPU", "Presel", "multichan"); 
	a->overlay2Files(rootfile, dolist[i], rootfile, dolist[i], "ELoPU", "EHiPU", "Presel", "multichan"); 

	a->overlay2Files(rootfile, dolist[i], rootfile, dolist[i], "BClosePV", "BFarPV", "Presel", "multichan"); 
	a->overlay2Files(rootfile, dolist[i], rootfile, dolist[i], "EClosePV", "EFarPV", "Presel", "multichan"); 
      }
      delete a;

    }
  }

  // -- muon id/trigger efficiency numbers
  if (mode & 32) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotEfficiencies *a = new plotEfficiencies(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    if (!doUseBDT) a->fDoUseBDT = false; 
    if (Mode > -1) {
      a->makeAll(Mode);
      delete a;
    } else {
      a->makeAll(1); 
      delete a;
      
      gROOT->Clear();  gROOT->DeleteAll();
      a = new plotEfficiencies(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
      if (!doUseBDT) a->fDoUseBDT = false; 
      a->makeAll(8); 
      delete a;

      gROOT->Clear();  gROOT->DeleteAll();
      a = new plotEfficiencies(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
      if (!doUseBDT) a->fDoUseBDT = false; 
      a->makeAll(16); 
      delete a;

      gROOT->Clear();  gROOT->DeleteAll();
      a = new plotEfficiencies(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
      if (!doUseBDT) a->fDoUseBDT = false; 
      a->makeAll(2); 
      delete a;

    }
  }

  // -- miscellaneous stuff
  if (mode & 64) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotMisc *a = new plotMisc(files.c_str(), dir.c_str(), cuts.c_str(), suffixMode);
    if (!doUseBDT) a->fDoUseBDT = false; 
    a->makeAll(); 
    delete a;
  }

  return 0;
}

