#include "bmm2Reader.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

#include "candAna.hh"
#include "candAnaMuMu.hh"
#include "candAnaBu2JpsiK.hh"
#include "candAnaBs2JpsiPhi.hh"

using namespace std;

// ----------------------------------------------------------------------
bmm2Reader::bmm2Reader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> bmm2Reader: constructor..." << endl;
  fVerbose = 0; 
}


// ----------------------------------------------------------------------
bmm2Reader::~bmm2Reader() {
  cout << "==> bmm2Reader: destructor..." << endl;
}


// ----------------------------------------------------------------------
void bmm2Reader::startAnalysis() {
  cout << "==> bmm2Reader: fVerbose = " << fVerbose << endl;
  fpJSON = new JSON(JSONFILE.c_str(), (fVerbose!=0?1:0)); 
}


// ----------------------------------------------------------------------
void bmm2Reader::eventProcessing() {
  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(0); 

  bool json = false;  
   
  if (fIsMC) {
    json = 1; 
  } else {
    json = fpJSON->good(fRun, fLS); 
    if (fVerbose > 1 && !json) {
      cout << "JSON = 0 for run = " << fRun << " and LS = " << fLS << endl;
    }
  }
  
  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    //    cout << "  calling " << lCandAnalysis[i]->fName << " analysis()" << endl;

    lCandAnalysis[i]->fJSON  = json; 
    lCandAnalysis[i]->evtAnalysis(fpEvt);
  }

}


// ----------------------------------------------------------------------
void bmm2Reader::bookHist() {
  fpHistFile->cd();
  TH1D *h;
  h = new TH1D("monEvents", "monEvents", 10, 0., 10.);

  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    //    cout << "  calling " << lCandAnalysis[i]->fName << " bookHist()" << endl;
    lCandAnalysis[i]->bookHist();
  }

}


// ----------------------------------------------------------------------
void bmm2Reader::readCuts(TString filename, int dump) {
  if (dump) cout << "==> bmm2Reader: Reading " << filename << " for classes setup" << endl;
  
  ifstream is(filename.Data());
  char buffer[1000]; 
  char className[200], cutFile[200]; 
  while (is.getline(buffer, 1000, '\n')) {
    sscanf(buffer, "%s %s", className, cutFile);

    // -- set up candidate analyzer classes
    if (!strcmp(className, "candAnaMuMu")) {
      candAna *a = new candAnaMuMu(this, "candAnaMuMu", cutFile); 
      a->fVerbose = fVerbose; 
      a->BLIND = BLIND; 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaBu2JpsiK")) {
      candAna *a = new candAnaBu2JpsiK(this, "candAnaBu2JpsiK", cutFile); 
      a->fVerbose = fVerbose; 
      lCandAnalysis.push_back(a); 
    }

    if (!strcmp(className, "candAnaBs2JpsiPhi")) {
      candAna *a = new candAnaBs2JpsiPhi(this, "candAnaBs2JpsiPhi", cutFile); 
      a->fVerbose = fVerbose; 
      lCandAnalysis.push_back(a); 
    }


    // -- all the rest ...
    if (!strcmp(className, "JSON")) {
      char json[1000];
      sscanf(buffer, "%s %s", className, json);
      JSONFILE = string(json);
      if (dump) cout << "JSON FILE:           " << JSONFILE << endl;
    }

    if (!strcmp(className, "ptSgMUID")) {
      char name[1000];
      sscanf(buffer, "%s %s", className, name);
      ptSgMUID = new PidTable(name); 
      if (dump) cout << "Seagulls MUID:           " << name << endl;
    }

    if (!strcmp(className, "ptCbMUID")) {
      char name[1000];
      sscanf(buffer, "%s %s", className, name);
      ptCbMUID = new PidTable(name); 
      if (dump) cout << "Cowboys MUID:           " << name << endl;
    }

    if (!strcmp(className, "ptSgMUT1")) {
      char name[1000];
      sscanf(buffer, "%s %s", className, name);
      ptSgMUT1 = new PidTable(name); 
      if (dump) cout << "Seagulls MUT1:           " << name << endl;
    }

    if (!strcmp(className, "ptCbMUT1")) {
      char name[1000];
      sscanf(buffer, "%s %s", className, name);
      ptCbMUT1 = new PidTable(name); 
      if (dump) cout << "Cowboys MUT1:           " << name << endl;
    }

    if (!strcmp(className, "ptSgMUT2")) {
      char name[1000];
      sscanf(buffer, "%s %s", className, name);
      ptSgMUT2 = new PidTable(name); 
      if (dump) cout << "Seagulls MUT2:           " << name << endl;
    }

    if (!strcmp(className, "ptCbMUT2")) {
      char name[1000];
      sscanf(buffer, "%s %s", className, name);
      ptCbMUT2 = new PidTable(name); 
      if (dump) cout << "Cowboys MUT2:           " << name << endl;
    }

  }

}


