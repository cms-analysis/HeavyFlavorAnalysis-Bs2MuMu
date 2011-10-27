#include "bmm2Reader.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

#include "candAna.hh"
#include "candAnaMuMu.hh"
#include "candAnaBu2JpsiKp.hh"
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
  cout << "==> bmm2Reader: setup PidTables and JSON file" << endl;

  fpMuonID = new PidTable("../macros/pidtables/110606/H2D_MuonIDEfficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 
  fpMuonTr1 = new PidTable("../macros/pidtables/110606/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 
  fpMuonTr2 = new PidTable("../macros/pidtables/110606/H2D_L3Efficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 

  fpJSON = new JSON(JSONFILE.c_str()); 
  
  cout << "==> bmm2Reader: start analysis" << endl;

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
    cout << "  calling " << lCandAnalysis[i]->fName << " bookHist()" << endl;
    lCandAnalysis[i]->bookHist();
  }

}


// ----------------------------------------------------------------------
void bmm2Reader::readCuts(TString filename, int dump) {
  cout << "==> bmm2Reader: readCuts" << endl;
  
  if (dump) cout << "==> bmm2Reader: Reading " << filename << " for classes setup" << endl;
  
  ifstream is(filename.Data());
  char buffer[1000]; 
  char className[200], cutFile[200]; 
  while (is.getline(buffer, 1000, '\n')) {
    sscanf(buffer, "%s %s", className, cutFile);

    if (!strcmp(className, "candAnaMuMu")) {
      candAna *a = new candAnaMuMu(this, "candAnaMuMu", cutFile); 
      a->fVerbose = fVerbose; 
      a->BLIND = BLIND; 
      lCandAnalysis.push_back(a); 
      if (dump) cout << "candAnaMuMu with          " << cutFile << endl;
    }

    if (!strcmp(className, "candAnaBu2JpsiKp")) {
      candAna *a = new candAnaBu2JpsiKp(this, "candAnaBu2JpsiKp", cutFile); 
      a->fVerbose = fVerbose; 
      lCandAnalysis.push_back(a); 
      if (dump) cout << "candAnaBu2JpsiKp with     " << cutFile << endl;
    }

    if (!strcmp(className, "candAnaBs2JpsiPhi")) {
      candAna *a = new candAnaBs2JpsiPhi(this, "candAnaBs2JpsiPhi", cutFile); 
      a->fVerbose = fVerbose; 
      lCandAnalysis.push_back(a); 
      if (dump) cout << "candAnaBs2JpsiPhi with     " << cutFile << endl;
    }

    if (!strcmp(className, "JSON")) {
      char json[1000];
      sscanf(buffer, "%s %s", className, json);
      JSONFILE = string(json);
      if (dump) cout << "JSON FILE:           " << JSONFILE << endl;
    }

  }

}


