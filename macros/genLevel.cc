#include "genLevel.hh"
#include "../interface/HFMasses.hh"

#include "TRandom.h"
#include <cmath>

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/genLevel.default.cuts
//           ./runTreeReaders -f test.root 
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
genLevel::genLevel(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> genLevel: constructor..." << endl;
}

// ----------------------------------------------------------------------
genLevel::~genLevel() {
  cout << "==> genLevel: destructor..." << endl;

}

// ----------------------------------------------------------------------
void genLevel::startAnalysis() {
  cout << "==> genLevel: Starting analysis..." << endl;
}


// ----------------------------------------------------------------------
void genLevel::eventProcessing() {

  TGenCand *pCand; 
  TH1D *h; 
  int muType; 
  double pt(0.), eta(0.); 
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);

    if (TMath::Abs(pCand->fID) == 13) {
      muType = muonType(pCand); 
      pt = pCand->fP.Perp();
      eta= pCand->fP.Eta();
      ((TH1D*)fpHistFile->Get("h1"))->Fill(muType); 

      if (muType > 0) ((TH1D*)fpHistFile->Get("h100"))->Fill(pt); 
      if (muType & 1) ((TH1D*)fpHistFile->Get("h101"))->Fill(pt); 
      if (muType & 2) ((TH1D*)fpHistFile->Get("h102"))->Fill(pt); 
      if (muType & 4) ((TH1D*)fpHistFile->Get("h104"))->Fill(pt); 
      if (muType & 8) ((TH1D*)fpHistFile->Get("h108"))->Fill(pt); 
	

      //      cout << muType << ": "; pCand->dump();
      //       if (muType >= 4) {
      // 	cout << "----------------------------------------------------------------------" << endl;
      // 	fpEvt->dumpGenBlock();
      // 	cout << "----------------------------------------------------------------------" << endl;
      //       }
    }

  }


  // -- initialize all variables
  initVariables(); 

}


// ----------------------------------------------------------------------
int genLevel::muonType(TGenCand *pCand) {
  int id(999), cnt(0); 
  int rest(0), ganz(0); 
  int momI = pCand->fMom1;
  TGenCand *pMom;
  bool foundT(false), foundC(false), foundB(false), light(false); 

  int result(1024); 

  while (momI > 1 && momI < fpEvt->nGenCands()) {
    pMom = fpEvt->getGenCand(momI); 
    momI = pMom->fMom1;

    id  = TMath::Abs(pMom->fID);

    rest =  id%1000;
    ganz =  id/1000;
    
    // -- hit a string
    if (rest == 92) {
      cout << "?????????????????????????????????????????????????????????????????????????" << endl;
      break;
    }

    // -- hit a quark
    if (rest == 5) {
      cout << "?????????????????????????????????????????????????????????????????????????" << endl;
      break;
    }

    if (ganz == 5) {
      foundB = 1; 
      break;
    }

    if (ganz == 4) {
      foundC = 1; 
      break;
    }
    
    if (rest == 15) {
      foundT = 1; 
    }

    if (15 < rest && rest < 300) {
      light = 1; 
    }

    if (rest > 399 && rest < 499) {
      foundC = 1; 
    }

    if (rest > 499 && rest < 599) {
      foundB = 1; 
      break;
    }
  }
  
  result = 0; 

  if (foundB) result += 1; 
  if (foundC) result += 2; 
  if (foundT) result += 4; 
  if (light)  result += 8; 
  return result; 
}

// ----------------------------------------------------------------------
void genLevel::initVariables() {

  

}


// ----------------------------------------------------------------------
void genLevel::fillHist() {


}

// ---------------------------------------------------------------------- 
void genLevel::bookHist() {
  cout << "==> genLevel: bookHist " << endl;
  
  TH1D *h; 
  h = new TH1D("h1", "mu type", 10, 0., 10.);

  h = new TH1D("h100", "mu pt 1", 100, 0., 20.);
  h = new TH1D("h101", "mu pt 1", 100, 0., 20.);
  h = new TH1D("h102", "mu pt 2", 100, 0., 20.);
  h = new TH1D("h104", "mu pt 4", 100, 0., 20.);
  h = new TH1D("h108", "mu pt 8", 100, 0., 20.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",     &fRun,     "run/I");

}

// ----------------------------------------------------------------------
void genLevel::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> genLevel: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> genLevel: Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  int ibin; 
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      if (dump) cout << "TYPE:           " << TYPE << endl;
    }

    if (!ok) cout << "==> genLevel: ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}
