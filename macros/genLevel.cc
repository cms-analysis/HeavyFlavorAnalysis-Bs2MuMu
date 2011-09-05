#include "genLevel.hh"
#include "../interface/HFMasses.hh"

#include "TRandom.h"
#include <cmath>
#include <map>

using std::cout;
using std::endl;
using std::string;
using std::map;

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

  printBdecays();
  //  bbbarCrossSection();
  // -- initialize all variables
  //  initVariables(); 

}



// ----------------------------------------------------------------------
void genLevel::endAnalysis() {

  double ftilde(0.11); 
  double f = 2.*ftilde*(1-ftilde) + ftilde*ftilde; 
  double n0 = ((TH1D*)fpHistFile->Get("e111"))->GetSumOfWeights(); 
  double n1 = ((TH1D*)fpHistFile->Get("e112"))->GetSumOfWeights();
  double c0 = ((TH1D*)fpHistFile->Get("h110"))->GetSumOfWeights(); 
  double c1 = ((TH1D*)fpHistFile->Get("h112"))->GetSumOfWeights();

  cout << "XSECTION:            " << XSECTION << endl;
  cout << "NTOTAL:              " << NTOTAL << endl;
  cout << "b->mu events:        " << n0 << endl;
  cout << "b->mu decays:        " << c0 << endl;
  cout << "b->mu events in acc: " << n1 << endl;
  cout << "b->mu decays in acc: " << c1 << endl;

  cout << "bbbar xsection:       " << n0/(f*NTOTAL)*XSECTION << endl;
  cout << "bbbar xsection:       " << c0/(2.*ftilde*NTOTAL)*XSECTION << endl;
  cout << "b-mu  xsection:       " << n0/(NTOTAL)*XSECTION << endl;
  cout << "b-mu in acc xsection: " << n1/(NTOTAL)*XSECTION << endl;
}

// ----------------------------------------------------------------------
void genLevel::printBdecays() {
  
  TGenCand *pCand, *pD; 
  int muType(0), evtType(0), nevt(0), bevt(0), bacc(0); 
  bool acc(false);
  double pt(0.), eta(0.); 
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);

    if (TMath::Abs(pCand->fID) == 531 || TMath::Abs(pCand->fID) == 511) {
      cout << "--- Event " << fEvent << " ------------------------------------------------------------------" << endl;
      pCand->dump();
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = fpEvt->getGenCand(iD);
	pD->dump();
	cout << "pT: " << pD->fP.Perp() << " eta: " << pD->fP.Eta() << endl;
      }
    }
  }


  TAnaCand *pC; 
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pC = fpEvt->getCand(iC);
    cout << "typ " << pC->fType << endl;
  }

  TAnaTrack *pT;       
  cout << "rec tracks: " << fpEvt->nRecTracks() << endl;
  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {
    pT = fpEvt->getRecTrack(it); 
    pT->dump();
  }


  TString a; 
  int ps(0); 
  bool result(false), wasRun(false), error(false); 

  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i]; 
    ps = fpEvt->fHLTPrescale[i]; 
    wasRun = fpEvt->fHLTWasRun[i]; 
    result = fpEvt->fHLTResult[i]; 
    error  = fpEvt->fHLTError[i]; 

    if (wasRun && result) {
      cout << a << endl;
    }      
  }


}


// ----------------------------------------------------------------------
void genLevel::bbbarCrossSection() {

  TGenCand *pCand; 
  int muType(0), evtType(0), nevt(0), bevt(0), bacc(0); 
  bool acc(false);
  double pt(0.), eta(0.); 
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);

    if (TMath::Abs(pCand->fID) == 13) {
      //      cout << "muon at " << iC << endl;
      muType = muonType(pCand); 
      pt = pCand->fP.Perp();
      eta= pCand->fP.Eta();

      evtType |= muType; 

      if (pt > 5 && eta < 2.1 && eta > -2.1) {
	acc = true; 
      } else {
	acc = false;
      }

      ((TH1D*)fpHistFile->Get("h1"))->Fill(muType); 
      if (muType > 0) ((TH1D*)fpHistFile->Get("h100"))->Fill(pt); 
      if (muType ==1) ((TH1D*)fpHistFile->Get("h110"))->Fill(pt); 
      if (muType ==1 && acc) ((TH1D*)fpHistFile->Get("h112"))->Fill(pt); 

      if (muType > 0)         ++nevt;
      if (muType == 1)        ++bevt;
      if (muType == 1 && acc) ++bacc;
      
    }

  }

  ((TH1D*)fpHistFile->Get("e1"))->Fill(evtType); 
  if (evtType > 0) ((TH1D*)fpHistFile->Get("e100"))->Fill(pt); 
  if (evtType ==1) ((TH1D*)fpHistFile->Get("e110"))->Fill(pt); 
  if (bevt    > 0) ((TH1D*)fpHistFile->Get("e111"))->Fill(pt); 
  if (bacc    > 0) ((TH1D*)fpHistFile->Get("e112"))->Fill(pt); 
  
}


// ----------------------------------------------------------------------
int genLevel::muonType(TGenCand *pCand) {
  int id(999); 
  int rest(0), ganz(0); 
  int momI = pCand->fMom1;
  TGenCand *pMom;
  bool foundT(false), foundC(false), foundB(false), light(false); 

  int result(1024); 

  string str("");
  str += Form(" %i", pCand->fID); 
  
  while (momI > 1 && momI < fpEvt->nGenCands()) {
    pMom = fpEvt->getGenCand(momI); 
    momI = pMom->fMom1;

    id  = TMath::Abs(pMom->fID);
    str += Form(" %i", pMom->fID);

    rest =  id%1000;
    ganz =  id/1000;
    
    // -- hit a string
    if (rest == 92) {
      break;
    }

    // -- hit a quark
    if (rest < 6) {
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

//   if (result > -1) {
//     cout << fEvent << " " << str << " --> " << result << endl;
//   }

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
  h = new TH1D("h1",   "mu type", 20, 0., 20.);
  h = new TH1D("h100", "mu pt 0", 100, 0., 20.);
  h = new TH1D("h110", "mu pt ==1", 100, 0., 20.);
  h = new TH1D("h112", "mu pt ==1", 100, 0., 20.);

  h = new TH1D("e1",  "evt type", 20, 0., 20.);
  h = new TH1D("e100", "evt pt 0", 100, 0., 20.);
  h = new TH1D("e110", "evt pt ==1", 100, 0., 20.);
  h = new TH1D("e111", "evt pt ==1", 100, 0., 20.);
  h = new TH1D("e112", "evt acc ==1", 100, 0., 20.);

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
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "NTOTAL")) {
      NTOTAL = int(CutValue); ok = 1;
      if (dump) cout << "NTOTAL:           " << NTOTAL << endl;
    }

    if (!strcmp(CutName, "XSECTION")) {
      XSECTION = double(CutValue); ok = 1;
      if (dump) cout << "XSECTION:           " << XSECTION << endl;
    }

    if (!ok) cout << "==> genLevel: ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}
