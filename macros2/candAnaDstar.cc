#include "candAnaDstar.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaDstar::candAnaDstar(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaDstar: name = " << name << ", reading cutsfile " << cutsFile << endl;

  readCuts(cutsFile, 1); 

}


// ----------------------------------------------------------------------
candAnaDstar::~candAnaDstar() {
  cout << "==> candAnaDstar: destructor..." << endl;
}

// ----------------------------------------------------------------------
void candAnaDstar::evtAnalysis(TAna01Event *evt) {
  fpEvt = evt; 
  
  TAnaCand *pCand(0), *pC(0);

  // -- loop over all seq vtx fit candidates for D*
  double mdstar(0.), mdz(0.), dm(0.); 
  bool tm(false); 
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    tm = false; 
    if (300054 == pCand->fType) {
      tm = truthMatch(pCand); 
      //       cout << " -> " << pCand->fType;
      //       cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
      //       cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
      //       cout << " mass: " << pCand->fMass << endl;
      mdstar = pCand->fMass; 
      // -- D0 
      if (pCand->fDau1 < 0) continue;
      pC = fpEvt->getCand(pCand->fDau1);
      cout << "  fDau1 " << pCand->fDau1;
      cout << " pC " << pC;
      cout << " type: " << pC->fType;
      cout << " sig: " << pC->fSig1 << " .. " << pC->fSig2;
      cout << ": dau " << pC->fDau1 << " .. " << pC->fDau2;
      cout << " mass: " << pC->fMass << endl;
      mdz = pC->fMass;
      dm = mdstar - mdz;
      if (tm) {
	((TH1D*)fHistDir->Get("mds"))->Fill(mdstar);
	((TH1D*)fHistDir->Get("mdz"))->Fill(mdz);
	((TH1D*)fHistDir->Get("dm"))->Fill(dm);
      } 
      ((TH1D*)fHistDir->Get("all_mds"))->Fill(mdstar);
      ((TH1D*)fHistDir->Get("all_mdz"))->Fill(mdz);
      ((TH1D*)fHistDir->Get("all_dm"))->Fill(dm);
    }
  }
  //  cout << fpCand->fType << " -> mass = " << fpCand->fMass << endl;
  
  
}

// ----------------------------------------------------------------------
// -- works ONLY for this dedicated decay mode with the seq vtx fitter!
bool candAnaDstar::truthMatch(TAnaCand *pCand) {

  // -- check slow pion
  TAnaTrack *pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex);
  if (fpEvt->getRecTrack(pT->fIndex)->fGenIndex < 0) return false; 
  TGenCand  *pG = fpEvt->getGenCand(fpEvt->getRecTrack(pT->fIndex)->fGenIndex); 
  if (0 == pG) return false;
  if (211 != TMath::Abs(pG->fID)) return false;
  pG = fpEvt->getGenCand(pG->fMom1); 
  if ((0 == pG) || 413 != pG->fID) return false;

  // -- check D0 
  if (pCand->fDau1 < 0) return false;
  TAnaCand *pC = fpEvt->getCand(pCand->fDau1);
  
  int type(0); 
  for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
    pT = fpEvt->getSigTrack(id);
    type = pT->fMCID;
    pT = fpEvt->getRecTrack(pT->fIndex);
    pG = fpEvt->getGenCand(pT->fGenIndex);
    //    cout << "dau cand sigtrack " << id << " with type = " << type << " and gen ID = " << pG->fID << endl;
    if (TMath::Abs(type) != TMath::Abs(pG->fID)) return false;
  }

  return true; 
}
  

// ----------------------------------------------------------------------
void candAnaDstar::candAnalysis() {

  if (0 == fpCand) return;

  //  candAna::candAnalysis();

}

// ----------------------------------------------------------------------
void candAnaDstar::moreBasicCuts() {
  cout << "   candAnaDstar: more basic cuts" << endl;
}


// ----------------------------------------------------------------------
void candAnaDstar::bookHist() {
  //  candAna::bookHist();
  cout << "==>candAnaDstar: bookHist" << endl;

  fHistDir->cd();

  TH1 *h = new TH1D("mds", "m(dstar)", 60, 1.8, 2.4);
  h = new TH1D("mdz", "m(d0)", 60, 1.5, 2.1);
  h = new TH1D("dm", "delta(m)", 60, 0.13, 0.16);

  h = new TH1D("all_mds", "m(dstar)", 60, 1.8, 2.4);
  h = new TH1D("all_mdz", "m(d0)", 60, 1.5, 2.1);
  h = new TH1D("all_dm", "delta(m)", 60, 0.13, 0.16);
}


// ----------------------------------------------------------------------
void candAnaDstar::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaDstar: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;
  int ok(0);

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

   
  }

}
