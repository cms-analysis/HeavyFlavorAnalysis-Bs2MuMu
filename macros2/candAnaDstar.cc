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

    if (300054 == pCand->fType) {
      tm = false; 
      tm = truthMatch(pCand); 
      if (tm) {
	fpEvt->dumpGenBlock(); 
	// 	cout << " -> " << pCand->fType;
	// 	cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
	// 	cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
	// 	cout << " mass: " << pCand->fMass << endl;

	for (int ih = 0; ih < fpEvt->nCands(); ++ih) {
	  pC = fpEvt->getCand(ih);
	  if (54 == pC->fType) {
	    cout << "DUMP HFTruthCandidate" << endl;
	    dumpHFTruthCand(pC); 
	  }
	}
	cout << "DUMP HFDstarCandidate" << endl;
	dumpHFDstarCand(pCand); 

      }
      mdstar = pCand->fMass; 

      // -- D0 
      if (pCand->fDau1 < 0) {
	if (pCand->fDau2 < 0) {
	  cout << "pCan->fDauX = -1!!!" << endl;
	}
	continue;
      }
      pC = fpEvt->getCand(pCand->fDau1);
      if (tm) {
	// 	cout << "  fDau1 " << pCand->fDau1;
	// 	cout << " pC " << pC;
	// 	cout << " type: " << pC->fType;
	// 	cout << " sig: " << pC->fSig1 << " .. " << pC->fSig2;
	// 	cout << ": dau " << pC->fDau1 << " .. " << pC->fDau2;
	// 	cout << " mass: " << pC->fMass << endl;
      }
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
void candAnaDstar::dumpHFTruthCand(TAnaCand *pC) {
  TAnaTrack *pT(0); 
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    pT->dump(); 
  }
}


// ----------------------------------------------------------------------
void candAnaDstar::dumpHFDstarCand(TAnaCand *pC) {
  TAnaTrack *pT(0); 
  // -- D0 daughters
  if (pC->fDau1 < 0) {
    cout << "XXXXXXXXX cannot get daughter cand of " << pC->fType << endl;
    return;
  }
  TAnaCand *pD = fpEvt->getCand(pC->fDau1); 
  for (int id = pD->fSig1; id <= pD->fSig2; ++id) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(id)->fIndex);
    pT->dump(); 
  }
  // -- slow pion
  pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pC->fSig1)->fIndex); 
  pT->dump(); 

}


// ----------------------------------------------------------------------
// -- works ONLY for this dedicated decay mode with the seq vtx fitter!
bool candAnaDstar::truthMatch(TAnaCand *pCand) {

  // -- check slow pion
  TAnaTrack *pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex);
  if (pT->fGenIndex < 0) return false; 
  TGenCand  *pG = fpEvt->getGenCand(fpEvt->getRecTrack(pT->fIndex)->fGenIndex); 
  if (0 == pG) return false;
  if (211 != TMath::Abs(pG->fID)) return false;
  pG = fpEvt->getGenCand(pG->fMom1); 
  if ((0 == pG) || 413 != pG->fID) return false;
  if (pG->fDau2 - pG->fDau1 > 1) return false; 
  if (fVerbose > 2) cout << "slow pion " << pT->fIndex << " OK, fGenIndex = " << pT->fGenIndex << endl;

  // -- check D0 
  if (pCand->fDau1 < 0) return false;
  TAnaCand *pC = fpEvt->getCand(pCand->fDau1);
  
  int type(0); 
  for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
    pT = fpEvt->getSigTrack(id);
    type = pT->fMCID;
    pT = fpEvt->getRecTrack(pT->fIndex);
    pG = fpEvt->getGenCand(pT->fGenIndex);
    if (fVerbose > 2) cout << "dau cand sigtrack " << id << " with type = " << type << " and gen ID = " << pG->fID << endl;
    if (TMath::Abs(type) != TMath::Abs(pG->fID)) return false;
  }

  if (pG->fDau2 - pG->fDau1 > 1) return false; 
  

  if (fVerbose > 0)   cout << "===> truth matching OK" << endl;
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

  TH1 *h = new TH1D("mds", "m(dstar)", 60, 1.9, 2.3);
  h = new TH1D("mdz", "m(d0)", 60, 1.8, 1.92);
  h = new TH1D("dm", "delta(m)", 60, 0.13, 0.16);

  h = new TH1D("all_mds", "m(dstar)", 60, 1.9, 2.3);
  h = new TH1D("all_mdz", "m(d0)", 60, 1.8, 1.92);
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
