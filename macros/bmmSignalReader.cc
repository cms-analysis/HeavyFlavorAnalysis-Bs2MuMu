#include "bmmSignalReader.hh"
#include "TRandom.h"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using std::string;
using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
bmmSignalReader::bmmSignalReader(TChain *tree, TString evtClassName): bmmReader(tree, evtClassName) {
  cout << "==> bmmSignalReader: constructor..." << endl;
}

// ----------------------------------------------------------------------
bmmSignalReader::~bmmSignalReader() {
  cout << "==> bmmSignalReader: destructor..." << endl;

}

// ----------------------------------------------------------------------
void bmmSignalReader::startAnalysis() {
  cout << "==> bmmSignalReader: Starting analysis ... " << (BLIND? "blinded": "NOT blinded") << endl;
}


// ----------------------------------------------------------------------
void bmmSignalReader::eventProcessing() {
  bmmReader::eventProcessing();
}


// ----------------------------------------------------------------------
void bmmSignalReader::initVariables() {
  bmmReader::initVariables();
 
}


// ----------------------------------------------------------------------
void bmmSignalReader::MCKinematics() {
  TGenCand *pC, *pM1, *pM2, *pB; 
  int nphotons(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (531 == TMath::Abs(pC->fID)) {
      pB = pC;
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenCand(id); 
	if (13 == TMath::Abs(pC->fID)) {
	  if (0 == pM1) {
	    pM1 = fpEvt->getGenCand(id); 
	  } else {
	    pM2 = fpEvt->getGenCand(id); 
	  }
	}
      }
      if (0 != pM1 && 0 != pM2) {
	goodMatch = true; 
	nphotons = pB->fDau2 - pB->fDau1 - 1; 
	break;
      }
    }
  }
  
  if (fVerbose > 4) {
    cout << "----------------------------------------------------------------------" << endl;
    if (goodMatch) {
      pB->dump(); 
      pM1->dump(); 
      pM2->dump();
    } else {
      cout << "no tm decay found" << endl;
    }
  }

  if (!goodMatch) {
    fGoodMCKinematics = false; 
    if (fVerbose > 2) cout << "No matched signal decay found" << endl;
    return;
  }

  // -- hard coded ?! FIXME
  if (pM1->fP.Perp() < 2.5) fGoodMCKinematics = false;  
  if (pM2->fP.Perp() < 2.5) fGoodMCKinematics = false;  
  if (TMath::Abs(pM1->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pM2->fP.Eta()) > 2.4) fGoodMCKinematics = false;  

}


// ----------------------------------------------------------------------
void bmmSignalReader::candidateSelection(int mode) {
  int nc0(fCands.size()), nc1(0);

  fpCand = 0; 
  ((TH1D*)fpHistFile->Get("bnc0"))->Fill(nc0); 

  if (0 == nc0) {
    nc1 = 0; 
    ((TH1D*)fpHistFile->Get("bnc1"))->Fill(nc1); 
    if (fVerbose > 0) cout << "bmmSignalReader> no candidate " << TYPE << " found" << endl;
    return; 
  }
  
  TAnaCand *pCand; 

  // -- Create list with candidates that fullfil basic cuts on their final state particles
  vector<int> clist; 
  if (nc0 > 1) {
    for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
      pCand = fCands[iC]; 
      
      TAnaTrack *pl1 = fpEvt->getSigTrack(pCand->fSig1); 
      TAnaTrack *pl2 = fpEvt->getSigTrack(pCand->fSig2); 

      if (pl1->fQ*pl2->fQ > 0) continue;

      if (false == fGoodTracks[iC])   continue; 
      if (false == fGoodTracksPT[iC]) continue; 

      if (false == fGoodMuonsID[iC])  continue; 
      if (false == fGoodMuonsPT[iC])  continue; 
      ++nc1; 
      clist.push_back(iC);
    }
  } 

  // -- Select best candidate
  if (1 == clist.size()) {
    fpCand = fCands[clist[0]];
    nc1 = 1; 
    if (fVerbose > 0) cout << "bmmSignalReader> found exactly one candidate " << TYPE 
			   << " passing the constituents selection" 
			   << " with sigTracks " << fpCand->fSig1 << " .. " << fpCand->fSig2
			   << endl;
  } else if (clist.size() > 1) {
    double drMin(99.), dr(0.); 
    for (unsigned int i = 0; i < clist.size(); ++i) {
      pCand = fCands[clist[i]]; 
      TAnaTrack *pl1 = fpEvt->getSigTrack(pCand->fSig1); 
      TAnaTrack *pl2 = fpEvt->getSigTrack(pCand->fSig1+1); 
      
      // -- closest in rphi
      if (0 == mode) {
	dr = pl1->fPlab.DeltaR(pl2->fPlab); 
	if (dr < drMin) {
	  fpCand = pCand; 
	  drMin = dr; 
	  if (fVerbose > 0) cout << "bmmSignalReader> found another better candidate " << TYPE 
				 << " with sigTracks " << fpCand->fSig1 << " .. " << fpCand->fSig2
				 << " with dr = " << dr 
				 << endl;
	}
      }
    }
  } else {
    fpCand = fCands[0];
    nc1 = 1; 
    if (fVerbose > 0) cout << "bmmSignalReader> found no candidate " << TYPE 
			   << " passing the constituents selection" 
			   << endl;
  }

  ((TH1D*)fpHistFile->Get("bnc1"))->Fill(nc1); 

  if (BLIND && fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX) fpCand = 0; 
  if (fpCand) fillCandidateVariables();
}


// ----------------------------------------------------------------------
int bmmSignalReader::tmCand(TAnaCand *pC) {
  
  int truth(0); 
  TAnaTrack *pT; 
  vector<TGenCand*> gCand; 
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    if (pT->fGenIndex < 0) return 0; 
    gCand.push_back(fpEvt->getGenCand(pT->fGenIndex)); 
  }

  TGenCand *pG, *pM; 
  int matched(0), genDaughters(0); 
  for (unsigned int i = 0; i < gCand.size(); ++i) {
    pG = gCand[i]; 
    if (pG->fMom1 > -1) {
      pM = fpEvt->getGenCand(pG->fMom1); 
    } else {
      return 0; 
    }
    if (531 != TMath::Abs(pM->fID)) {
      return 0;  
    } else {
      ++matched; 
    }
  }
  if (matched == 2) {
    genDaughters = pM->fDau2 - pM->fDau1 +1; 
  }
  
  if (genDaughters == matched) {
    return 1; 
  } else {
    return 2;
  }

}



// ----------------------------------------------------------------------
void bmmSignalReader::fillHist() {
  bmmReader::fillHist(); 
}

// ---------------------------------------------------------------------- 
void bmmSignalReader::bookHist() {
  bmmReader::bookHist();
  cout << "==> bmmReader: bookHist " << endl;
  TH1D *h1; 
  h1 = new TH1D("h100", "", 100, 0., 100.);
  fTree->Branch("m2",      &fCandM,  "m2/D");
}

// ----------------------------------------------------------------------
void bmmSignalReader::readCuts(TString filename, int dump) {
  bmmReader::readCuts(filename, dump); 

  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> bmmReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> bmmReader: Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  TH1D *hcuts = (TH1D*)gFile->Get("hcuts"); 
  if (0 == hcuts) {
    new TH1D("hcuts", "", 1000, 0., 1000.);
    hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  }
  int ibin; 
  string cstring; 
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
    
    ok = checkCut(CutName, hcuts); 
    if (!ok) cout << "==> bmmSignalReader: error? nothing done with " << CutName << "!!" << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}

