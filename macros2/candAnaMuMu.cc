#include "candAnaMuMu.hh"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaMuMu::candAnaMuMu(bmm2Reader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candMuMuAna: constructor..." << endl;
  readCuts(cutsFile, 1); 
}


// ----------------------------------------------------------------------
candAnaMuMu::~candAnaMuMu() {
  cout << "==> candMuMuAna: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaMuMu::candAnalysis() {
  candAna::candAnalysis();
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(2); 

}

// ----------------------------------------------------------------------
void candAnaMuMu::processType() {

}


// ----------------------------------------------------------------------
void candAnaMuMu::genMatch() {
  fGenM1Tmi = fGenM2Tmi = -1; 
  fNGenPhotons = 0; 

  int id1(13), id2(13); 
  if (1000095 == TYPE) {id1 = 211; id2 = 13;}
  if (1000086 == TYPE) {id1 = 321; id2 = 13;}

  TGenCand *pC(0), *pM1(0), *pM2(0), *pB(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (TRUTHCAND == TMath::Abs(pC->fID)) {
      pM1 = pM2 = 0; 
      pB = pC;
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenCand(id); 
	if (id1 == TMath::Abs(pC->fID) || id2 == TMath::Abs(pC->fID)) {
	  if (0 == pM1) {
	    pM1 = fpEvt->getGenCand(id); 
	  } else {
	    pM2 = fpEvt->getGenCand(id); 
	  }
	}
      }
      if (0 != pM1 && 0 != pM2) {
	goodMatch = true; 
	fNGenPhotons = pB->fDau2 - pB->fDau1 - 1; 
	if (fVerbose > 10) {
	  cout << "found gen match for B gen idx = " << pB->fNumber << endl;
	}
	break;
      }
    }
  }

  
  fGenBTmi = -1; 
  if (goodMatch) {
    fGenBTmi = pB->fNumber; 
    double m = pB->fP.Mag();
    double p = pB->fP.P();
    // Meson pointer
    TGenCand *pM = fpEvt->getGenCand(pB->fMom1); 
    // the meson is the original except if it oscillated
    if (531 != TMath::Abs(pM->fID)) pM = pB;
    double x = (pM1->fV - pM->fV).Mag(); 
    fGenLifeTime = x*m/p/TMath::Ccgs();
    //     if (fGenLifeTime < 0.0004) {
    //       cout << "idx = " << pB->fNumber << " t: " << fGenLifeTime << " p: " << p  << " x: " << x 
    // 	   << " ID: " << pB->fID 
    // 	   << " m: " << pM1->fV.X() << ", " << pM1->fV.Y()  << ", " << pM1->fV.Z()
    // 	   << " B: " << pB->fV.X() << ", " << pB->fV.Y()  << ", " << pB->fV.Z()
    // 	   << endl;
    //       fpEvt->dumpGenBlock();
    //     }

    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fGenM1Tmi = pM1->fNumber; 
      fGenM2Tmi = pM2->fNumber; 
    } else {
      fGenM1Tmi = pM2->fNumber; 
      fGenM2Tmi = pM1->fNumber; 
    }
  } else {
    fGenM1Tmi = -1; 
    fGenM2Tmi = -1; 
  }
  
  if (fVerbose > 10) {
    cout << "fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fGenM2Tmi = " << fGenM2Tmi << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaMuMu::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = -1; 
  TAnaTrack *pT(0);
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    pT = fpEvt->getRecTrack(i); 
    if (pT->fGenIndex < 0) continue;

    // -- muon 1
    if (pT->fGenIndex == fGenM1Tmi) {
      fRecM1Tmi = i; 
    }

    // -- muon 2
    if (pT->fGenIndex == fGenM2Tmi) {
      fRecM2Tmi = i; 
    }

    // -- skip rest if both matches found
    if (fRecM1Tmi > -1 && fRecM2Tmi > -1) break;
  }


  if (fVerbose > 10) {
    cout << "fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaMuMu::candMatch() {
  fCandTmi = -1;   
  int idx(-1); 
  int d1Matched(0), d2Matched(0); 
  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC); 
    if (TYPE != pCand->fType) continue;
    
    d1Matched = d2Matched = 0; 
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex; 
      if (fVerbose > 10) {
	cout << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << endl;
      }
      if (idx == fRecM1Tmi) {
	d1Matched = 1; 
      }
      if (idx == fRecM2Tmi) {
	d2Matched = 1; 
      }
    }
    
    if (d1Matched && d2Matched) {
      fCandTmi = iC;
      break;
    }
  }
  if (fVerbose > 10) {
    cout << "fCandTmi = " << fCandTmi << endl;
  }
}


// ----------------------------------------------------------------------
void candAnaMuMu::efficiencyCalculation() {
  fGoodEffCand = false;

  // -- gen level 
  TGenCand *pB(0), *pM1(0), *pM2(0); 
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }
  pB  = fpEvt->getGenCand(fGenBTmi); 
  pM1 = fpEvt->getGenCand(fGenM1Tmi); 
  pM2 = fpEvt->getGenCand(fGenM2Tmi); 

  // -- reco level
  TAnaTrack *prM1(0), *prM2(0); 
  int m1Matched(0), m2Matched(0), m1ID(0), m2ID(0), m1GT(0), m2GT(0);
  if (fRecM1Tmi > -1) {
    m1Matched = 1; 
    prM1 = fpEvt->getRecTrack(fRecM1Tmi); 
    if (goodMuon(prM1)) m1ID = 1; 
    if (TRACKQUALITY > 0 && (0 == (prM1->fTrackQuality & TRACKQUALITY))) {
      m1GT = 0; 
    } else {
      m1GT = 1;
    }
  }

  if (fRecM2Tmi > -1) {
    m2Matched = 1; 
    prM2 = fpEvt->getRecTrack(fRecM2Tmi); 
    if (goodMuon(prM2)) m2ID = 1; 
    if (TRACKQUALITY > 0 && (0 == (prM2->fTrackQuality & TRACKQUALITY))) {
      m2GT = 0; 
    } else {
      m2GT = 1;
    }
  } 

  // -- cand level 
  TAnaCand *pCand(0);
  if (fCandTmi > -1) {
    pCand = fpEvt->getCand(fCandTmi);
  }

  // -- EffTree filling for all events with a signal decay
  fETgpt   = pB->fP.Perp(); 
  fETgeta  = pB->fP.Eta(); 
  fETg1pt  = pM1->fP.Perp(); 
  fETg1eta = pM1->fP.Eta(); 
  fETg2pt  = pM2->fP.Perp(); 
  fETg2eta = pM2->fP.Eta(); 
  if (m1Matched) {
    fETm1pt  = prM1->fPlab.Perp(); 
    fETm1eta = prM1->fPlab.Eta(); 
    fETm1q   = prM1->fQ;
    fETm1gt  = (m1GT>0?true:false); 
    fETm1id  = (m1ID>0?true:false);
  } else {
    fETm1pt  = -99.; 
    fETm1eta = -99.; 
    fETm1q   = -99;
    fETm1gt  = false; 
    fETm1id  = false;
  }
  if (m2Matched) {
    fETm2pt  = prM2->fPlab.Perp(); 
    fETm2eta = prM2->fPlab.Eta(); 
    fETm2q   = prM2->fQ;
    fETm2gt  = (m2GT>0?true:false); 
    fETm2id  = (m2ID>0?true:false);
  } else {
    fETm2pt  = -99.; 
    fETm2eta = -99.; 
    fETm2q   = -99;
    fETm2gt  = false; 
    fETm2id  = false;
  }
  if (pCand) {
    fETcandMass = pCand->fMass; 
  } else {
    fETcandMass = -99.;
  }
  fEffTree->Fill(); 

}


// ----------------------------------------------------------------------
void candAnaMuMu::bookHist() {
  cout << "==>candAnaMuMu: bookHist" << endl;
  candAna::bookHist();

}


// ----------------------------------------------------------------------
void candAnaMuMu::readCuts(string filename, int dump) {

  candAna::readCuts(filename, dump); 

}
