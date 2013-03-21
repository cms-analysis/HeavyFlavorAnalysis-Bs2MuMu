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

  // -- modifications for rare backgrounds
  if (1000082 == TYPE) {id1 = 321; id2 = 321;}
  if (1000083 == TYPE) {id1 = 321; id2 = 211;}
  if (1000084 == TYPE) {id1 = 211; id2 = 211;}
  if (1000085 == TYPE) {id1 = 211; id2 = 13;}
  if (1000086 == TYPE) {id1 = 321; id2 = 13;}

  if (1000091 == TYPE) {id1 = 211; id2 = 211;}
  if (1000092 == TYPE) {id1 = 321; id2 = 211;}
  if (1000093 == TYPE) {id1 = 321; id2 = 321;}
  if (1000095 == TYPE) {id1 = 211; id2 = 13;}
  if (1000098 == TYPE) {id1 = 211; id2 = 211;}

  if (1000060 == TYPE) {id1 = 2212; id2 = 211;}
  if (1000061 == TYPE) {id1 = 2212; id2 = 321;}
  if (1000062 == TYPE) {id1 = 2212; id2 = 13;}

  TGenCand *pC(0), *pM1(0), *pM2(0), *pB(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenT(); ++i) {
    pC = fpEvt->getGenT(i); 
    if (TRUTHCAND == TMath::Abs(pC->fID)) {
      pM1 = pM2 = 0; 
      pB = pC;
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenTWithIndex(id); 
	if (id1 == TMath::Abs(pC->fID) || id2 == TMath::Abs(pC->fID)) {
	  if (0 == pM1) {
	    pM1 = fpEvt->getGenTWithIndex(id); 
	  } else {
	    pM2 = fpEvt->getGenTWithIndex(id); 
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
    //    TGenCand *pM = fpEvt->getGenCand(pB->fMom1); 
    TGenCand *pM = fpEvt->getGenTWithIndex(pB->fMom1); 
    // the meson is the original except if it oscillated
    if (531 != TMath::Abs(pM->fID)) pM = pB;
    double x = (pM1->fV - pM->fV).Mag(); 
    fGenLifeTime = x*m/p/TMath::Ccgs();

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

//   // -- check that only one reco track is matched to each gen particle
//   //    else skip the *event*!
//   if (fGenM1Tmi > -1 && fGenM2Tmi > -1) {
//     TSimpleTrack *pT(0);
//     int cntM1(0), cntM2(0); 
//     for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
//       pT = fpEvt->getSimpleTrack(i); 
//       if (fGenM1Tmi == pT->getGenIndex()) ++cntM1;
//       if (fGenM2Tmi == pT->getGenIndex()) ++cntM2;
//     }
    
//     static int cntBadEvents = 0; 
//     if (cntM1 > 1 || cntM2 > 1) {
//       cout << "BAD BAD event: multiple reco tracks matched to the same gen particle! " << ++cntBadEvents 
// 	   << ": " << cntM1 << " .. " << cntM2
// 	   << " (gen-reco matches) " << endl;
//       fBadEvent = true; 
//     }
//   }

}


// ----------------------------------------------------------------------
void candAnaMuMu::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = -1; 
  TSimpleTrack *pT(0);
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i); 
    if (pT->getGenIndex() < 0) continue;

    // -- muon 1
    if (pT->getGenIndex() == fGenM1Tmi) {
      fRecM1Tmi = i; 
    }

    // -- muon 2
    if (pT->getGenIndex() == fGenM2Tmi) {
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
  if (fRecM1Tmi < 0 ||  fRecM2Tmi < 0) return; 
  
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC); 
    if (TYPE != pCand->fType)  continue;
    
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
  pB  = fpEvt->getGenTWithIndex(fGenBTmi); 
  pM1 = fpEvt->getGenTWithIndex(fGenM1Tmi); 
  pM2 = fpEvt->getGenTWithIndex(fGenM2Tmi); 

  // -- reco level
  TSimpleTrack *prM1(0), *prM2(0); 
  double bla(0); 
  int m1Matched(0), m2Matched(0), m1ID(0), m1tmID(0), m1mvaID(0), m2ID(0), m2tmID(0), m2mvaID(0), m1GT(0), m2GT(0);
  if (fRecM1Tmi > -1) {
    m1Matched = 1; 
    prM1 = fpEvt->getSimpleTrack(fRecM1Tmi); 
    if (tightMuon(prM1)) m1tmID = 1; 
    if (mvaMuon(prM1, bla)) m1mvaID = 1; 
    if (prM1->getHighPurity()) {
      m1GT = 1; 
    } else {
      m1GT = 0;
    }
  }

  if (fRecM2Tmi > -1) {
    m2Matched = 1; 
    prM2 = fpEvt->getSimpleTrack(fRecM2Tmi); 
    if (tightMuon(prM2)) m2tmID = 1; 
    if (mvaMuon(prM2, bla)) m2mvaID = 1; 
    if (prM2->getHighPurity()) {
      m2GT = 1; 
    } else {
      m2GT = 0;
    }
  } 

  m1ID = m1tmID; 
  m2ID = m2tmID; 

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
    fETm1pt    = prM1->getP().Perp(); 
    fETm1eta   = prM1->getP().Eta(); 
    fETm1q     = prM1->getCharge();
    fETm1gt    = (m1GT>0?true:false); 
    fETm1id    = (m1ID>0?true:false);
    fETm1tmid  = (m1tmID>0?true:false);
    fETm1mvaid = (m1mvaID>0?true:false);
  } else {
    fETm1pt  = -99.; 
    fETm1eta = -99.; 
    fETm1q   = -99;
    fETm1gt  = false; 
    fETm1id  = false;
    fETm1tmid  = false;
    fETm1mvaid = false;
  }
  if (m2Matched) {
    fETm2pt    = prM2->getP().Perp(); 
    fETm2eta   = prM2->getP().Eta(); 
    fETm2q     = prM2->getCharge();
    fETm2gt    = (m2GT>0?true:false); 
    fETm2id    = (m2ID>0?true:false);
    fETm2tmid  = (m2tmID>0?true:false);
    fETm2mvaid = (m2mvaID>0?true:false);
  } else {
    fETm2pt  = -99.; 
    fETm2eta = -99.; 
    fETm2q   = -99;
    fETm2gt  = false; 
    fETm2id  = false;
    fETm2tmid  = false;
    fETm2mvaid = false;
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
