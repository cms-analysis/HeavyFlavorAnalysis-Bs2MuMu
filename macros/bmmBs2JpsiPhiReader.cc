#include "bmmBs2JpsiPhiReader.hh"
#include "TRandom.h"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using std::string;
using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
bmmBs2JpsiPhiReader::bmmBs2JpsiPhiReader(TChain *tree, TString evtClassName): bmmReader(tree, evtClassName) {
  cout << "==> bmmBs2JpsiPhiReader: constructor..." << endl;
}


// ----------------------------------------------------------------------
bmmBs2JpsiPhiReader::~bmmBs2JpsiPhiReader() {
  cout << "==> bmmBs2JpsiPhiReader: destructor..." << endl;

}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::moreBasicCuts() {
  cout << "   bmmBs2JpsiPhiReader: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodJpsiMass", "m(J/psi)", fGoodJpsiMass); 
  fAnaCuts.addCut("fGoodDeltaR", "Delta R(KK)", fGoodDeltaR); 
  fAnaCuts.addCut("fGoodMKK", "m(KK) [GeV]", fGoodMKK); 
}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::startAnalysis() {
  bmmReader::startAnalysis();
  cout << "==> bmmBs2JpsiPhiReader: Summary of analysis cuts:" << endl;
  fAnaCuts.dumpAll(); 
}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::eventProcessing() {
  bmmReader::eventProcessing();
}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::initVariables() {
  bmmReader::initVariables();
  fDeltaR = fMKK = 99.;
  fGoodDeltaR = fGoodMKK = false;
  fKa1PtNrf = fKa2PtNrf = -99.;
  fKa1PtGen = fKa2PtGen = fKa1EtaGen = fKa2EtaGen = -99.;

}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::insertCand(TAnaCand* pCand) {
  //  cout << "bmmBs2JpsiPhiReader::insertCand" << endl;
  bmmReader::insertCand(pCand);
}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::efficiencyCalculation() {

  fGoodEffCand = false;

  ((TH1D*)fpHistFile->Get("efficiency"))->Fill(0); 
  ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(1, "all events"); 


  // -- gen level 
  TGenCand *pM1(0), *pM2(0), *pK1(0), *pK2(0); 
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi || -1 == fGenK1Tmi || -1 == fGenK2Tmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }

  pM1 = fpEvt->getGenCand(fGenM1Tmi); 
  pM2 = fpEvt->getGenCand(fGenM2Tmi); 
  pK1 = fpEvt->getGenCand(fGenK1Tmi); 
  pK2 = fpEvt->getGenCand(fGenK2Tmi); 

  ((TH1D*)fpHistFile->Get("efficiency"))->Fill(1); 
  ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(2, "gen signal decays"); 


  // -- reco level
  TAnaTrack *prM1(0), *prM2(0), *prK1(0), *prK2(0); 
  int m1Matched(0), m2Matched(0), k1Matched(0), k2Matched(0), m1ID(0), m2ID(0), m1GT(0), m2GT(0), k1GT(0), k2GT(0);
  if (fRecM1Tmi > -1) {
    m1Matched = 1; 
    prM1 = fpEvt->getRecTrack(fRecM1Tmi); 
    if (muonID(prM1)) m1ID = 1; 
    if (TRACKQUALITY > 0 && (0 == (prM1->fTrackQuality & TRACKQUALITY))) {
      m1GT = 0; 
    } else {
      m1GT = 1;
    }
  }

  if (fRecM2Tmi > -1) {
    m2Matched = 1; 
    prM2 = fpEvt->getRecTrack(fRecM2Tmi); 
    if (muonID(prM2)) m2ID = 1; 
    if (TRACKQUALITY > 0 && (0 == (prM2->fTrackQuality & TRACKQUALITY))) {
      m2GT = 0; 
    } else {
      m2GT = 1;
    }
  } 

  if (fRecK1Tmi > -1) {
    k1Matched = 1; 
    prK1 = fpEvt->getRecTrack(fRecK1Tmi); 
    if (TRACKQUALITY > 0 && (0 == (prK1->fTrackQuality & TRACKQUALITY))) {
      k1GT = 0; 
    } else {
      k1GT = 1;
    }
  } 

  if (fRecK2Tmi > -1) {
    k2Matched = 1; 
    prK2 = fpEvt->getRecTrack(fRecK2Tmi); 
    if (TRACKQUALITY > 0 && (0 == (prK2->fTrackQuality & TRACKQUALITY))) {
      k2GT = 0; 
    } else {
      k2GT = 1;
    }
  } 

  // -- cand level 
  TAnaCand *pCand(0);
  if (fCandTmi > -1) {
    pCand = fCands[fCandTmi];
  }


  // -- EffTree filling for all events with a signal decay
  fETg1pt  = pM1->fP.Perp(); 
  fETg1eta = pM1->fP.Eta(); 
  fETg2pt  = pM2->fP.Perp(); 
  fETg2eta = pM2->fP.Eta(); 
  fETg3pt  = pK1->fP.Perp(); 
  fETg3eta = pK1->fP.Eta(); 
  fETg4pt  = pK2->fP.Perp(); 
  fETg4eta = pK2->fP.Eta(); 
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
  if (k1Matched) {
    fETk1pt  = prK1->fPlab.Perp(); 
    fETk1eta = prK1->fPlab.Eta(); 
    fETk1q   = prK1->fQ;
    fETk1gt  = (k1GT>0?true:false); 
  } else {
    fETk1pt  = -99.; 
    fETk1eta = -99.; 
    fETk1q   = -99;
    fETk1gt  = false; 
  }
  if (k2Matched) {
    fETk2pt  = prK2->fPlab.Perp(); 
    fETk2eta = prK2->fPlab.Eta(); 
    fETk2q   = prK2->fQ;
    fETk2gt  = (k2GT>0?true:false); 
  } else {
    fETk2pt  = -99.; 
    fETk2eta = -99.; 
    fETk2q   = -99;
    fETk2gt  = false; 
  }
  if (pCand) {
    fETcandMass = pCand->fMass; 
  } else {
    fETcandMass = -99.;
  }
  fEffTree->Fill(); 


  // -- results...
  if ((TMath::Abs(pM1->fP.Eta()) < 2.5)
      && (TMath::Abs(pM2->fP.Eta()) < 2.5) 
      && (TMath::Abs(pK1->fP.Eta()) < 2.5)
      && (TMath::Abs(pK2->fP.Eta()) < 2.5)
      ) {
    ((TH1D*)fpHistFile->Get("efficiency"))->Fill(2); 
    ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(3, "+ eta cuts"); 
    
    if ((pM1->fP.Perp() > 1.0)  
	&& (pM2->fP.Perp() > 1.0)
	&& (pK1->fP.Perp() > 0.4)
	&& (pK2->fP.Perp() > 0.4)
	) {
      ((TH1D*)fpHistFile->Get("efficiency"))->Fill(3); 
      ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(4, "+ pT cuts"); 

      if (m1Matched && m2Matched && k1Matched && k2Matched
	  && (prM1->fPlab.Perp() > 3.0) && (TMath::Abs(prM1->fPlab.Eta()) < 2.4)
	  && (prM2->fPlab.Perp() > 3.0) && (TMath::Abs(prM2->fPlab.Eta()) < 2.4)
	  && (prK1->fPlab.Perp() > 0.5) && (TMath::Abs(prK1->fPlab.Eta()) < 2.4)
	  && (prK2->fPlab.Perp() > 0.5) && (TMath::Abs(prK2->fPlab.Eta()) < 2.4)
	  && (prM1->fQ*prM2->fQ < 0)
	  && m1GT && m2GT && k1GT && k2GT
	  ) {
	((TH1D*)fpHistFile->Get("efficiency"))->Fill(4); 
	((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(5, "+ reco tracks/cuts"); 
	
	if (m1ID && m2ID) {
	  ((TH1D*)fpHistFile->Get("efficiency"))->Fill(5); 
	  ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(6, "+ muon ID"); 
	  
	  if (fGoodHLT) {
	    ((TH1D*)fpHistFile->Get("efficiency"))->Fill(6); 
	    ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(7, "+ trigger"); 
	    
	    if (pCand) {
	      ((TH1D*)fpHistFile->Get("efficiency"))->Fill(7); 
	      ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(8, "+ candidate"); 
	      fGoodEffCand = true;
	      ((TH1D*)fpHistFile->Get("effMass"))->Fill(pCand->fMass);
	    } 
	  }
	} 
      } 
    } 
  } 

}



// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::genMatch() {

  fGenM1Tmi = fGenM2Tmi = fGenK1Tmi = -1; 
  fNGenPhotons = 0; 

  TGenCand *pC(0), *pB(0), *pPsi(0), *pPhi(0), *pM1(0), *pM2(0), *pK1(0), *pK2(0); 
  int npsi(0), nphi(0), nb(0), ngamma(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (531 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = pB->fDau2 - pB->fDau1 + 1; 
      if (nb > 2) continue; // skip B decays where more than J/psi and phi came from B
      ngamma = 0; 
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenCand(id); 
	if (22 == TMath::Abs(pC->fID)) ++ngamma;
	if (443 == TMath::Abs(pC->fID)) {
	  pPsi = pC; 
	  npsi = pPsi->fDau2 - pPsi->fDau1 + 1; 
	  pM1 = pM2 = 0;
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenCand(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenCand(idd); 
	      } else {
		pM2 = fpEvt->getGenCand(idd); 
	      }
	    }
	  }
	} else if (333 == TMath::Abs(pC->fID)) {
	  pPhi = fpEvt->getGenCand(id); 
	  nphi = pPhi->fDau2 - pPhi->fDau1 + 1; 
	  pK1 = pK2 = 0;
	  for (int idd = pPhi->fDau1; idd <= pPhi->fDau2; ++idd) {
	    pC = fpEvt->getGenCand(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (321 == TMath::Abs(pC->fID)) {
	      if (0 == pK1) {
		pK1 = fpEvt->getGenCand(idd); 
	      } else {
		pK2 = fpEvt->getGenCand(idd); 
	      }
	    }
	  }
	}
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK1 && 0 != pK2 
	  && (pPsi->fMom1 == pPhi->fMom1)
	  ) {
	goodMatch = true; 
	fNGenPhotons = ngamma;
	break;
      }
    }
  }

  if (!goodMatch) {
    fGoodMCKinematics = false; 
    if (fVerbose > 2) cout << "No matched signal decay found" << endl;
    return;
  }

  if (goodMatch) {
    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fGenM1Tmi = pM1->fNumber; 
      fGenM2Tmi = pM2->fNumber; 
    } else {
      fGenM1Tmi = pM2->fNumber; 
      fGenM2Tmi = pM1->fNumber; 
    }
    if (pK1->fP.Perp() > pK2->fP.Perp()) {
      fGenK1Tmi = pK1->fNumber; 
      fGenK2Tmi = pK2->fNumber; 
    } else {
      fGenK1Tmi = pK2->fNumber; 
      fGenK2Tmi = pK1->fNumber; 
    }
  } else {
    fGenM1Tmi = -1; 
    fGenM2Tmi = -1; 
    fGenK1Tmi = -1;  
    fGenK2Tmi = -1;  
  }

}

// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::recoMatch() {
  fRecM1Tmi = fRecM2Tmi = fRecK1Tmi = fRecK2Tmi =-1; 
  TAnaTrack *pT(0);
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    pT = fpEvt->getRecTrack(i); 
    if (pT->fGenIndex < 0) continue;

    // -- muon 1
    if (fGenM1Tmi > -1 && pT->fGenIndex == fGenM1Tmi) {
      fRecM1Tmi = i; 
    }

    // -- muon 2
    if (fGenM2Tmi > -1 && pT->fGenIndex == fGenM2Tmi) {
      fRecM2Tmi = i; 
    }

    // -- kaon 1
    if (fGenK1Tmi > -1 && pT->fGenIndex == fGenK1Tmi) {
      fRecK1Tmi = i; 
    }

    // -- kaon 2
    if (fGenK2Tmi > -1 && pT->fGenIndex == fGenK2Tmi) {
      fRecK2Tmi = i; 
    }

    // -- skip rest if all matches found
    if (fRecM1Tmi > -1 && fRecM2Tmi > -1 && fRecK1Tmi > -1 && fRecK2Tmi > -1) break;
  }


  if (fVerbose > 10) {
    cout << "fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi << endl;
    cout << "fRecK1Tmi = " << fRecK1Tmi << " matched to fGenK1Tmi = " << fGenK1Tmi << endl;
    cout << "fRecK2Tmi = " << fRecK2Tmi << " matched to fGenK2Tmi = " << fGenK2Tmi << endl;
  }

}

// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::candMatch() {
  fCandTmi = -1;   
  int idx(-1), type(-1); 
  int d1Matched(0), d2Matched(0), d3Matched(0), d4Matched(0); 
  TAnaCand *pCand(0);
  for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
    pCand = fCands[iC]; 
    
    d1Matched = d2Matched = d3Matched = d4Matched = 0; 
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex; 
      type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
      //       if (fGenM1Tmi > -1) cout << "  --> " << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi << endl;
      if (fVerbose > 10) {
	cout << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi << endl;
      }
      if (fRecM1Tmi > -1 && type == 13 && idx == fRecM1Tmi) {
	d1Matched = 1; 
      }
      if (fRecM2Tmi > -1 && type == 13 && idx == fRecM2Tmi) {
	d2Matched = 1; 
      }
      if (fRecK1Tmi > -1 && type == 321 && idx == fRecK1Tmi) {
	d3Matched = 1; 
      }
      if (fRecK2Tmi > -1 && type == 321 && idx == fRecK2Tmi) {
	d4Matched = 1; 
      }
    }
    
    if (d1Matched && d2Matched && d3Matched && d4Matched) {
      fCandTmi = iC;
      break;
    }
  }
  if (fVerbose > 10) {
    cout << "fCandTmi = " << fCandTmi << " matched to rec tracks " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi  << " " << fRecK2Tmi 
	 << endl;
  }

}


// ----------------------------------------------------------------------
int bmmBs2JpsiPhiReader::tmCand(TAnaCand *pC) {
  TAnaCand *pCand(0);
  for (int iC = 0; iC < static_cast<int>(fCands.size()); ++iC) {
    pCand = fCands[iC]; 
    if (pCand == pC) {
      if (iC == fCandTmi) {
	if (fNGenPhotons) {
	  return 2; 
	} else {
	  return 1;
	}
      } else {
	return 0; 
      }
    }
  }
  return 0;
}



// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::MCKinematics() {
  ((TH1D*)fpHistFile->Get("genStudy"))->Fill(1); 
  fGoodMCKinematics = true; 
  TGenCand *pC(0), *pB(0), *pPsi(0), *pPhi(0), *pM1(0), *pM2(0), *pK1(0), *pK2(0); 
  int npsi(0), nphi(0), nb(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (531 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = pB->fDau2 - pB->fDau1 + 1; 
      if (nb > 2) continue; // skip B decays where more than J/psi and phi came from B
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenCand(id); 
	if (443 == TMath::Abs(pC->fID)) {
	  pPsi = pC; 
	  npsi = pPsi->fDau2 - pPsi->fDau1 + 1; 
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenCand(idd); 
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenCand(idd); 
	      } else {
		pM2 = fpEvt->getGenCand(idd); 
	      }
	    }
	  }
	} else if (333 == TMath::Abs(pC->fID)) {
	  pPhi = fpEvt->getGenCand(id); 
	  nphi = pPhi->fDau2 - pPhi->fDau1 + 1; 
	  for (int idd = pPhi->fDau1; idd <= pPhi->fDau2; ++idd) {
	    pC = fpEvt->getGenCand(idd); 
	    if (321 == TMath::Abs(pC->fID)) {
	      if (0 == pK1) {
		pK1 = fpEvt->getGenCand(idd); 
	      } else {
		pK2 = fpEvt->getGenCand(idd); 
	      }
	    }
	  }
	}
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK1 && 0 != pK2) {
	goodMatch = true; 
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
      pK1->dump();
      pK2->dump();
    } else {
      cout << "no tm decay found" << endl;
    }
  }

  if (!goodMatch) {
    fGoodMCKinematics = false; 
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    //    fpEvt->dumpGenBlock(); 
    return;
  } else {
    if (fVerbose > 2 ) cout << "--------------------> Did find a matched signal decay" << endl;
    //    fpEvt->dumpGenBlock(); 
  }


  // ----------------------------------------------------------------------
  // -- Acceptance: 
  TAnaTrack *pT(0), *prM1(0), *prM2(0); 
  int m1Matched(0), m2Matched(0), k1Matched(0), k2Matched(0);
  int m1Acc(0), m2Acc(0), k1Acc(0), k2Acc(0);
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    pT = fpEvt->getRecTrack(i); 
    if (pT->fGenIndex == pM1->fNumber) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(11); 
      prM1 = pT; 
      m1Matched = 1; 
      if ((pT->fPlab.Perp() > 1.) && (TMath::Abs(pT->fPlab.Eta()) < 2.4)) {
	m1Acc = 1; 
	((TH1D*)fpHistFile->Get("acceptance"))->Fill(21); 
      }
    }
    if (pT->fGenIndex == pM2->fNumber) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(12); 
      prM2 = pT; 
      m2Matched = 1; 
      if ((pT->fPlab.Perp() > 1.) && (TMath::Abs(pT->fPlab.Eta()) < 2.4)) {
	m2Acc = 1; 
	((TH1D*)fpHistFile->Get("acceptance"))->Fill(22); 
      }
    }
    if (pT->fGenIndex == pK1->fNumber) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(13); 
      k1Matched = 1; 
      if ((pT->fPlab.Perp() > 0.5) && (TMath::Abs(pT->fPlab.Eta()) < 2.4)) {
	k1Acc = 1; 
	((TH1D*)fpHistFile->Get("acceptance"))->Fill(23); 
      }
    }
    if (pT->fGenIndex == pK2->fNumber) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(14); 
      k2Matched = 1; 
      if ((pT->fPlab.Perp() > 0.5) && (TMath::Abs(pT->fPlab.Eta()) < 2.4)) {
	k2Acc = 1; 
	((TH1D*)fpHistFile->Get("acceptance"))->Fill(24); 
      }
    }
  }

  ((TH1D*)fpHistFile->Get("acceptance"))->Fill(30); 
  ((TH1D*)fpHistFile->Get("acceptance"))->Fill(1); // denominator
  if (m1Matched && m1Acc && m2Matched && m2Acc) {
    ((TH1D*)fpHistFile->Get("acceptance"))->Fill(31); 
    if (k1Matched && k1Acc && k2Matched && k2Acc) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(32); 
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(2);  // numerator
    }
  }



  // ----------------------------------------------------------------------
  // -- preselection efficiency
  TAnaCand *pCand; 
  ((TH1D*)fpHistFile->Get("presel"))->Fill(10); 
  if (m1Matched && m1Acc && m2Matched && m2Acc && k1Matched && k1Acc && k2Matched && k2Acc) {
    ((TH1D*)fpHistFile->Get("presel"))->Fill(11); 
    if (muonID(prM1) && muonID(prM2)) {
      ((TH1D*)fpHistFile->Get("presel"))->Fill(12); 
      ((TH1D*)fpHistFile->Get("presel"))->Fill(1);  // denominator
    }

    int tm(0);
    for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
      pCand = fCands[iC]; 

      tm = 0; 
      for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
	pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
	if (pT->fGenIndex == pM1->fNumber) {
	  prM1 = pT; 
	  ++tm;
	}
	if (pT->fGenIndex == pM2->fNumber) {
	  prM2 = pT; 
	  ++tm;
	}
	if (pT->fGenIndex == pK1->fNumber) {
	  ++tm;
	}
	if (pT->fGenIndex == pK2->fNumber) {
	  ++tm;
	}
      }

      // -- all signal tracks must be matched
      if (tm != (pCand->fSig2-pCand->fSig1+1)) {
	//	cout << "Did not find matched candidate " << pCand->fType << " at " << iC << endl;
	continue;
      }
      //      cout << "Found matched candidate "  << pCand->fType << " at " << iC << endl;
      ((TH1D*)fpHistFile->Get("presel"))->Fill(20); 
      if (muonID(prM1) && muonID(prM2)) {
	((TH1D*)fpHistFile->Get("presel"))->Fill(21); 
	((TH1D*)fpHistFile->Get("presel"))->Fill(2);  // numerator
      }
      break;
    }
  }

  // -- hard coded ?! FIXME
  if (pM1->fP.Perp() < 3.0) fGoodMCKinematics = false;  
  if (pM2->fP.Perp() < 3.0) fGoodMCKinematics = false;  
  if (pK1->fP.Perp() < 0.5) fGoodMCKinematics = false;  
  if (pK2->fP.Perp() < 0.5) fGoodMCKinematics = false;  
  if (TMath::Abs(pM1->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pM2->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pK1->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pK2->fP.Eta()) > 2.4) fGoodMCKinematics = false;  

}



// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::candidateSelection(int mode) {
  int nc0(fCands.size()), nc1(0);

  fpCand = 0; 
  fCandIdx = -1; 
  ((TH1D*)fpHistFile->Get("bnc0"))->Fill(nc0); 
  
  if (0 == nc0) {
    nc1 = 0; 
    fCandIdx = -1; 
    ((TH1D*)fpHistFile->Get("bnc1"))->Fill(nc1); 
    if (fVerbose > 0) cout << "bmmBs2JpsiPhiReader> no candidate " << TYPE << " found" << endl;
    return;
  }

  TAnaCand *pCand; 

  // -- Create list with candidates that fullfil basic cuts on their final state particles
  vector<int> clist; 
  if (nc0 > 1) {
    for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
      pCand = fCands[iC]; 
      
      TAnaTrack *ps(0), *pl1(0), *pl2(0);
      for (int it = pCand->fSig1; it <= pCand->fSig2; ++it) {
	ps = fpEvt->getSigTrack(it); 
	if (TMath::Abs(ps->fMCID) != 13) continue;
	if (0 == pl1) {
	  pl1 = ps;
	} else {
	  pl2 = ps;
	}
      }

      if (pl1->fQ*pl2->fQ > 0) continue;

      if (false == fvGoodTracks[iC])   continue; 
      if (false == fvGoodTracksPt[iC]) continue; 
      if (false == fvGoodTracksEta[iC]) continue; 

      if (false == fvGoodMuonsID[iC])  continue; 
      if (false == fvGoodMuonsPt[iC])  continue; 
      if (false == fvGoodMuonsEta[iC])  continue; 
      fvGoodCand[iC] = true; 
      ++nc1; 
      clist.push_back(iC);
    }
  } 

  // -- 'best' candidate selection
  if (1 == clist.size()) {
    // -- ONE candidate fulfilled the basic requirements. Take it.
    pCand = fCands[clist[0]]; 
    if (BLIND && pCand->fMass > SIGBOXMIN && pCand->fMass < SIGBOXMAX) {
      fvGoodCand[clist[0]] = false; 
      fpCand = 0; 
      fCandIdx = -1;
    } else{
      fpCand = fCands[clist[0]];
      fCandIdx = clist[0]; 
      nc1 = 1; 
      if (fVerbose > 0) cout << "bmmBs2JpsiPhiReader> found exactly one candidate " << TYPE 
			     << " passing the constituents selection at "  << clist[0]
			     << endl;
    }
  } else if (clist.size() > 1) {
    // -- SEVERAL candidate fulfilled the basic requirements. Choose one. 
    double drMin(99.), dr(0.); 
    for (unsigned int i = 0; i < clist.size(); ++i) {
      pCand = fCands[clist[i]]; 
      if (BLIND && pCand->fMass > SIGBOXMIN && pCand->fMass < SIGBOXMAX) {
	fvGoodCand[clist[i]] = false; 
	continue; 
      }
      TAnaTrack *pl1 = fpEvt->getSigTrack(pCand->fSig1); 
      TAnaTrack *pl2 = fpEvt->getSigTrack(pCand->fSig1+1); 
      
      // -- if all candidates are taken, then fpCand at this point is pointing to a random cand
      if (mode < 0) {
	if (fVerbose > 0) cout << "bmmBs2JpsiPhiReader> changing to candidate at " << clist[i] << endl;
	fpCand = pCand; 
	fCandIdx = clist[i]; 
      }
      // -- closest in rphi. This is pointless for this mode ...
      if (0 == mode) {
	dr = pl1->fPlab.DeltaR(pl2->fPlab); 
	if (dr < drMin) {
	  fpCand = pCand; 
	  fCandIdx = clist[i]; 
	  fvGoodCand[clist[i]] = true; 
	  drMin = dr; 
	  if (fVerbose > 0) cout << "bmmBs2JpsiPhiReader> found another better candidate " << TYPE 
				 << " with dr = " << dr << " at " << clist[i]
				 << endl;
	}
      }
      // -- higher pT
      if (1 == mode) {
	dr = pCand->fPlab.Perp(); 
	if (dr > drMin) {
	  fpCand = pCand; 
	  fCandIdx = clist[i]; 
	  fvGoodCand[fCandIdx] = true; 
	  drMin = dr; 
	  if (fVerbose > 0) cout << "bmmBs2JpsiPhiReader> found another better candidate " << TYPE << " with pT = " << dr << endl;
	}
      } 
    }
  } else {
    // -- NO candidate fulfilled the basic requirements. Then take the first one
    pCand = fCands[0]; 
    if (BLIND && pCand->fMass > SIGBOXMIN && pCand->fMass < SIGBOXMAX) {
      fvGoodCand[0] = false; 
      fpCand = 0;
      fCandIdx = -1; 
    } else {
      fpCand = fCands[0];
      fCandIdx = 0; 
      fvGoodCand[0] = true; 
      nc1 = 1; 
      if (fVerbose > 0) cout << "bmmBs2JpsiPhiReader> found no candidate " << TYPE 
			     << " passing the constituents selection, filling 0"
			     << " with signal tracks " << fpCand->fSig1 << ".." << fpCand->fSig2
			     << endl;
    }
  }
  ((TH1D*)fpHistFile->Get("bnc1"))->Fill(nc1); 

  //  if (fpCand) fillCandidateVariables();
}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::fillCandidateVariables() {

  if (0 == fpCand) return;

  // -- Check for J/psi mass
  TAnaCand *pD = 0; 
  fGoodJpsiMass = false; 
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i); 
    //    cout << "i = " << i << " pD = " << pD << endl;
    if (pD->fType == JPSITYPE) {
      if ((JPSIMASSLO < pD->fMass) && (pD->fMass < JPSIMASSHI)) fGoodJpsiMass = true;
      fJpsiMass = pD->fMass;
      //       cout << "type = " << pD->fType 
      //        	   << " with mass = " << pD->fMass 
      // 	   << " fGoodJpsiMass = " << fGoodJpsiMass 
      // 	   << endl;
    }
  }

  // -- Get Kaons
  TAnaTrack *p0; 
  TAnaTrack *p1(0), *ps1(0);
  TAnaTrack *p2(0), *ps2(0); 
  
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    if (TMath::Abs(p0->fMCID) != 321) continue;
    if (0 == p1) {
      p1 = p0; 
    } else {
      p2 = p0; 
    }
  }
 
  // -- switch to RecTracks!
  ps1 = p1; 
  p1 = fpEvt->getRecTrack(p1->fIndex);
  ps2 = p2; 
  p2 = fpEvt->getRecTrack(p2->fIndex);


  fKa1Pt        = p1->fPlab.Perp(); 
  fKa1Eta       = p1->fPlab.Eta(); 
  fKa1Phi       = p1->fPlab.Phi(); 
  fKa1TkQuality = p1->fTrackQuality & TRACKQUALITY;
  fKa1PtNrf     = ps1->fPlab.Perp();
  fKa1EtaNrf    = ps1->fPlab.Eta();

  fKa2Pt        = p2->fPlab.Perp(); 
  fKa2Eta       = p2->fPlab.Eta(); 
  fKa2Phi       = p2->fPlab.Phi(); 
  fKa2TkQuality = p2->fTrackQuality & TRACKQUALITY;
  fKa2PtNrf     = ps2->fPlab.Perp();
  fKa2EtaNrf    = ps2->fPlab.Eta();

  if (tmCand(fpCand)) {
    TGenCand *pg1 = fpEvt->getGenCand(fpEvt->getRecTrack(p1->fIndex)->fGenIndex);
    fKa1PtGen     = pg1->fP.Perp();
    fKa1EtaGen    = pg1->fP.Eta();
    TGenCand *pg2 = fpEvt->getGenCand(fpEvt->getRecTrack(p2->fIndex)->fGenIndex);
    fKa2PtGen     = pg2->fP.Perp();
    fKa2EtaGen    = pg2->fP.Eta();
  } else {
    fKa1PtGen     = -99.;
    fKa1EtaGen    = -99.;
    fKa2PtGen     = -99.;
    fKa2EtaGen    = -99.;
  }

  fDeltaR  = p1->fPlab.DeltaR(p2->fPlab); 

  TLorentzVector ka1, ka2; 
  ka1.SetPtEtaPhiM(fKa1Pt, fKa1Eta, fKa1Phi, MKAON); 
  ka2.SetPtEtaPhiM(fKa2Pt, fKa2Eta, fKa2Phi, MKAON); 

  TLorentzVector phiCand = ka1 + ka2; 
  fMKK     = phiCand.M();

  fGoodDeltaR = (fDeltaR < DELTAR);
  fGoodMKK    = ((MKKLO < fMKK ) && (fMKK < MKKHI)); 
  
  bmmReader::fillCandidateVariables();


}

// ----------------------------------------------------------------------
int bmmBs2JpsiPhiReader::tmCand2(TAnaCand *pC) {
  
  TAnaTrack *pT; 
  // -- Get gen tracks of this cand (if available)
  vector<TGenCand*> gCand; 
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    if (pT->fGenIndex < 0) return 0; 
    gCand.push_back(fpEvt->getGenCand(pT->fGenIndex)); 
  }

  TGenCand *pG(0), *pM1(0), *pM2(0), *pPsi(0), *pK1(0), *pK2(0), *pPhi(0), *pB(0); 
  for (unsigned int i = 0; i < gCand.size(); ++i) {
    pG = gCand[i]; 
    if (0 == pM1 && 13 == TMath::Abs(pG->fID)) {
      pM1 = pG; 
    } else if (0 == pM2 && 13 == TMath::Abs(pG->fID)) {
      pM2 = pG; 
    } else if (0 == pK1 && 321 == TMath::Abs(pG->fID)) {
      pK1 = pG; 
    } else if (0 == pK2 && 321 == TMath::Abs(pG->fID)) {
      pK2 = pG; 
    }
  }

  if (0 == pM1 || 0 == pM2 || 0 == pK1 || 0 == pK2) {
    if (fVerbose > 5) {
      if (0 == pM1) cout << "could not match m1" << endl;
      if (0 == pM2) cout << "could not match m2" << endl;
      if (0 == pK1) cout << "could not match k1" << endl;
      if (0 == pK2) cout << "could not match k2" << endl;
    }
    return 0; 
  }

  // -- check that the two muons come from the same J/psi
  if (pM1->fMom1 != pM2->fMom1) {
    if (fVerbose > 5) cout << "m1 has different mother than m2" << endl;
    return 0; 
  }
  pPsi = fpEvt->getGenCand(pM1->fMom1);
  if (TMath::Abs(pPsi->fID) != 443) {
    if (fVerbose > 5) cout << "m1 has mother != 443" << endl;
    return 0;
  }
  
  // -- check that the two kaons come from the same phi
  if (pK1->fMom1 != pK2->fMom1) {
    if (fVerbose > 5) cout << "k1 has different mother than k2" << endl;
    return 0; 
  }

  pPhi = fpEvt->getGenCand(pK1->fMom1);
  if (TMath::Abs(pPhi->fID) != 333) {
    if (fVerbose > 5) cout << "k1 has mother != 333" << endl;
    return 0;
  }

  // -- check that the phi and psi have the same mother
  if (pPsi->fMom1 != pPhi->fMom1) {
    if (fVerbose > 5) cout << "phi has mother != psi" << endl;
    return 0; 
  }

  pB = fpEvt->getGenCand(pPsi->fMom1); 
  if (TMath::Abs(pB->fID) != 531) {
    if (fVerbose > 5) cout << "B is != 531" << endl;
    return 0;
  }

  // -- Now check for the number of generator daughters
  int nphotons(0), ngendaughters(0);
  for (int i = pB->fDau1; i <= pB->fDau2; ++i) {
    if (22 == fpEvt->getGenCand(i)->fID) {
      ++nphotons;
      continue;
    }
    if (443 == fpEvt->getGenCand(i)->fID) continue;
    if (333 == fpEvt->getGenCand(i)->fID) continue;
    ++ngendaughters;
  }

  for (int i = pPsi->fDau1; i <= pPsi->fDau2; ++i) {
    if (22 == fpEvt->getGenCand(i)->fID) {
      ++nphotons;
      continue;
    }
    ++ngendaughters;
  }

  for (int i = pPhi->fDau1; i <= pPhi->fDau2; ++i) {
    if (22 == fpEvt->getGenCand(i)->fID) {
      ++nphotons;
      continue;
    }
    ++ngendaughters;
  }

  if (ngendaughters == 4) {
    if (nphotons > 0) {
      return 2;
    } else {
      if (fVerbose > 5) {
	cout << "tmCand: " << pB->fNumber << " " << pM1->fNumber << " " << pM2->fNumber << " " << pK1->fNumber << " " << pK2->fNumber << endl;
	cout << "        ";
	for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
	  pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
	  cout << pT->fIndex << " (" << pT->fPlab.Perp() << " " << pT->fPlab.Eta() << " " 
	       << fpEvt->getSigTrack(i)->fMCID
	       << ") " ;
	}
	cout << endl;
      }
      return 1; 
    }
  }
  
  return 0; 
}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::fillCandidateHistograms() {
  //  cout << "bmmBs2JpsiPhiReader::fillCandidateHistograms()" << endl;
  //   fAnaCuts.update(); 
  //   fAnaCuts.dumpAll(); 

  fpDeltaR->fill(fDeltaR, fCandM); 
  fpMKK->fill(fMKK, fCandM); 

  fpMpsi->fill(fJpsiMass, fCandM); 
  fpTracksPt->fill(fKa1Pt, fCandM);
  fpTracksPt->fill(fKa2Pt, fCandM);
  fpTracksEta->fill(fKa1Eta, fCandM);
  fpTracksEta->fill(fKa2Eta, fCandM);

  bmmReader::fillCandidateHistograms(); 
}


// ---------------------------------------------------------------------- 
void bmmBs2JpsiPhiReader::bookHist() {
  bmmReader::bookHist();
  cout << "==> bmmBs2JpsiPhiReader: bookHist " << endl;

  // -- Additional analysis distributions
  fpMpsi     = bookDistribution("mpsi", "m(J/#psi) [GeV]", "fGoodJpsiMass", 40, 2.8, 3.4);           
  fpDeltaR   = bookDistribution("deltar", "#Delta R", "fGoodDeltaR", 50, 0., 1.);           
  fpMKK      = bookDistribution("mkk", "m(KK) [GeV]", "fGoodMKK", 50, 0.95, 1.15);           

  // -- Additional reduced tree variables
  fTree->Branch("mpsi",  &fJpsiMass, "mpsi/D");
  fTree->Branch("mkk",   &fMKK,      "mkk/D");
  fTree->Branch("dr",    &fDeltaR,   "dr/D");

  fTree->Branch("k1pt",  &fKa1Pt,    "k1pt/D");
  fTree->Branch("k1eta", &fKa1Eta,   "k1eta/D");
  fTree->Branch("k1phi", &fKa1Phi,   "k1phi/D");
  fTree->Branch("k1gt",  &fKa1TkQuality,"k1gt/I");
  fTree->Branch("k2pt",  &fKa2Pt,    "k2pt/D");
  fTree->Branch("k2eta", &fKa2Eta,   "k2eta/D");
  fTree->Branch("k2phi", &fKa2Phi,   "k2phi/D");
  fTree->Branch("k2gt",  &fKa2TkQuality,"k2gt/I");

  fTree->Branch("t3pt",  &fKa1PtNrf, "t3pt/D");
  fTree->Branch("t3eta", &fKa1EtaNrf,"t3eta/D");

  fTree->Branch("t4pt",  &fKa2PtNrf, "t4pt/D");
  fTree->Branch("t4eta", &fKa2EtaNrf,"t4eta/D");

  fTree->Branch("g3pt", &fKa1PtGen,  "g3pt/D");
  fTree->Branch("g3eta",&fKa1EtaGen, "g3eta/D");
  fTree->Branch("g4pt", &fKa2PtGen,  "g4pt/D");
  fTree->Branch("g4eta",&fKa2EtaGen, "g4eta/D");

  // -- Additional effTree variables
  fEffTree->Branch("k1pt",   &fETk1pt,            "k1pt/F");
  fEffTree->Branch("g3pt",   &fETg3pt,            "g3pt/F");
  fEffTree->Branch("k1eta",  &fETk1eta,           "k1eta/F");
  fEffTree->Branch("g3eta",  &fETg3eta,           "g3eta/F");
  fEffTree->Branch("k1q",    &fETk1q,             "k1q/I");
  fEffTree->Branch("k1gt",   &fETk1gt,            "k1gt/O");

  fEffTree->Branch("k2pt",   &fETk2pt,            "k2pt/F");
  fEffTree->Branch("g4pt",   &fETg4pt,            "g4pt/F");
  fEffTree->Branch("k2eta",  &fETk2eta,           "k2eta/F");
  fEffTree->Branch("g4eta",  &fETg4eta,           "g4eta/F");
  fEffTree->Branch("k2q",    &fETk2q,             "k2q/I");
  fEffTree->Branch("k2gt",   &fETk2gt,            "k2gt/O");

}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::readCuts(TString filename, int dump) {
  bmmReader::readCuts(filename, dump); 

  fCutFile = filename;
  if (dump) cout << "==> bmmBs2JpsiPhiReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(string(fCutFile.Data()), cutLines);

  char  buffer[200];
  char CutName[100];
  float CutValue;
  int ok(0);

  TH1D *hcuts = (TH1D*)gFile->Get("hcuts"); 
  if (0 == hcuts) {
    new TH1D("hcuts", "", 1000, 0., 1000.);
    hcuts->GetXaxis()->SetBinLabel(1, fCutFile.Data());
  }
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 

    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
    
    if (!strcmp(CutName, "MKKLO")) {
      MKKLO = CutValue; ok = 1;
      if (dump) cout << "MKKLO:           " << MKKLO << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKKLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{min}(KK) :: %3.1f", CutName, MKKLO));
    }

    if (!strcmp(CutName, "MKKHI")) {
      MKKHI = CutValue; ok = 1;
      if (dump) cout << "MKKHI:           " << MKKHI << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKKHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{max}(KK) :: %3.1f", CutName, MKKHI));
    }

    if (!strcmp(CutName, "DELTAR")) {
      DELTAR = CutValue; ok = 1;
      if (dump) cout << "DELTAR:           " << DELTAR << endl;
      ibin = 301;
      hcuts->SetBinContent(ibin, DELTAR);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #Delta R(KK) :: %3.1f", CutName, DELTAR));
    }
    
    ok = checkCut(CutName, hcuts); 
    if (!ok) cout << "==> bmmBs2JpsiPhiReader: error? nothing done with " << CutName << "!!" << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}


