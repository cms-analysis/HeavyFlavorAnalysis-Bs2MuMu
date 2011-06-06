#include "bmmNormalizationReader.hh"
#include "TRandom.h"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using std::string;
using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
bmmNormalizationReader::bmmNormalizationReader(TChain *tree, TString evtClassName): bmmReader(tree, evtClassName) {
  cout << "==> bmmNormalizationReader: constructor..." << endl;
}

// ----------------------------------------------------------------------
bmmNormalizationReader::~bmmNormalizationReader() {
  cout << "==> bmmNormalizationReader: destructor..." << endl;

}

// ----------------------------------------------------------------------
void bmmNormalizationReader::moreBasicCuts() {
  cout << "   bmmNormalizationReader: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodJpsiMass", "m(J/psi)", fGoodJpsiMass); 
}

// ----------------------------------------------------------------------
void bmmNormalizationReader::startAnalysis() {
  bmmReader::startAnalysis();
  cout << "==> bmmNormalizationReader: Summary of analysis cuts:" << endl;
  fAnaCuts.dumpAll(); 
}


// ----------------------------------------------------------------------
void bmmNormalizationReader::initVariables() {
  bmmReader::initVariables();
 
}

// ----------------------------------------------------------------------
void bmmNormalizationReader::eventProcessing() {
  bmmReader::eventProcessing();

}

// ----------------------------------------------------------------------
void bmmNormalizationReader::MCKinematics() {
  fGoodMCKinematics = true; 
  TGenCand *pC(0), *pM1(0), *pM2(0), *pK(0), *pB(0), *pPsi(0); 
  int nphotons(0), npsi(0), nb(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (521 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = pB->fDau2 - pB->fDau1 + 1; 
      if (nb > 2) continue; // skip B decays where more than J/psi and kaon came from B

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
	} else if (321 == TMath::Abs(pC->fID)) {
	  pK = fpEvt->getGenCand(id); 
	}
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK) {
	goodMatch = true; 
	nphotons = 3 - nb - npsi;
	break;
      }
    }
  }
  
  if (fVerbose > 4) {
    if (goodMatch) {
      pB->dump(); 
      pM1->dump(); 
      pM2->dump();
      pK->dump();
    } else {
      cout << "no tm decay found" << endl;
    }
  }

  if (!goodMatch) {
    fGoodMCKinematics = false; 
    if (fVerbose > 2) cout << "No matched signal decay found" << endl;
    return;
  }


  // ----------------------------------------------------------------------
  // -- Acceptance: 
  TAnaTrack *pT(0), *prM1(0), *prM2(0), *prK(0); 
  int m1Matched(0), m2Matched(0), KMatched(0);
  int m1Acc(0), m2Acc(0), KAcc(0);
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
    if (pT->fGenIndex == pK->fNumber) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(13); 
      prK = pT; 
      KMatched = 1; 
      if ((pT->fPlab.Perp() > 0.5) && (TMath::Abs(pT->fPlab.Eta()) < 2.4)) {
	KAcc = 1; 
	((TH1D*)fpHistFile->Get("acceptance"))->Fill(23); 
      }
    }
  }

  ((TH1D*)fpHistFile->Get("acceptance"))->Fill(30); 
  ((TH1D*)fpHistFile->Get("acceptance"))->Fill(1); // denominator
  if (m1Matched && m1Acc && m2Matched && m2Acc) {
    ((TH1D*)fpHistFile->Get("acceptance"))->Fill(31); 
    if (KMatched && KAcc) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(32); 
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(2);  // numerator
    }
  }

  // ----------------------------------------------------------------------
  // -- preselection efficiency
  TAnaCand *pCand; 
  ((TH1D*)fpHistFile->Get("presel"))->Fill(10); 
  if (m1Matched && m1Acc && m2Matched && m2Acc && KMatched && KAcc) {
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
	if (pT->fGenIndex == pK->fNumber) {
	  ++tm;
	}
      }

      // -- all signal tracks must be matched
      if (tm != (pCand->fSig2-pCand->fSig1+1)) {
	//	cout << "Did not find matched candidate " << pCand->fType << " at " << iC << endl;
	continue;
      }
      ((TH1D*)fpHistFile->Get("presel"))->Fill(20); 
      if (muonID(prM1) && muonID(prM2)) {
	((TH1D*)fpHistFile->Get("presel"))->Fill(21); 
	((TH1D*)fpHistFile->Get("presel"))->Fill(2);  // numerator
      }
      break;
    }
  }


  // -- hard coded ?! FIXME
  if (pM1->fP.Perp() < 1.0) fGoodMCKinematics = false;  
  if (pM2->fP.Perp() < 1.0) fGoodMCKinematics = false;  
  if (pK->fP.Perp()  < 0.5) fGoodMCKinematics = false;  
  if (TMath::Abs(pM1->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pM2->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pK->fP.Eta()) > 2.4) fGoodMCKinematics = false;  


  if (fEvt == 1382086) {
    fpEvt->dumpGenBlock();
  }

}



// ----------------------------------------------------------------------
void bmmNormalizationReader::efficiencyCalculation() {
  fGoodEffCand = false;

  ((TH1D*)fpHistFile->Get("efficiency"))->Fill(0); 
  ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(1, "all events"); 


  // -- gen level 
  TGenCand *pB(0), *pM1(0), *pM2(0), *pK(0); 
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi || -1 == fGenK1Tmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }
  pB  = fpEvt->getGenCand(fGenBTmi); 
  pM1 = fpEvt->getGenCand(fGenM1Tmi); 
  pM2 = fpEvt->getGenCand(fGenM2Tmi); 
  pK  = fpEvt->getGenCand(fGenK1Tmi); 

  ((TH1D*)fpHistFile->Get("efficiency"))->Fill(1); 
  ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(2, "gen signal decays"); 


  // -- reco level
  TAnaTrack *prM1(0), *prM2(0), *prK(0); 
  int m1Matched(0), m2Matched(0), kMatched(0), m1ID(0), m2ID(0), m1GT(0), m2GT(0), kGT(0);
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
    kMatched = 1; 
    prK = fpEvt->getRecTrack(fRecK1Tmi); 
    if (TRACKQUALITY > 0 && (0 == (prK->fTrackQuality & TRACKQUALITY))) {
      kGT = 0; 
    } else {
      kGT = 1;
    }
  } 

  // -- cand level 
  TAnaCand *pCand(0);
  if (fCandTmi > -1) {
    pCand = fCands[fCandTmi];
  }
    
  // -- EffTree filling for all events with a signal decay
  fETgpt   = pB->fP.Perp(); 
  fETgeta  = pB->fP.Eta(); 
  fETg1pt  = pM1->fP.Perp(); 
  fETg1eta = pM1->fP.Eta(); 
  fETg2pt  = pM2->fP.Perp(); 
  fETg2eta = pM2->fP.Eta(); 
  fETg3pt  = pK->fP.Perp(); 
  fETg3eta = pK->fP.Eta(); 
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
  if (kMatched) {
    fETk1pt  = prK->fPlab.Perp(); 
    fETk1eta = prK->fPlab.Eta(); 
    fETk1q   = prK->fQ;
    fETk1gt  = (kGT>0?true:false); 
  } else {
    fETk1pt  = -99.; 
    fETk1eta = -99.; 
    fETk1q   = -99;
    fETk1gt  = false; 
  }
  if (pCand) {
    fETcandMass = pCand->fMass; 
  } else {
    fETcandMass = -99.;
  }
  fEffTree->Fill(); 


  // -- results...
  if ((TMath::Abs(pM1->fP.Eta()) < 2.5) && (TMath::Abs(pM2->fP.Eta()) < 2.5) && (TMath::Abs(pK->fP.Eta()) < 2.5)) {
    ((TH1D*)fpHistFile->Get("efficiency"))->Fill(2); 
    ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(3, "+ eta cuts"); 
    
    if ((pM1->fP.Perp() > 1.) && (pM2->fP.Perp() > 1.) && (pK->fP.Perp() > 0.4)) {
      ((TH1D*)fpHistFile->Get("efficiency"))->Fill(3); 
      ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(4, "+ pT cuts"); 
      
      if (m1Matched && m2Matched && kMatched
	  && (prM1->fPlab.Perp() > 3.0) && (prM2->fPlab.Perp() > 3.0) && (prK->fPlab.Perp() > 0.5)
	  && (TMath::Abs(prM1->fPlab.Eta()) < 2.4) && (TMath::Abs(prM2->fPlab.Eta()) < 2.4) && (TMath::Abs(prK->fPlab.Eta()) < 2.4)
	  && (prM1->fQ*prM2->fQ < 0)
	  && m1GT && m2GT && kGT
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
int bmmNormalizationReader::partialReco(TAnaCand *pCand) {
  struct bcounter{
    int mcid;
    int index;
    int matches;
  };  

  //  cout << "evt: " << fEvt << " event: " << fEvent << endl;

  bcounter bco[] = {{511, -1, 0}, {511, -1, 0}, {521, -1, 0}, {521, -1, 0}, {531, -1, 0}, {531, -1, 0}, {5122, -1, 0}, {5122, -1, 0}};

  int idx(0), type(0), bIdx(0), bId(0), nsame(0), nother(0), pr(0); 
  TAnaTrack *pT(0); 
  TGenCand *pG(0); 
  for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
    idx = fpEvt->getSigTrack(i)->fIndex; 
    type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
    pT = fpEvt->getRecTrack(idx);
    if (pT->fGenIndex <0) continue;
    pG = fpEvt->getGenCand(pT->fGenIndex);
    bIdx = fromB(pG); 
    bId = 0; 
    if (bIdx > -1) {
      bId =  TMath::Abs(fpEvt->getGenCand(bIdx)->fID);
    } else {
      continue;
    }

    for (int ib = 0 ; ib < 8; ++ib) {
      if (bco[ib].mcid != bId) continue; 
      if (bco[ib].index == bIdx) {
	//cout << "adding old to bco[" << ib << "]: "  << bId << " at " << bIdx << " for track " << i << endl;
	bco[ib].matches++;
	break;
      } else {
	//cout << "not adding to bco[" << ib << "].index = " << bco[ib].index << " bId = " << bId << " at " << bIdx << " for track " << i << endl;
      }

      if (bco[ib].index < 0) {
	//cout << "adding new to bco[" << ib << "]: " << bId << " at " << bIdx << " for track " << i << endl;
	bco[ib].index = bIdx; 
	bco[ib].matches++;
	break;
      }
    }
  }

  for (int i = 0 ; i < 8; ++i) {
    if (bco[i].mcid == 511 && bco[i].matches == pCand->fSig2 - pCand->fSig1 + 1) {
      //cout << "-> matched " << bco[i].matches << " tracks to a Bz at " << bco[i].index << endl;
      return 511;
    }

    if (bco[i].mcid == 521 && bco[i].matches == pCand->fSig2 - pCand->fSig1 + 1) {
      //cout << "-> matched " << bco[i].matches << " tracks to a B+ at " << bco[i].index << endl;
      return 521;
    }

    if (bco[i].mcid == 531 && bco[i].matches == pCand->fSig2 - pCand->fSig1 + 1) {
      //cout << "-> matched " << bco[i].matches << " tracks to a Bs at " << bco[i].index << endl;
      return 531;
    }

    if (bco[i].mcid == 5122 && bco[i].matches == pCand->fSig2 - pCand->fSig1 + 1) {
      //cout << "-> matched " << bco[i].matches << " tracks to a LambdaB at " << bco[i].index << endl;
      return 5122;
    }
    
  }

  return 0; 
}


// ----------------------------------------------------------------------
int bmmNormalizationReader::fromB(TGenCand *pCand) {
  TGenCand *pMom(0);
  int idx = -1; 
  
  //  cout << pCand << endl;
  if (pCand && pCand->fMom1 > -1) {
    //    cout << "fromB: " << pCand->fMom1 << endl;
    pMom = fpEvt->getGenCand(pCand->fMom1);
  } else {
    return -1;
  } 

  while (pMom) {
    idx = pMom->fNumber; 
    if (521  == TMath::Abs(pMom->fID)) break;
    if (511  == TMath::Abs(pMom->fID)) break;
    if (531  == TMath::Abs(pMom->fID)) break;
    if (5122 == TMath::Abs(pMom->fID)) break;
    if (pMom->fMom1 > -1) {
      pMom = fpEvt->getGenCand(pMom->fMom1);
    }
    if (pMom->fNumber < 10) {
      idx = -1; 
      break;
    }
  }

  return idx; 
}

// ----------------------------------------------------------------------
void bmmNormalizationReader::genMatch() {

  fGenM1Tmi = fGenM2Tmi = fGenK1Tmi = -1; 
  fNGenPhotons = 0; 

  TGenCand *pC(0), *pM1(0), *pM2(0), *pK(0), *pB(0), *pPsi(0); 
  int nb(0), npsi(0), nphotons(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (521 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = 0; 
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenCand(id); 
	if (443 == TMath::Abs(pC->fID)) {
	  ++nb;
	  pPsi = pC; 
	  npsi = pPsi->fDau2 - pPsi->fDau1 + 1; 
	  pM1 = pM2 = 0; 
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenCand(idd); 
	    if (22 == TMath::Abs(pC->fID)) ++nphotons;
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenCand(idd); 
	      } else {
		pM2 = fpEvt->getGenCand(idd); 
	      }
	    }
	  }
	} else if (321 == TMath::Abs(pC->fID)) {
	  ++nb;
	  pK = fpEvt->getGenCand(id); 
	} else 	if (22 == TMath::Abs(pC->fID)) {
	  ++nphotons;
	} else {
	  ++nb;
	}
      }
      if (nb > 2) {
	pM1 = pM2 = pK = pPsi = 0; 
	continue; // skip B decays where more than J/psi and kaon came from B
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK && (pPsi->fMom1 == pK->fMom1)) {
	goodMatch = true; 
	fNGenPhotons = nphotons;
	break;
      }
    }
  }

  if (goodMatch) {
    fGenBTmi = pB->fNumber; 
    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fGenM1Tmi = pM1->fNumber; 
      fGenM2Tmi = pM2->fNumber; 
    } else {
      fGenM1Tmi = pM2->fNumber; 
      fGenM2Tmi = pM1->fNumber; 
    }
    fGenK1Tmi = pK->fNumber; 
  } else {
    fGenM1Tmi = -1; 
    fGenM2Tmi = -1; 
    fGenK1Tmi = -1;  
  }
  
  if (fVerbose > 10) {
    cout << "fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fGenM2Tmi = " << fGenM2Tmi << endl;
    cout << "fGenK1Tmi = " << fGenK1Tmi << endl;
  }

  if (!goodMatch && evtFoundInCN(fEvt)) {
    cout << "-> Event: " << fEvt << endl;
  }


//   if (fGenM1Tmi > -1) {
//     cout << "----------------------------------------------------------------------" << endl;
//     cout << "fGenM1Tmi = " << fGenM1Tmi << endl;
//     cout << "fGenM2Tmi = " << fGenM2Tmi << endl;
//     cout << "fGenK1Tmi = " << fGenK1Tmi << endl;
//   }
}

// ----------------------------------------------------------------------
void bmmNormalizationReader::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = fRecK1Tmi = -1; 
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

    // -- kaon
    if (fGenK1Tmi > -1 && pT->fGenIndex == fGenK1Tmi) {
      fRecK1Tmi = i; 
    }

    // -- skip rest if all matches found
    if (fRecM1Tmi > -1 && fRecM2Tmi > -1 && fRecK1Tmi > -1) break;
  }


  if (fVerbose > 10) {
    cout << "fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi << endl;
    cout << "fRecK1Tmi = " << fRecK1Tmi << " matched to fGenK1Tmi = " << fGenK1Tmi << endl;
  }

//   if (fRecM1Tmi > -1) {
//     cout << "fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi << endl;
//     cout << "fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi << endl;
//     cout << "fRecK1Tmi = " << fRecK1Tmi << " matched to fGenK1Tmi = " << fGenK1Tmi << endl;
//   }
  
}

// ----------------------------------------------------------------------
void bmmNormalizationReader::candMatch() {
  fCandTmi = -1;   
  int idx(-1), type(-1); 
  int d1Matched(0), d2Matched(0), d3Matched(0); 
  TAnaCand *pCand(0);
  for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
    pCand = fCands[iC]; 
    
    d1Matched = d2Matched = d3Matched = 0; 
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex; 
      type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
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
    }
    
    if (d1Matched && d2Matched && d3Matched) {
      fCandTmi = iC;
      break;
    }
  }
  if (fVerbose > 10) {
    cout << "fCandTmi = " << fCandTmi << " matched to rec tracks " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi << endl;
  }

//   if (fGenM1Tmi > -1) cout << "fCandTmi = " << fCandTmi << " matched to rec tracks " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi << endl;

}


// ----------------------------------------------------------------------
int bmmNormalizationReader::tmCand(TAnaCand *pC) {
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
void bmmNormalizationReader::candidateSelection(int mode) {
  int nc0(fCands.size()), nc1(0);

  fpCand = 0; 
  ((TH1D*)fpHistFile->Get("bnc0"))->Fill(nc0); 
  
  if (0 == nc0) {
    nc1 = 0; 
    fCandIdx = -1; 
    ((TH1D*)fpHistFile->Get("bnc1"))->Fill(nc1); 
    if (fVerbose > 0) cout << "bmmNormalizationReader> no candidate " << TYPE << " found" << endl;
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
    fpCand = fCands[clist[0]];
    fCandIdx = clist[0]; 
    nc1 = 1; 
    if (fVerbose > 0) cout << "bmmNormalizationReader> found exactly one candidate " << TYPE << " passing the constituents selection" << endl;
  } else if (clist.size() > 1) {
    // -- SEVERAL candidate fulfilled the basic requirements. Choose one. 
    double drMin(99.), dr(0.); 
    for (unsigned int i = 0; i < clist.size(); ++i) {
      pCand = fCands[clist[i]]; 
      if (BLIND && pCand->fMass > SIGBOXMIN && pCand->fMass < SIGBOXMAX) continue; 
      TAnaTrack *pl1 = fpEvt->getSigTrack(pCand->fSig1); 
      TAnaTrack *pl2 = fpEvt->getSigTrack(pCand->fSig1+1); 
      
      // -- if all candidates are taken, then fpCand at this point is pointing to a random cand
      fpCand = pCand; 
      fCandIdx = clist[i]; 
      // -- closest in rphi. This is pointless for this mode ...
      if (0 == mode) {
	dr = pl1->fPlab.DeltaR(pl2->fPlab); 
	if (dr < drMin) {
	  fpCand = pCand; 
	  fCandIdx = clist[i]; 
	  drMin = dr; 
	  if (fVerbose > 0) cout << "bmmNormalizationReader> found another better candidate " << TYPE << " with dr = " << dr << endl;
	}
      }
      // -- higher pT
      if (1 == mode) {
	dr = pCand->fPlab.Perp(); 
	if (dr > drMin) {
	  fpCand = pCand; 
	  fCandIdx = clist[i]; 
	  drMin = dr; 
	  if (fVerbose > 0) cout << "bmmNormalizationReader> found another better candidate " << TYPE << " with pT = " << dr << endl;
	}
      }
    }
  } else {
    // -- NO candidate fulfilled the basic requirements. Then take the first one. 
    fpCand = fCands[0];
    fCandIdx = 0; 
    fvGoodCand[0] = true; 
    nc1 = 1; 
    if (fVerbose > 0) cout << "bmmNormalizationReader> found no candidate " << TYPE << " passing the constituents selection" << endl;
  }

  ((TH1D*)fpHistFile->Get("bnc1"))->Fill(nc1); 

  //  if (fpCand) fillCandidateVariables();
}


// ----------------------------------------------------------------------
void bmmNormalizationReader::fillCandidateVariables() {
  if (0 == fpCand) return;

  TAnaTrack *p0, *pk(0), *pks(0); 
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    if (321 == TMath::Abs(p0->fMCID)) {
      pks = p0; 
    }
  }

  pk = fpEvt->getRecTrack(pks->fIndex);

  fKaonPt        = pk->fPlab.Perp(); 
  fKaonEta       = pk->fPlab.Eta();  
  fKaonPhi       = pk->fPlab.Phi(); 
  fKaonTkQuality = pk->fTrackQuality & TRACKQUALITY;
  fKaonPtNrf     = pks->fPlab.Perp();
  fKaonEtaNrf    = pks->fPlab.Eta();

  
  if (fIsMC) {
    fGenBpartial = partialReco(fpCand); 
    if (fpCand->fMass > 5.5 && fGenBpartial == 511) fpEvt->dumpGenBlock();
  }

  if (tmCand(fpCand)) {
    TGenCand *pg1 = fpEvt->getGenCand(fpEvt->getRecTrack(pk->fIndex)->fGenIndex);
    fKPtGen     = pg1->fP.Perp();
    fKEtaGen    = pg1->fP.Eta();
  } else {
    fKPtGen     = -99.;
    fKEtaGen    = -99.;
  }
  

  // -- Check for J/psi mass
  TAnaCand *pD = 0; 
  fGoodJpsiMass = false; 
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i); 
    if (pD->fType == JPSITYPE) {
      if ((JPSIMASSLO < pD->fMass) && (pD->fMass < JPSIMASSHI)) fGoodJpsiMass = true;
      fJpsiMass = pD->fMass;
      break;
      //       cout << "type = " << pD->fType 
      // 	   << " with mass = " << pD->fMass 
      // 	   << " fGoodJpsiMass = " << fGoodJpsiMass 
      // 	   << endl;
    }
  }
  
  bmmReader::fillCandidateVariables();

  fPreselection = fPreselection && fGoodJpsiMass;
}

// ----------------------------------------------------------------------
int bmmNormalizationReader::tmCand2(TAnaCand *pC) {
  
  TAnaTrack *pT; 
  vector<TGenCand*> gCand; 
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    if (pT->fGenIndex < 0) return 0; 
    gCand.push_back(fpEvt->getGenCand(pT->fGenIndex)); 
  }

  TGenCand *pG(0), *pM1(0), *pM2(0), *pPsi(0), *pK(0), *pB(0); 
  for (unsigned int i = 0; i < gCand.size(); ++i) {
    pG = gCand[i]; 
    if (0 == pM1 && 13 == TMath::Abs(pG->fID)) {
      pM1 = pG; 
    } else if (0 == pM2 && 13 == TMath::Abs(pG->fID)) {
      pM2 = pG; 
    }
    if (0 == pK && 321 == TMath::Abs(pG->fID)) {
      pK = pG; 
    }
  }

  if (0 == pM1 || 0 == pM2 || 0 == pK) return 0; 

  // -- check that the two muons come from the same J/psi
  if (pM1->fMom1 != pM2->fMom1) return 0; 
  pPsi = fpEvt->getGenCand(pM1->fMom1);
  if (TMath::Abs(pPsi->fID) != 443) return 0;
  
  // -- check that the kaon has the same mother like the J/psi
  if (pK->fMom1 != pPsi->fMom1) return 0; 
  pB = fpEvt->getGenCand(pK->fMom1); 

  // -- Now check for the number of generator daughters
  int nphotons(0), ngendaughters(0);
  for (int i = pB->fDau1; i <= pB->fDau2; ++i) {
    if (22 == fpEvt->getGenCand(i)->fID) {
      ++nphotons;
      continue;
    }
    if (443 == fpEvt->getGenCand(i)->fID) continue;
    ++ngendaughters;
  }

  for (int i = pPsi->fDau1; i <= pPsi->fDau2; ++i) {
    if (22 == fpEvt->getGenCand(i)->fID) {
      ++nphotons;
      continue;
    }
    ++ngendaughters;
  }

  if (3 == ngendaughters) {
    if (nphotons > 0) {
      return 2;
    } else {
      return 1;
    }
  }

  return -1; 
}


// ----------------------------------------------------------------------
void bmmNormalizationReader::fillCandidateHistograms() {

  fpMpsi->fill(fJpsiMass, fCandM); 
  fpTracksPt->fill(fKaonPt, fCandM);
  fpTracksEta->fill(fKaonEta, fCandM);

  bmmReader::fillCandidateHistograms(); 

}

// ---------------------------------------------------------------------- 
void bmmNormalizationReader::bookHist() {
  bmmReader::bookHist();

  // -- Additional analysis distributions
  fpMpsi   = bookDistribution("mpsi", "m(J/#psi) [GeV]", "fGoodJpsiMass", 40, 2.8, 3.4);           

  // -- Additional reduced tree variables
  fTree->Branch("mpsi", &fJpsiMass,  "mpsi/D");
  fTree->Branch("kpt",  &fKaonPt,    "kpt/D");
  fTree->Branch("keta", &fKaonEta,   "keta/D");
  fTree->Branch("kphi", &fKaonPhi,   "kphi/D");
  fTree->Branch("kgt",  &fKaonTkQuality,"kgt/I");
  fTree->Branch("t3pt", &fKaonPtNrf, "t3pt/D");
  fTree->Branch("t3eta",&fKaonEtaNrf,"t3eta/D");
  fTree->Branch("g3pt", &fKPtGen,    "g3pt/D");
  fTree->Branch("g3eta",&fKEtaGen,   "g3eta/D");

  // -- Additional effTree variables
  fEffTree->Branch("k1pt",   &fETk1pt,            "k1pt/F");
  fEffTree->Branch("g3pt",   &fETg3pt,            "g3pt/F");
  fEffTree->Branch("k1eta",  &fETk1eta,           "k1eta/F");
  fEffTree->Branch("g3eta",  &fETg3eta,           "g3eta/F");
  fEffTree->Branch("k1q",    &fETk1q,             "k1q/I");
  fEffTree->Branch("k1gt",   &fETk1gt,            "k1gt/O");
}


// ----------------------------------------------------------------------
void bmmNormalizationReader::readCuts(TString filename, int dump) {
  bmmReader::readCuts(filename, dump); 

  fCutFile = filename;
  if (dump) cout << "==> bmmNormalizationReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
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
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 

    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    ok = checkCut(CutName, hcuts); 
    if (!ok) cout << "==> bmmNormalizationReader: error? nothing done with " << CutName << "!!" << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}
