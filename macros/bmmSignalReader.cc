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
  bmmReader::startAnalysis();
  cout << "==> bmmSignalReader: Summary of analysis cuts:" << endl;
  fAnaCuts.dumpAll(); 
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
  ((TH1D*)fpHistFile->Get("genStudy"))->Fill(1); 
  fGoodMCKinematics = true; 
  TGenCand *pC(0), *pM1(0), *pM2(0), *pB(0); 
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
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }


  // ----------------------------------------------------------------------
  // -- generator-level acceptance numbers
  if ((pM1->fP.Perp() > 0.) && (pM2->fP.Perp() > 0.)) {
    ((TH1D*)fpHistFile->Get("genStudy"))->Fill(1); 
  }
  if ((pM1->fP.Perp() > 1.) && (pM2->fP.Perp() > 1.)) {
    ((TH1D*)fpHistFile->Get("genStudy"))->Fill(2); 
  }
  if ((pM1->fP.Perp() > 1.) && (pM2->fP.Perp() > 1.)
      && (TMath::Abs(pM1->fP.Eta()) < 2.4) && (TMath::Abs(pM2->fP.Eta()) < 2.4)) {
    ((TH1D*)fpHistFile->Get("genStudy"))->Fill(3); 
  }

  // ----------------------------------------------------------------------
  // -- Acceptance: 
  TAnaTrack *pT(0), *prM1(0), *prM2(0); 
  int m1Matched(0), m2Matched(0);
  int m1Acc(0), m2Acc(0);
  int m1ID(0), m2ID(0);
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    pT = fpEvt->getRecTrack(i); 
    if (pT->fGenIndex == pM1->fNumber) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(11); 
      prM1 = pT; 
      m1Matched = 1; 
      if ((pT->fPlab.Perp() > 1.) && (TMath::Abs(pT->fPlab.Eta()) < 2.4)) {
	m1Acc = 1; 
	((TH1D*)fpHistFile->Get("acceptance"))->Fill(21); 
	if (muonID(prM1)) {
	  m1ID = 1; 
	}
      }
    }
    if (pT->fGenIndex == pM2->fNumber) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(12); 
      prM2 = pT; 
      m2Matched = 1; 
      if ((pT->fPlab.Perp() > 1.) && (TMath::Abs(pT->fPlab.Eta()) < 2.4)) {
	m2Acc = 1; 
	((TH1D*)fpHistFile->Get("acceptance"))->Fill(22); 
	if (muonID(prM2)) {
	  m2ID = 1; 
	}
      }
    }
  }

  ((TH1D*)fpHistFile->Get("acceptance"))->Fill(30); 
  ((TH1D*)fpHistFile->Get("acceptance"))->Fill(1); // denominator
  if (m1Matched && m1Acc && m2Matched && m2Acc) {
    ((TH1D*)fpHistFile->Get("acceptance"))->Fill(31); 
    ((TH1D*)fpHistFile->Get("acceptance"))->Fill(2);  // numerator
    if (m1ID && m2ID) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(3); 
    }
    if (prM2->fPlab.Perp() > 3 && prM1->fPlab.Perp() > 3) {
      ((TH1D*)fpHistFile->Get("acceptance"))->Fill(5); 
      if (m1ID && m2ID) {
	((TH1D*)fpHistFile->Get("acceptance"))->Fill(6); 
      }
    }
  }
  
  

  // ----------------------------------------------------------------------
  // -- preselection efficiency
  TAnaCand *pCand; 
  ((TH1D*)fpHistFile->Get("presel"))->Fill(10); 
  if (m1Matched && m1Acc && m2Matched && m2Acc) {
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
  if (pM1->fP.Perp() < 1.0) fGoodMCKinematics = false;  
  if (pM2->fP.Perp() < 1.0) fGoodMCKinematics = false;  
  if (TMath::Abs(pM1->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pM2->fP.Eta()) > 2.4) fGoodMCKinematics = false;  

}


// ----------------------------------------------------------------------
void bmmSignalReader::efficiencyCalculation() {
  fGoodEffCand = false;

  ((TH1D*)fpHistFile->Get("efficiency"))->Fill(0); 
  ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(1, "all events"); 

  // -- gen level 
  TGenCand *pM1(0), *pM2(0); 
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }

  pM1 = fpEvt->getGenCand(fGenM1Tmi); 
  pM2 = fpEvt->getGenCand(fGenM2Tmi); 

  ((TH1D*)fpHistFile->Get("efficiency"))->Fill(1); 
  ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(2, "gen signal decays"); 


  // -- reco level
  TAnaTrack *prM1(0), *prM2(0); 
  int m1Matched(0), m2Matched(0), m1ID(0), m2ID(0), m1GT(0), m2GT(0);
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

  // -- results...
  if ((TMath::Abs(pM1->fP.Eta()) < 2.5) && (TMath::Abs(pM2->fP.Eta()) < 2.5)) {
    ((TH1D*)fpHistFile->Get("efficiency"))->Fill(2); 
    ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(3, "+ eta cuts"); 
    
    if ((pM1->fP.Perp() > 1.) && (pM2->fP.Perp() > 1.)) {
      ((TH1D*)fpHistFile->Get("efficiency"))->Fill(3); 
      ((TH1D*)fpHistFile->Get("efficiency"))->GetXaxis()->SetBinLabel(4, "+ pT cuts"); 
      
      if (m1Matched && m2Matched
	  && (prM1->fPlab.Perp() > 3.0) && (prM2->fPlab.Perp() > 3.0)
	  && (TMath::Abs(prM1->fPlab.Eta()) < 2.4) && (TMath::Abs(prM2->fPlab.Eta()) < 2.4)
	  && (prM1->fQ*prM2->fQ < 0)
	  && m1GT && m2GT
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
void bmmSignalReader::candidateSelection(int mode) {
  int nc0(fCands.size()), nc1(0);

  fpCand = 0; 
  fCandIdx = -1; 
  ((TH1D*)fpHistFile->Get("bnc0"))->Fill(nc0); 
  
  if (0 == nc0) {
    nc1 = 0; 
    fCandIdx = -1; 
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
      if (fVerbose > 0) cout << "bmmSignalReader> found exactly one candidate " << TYPE 
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
	if (fVerbose > 0) cout << "bmmSignalReader> changing to candidate at " << clist[i] << endl;
	fpCand = pCand; 
	fCandIdx = clist[i]; 
      }
      // -- closest in rphi
      if (0 == mode) {
	dr = pl1->fPlab.DeltaR(pl2->fPlab); 
	if (dr < drMin) {
	  fpCand = pCand; 
	  fCandIdx = clist[i]; 
	  fvGoodCand[clist[i]] = true; 
	  drMin = dr; 
	  if (fVerbose > 0) cout << "bmmSignalReader> found another better candidate " << TYPE 
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
	  if (fVerbose > 0) cout << "bmmSignalReader> found another better candidate " << TYPE << " with pT = " << dr << endl;
	}
      } 
    }
  } else {
    // -- NO candidate fulfilled the basic requirements. Then take the first one
    pCand = fCands[0]; 
    if (BLIND && pCand->fMass > SIGBOXMIN && pCand->fMass < SIGBOXMAX) {
      if (fVerbose > 0) cout << "bmmSignalReader> found one candidate " << TYPE 
			     << " passing none of the constituents selection, filling none"
			     << " with mass  " << pCand->fMass 
			     << endl;
      fvGoodCand[0] = false; 
      fpCand = 0;
      fCandIdx = -1; 
    } else {
      fpCand = fCands[0];
      fCandIdx = 0; 
      fvGoodCand[0] = true; 
      nc1 = 1; 
      if (fVerbose > 0) cout << "bmmSignalReader> found no candidate " << TYPE 
			     << " passing the constituents selection, filling 0"
			     << " with signal tracks " << fpCand->fSig1 << ".." << fpCand->fSig2
			     << endl;
    }
  }
  ((TH1D*)fpHistFile->Get("bnc1"))->Fill(nc1); 

}


// ----------------------------------------------------------------------
void bmmSignalReader::genMatch() {

  fGenM1Tmi = fGenM2Tmi = -1; 
  fNGenPhotons = 0; 

  TGenCand *pC(0), *pM1(0), *pM2(0), *pB(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (531 == TMath::Abs(pC->fID)) {
      pM1 = pM2 = 0; 
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
	fNGenPhotons = pB->fDau2 - pB->fDau1 - 1; 
	if (fVerbose > 10) {
	  cout << "found gen match for B gen idx = " << pB->fNumber << endl;
	}
	break;
      }
    }
  }

  
  if (goodMatch) {
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
void bmmSignalReader::recoMatch() {

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
void bmmSignalReader::candMatch() {
  fCandTmi = -1;   
  int idx(-1); 
  int d1Matched(0), d2Matched(0); 
  TAnaCand *pCand(0);
  TAnaTrack *pT(0); 
  for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
    pCand = fCands[iC]; 
    
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
int bmmSignalReader::tmCand(TAnaCand *pC) {
  TAnaCand *pCand(0);
  for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
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
int bmmSignalReader::tmCand2(TAnaCand *pC) {
  // this is just not the correct way how to truth math!!
  TAnaTrack *pT; 
  vector<TGenCand*> gCand; 
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    if (pT->fGenIndex < 0) return 0; 
    gCand.push_back(fpEvt->getGenCand(pT->fGenIndex)); 
  }

  TGenCand *pG, *pM, *pB(0); 
  int matched(0), genDaughters(0); 
  for (unsigned int i = 0; i < gCand.size(); ++i) {
    pG = gCand[i]; 
    if (fVerbose > 10) cout << pG->fID << endl;
    if (pG->fMom1 > -1) {
      pM = fpEvt->getGenCand(pG->fMom1); 
      cout << pG->fNumber << "  " << pG->fID << "     " << pM->fNumber << "  " << pM->fID << endl;
    } else {
      return 0; 
    }
    if (531 != TMath::Abs(pM->fID) ) {
      return 0;  
    } else {
      cout << "    with mother 531 at " << pM->fNumber << endl;
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
  cout << "==> bmmSignalReader: bookHist " << endl;
  TH1D *h1; 
  h1 = new TH1D("h100", "", 100, 0., 100.);
  fTree->Branch("m2",      &fCandM,  "m2/D");
}


// ----------------------------------------------------------------------
void bmmSignalReader::readCuts(TString filename, int dump) {
  bmmReader::readCuts(filename, dump); 

  fCutFile = filename;
  if (dump) cout << "==> bmmSignalReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
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
    if (!ok) cout << "==> bmmSignalReader: error? nothing done with " << CutName << "!!" << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}


