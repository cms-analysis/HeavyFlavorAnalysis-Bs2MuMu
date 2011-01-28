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
      nb = pB->fDau2 - pB->fDau1; 
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
    cout << "----------------------------------------------------------------------" << endl;
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
  TAnaTrack *pT(0), *prM1(0), *prM2(0); 
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
	  cout << "prM1 = " << prM1 << " pT = " << pT << endl;
	  prM1 = pT; 
	  ++tm;
	}
	if (pT->fGenIndex == pM2->fNumber) {
	  cout << "prM2 = " << prM2 << " pT = " << pT << endl;
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
  if (pK->fP.Perp()  < 0.5) fGoodMCKinematics = false;  
  if (TMath::Abs(pM1->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pM2->fP.Eta()) > 2.4) fGoodMCKinematics = false;  
  if (TMath::Abs(pK->fP.Eta()) > 2.4) fGoodMCKinematics = false;  




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
      
      TAnaTrack *pl1 = fpEvt->getSigTrack(pCand->fSig1); 
      TAnaTrack *pl2 = fpEvt->getSigTrack(pCand->fSig1+1); 

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

  TAnaTrack *p1 = fpEvt->getSigTrack(fpCand->fSig2); 

  fKaonPt  = p1->fPlab.Perp(); 
  fKaonEta = p1->fPlab.Eta();  
  fKaonPhi = p1->fPlab.Phi(); 

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
}

// ----------------------------------------------------------------------
int bmmNormalizationReader::tmCand(TAnaCand *pC) {
  
  TAnaTrack *pT; 
  vector<TGenCand*> gCand; 
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    if (pT->fGenIndex < 0) return 0; 
    gCand.push_back(fpEvt->getGenCand(pT->fGenIndex)); 
  }

  TGenCand *pG(0), *pM1(0), *pM2(0), *pPsi(0), *pK(0), *pB; 
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
  // FIXME: ancestor instead of mother!
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


  if (ngendaughters == 3) {
    return 1;
  } else {
    if (nphotons > 0) {
      return 3;
    } else {
      return 2; 
    }
  }

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
  fTree->Branch("mpsi", &fJpsiMass, "mpsi/D");
  fTree->Branch("kpt",  &fKaonPt,   "kpt/D");
  fTree->Branch("keta", &fKaonEta,  "keta/D");
  fTree->Branch("kphi", &fKaonPhi,  "kphi/D");

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
