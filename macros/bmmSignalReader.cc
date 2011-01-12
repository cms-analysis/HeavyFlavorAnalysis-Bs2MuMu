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
    if (fVerbose > -1 ) cout << "--------------------> No matched signal decay found" << endl;
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
  TAnaTrack *pT, *prM1, *prM2; 
  int m1Matched(0), m2Matched(0);
  int m1Acc(0), m2Acc(0);
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
  }

  ((TH1D*)fpHistFile->Get("acceptance"))->Fill(30); 
  if (m1Matched && m1Acc && m2Matched && m2Acc) {
    ((TH1D*)fpHistFile->Get("acceptance"))->Fill(31); 
  }
  

  // ----------------------------------------------------------------------
  // -- preselection efficiency
  TAnaCand *pCand; 
  ((TH1D*)fpHistFile->Get("presel"))->Fill(10); 
  if (m1Matched && m1Acc && m2Matched && m2Acc) {
    ((TH1D*)fpHistFile->Get("presel"))->Fill(11); 
    if (muonID(prM1) && muonID(prM2)) {
      ((TH1D*)fpHistFile->Get("presel"))->Fill(12); 
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

  //  if (fpCand) fillCandidateVariables();
}


// ----------------------------------------------------------------------
int bmmSignalReader::tmCand(TAnaCand *pC) {
  
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


