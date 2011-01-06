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
}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::insertCand(TAnaCand* pCand) {
  //  cout << "bmmBs2JpsiPhiReader::insertCand" << endl;
  bmmReader::insertCand(pCand);
}


// ----------------------------------------------------------------------
void bmmBs2JpsiPhiReader::MCKinematics() {
  ((TH1D*)fpHistFile->Get("genStudy"))->Fill(1); 
  fGoodMCKinematics = true; 
  TGenCand *pC(0), *pB(0), *pPsi(0), *pPhi(0), *pM1(0), *pM2(0), *pK1(0), *pK2(0); 
  int npsi(0), nphi(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pC = fpEvt->getGenCand(i); 
    if (531 == TMath::Abs(pC->fID)) {
      pB = pC;
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

  // -- hard coded ?! FIXME
  if (pM1->fP.Perp() < 2.0) fGoodMCKinematics = false;  
  if (pM2->fP.Perp() < 2.0) fGoodMCKinematics = false;  
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

  // -- Get Kaons
  TAnaTrack *p0; 
  TAnaTrack *p1 = 0;
  TAnaTrack *p2 = 0; fpEvt->getSigTrack(fpCand->fSig1+1); 
  
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    if (TMath::Abs(p0->fMCID) != 321) continue;
    if (0 == p1) {
      p1 = p0; 
    } else {
      p2 = p0; 
    }
  }

  fKa1Pt  = p1->fPlab.Perp(); 
  fKa1Eta = p1->fPlab.Eta(); 
  fKa1Phi = p1->fPlab.Phi(); 

  fKa2Pt  = p2->fPlab.Perp(); 
  fKa2Eta = p2->fPlab.Eta(); 
  fKa2Phi = p2->fPlab.Phi(); 

  fDeltaR  = p1->fPlab.DeltaR(p2->fPlab); 

  TLorentzVector ka1, ka2; 
  ka1.SetPtEtaPhiM(fKa1Pt, fKa1Eta, fKa1Phi, MKAON); 
  ka2.SetPtEtaPhiM(fKa2Pt, fKa2Eta, fKa2Phi, MKAON); 

  TLorentzVector phiCand = ka1 + ka2; 
  fMKK     = phiCand.M();

  fGoodDeltaR = (fDeltaR < DELTAR);
  fGoodMKK    = (fMKK < MKKHI); 
  
  bmmReader::fillCandidateVariables();


}

// ----------------------------------------------------------------------
int bmmBs2JpsiPhiReader::tmCand(TAnaCand *pC) {
  
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
void bmmBs2JpsiPhiReader::fillCandidateHistograms() {
  //  cout << "bmmBs2JpsiPhiReader::fillCandidateHistograms()" << endl;
  //   fAnaCuts.update(); 
  //   fAnaCuts.dumpAll(); 

  fpDeltaR->fill(fDeltaR, fCandM); 
  fpMKK->fill(fMKK, fCandM); 

  bmmReader::fillCandidateHistograms(); 
}


// ---------------------------------------------------------------------- 
void bmmBs2JpsiPhiReader::bookHist() {
  bmmReader::bookHist();
  cout << "==> bmmBs2JpsiPhiReader: bookHist " << endl;

  fpDeltaR = bookDistribution("deltar", "#Delta R", "fGoodDeltaR", 50, 0., 1.);           
  fpMKK    = bookDistribution("mkk", "m(KK) [GeV]", "fGoodMKK", 50, 0.95, 1.15);           

  TH1D *h1; 
  h1 = new TH1D("h100", "", 100, 0., 100.);
  fTree->Branch("m2",      &fCandM,  "m2/D");
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

