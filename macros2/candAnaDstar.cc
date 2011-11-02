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
  
  TAnaCand *pCand(0), *pC(0), *pC1(0);
  TAnaTrack *pK, *pPi, *pPis; 
  //  cout << "Evt: " << fEvt << " ----------------------------------------------------------------------" << endl;
  // -- loop over all seq vtx fit candidates for D*
  double mdstar(0.), mdz(0.), dm(0.), fls3d(0.), flsxy(0.), prob(0.), chi2(0.), alpha(0.), pt(0.); 
  bool tm(false); 
  int ncand(0); 
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (300054 == pCand->fType) ++ncand; 
  }

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    //     if (54 == pCand->fType) {
    //       cout << "DUMP HFTruthCandidate with mass = " << pCand->fMass << endl;
    //       dumpHFTruthCand(pCand); 
    //     }

    if (300054 == pCand->fType) {
      tm = false; 
      tm = truthMatch(pCand); 
      //      if (tm) {
	// 	cout << " -> " << pCand->fType;
	// 	cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
	// 	cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
	// 	cout << " mass: " << pCand->fMass << endl;

	// 	cout << "DUMP HFDstarCandidate with mass = " << pCand->fMass << endl;
	// 	dumpHFDstarCand(pCand); 
	// 	fpEvt->dumpGenBlock(); 
      //      }
      mdstar = pCand->fMass; 

      // -- D0 
      if (pCand->fDau1 < 0) {
	if (pCand->fDau2 < 0) {
	  cout << "pCand->fDauX = -1!!! " << pCand->fType << endl;
	  pCand->dump();
	}
	continue;
      }
      pC = fpEvt->getCand(pCand->fDau1);
      pK = 0;
      pPi =0; 

      for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
	if (211 == fpEvt->getSigTrack(id)->fMCID) {
	  pPi = fpEvt->getRecTrack(fpEvt->getSigTrack(id)->fIndex);
	} else {
	  pK = fpEvt->getRecTrack(fpEvt->getSigTrack(id)->fIndex);
	}
      }
      pPis = fpEvt->getRecTrack(fpEvt->getSigTrack(pC->fSig1)->fIndex); 
	
      mdz = pC->fMass;
      dm = mdstar - mdz;

      TAnaVertex sv = pCand->fVtx;
      int pvidx = (pCand->fPvIdx > -1? pCand->fPvIdx : 0); 
      TVector3 svpv(sv.fPoint - fpEvt->getPV(pvidx)->fPoint); 
      fls3d = sv.fD3d/sv.fD3dE; 
      flsxy = sv.fDxy/sv.fDxyE; 
      prob  = sv.fProb;
      chi2  = sv.fChi2;
      pt    = pCand->fPlab.Perp(); 
      alpha = svpv.Angle(pCand->fPlab);


      if (prob < 0.05) continue;
      if (pt < 3) continue;
      if (alpha > 0.05) continue;

      if (tm) {
	((TH1D*)fHistDir->Get("mds"))->Fill(mdstar);
	((TH1D*)fHistDir->Get("mdz"))->Fill(mdz);
	((TH1D*)fHistDir->Get("dm"))->Fill(dm);
	((TH2D*)fHistDir->Get("h2"))->Fill(mdz, dm);
	((TH1D*)fHistDir->Get("fls3d"))->Fill(fls3d);
	((TH1D*)fHistDir->Get("flsxy"))->Fill(flsxy);
	((TH1D*)fHistDir->Get("prob"))->Fill(prob);
	((TH1D*)fHistDir->Get("chi2"))->Fill(chi2);
	((TH1D*)fHistDir->Get("alpha"))->Fill(alpha);
	((TH1D*)fHistDir->Get("pt"))->Fill(pt);
	
      } 

      ((TH1D*)fHistDir->Get("all_mds"))->Fill(mdstar);
      ((TH1D*)fHistDir->Get("all_mdz"))->Fill(mdz);
      ((TH2D*)fHistDir->Get("all_h2"))->Fill(mdz, dm);
      ((TH1D*)fHistDir->Get("all_dm"))->Fill(dm);
      ((TH1D*)fHistDir->Get("all_fls3d"))->Fill(fls3d);
      ((TH1D*)fHistDir->Get("all_flsxy"))->Fill(flsxy);
      ((TH1D*)fHistDir->Get("all_prob"))->Fill(prob);
      ((TH1D*)fHistDir->Get("all_chi2"))->Fill(chi2);
      ((TH1D*)fHistDir->Get("all_alpha"))->Fill(alpha);
      ((TH1D*)fHistDir->Get("all_pt"))->Fill(pt);

      if (1) {
	if (dm < 0.147 && dm > 0.145 && !tm) {
	  cout << "%%%%%%%%  ??" << endl;
	  tm = truthMatch(pCand, 1); 
	  cout << "  %%%" << endl;
	  dumpHFDstarCand(pCand); 
	  cout << "  %%%" << endl;
	  for (int i = 0; i < fpEvt->nCands(); ++i) {
	    pC1 = fpEvt->getCand(i);
	    if (300054 == pC1->fType) {
	      dumpHFDstarCand(pC1); 
	    }
	  }
	}
      }
    }

  }

  ((TH1D*)fHistDir->Get("all_ncand"))->Fill(ncand);

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
  cout << "HFDstarCand: idx = " << pC->fIndex << " type = " << pC->fType
       << " m* = " << pC->fMass << " m0 = " << pD->fMass << " dm = " << pC->fMass-pD->fMass << endl;
  //  if (fVerbose > -1) cout << "   mdz = " << pD->fMass << endl;
  for (int id = pD->fSig1; id <= pD->fSig2; ++id) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(id)->fIndex);
    cout << fpEvt->getSigTrack(id)->fMCID << " " ; 
    pT->dump(); 
  }
  // -- slow pion
  pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pC->fSig1)->fIndex); 
  cout << fpEvt->getSigTrack(pC->fSig1)->fMCID << " " ; 
  pT->dump(); 

}


// ----------------------------------------------------------------------
// -- works ONLY for this dedicated decay mode with the seq vtx fitter!
bool candAnaDstar::truthMatch(TAnaCand *pCand, int verbose) {

  // -- check slow pion
  TAnaTrack *pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex);
  if (pT->fGenIndex < 0) {
    if (verbose > 0) cout << "pT->fGenIndex < 0" << endl;
    return false; 
  }
  TGenCand  *pG = fpEvt->getGenCand(fpEvt->getRecTrack(pT->fIndex)->fGenIndex); 
  if (0 == pG) {
    if (verbose > 0) cout << "0 == pG" << endl;
    return false;
  }
  if (211 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "211 != TMath::Abs(pG->fID)" << endl;
    return false;
  }
  int moSlowPion = pG->fMom1;
  pG = fpEvt->getGenCand(pG->fMom1); 
  if ((0 == pG) || 413 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "(0 == pG) || 413 != pG->fID, pG->fID = " << pG->fID << " moSlowPion = " << moSlowPion << endl;
    return false;
  }
  if (pG->fDau2 - pG->fDau1 > 1) {
    if (verbose > 0) cout << "pG->fDau2 - pG->fDau1 > 1" << endl;
    return false; 
  }
  if (verbose > 0) cout << "slow pion " << pT->fIndex << " OK, fGenIndex = " << pT->fGenIndex << endl;


  // -- check D0 
  if (pCand->fDau1 < 0) {
    if (verbose > 0) cout << "no pCand->fDau1" << endl;
    return false;
  }
  TAnaCand *pC = fpEvt->getCand(pCand->fDau1);
  
  // -- check D0 daughters
  int type(0), moIdx(-1); 
  for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
    pT = fpEvt->getSigTrack(id);
    type = pT->fMCID;
    pT = fpEvt->getRecTrack(pT->fIndex);
    if (pT->fGenIndex < 0) {
      if (verbose > 0) cout << "no pT->fGenIndex" << endl;
      return false;
    }
    pG = fpEvt->getGenCand(pT->fGenIndex);
    if (moIdx < 0) {
      moIdx = pG->fMom1;
    } else {
      if (moIdx != pG->fMom1) {
	if (verbose > 0) cout << "moIdx != pG->fMom1" << endl;
	return false;
      }
    }
    if (verbose > 0) 
      cout << "dau cand sigtrack " << id 
 	   << " with type = " << type 
 	   << " and gen ID = " << pG->fID 
 	   << " at gen idx = " << pT->fGenIndex 
 	   << endl;
    if (TMath::Abs(type) != TMath::Abs(pG->fID)) {
      if (verbose > 0) cout << "TMath::Abs(type) != TMath::Abs(pG->fID), type = " << type << " pG->fID = " << pG->fID 
			    << " track " << pT->fIndex << endl;
      return false;
    }
  }

  // -- Get gen-level D0
  if (pG->fMom1 < 0) {
    if (verbose > 0) cout << "pG->fMom1 < 0" << endl;
    return false; 
  }
  pG = fpEvt->getGenCand(pG->fMom1);
  if (pG->fMom1 != moSlowPion) {
    if (verbose > 0) cout << "pG->fMom1 != moSlowPion" << endl;
    return false; 
  }
  if (pG->fDau2 - pG->fDau1 > 1) {
    if (verbose > 0) cout << "pG->fDau2 - pG->fDau1 > 1" << endl;
    return false; 
  }

  if (verbose > 0)   cout << "===> truth matching OK" << endl;
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

  TH1 *h = new TH1D("mds", "m(dstar)", 60, 1.9, 2.2);
  h = new TH1D("mdz", "m(d0)", 60, 1.8, 1.92);
  h = new TH1D("dm", "delta(m)", 60, 0.13, 0.16);
  h = new TH1D("ncand", "ncand", 200, 0., 200);
  h = new TH1D("fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("chi2", "chi2", 100, 0., 10.);
  h = new TH1D("alpha", "alpha", 50, 0., 0.2);
  h = new TH1D("pt", "pT", 50, 0., 25);
  h = new TH1D("ptK",   "pT", 50, 0., 10);
  h = new TH1D("ptPi",  "pT", 50, 0., 10);
  h = new TH1D("ptPis", "pT", 50, 0., 10);
  TH2D *h2 = new TH2D("h2", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);

  h = new TH1D("all_mds", "m(dstar)", 60, 1.9, 2.2);
  h = new TH1D("all_mdz", "m(d0)", 60, 1.8, 1.92);
  h = new TH1D("all_dm", "delta(m)", 60, 0.13, 0.16);
  h = new TH1D("all_ncand", "ncand", 200, 0., 200);
  h = new TH1D("all_fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("all_flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("all_prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("all_chi2", "chi2", 100, 0., 10.);
  h = new TH1D("all_alpha", "alpha", 50, 0., 0.2);
  h = new TH1D("all_pt", "pT", 50, 0., 25);
  h2 = new TH2D("all_h2", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);

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
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

   
  }

}
