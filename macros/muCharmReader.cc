#include "muCharmReader.hh"
#include "../interface/HFMasses.hh"

#include "TRandom.h"
#include <cmath>

using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/muCharmReader.default.cuts
//           ./runTreeReaders -f test.root 
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
muCharmReader::muCharmReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> muCharmReader: constructor..." << endl;
  fpJSON = new JSON("/shome/ursl/json/1.json"); 
}

// ----------------------------------------------------------------------
muCharmReader::~muCharmReader() {
  cout << "==> muCharmReader: destructor..." << endl;

}

// ----------------------------------------------------------------------
void muCharmReader::startAnalysis() {
  cout << "==> muCharmReader: Starting analysis..." << endl;
}


// ----------------------------------------------------------------------
void muCharmReader::eventProcessing() {

  cout << fRun << " " << fLS << ": " << fpJSON->good(fRun, fLS) << endl;

  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 
  ((TH1D*)fpHistFile->Get("h2"))->Fill(fpEvt->nCands()); 

  //   for (int i = 0; i < 100; ++i) {
  //     cout << fpEvt->fHLTNames[i] << endl;
  //   }
  //   return;


  TAnaCand *pCand; 
  TH1D *h; 
  int n1313(0), n1300(0); 
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    ((TH1D*)fpHistFile->Get("h3"))->Fill(pCand->fType); 
    if (pCand->fType < 100) {
      //      cout << Form("%6i %i", iC, pCand->fType) << endl;
    }

    if (h = (TH1D*)fpHistFile->Get(Form("m%d", pCand->fType))) {
      if (pCand->fType > 20000 && pCand->fType < 21000) {
	h->Fill(pCand->fVar1); 
      } else {
	h->Fill(pCand->fMass); 
      }
    } else {
      //      cout << "Unknown candidate " << pCand->fType << endl;
    }

    if (100521 == pCand->fType) doBplus(pCand); 
    if (20010 == pCand->fType) doDzero(pCand); 
    if (1300 == pCand->fType) {
      ++n1300;
    }
    if (1313 == pCand->fType) {
      ++n1313;
      doJpsi(pCand); 
      doUpsilon(pCand); 
    }
  }
  
  ((TH1D*)fpHistFile->Get("h100"))->Fill(n1313); 
  ((TH1D*)fpHistFile->Get("h101"))->Fill(n1300); 

  return;

  // -- initialize all variables
  initVariables(); 

  if (!goodRun()) {
    //cout << "not a good run: " << fRun << endl;
    return; 
  } 

  // -- track selection for all candidates
  trackSelection(); 

  // -- Select a candidate
  candidateSelection(0); 

//   if (0 != fpCand) {
//     fillHist(); 
//   }

}


// ----------------------------------------------------------------------
void muCharmReader::initVariables() {

  

}



// ----------------------------------------------------------------------
void muCharmReader::doBplus(TAnaCand *pCand ) {
  
  ((TH1D*)fpHistFile->Get("m100521h0"))->Fill(pCand->fMass);

  if (pCand->fVtx.fChi2 > 5.) return;
  ((TH1D*)fpHistFile->Get("m100521h1"))->Fill(pCand->fMass);

  ((TH1D*)fpHistFile->Get("m100521h10"))->Fill(pCand->fPlab.Perp());
  
  if (pCand->fPlab.Perp() < 3.) return;
  ((TH1D*)fpHistFile->Get("m100521h2"))->Fill(pCand->fMass);

  if (pCand->fVtx.fDxy/pCand->fVtx.fDxyE < 1) return; 
  ((TH1D*)fpHistFile->Get("m100521h3"))->Fill(pCand->fMass);


}


// ----------------------------------------------------------------------
void muCharmReader::doDzero(TAnaCand *pCand ) {

  TAnaTrack *pt1, *pt2, *pt3; 
  //  cout << pCand->fSig1 << " .. " << pCand->fSig2 << endl;
  pt1 = fpEvt->getSigTrack(pCand->fSig1); 
  pt2 = fpEvt->getSigTrack(pCand->fSig1+1); 
  pt3 = fpEvt->getSigTrack(pCand->fSig2); 
  
  ((TH1D*)fpHistFile->Get("m20010h0"))->Fill(pCand->fVar1);
  if (pCand->fMaxDoca > 0.01) return;
  ((TH1D*)fpHistFile->Get("m20010h1"))->Fill(pCand->fVar1);

  if (pCand->fPlab.Perp() < 4) return;
  ((TH1D*)fpHistFile->Get("m20010h2"))->Fill(pCand->fVar1);
}


// ----------------------------------------------------------------------
void muCharmReader::doJpsi(TAnaCand *pCand ) {

  TAnaTrack *pt1, *pt2; 
  //  cout << pCand->fSig1 << " .. " << pCand->fSig2 << endl;
  pt1 = fpEvt->getSigTrack(pCand->fSig1); 
  pt2 = fpEvt->getSigTrack(pCand->fSig2); 

  fMass    = pCand->fMass;
  fDocaMax = pCand->fMaxDoca;
  fPt      = pCand->fPlab.Perp();
  fChi2    = pCand->fVtx.fChi2;
  
  fm1pt    = pt1->fPlab.Perp(); 
  fm2pt    = pt2->fPlab.Perp(); 

  if (pt1->fMuID == -1) fm1m = 0; else  fm1m     = pt1->fMuID;
  if (pt2->fMuID == -1) fm2m = 0; else  fm2m     = pt2->fMuID;

  fTree->Fill();
  
  ((TH1D*)fpHistFile->Get("m1300h0"))->Fill(pCand->fMass);

  if (fDocaMax > 0.01) return;
  ((TH1D*)fpHistFile->Get("m1300h1"))->Fill(pCand->fMass);

  if ((pt2->fMuID & 4) == 0) return;
  ((TH1D*)fpHistFile->Get("m1300h2"))->Fill(pCand->fMass);
  

}


// ----------------------------------------------------------------------
void muCharmReader::doUpsilon(TAnaCand *pCand ) {

  TAnaTrack *pt1, *pt2; 
  //  cout << pCand->fSig1 << " .. " << pCand->fSig2 << endl;
  pt1 = fpEvt->getSigTrack(pCand->fSig1); 
  pt2 = fpEvt->getSigTrack(pCand->fSig2); 

  fMass    = pCand->fMass;
  fDocaMax = pCand->fMaxDoca;
  fPt      = pCand->fPlab.Perp();
  fChi2    = pCand->fVtx.fChi2;
  
  fm1pt    = pt1->fPlab.Perp(); 
  fm2pt    = pt2->fPlab.Perp(); 

  if (pt1->fMuID == -1) fm1m = 0; else  fm1m     = pt1->fMuID;
  if (pt2->fMuID == -1) fm2m = 0; else  fm2m     = pt2->fMuID;

  fTree->Fill();
  
  ((TH1D*)fpHistFile->Get("ups1313h0"))->Fill(pCand->fMass);

  if (fm1m & 6 == 6) {
    ((TH1D*)fpHistFile->Get("ups1313h1"))->Fill(pCand->fMass);
    if (fm2m & 6 == 6) {
      ((TH1D*)fpHistFile->Get("ups1313h2"))->Fill(pCand->fMass);
    }
  }; 

}

// ----------------------------------------------------------------------
void muCharmReader::candidateSelection(int mode) {

  TAnaCand *pCand(0);
  vector<int> lCands;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (MCTYPE == pCand->fType) {
      fillTMCand(pCand, 50);
    }

    if (TYPE != pCand->fType) continue;

    // -- check that candidate fulfilled track selection
    TAnaTrack *pc1 = fpEvt->getSigTrack(pCand->fSig1); 
    TAnaTrack *pc2 = fpEvt->getSigTrack(pCand->fSig1+1); 
    if (0 == pc1->fInt1) continue;
    if (0 == pc2->fInt1) continue;

    if (pCand->fVtx.fChi2 > VTXCHI2) continue; 
    if (pCand->fPlab.Perp() < CHARMPTLO) continue; 
    
    lCands.push_back(iC); 
  }

  int nc(lCands.size());

  ((TH1D*)fpHistFile->Get("h2"))->Fill(nc); 
  if (0 == nc) return; 

  int best(0); 
  if (nc > 1) {
    double ptMax(-1), pt(0.); 
    for (unsigned int iC = 0; iC < lCands.size(); ++iC) {
      pCand = fpEvt->getCand(lCands[iC]); 

      pt = pCand->fPlab.Perp(); 

      if (0 == mode) {
	if (pt > ptMax) {
	  ptMax = pt; 
	  best = lCands[iC]; 
	}
      }
    }
  }

  if (best > -1) {
    fpCand1 = fpEvt->getCand(best); 

  //     TAnaTrack *pc1 = fpEvt->getSigTrack(pCand->fSig1); 
  //     TAnaTrack *pc2 = fpEvt->getSigTrack(pCand->fSig1+1); 
  
  //     TLorentzVector a1, a2, a0; 
  //     a1.SetVectM(pc1->fPlab, MKAON); 
  //     a2.SetVectM(pc2->fPlab, MPION); 
  //     a0 = a1 + a2; 
  
  //     fCandMass = a0.M();

  }
}


// ----------------------------------------------------------------------
void muCharmReader::fillTMCand(TAnaCand *pCand, int type) {
  
  // -- fill simple TM mass
  ((TH1D*)fpHistFile->Get(Form("h%i", 1000+type)))->Fill(pCand->fMass); 

   // -- now try to find an additional muon
   TAnaTrack *pc1 = fpEvt->getSigTrack(pCand->fSig1); 
   TAnaTrack *pc2 = fpEvt->getSigTrack(pCand->fSig1+1); 

   int pgI = pc1->fGenIndex;
   TGenCand *pB;
   int foundB(0); 
   if (pgI < fpEvt->nGenCands()) {
     TGenCand *pg = fpEvt->getGenCand(pgI); 
     int momI = pg->fMom1; 
     TGenCand *pD0 = fpEvt->getGenCand(momI); 

     pB = pD0;
     int id(999), cnt(0); 
     while (id > 100) {
       momI = pB->fMom1;
       pB = fpEvt->getGenCand(momI); 
       id  = TMath::Abs(pB->fID%1000);
       ++cnt;
       if (id > 499 && id < 599) {
	 foundB = 1; 
	 break;
       }
     }

     if (1 == foundB) {
       cout << "========> Found a B" << endl;
       pB->dump();
       TGenCand *pG; 
       for (int ig = pB->fDau1; ig <= pB->fDau2; ++ig) {
	 pG = fpEvt->getGenCand(ig); 
	 pG->dump(); 
	 if (13 == TMath::Abs(pG->fID)) {
	   cout << "++++++++++++++++++++++++++++ with a MUON" << endl;
	 }
       }
     }
      
   }

}


// ----------------------------------------------------------------------
void muCharmReader::MCKinematics() {

}

// ----------------------------------------------------------------------
void muCharmReader::L1TSelection() {

}

// ----------------------------------------------------------------------
void muCharmReader::HLTSelection() {

}

// ----------------------------------------------------------------------
void muCharmReader::trackSelection() {

  TAnaCand *pCand;
  TAnaTrack *pt, *ps[2]; 

  TLorentzVector pb, pm1, pm2; 

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (TYPE != pCand->fType) continue;
    
    // -- Get the 2 signal tracks
    ps[0] = fpEvt->getSigTrack(pCand->fSig1); 
    ps[1] = fpEvt->getSigTrack(pCand->fSig2); 

    for (int i = 0; i < 2; ++i) {

      // -- Get the corresponding RecTrack
      pt = fpEvt->getRecTrack(ps[i]->fIndex); 

      // -- Use spare variables in the RecTrack as bookmark whether it passed the signal selection
      ps[i]->fInt1 = 1; 
      
      if (0 == i && pt->fPlab.Perp() < KAPTLO) ps[i]->fInt1 = 0; 
      if (1 == i && pt->fPlab.Perp() < PIPTLO) ps[i]->fInt1 = 0; 
    }

  }

}


// ----------------------------------------------------------------------
void muCharmReader::muonSelection() {
 
}


// ----------------------------------------------------------------------
void muCharmReader::fillHist() {

  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 

  if (0 != fpCand1) {
    ((TH1D*)fpHistFile->Get("h10"))->Fill(fpCand1->fPlab.Perp()); 
    ((TH1D*)fpHistFile->Get("h11"))->Fill(fpCand1->fMass); 
    ((TH1D*)fpHistFile->Get("h12"))->Fill(fpCand1->fVtx.fChi2); 

    fTree->Fill(); 
  }

}

// ---------------------------------------------------------------------- 
void muCharmReader::bookHist() {
  cout << "==> muCharmReader: bookHist " << endl;
  
  TH1D *h; 
  h = new TH1D("h1", "Ntrk", 500, 0., 1000.);
  h = new TH1D("h2", "NCand", 20, 0., 20.);
  h = new TH1D("h3", "cand ID", 1000100, -100., 1000000.);
  h = new TH1D("h100", "1313 multiplicity", 20, 0., 20.);
  h = new TH1D("h101", "1300 multiplicity", 100, 0., 100.);
  h = new TH1D("h10", "pT", 40, 0., 20.);
  h = new TH1D("h11", "mass", 50, 1.6, 2.1);
  h = new TH1D("h12", "chi2", 50, 0., 10.);

  h = new TH1D("h1050", "TM D0->Kpi", 50, 1.6, 2.1);
  h = new TH1D("h2050", "TM mu D0->Kpi", 50, 1.6, 2.1);

  h = new TH1D("m1300", "mass 1300", 40, 2.7, 3.5);
  h = new TH1D("m1313", "mass 1313", 40, 2.7, 3.5);
  h = new TH1D("m20010", "mass 20010", 70, 1.7, 2.4);
  h = new TH1D("m20020", "mass 20020", 70, 1.7, 2.4);
  h = new TH1D("m20030", "mass 20030", 70, 1.7, 2.4);
  h = new TH1D("m20040", "mass 20040", 70, 1.7, 2.4);
  h = new TH1D("m20050", "mass 20050", 70, 1.7, 2.4);
  h = new TH1D("m20060", "mass 20060", 70, 1.7, 2.4);

  h = new TH1D("ups1313h0", "mass 1313", 40, 9.0, 11.0);
  h = new TH1D("ups1313h1", "mass 1313", 40, 9.0, 11.0);
  h = new TH1D("ups1313h2", "mass 1313", 40, 9.0, 11.0);


  h = new TH1D("m443", "mass 443", 40, 2.7, 3.5);

  h = new TH1D("m100521", "mass 100521", 60, 5.0, 5.6);
  h = new TH1D("m100511", "mass 100511", 60, 5.0, 5.6);
  h = new TH1D("m100531", "mass 100531", 60, 5.0, 5.6);

  h = new TH1D("m100521h0",  "mass 100521", 50, 5.0, 5.5);
  h = new TH1D("m100521h1",  "mass 100521", 60, 5.0, 5.6);
  h = new TH1D("m100521h2",  "mass 100521", 60, 5.0, 5.6);
  h = new TH1D("m100521h3",  "mass 100521", 60, 5.0, 5.6);

  h = new TH1D("m100521h10",  "pT 100521", 40, 0., 10.);

  h = new TH1D("m20010h0",  "mass 100521", 60, 1.7, 2.0);
  h = new TH1D("m20010h1",  "mass 100521", 60, 1.7, 2.0);
  h = new TH1D("m20010h2",  "mass 100521", 60, 1.7, 2.0);
  h = new TH1D("m20010h3",  "mass 100521", 60, 1.7, 2.0);

  h = new TH1D("m1300h0", "mass 1300", 40, 2.7, 3.5);
  h = new TH1D("m1300h1", "mass 1300", 40, 2.7, 3.5);
  h = new TH1D("m1300h2", "mass 1300", 40, 2.7, 3.5);
  h = new TH1D("m1300h3", "mass 1300", 40, 2.7, 3.5);
  h = new TH1D("m1300h4", "mass 1300", 40, 2.7, 3.5);
  h = new TH1D("m1300h5", "mass 1300", 40, 2.7, 3.5);


  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",     &fRun,     "run/I");
  fTree->Branch("m",       &fMass,    "m/D");
  fTree->Branch("pt",      &fPt,      "pt/D");
  fTree->Branch("docamax", &fDocaMax, "docamax/D");
  fTree->Branch("chi2",    &fChi2,    "chi2/D");
  fTree->Branch("m1m",     &fm1m,     "m1m/I");
  fTree->Branch("m1pt",    &fm1pt,    "m1pt/D");
  fTree->Branch("m2m",     &fm2m,     "m2m/I");
  fTree->Branch("m2pt",    &fm2pt,    "m2pt/D");

}

// ----------------------------------------------------------------------
void muCharmReader::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> muCharmReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> muCharmReader: Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  int ibin; 
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      if (dump) cout << "TYPE:           " << TYPE << endl;
    }

    if (!strcmp(CutName, "MCTYPE")) {
      MCTYPE = int(CutValue); ok = 1;
      if (dump) cout << "MCTYPE:           " << MCTYPE << endl;
    }

    if (!strcmp(CutName, "CHARMPTLO")) {
      CHARMPTLO = CutValue; ok = 1;
      if (dump) cout << "CHARMPTLO:           " << CHARMPTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, CHARMPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(B_{s}) [GeV]");
    }

    if (!strcmp(CutName, "VTXCHI2")) {
      VTXCHI2 = CutValue; ok = 1;
      if (dump) cout << "VTXCHI2:           " << VTXCHI2 << " " << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, VTXCHI2);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#chi^{2}");
    }

    if (!strcmp(CutName, "CHARMETALO")) {
      CHARMETALO = CutValue; ok = 1;
      if (dump) cout << "CHARMETALO:           " << CHARMETALO << endl;
      ibin = 12;
      hcuts->SetBinContent(ibin, CHARMETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{min}(B_{s})");
    }

    if (!strcmp(CutName, "CHARMETAHI")) {
      CHARMETAHI = CutValue; ok = 1;
      if (dump) cout << "CHARMETAHI:           " << CHARMETAHI << endl;
      ibin = 13;
      hcuts->SetBinContent(ibin, CHARMETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{max}(B_{s})");
    }

    if (!strcmp(CutName, "KAPTLO")) {
      KAPTLO = CutValue; ok = 1;
      if (dump) cout << "KAPTLO:           " << KAPTLO << " GeV" << endl;
      ibin = 21;
      hcuts->SetBinContent(ibin, KAPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(K)");
    }

    if (!strcmp(CutName, "PIPTLO")) {
      PIPTLO = CutValue; ok = 1;
      if (dump) cout << "PIPTLO:           " << PIPTLO << " GeV" << endl;
      ibin = 21;
      hcuts->SetBinContent(ibin, PIPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(#pi)");
    }
    if (!strcmp(CutName, "MUPTLO")) {
      MUPTLO = CutValue; ok = 1;
      if (dump) cout << "MUPTLO:           " << MUPTLO << " GeV" << endl;
      ibin = 21;
      hcuts->SetBinContent(ibin, MUPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(#mu)");
    }

    if (!strcmp(CutName, "MUPTHI")) {
      MUPTHI = CutValue; ok = 1;
      if (dump) cout << "MUPTHI:           " << MUPTHI << " GeV" << endl;
      ibin = 22;
      hcuts->SetBinContent(ibin, MUPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{max}(#mu)");
    }

    if (!strcmp(CutName, "MUETALO")) {
      MUETALO = CutValue; ok = 1;
      if (dump) cout << "MUETALO:           " << MUETALO << endl;
      ibin = 23;
      hcuts->SetBinContent(ibin, MUETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{min}(#mu)");
    }

    if (!strcmp(CutName, "MUETAHI")) {
      MUETAHI = CutValue; ok = 1;
      if (dump) cout << "MUETAHI:           " << MUETAHI << endl;
      ibin = 24;
      hcuts->SetBinContent(ibin, MUETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{max}(#mu)");
    }


    if (!ok) cout << "==> muCharmReader: ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}
