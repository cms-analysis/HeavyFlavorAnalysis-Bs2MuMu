#include "muCharmReader.hh"
#include "../interface/HFMasses.hh"

#include "TRandom.h"
#include <cmath>

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/muCharmReader.default.cuts
//           ./runTreeReaders -f test.root 
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
muCharmReader::muCharmReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> muCharmReader: constructor..." << endl;
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

  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 

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

  if (0 != fpCand) {
    fillHist(); 
  }

}


// ----------------------------------------------------------------------
void muCharmReader::initVariables() {

  
  fGoodMCKinematics = fGoodL1 = fGoodHLT = fGoodEvent = false; 
  fGoodMuonsID = fGoodMuonsPT = false; 
  fGoodTracks = fGoodTracksPT = false; 
  fGoodCandPT = false; 

}


// ----------------------------------------------------------------------
void muCharmReader::candidateSelection(int mode) {

  fCandPt = fCandMass = -1.; 
  fpCand = 0; 
  
  TAnaCand *pCand(0);
  vector<int> lCands;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (50 == pCand->fType) {
      fillTMCand(pCand, 50);
    }
    if (TYPE != pCand->fType) continue;

    lCands.push_back(iC); 
  }

  int nc(lCands.size());

  ((TH1D*)fpHistFile->Get("h2"))->Fill(nc); 
  if (0 == nc) return; 

  int best(0); 
  if (nc > 1) {
    double chi2Min(99.), chi2(0.); 
    for (unsigned int iC = 0; iC < lCands.size(); ++iC) {
      pCand = fpEvt->getCand(lCands[iC]); 

      chi2 = pCand->fVtx.fChi2; 

      if (0 == mode) {
	if (chi2 < chi2Min) {
	  chi2Min = chi2; 
	  best = lCands[iC]; 
	}
      }
    }
  }


  if (best > -1) {
    fpCand = fpEvt->getCand(best); 
    fCandPt   = pCand->fPlab.Perp();
    
    TAnaTrack *pc1 = fpEvt->getSigTrack(pCand->fSig1); 
    TAnaTrack *pc2 = fpEvt->getSigTrack(pCand->fSig1+1); 
    
    TLorentzVector a1, a2, a0; 
    a1.SetVectM(pc1->fPlab, MKAON); 
    a2.SetVectM(pc2->fPlab, MPION); 
    a0 = a1 + a2; 

    fCandMass = a0.M();
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
      ps[i]->fInt2 = 1; 
      
      //      if (pt->fTrackQuality < 0)     ps[i]->fInt1 = 0; 
      if (pt->fPlab.Perp() < MUPTLO) ps[i]->fInt2 = 0; 

    }

  }

}


// ----------------------------------------------------------------------
void muCharmReader::muonSelection() {
 
}


// ----------------------------------------------------------------------
void muCharmReader::fillHist() {

  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 

  if (0 != fpCand) {
    ((TH1D*)fpHistFile->Get("h10"))->Fill(fCandPt); 
    ((TH1D*)fpHistFile->Get("h11"))->Fill(fCandMass); 
    ((TH1D*)fpHistFile->Get("h12"))->Fill(fpCand->fVtx.fChi2); 

    fTree->Fill(); 
  }

}

// ---------------------------------------------------------------------- 
void muCharmReader::bookHist() {
  cout << "==> muCharmReader: bookHist " << endl;
  
  TH1D *h; 
  h = new TH1D("h1", "Ntrk", 200, 0., 200.);
  h = new TH1D("h2", "NCand", 20, 0., 20.);
  h = new TH1D("h10", "pT", 40, 0., 20.);
  h = new TH1D("h11", "mass", 50, 1.6, 2.1);
  h = new TH1D("h12", "chi2", 50, 0., 10.);

  h = new TH1D("h1050", "TM D0->Kpi", 50, 1.6, 2.1);
  h = new TH1D("h2050", "TM mu D0->Kpi", 50, 1.6, 2.1);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun, "run/I");
  fTree->Branch("pt",     &fCandPt, "pt/D");
  fTree->Branch("m",      &fCandMass, "pt/D");

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

    if (!strcmp(CutName, "CHARMPTLO")) {
      CHARMPTLO = CutValue; ok = 1;
      if (dump) cout << "CHARMPTLO:           " << CHARMPTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, CHARMPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(B_{s}) [GeV]");
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
