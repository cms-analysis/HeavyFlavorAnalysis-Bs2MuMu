#include "bmmReader.hh"
#include "TRandom.h"
#include <cmath>

#include "../interface/HFMasses.hh"

using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/bmmReader.default.cuts
//           ./runTreeReaders -f test.root 
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
bmmReader::bmmReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> bmmReader: constructor..." << endl;
}

// ----------------------------------------------------------------------
bmmReader::~bmmReader() {
  cout << "==> bmmReader: destructor..." << endl;

}

// ----------------------------------------------------------------------
void bmmReader::startAnalysis() {
  cout << "==> bmmReader: Starting analysis..." << endl;
}


// ----------------------------------------------------------------------
void bmmReader::eventProcessing() {

  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 

  // -- initialize all variables
  initVariables(); 

  // -- track selection for all candidates
  pvStudy(); 

  // -- track selection for all candidates
  trackSelection(); 

  // -- Select a candidate
  candidateSelection(0); 

  if (0 != fpCand) {
    fillHist(); 
  }

}


// ----------------------------------------------------------------------
void bmmReader::initVariables() {
  
  fGoodMCKinematics = fGoodL1 = fGoodHLT = fGoodEvent = false; 
  fGoodMuonsID = fGoodMuonsPT = false; 
  fGoodTracks = fGoodTracksPT = false; 
  fGoodCandPT = false; 

}


// ----------------------------------------------------------------------
void bmmReader::pvStudy() {

  TAnaVertex *pV; 
  cout << "Found " << fpEvt->nPV() << " PVs, best PV at "  << fpEvt->fBestPV
       << endl;
  for (int i = 0; i < fpEvt->nPV(); ++i) {
    pV = fpEvt->getPV(i); 
    cout << "  " << i << " -> with " << pV->getNtracks() << " tracks: ";
    for (int p = 0; p < pV->getNtracks(); ++p) {
      cout << pV->getTrack(p) << " "; 
    }
    cout << endl;
  }

}


// ----------------------------------------------------------------------
void bmmReader::candidateSelection(int mode) {

  fCandPt = fCandMass = -1.; 
  fpCand = 0; 
  
  TAnaCand *pCand(0);
  vector<int> lCands;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (TYPE != pCand->fType) continue;
    lCands.push_back(iC); 
  }

  int nc(lCands.size());

  ((TH1D*)fpHistFile->Get("h2"))->Fill(nc); 
  if (0 == nc) return; 

  int best(0); 
  if (nc > 1) {
    double drMin(99.), dr(0.); 
    for (unsigned int iC = 0; iC < lCands.size(); ++iC) {
      pCand = fpEvt->getCand(lCands[iC]); 

      TAnaTrack *pl1 = fpEvt->getSigTrack(pCand->fSig1); 
      TAnaTrack *pl2 = fpEvt->getSigTrack(pCand->fSig2); 


      // -- FIXME: SHOULD THESE CUTS BE APPLIED??
      // -- Check that candidate tracks passed track selection
      if (0 == pl1->fInt1 || 0 == pl1->fInt2) continue; 
      if (0 == pl2->fInt1 || 0 == pl2->fInt2) continue; 

      if (0 == mode) {
	dr = pl1->fPlab.DeltaR(pl2->fPlab); 
	if (dr < drMin) {
	  best = lCands[iC]; 
	  drMin = dr; 
	}
      }
    }
  }


  if (best > -1) {
    fpCand = fpEvt->getCand(best); 
    fCandPt   = pCand->fPlab.Perp();
    fCandMass = pCand->fMass;
  }


}

// ----------------------------------------------------------------------
void bmmReader::MCKinematics() {

}

// ----------------------------------------------------------------------
void bmmReader::L1TSelection() {

}

// ----------------------------------------------------------------------
void bmmReader::HLTSelection() {

}

// ----------------------------------------------------------------------
void bmmReader::trackSelection() {

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
void bmmReader::muonSelection() {
 
}


// ----------------------------------------------------------------------
void bmmReader::fillHist() {

  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 

  if (0 != fpCand) {
    ((TH1D*)fpHistFile->Get("h10"))->Fill(fCandPt); 
    ((TH1D*)fpHistFile->Get("h11"))->Fill(fCandMass); 

    fTree->Fill(); 
  }

}

// ---------------------------------------------------------------------- 
void bmmReader::bookHist() {
  cout << "==> bmmReader: bookHist " << endl;
  
  TH1D *h; 
  h = new TH1D("h1", "Ntrk", 200, 0., 200.);
  h = new TH1D("h2", "NCand", 20, 0., 20.);
  h = new TH1D("h10", "pT", 40, 0., 20.);
  h = new TH1D("h11", "mass", 60, 4.5, 6.0);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun, "run/I");
  fTree->Branch("pt",     &fCandPt, "pt/D");
  fTree->Branch("m",      &fCandMass, "pt/D");

}

// ----------------------------------------------------------------------
void bmmReader::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> bmmReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> bmmReader: Cut file  " << fCutFile.Data() << endl;
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

    if (!strcmp(CutName, "BSPTLO")) {
      BSPTLO = CutValue; ok = 1;
      if (dump) cout << "BSPTLO:           " << BSPTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, BSPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(B_{s}) [GeV]");
    }

    if (!strcmp(CutName, "BSETALO")) {
      BSETALO = CutValue; ok = 1;
      if (dump) cout << "BSETALO:           " << BSETALO << endl;
      ibin = 12;
      hcuts->SetBinContent(ibin, BSETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{min}(B_{s})");
    }

    if (!strcmp(CutName, "BSETAHI")) {
      BSETAHI = CutValue; ok = 1;
      if (dump) cout << "BSETAHI:           " << BSETAHI << endl;
      ibin = 13;
      hcuts->SetBinContent(ibin, BSETAHI);
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


    if (!ok) cout << "==> bmmReader: ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}
