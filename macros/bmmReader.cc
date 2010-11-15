#include "bmmReader.hh"
#include "TRandom.h"
#include <cmath>
#include <string>

#include "../interface/HFMasses.hh"

using std::string;
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

  if (fVerbose > 4) cout << "event: " << fEvent << endl;

  // -- initialize all variables
  initVariables(); 

  L1TSelection();
  HLTSelection();
  trackSelection(); 
  muonSelection();
  candidateSelection(0); 
  
  fillHist();

}


// ----------------------------------------------------------------------
void bmmReader::initVariables() {

  // -- Fill candidates vector
  fCands.clear();
  fGoodMCKinematics = true; 
  fGoodL1 = fGoodHLT = true; 

  fGoodMuonsID.clear();
  fGoodMuonsPT.clear();
  fGoodTracks.clear();
  fGoodTracksPT.clear();

  fGoodCandPT.clear();

  fGoodEvent = true;
  
  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (TYPE != pCand->fType) continue;
    insertCand(pCand);
  }


  // -- clear reduced tree variables
  fCandPt = fCandM = -1.;

}


// ----------------------------------------------------------------------
void bmmReader::insertCand(TAnaCand* pCand) {
    fCands.push_back(pCand); 
    fGoodMuonsID.push_back(true); 
    fGoodMuonsPT.push_back(true);
    fGoodTracks.push_back(true);
    fGoodTracksPT.push_back(true);
    fGoodCandPT.push_back(true);
}


// ----------------------------------------------------------------------
void bmmReader::pvStudy() {

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
  
  TAnaTrack *pt, *ps; 
  
  for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
    for (int it = fCands[iC]->fSig1; it <= fCands[iC]->fSig2; ++it) {
      ps = fpEvt->getSigTrack(it); 
      ps->fInt1 = 1; 
      ps->fInt2 = 1; 
      pt = fpEvt->getRecTrack(ps->fIndex); 

      // FIXME: replace (pt->fTrackQuality < 0) with something meaningful!
      if (TRACKQUALITY > 0 && (pt->fTrackQuality < 0)) {
	if (fVerbose > 1) cout << "track " << ps->fIndex << " failed track quality: " << pt->fTrackQuality << endl;
	fGoodTracks[iC] = false; 
      }

      if (TMath::Abs(pt->fTip) > TRACKTIP) {
	if (fVerbose > 1) cout << "track " << ps->fIndex << " failed tip: " << pt->fTip << endl;
	fGoodTracks[iC] = false; 
      }
      
      if (TMath::Abs(pt->fLip) > TRACKLIP) { 
	if (fVerbose > 1) cout << "track " << ps->fIndex << " failed lip: " << pt->fLip << endl;          
	fGoodTracks[iC] = false; 
      }

      if (pt->fPlab.Perp() < TRACKPTLO) {
	if (fVerbose > 1) cout << "track " << ps->fIndex << " failed pt: " << pt->fPlab.Perp() << endl;          
	fGoodTracksPT[iC] = false; 
      }
      
      if (pt->fPlab.Perp() > TRACKPTHI) {
	if (fVerbose > 1) cout << "track " << ps->fIndex << " failed pt: " << pt->fPlab.Perp() << endl;          
	fGoodTracksPT[iC] = false; 
      }
    }
  }

}


// ----------------------------------------------------------------------
void bmmReader::muonSelection() {

  TAnaTrack *pt, *ps; 
  
  for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
    for (int it = fCands[iC]->fSig1; it <= fCands[iC]->fSig2; ++it) {
      ps = fpEvt->getSigTrack(it); 
      pt = fpEvt->getRecTrack(ps->fIndex); 
      if (0 == (pt->fMuID & MUID)) {
	if (fVerbose > 1) cout << "muon " << ps->fIndex << " failed MUID: " << pt->fMuID << endl;          
	fGoodMuonsID[iC] = false; 
      }
      if (pt->fPlab.Perp() < MUPTLO) {
	if (fVerbose > 1) cout << "muon " << ps->fIndex << " failed MUPTLO: " << pt->fPlab.Perp() << endl;          
	fGoodMuonsPT[iC] = false; 
      }
      if (pt->fPlab.Perp() > MUPTHI) {
	if (fVerbose > 1) cout << "muon " << ps->fIndex << " failed MUPTHI: " << pt->fPlab.Perp() << endl;          
	fGoodMuonsPT[iC] = false; 
      }
    }
  }                                     
}

// ----------------------------------------------------------------------
void bmmReader::candidateSelection(int mode) {
  cout << "==>bmmReader::candidateSelection: nothing implemented" << endl;
}

// ----------------------------------------------------------------------
void bmmReader::fillCandidateVariables() {
  fCandPt = fpCand->fPlab.Perp();
  fCandM  = fpCand->fMass;
}


// ----------------------------------------------------------------------
void bmmReader::fillHist() {

  ((TH1D*)fpHistFile->Get("b1"))->Fill(fpEvt->nRecTracks()); 

  // -- only candidate histograms below
  if (0 == fpCand) return;

  fTree->Fill(); 

}

// ---------------------------------------------------------------------- 
void bmmReader::bookHist() {
  cout << "==> bmmReader: bookHist " << endl;
  
  TH1D *h; 
  h = new TH1D("b1", "Ntrk", 200, 0., 200.);
  h = new TH1D("bnc0", "NCand before selection", 20, 0., 20.);
  h = new TH1D("bnc1", "NCand after selection", 20, 0., 20.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun, "run/I");
  fTree->Branch("pt",     &fCandPt, "pt/D");
  fTree->Branch("m",      &fCandM,  "m/D");
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
  string cstring; 
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      if (dump) cout << "TYPE:           " << TYPE << endl;
      if (1313 == TYPE) cstring = "#mu^{+}#mu^{-}";
      if (200521 == TYPE) cstring = "#mu^{+}#mu^{-}K^{+}";
      ibin = 1;
      hcuts->SetBinContent(ibin, TYPE);
    }

    if (!strcmp(CutName, "CANDPTLO")) {
      CANDPTLO = CutValue; ok = 1;
      if (dump) cout << "CANDPTLO:           " << CANDPTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, CANDPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("p_{T}^{min}(%s) [GeV]", cstring.c_str()));
    }

    if (!strcmp(CutName, "CANDETALO")) {
      CANDETALO = CutValue; ok = 1;
      if (dump) cout << "CANDETALO:           " << CANDETALO << endl;
      ibin = 12;
      hcuts->SetBinContent(ibin, CANDETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("#eta^{min}(%s)", cstring.c_str()));
    }

    if (!strcmp(CutName, "CANDETAHI")) {
      CANDETAHI = CutValue; ok = 1;
      if (dump) cout << "CANDETAHI:           " << CANDETAHI << endl;
      ibin = 13;
      hcuts->SetBinContent(ibin, CANDETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("#eta^{max}(%s)", cstring.c_str()));
    }

    if (!strcmp(CutName, "SIGBOXMIN")) {
      SIGBOXMIN = CutValue; ok = 1;
      if (dump) cout << "SIGBOXMIN:           " << SIGBOXMIN << endl;
      ibin = 14;
      hcuts->SetBinContent(ibin, SIGBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, "SIGBOXMIN");
    }

    if (!strcmp(CutName, "SIGBOXMAX")) {
      SIGBOXMAX = CutValue; ok = 1;
      if (dump) cout << "SIGBOXMAX:           " << SIGBOXMAX << endl;
      ibin = 15;
      hcuts->SetBinContent(ibin, SIGBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, "SIGBOXMAX");
    }

    if (!strcmp(CutName, "BGLBOXMIN")) {
      BGLBOXMIN = CutValue; ok = 1;
      if (dump) cout << "BGLBOXMIN:           " << BGLBOXMIN << endl;
      ibin = 16;
      hcuts->SetBinContent(ibin, BGLBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, "BGLBOXMIN");
    }

    if (!strcmp(CutName, "BGLBOXMAX")) {
      BGLBOXMAX = CutValue; ok = 1;
      if (dump) cout << "BGLBOXMAX:           " << BGLBOXMAX << endl;
      ibin = 17;
      hcuts->SetBinContent(ibin, BGLBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, "BGLBOXMAX");
    }

    if (!strcmp(CutName, "BGHBOXMIN")) {
      BGHBOXMIN = CutValue; ok = 1;
      if (dump) cout << "BGHBOXMIN:           " << BGHBOXMIN << endl;
      ibin = 18;
      hcuts->SetBinContent(ibin, BGHBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, "BGHBOXMIN");
    }

    if (!strcmp(CutName, "BGHBOXMAX")) {
      BGHBOXMAX = CutValue; ok = 1;
      if (dump) cout << "BGHBOXMAX:           " << BGHBOXMAX << endl;
      ibin = 19;
      hcuts->SetBinContent(ibin, BGHBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, "BGHBOXMAX");
    }

    if (!strcmp(CutName, "TRACKQUALITY")) {
      TRACKQUALITY = CutValue; ok = 1;
      if (dump) cout << "TRACKQUALITY:           " << TRACKQUALITY << " " << endl;
      ibin = 20;
      hcuts->SetBinContent(ibin, TRACKQUALITY);
      hcuts->GetXaxis()->SetBinLabel(ibin, "track quality");
    }

    if (!strcmp(CutName, "TRACKPTLO")) {
      TRACKPTLO = CutValue; ok = 1;
      if (dump) cout << "TRACKPTLO:           " << TRACKPTLO << " GeV" << endl;
      ibin = 21;
      hcuts->SetBinContent(ibin, TRACKPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(track)");
    }

    if (!strcmp(CutName, "TRACKPTHI")) {
      TRACKPTHI = CutValue; ok = 1;
      if (dump) cout << "TRACKPTHI:           " << TRACKPTHI << " GeV" << endl;
      ibin = 22;
      hcuts->SetBinContent(ibin, TRACKPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{max}(track)");
    }

    if (!strcmp(CutName, "TRACKTIP")) {
      TRACKTIP = CutValue; ok = 1;
      if (dump) cout << "TRACKTIP:           " << TRACKTIP << " cm" << endl;
      ibin = 23;
      hcuts->SetBinContent(ibin, TRACKTIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, "doca_{xy}(track)");
    }

    if (!strcmp(CutName, "TRACKLIP")) {
      TRACKLIP = CutValue; ok = 1;
      if (dump) cout << "TRACKLIP:           " << TRACKLIP << " cm" << endl;
      ibin = 24;
      hcuts->SetBinContent(ibin, TRACKLIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, "doca_{z}(track)");
    }




    if (!strcmp(CutName, "MUID")) {
      MUID = int(CutValue); ok = 1;
      if (dump) cout << "MUID:           " << MUID << endl;
      ibin = 30;
      hcuts->SetBinContent(ibin, MUID);
      hcuts->GetXaxis()->SetBinLabel(ibin, "MuID");
    }

    if (!strcmp(CutName, "MUPTLO")) {
      MUPTLO = CutValue; ok = 1;
      if (dump) cout << "MUPTLO:           " << MUPTLO << " GeV" << endl;
      ibin = 31;
      hcuts->SetBinContent(ibin, MUPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(#mu)");
    }

    if (!strcmp(CutName, "MUPTHI")) {
      MUPTHI = CutValue; ok = 1;
      if (dump) cout << "MUPTHI:           " << MUPTHI << " GeV" << endl;
      ibin = 32;
      hcuts->SetBinContent(ibin, MUPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{max}(#mu)");
    }

    if (!strcmp(CutName, "MUETALO")) {
      MUETALO = CutValue; ok = 1;
      if (dump) cout << "MUETALO:           " << MUETALO << endl;
      ibin = 33;
      hcuts->SetBinContent(ibin, MUETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{min}(#mu)");
    }

    if (!strcmp(CutName, "MUETAHI")) {
      MUETAHI = CutValue; ok = 1;
      if (dump) cout << "MUETAHI:           " << MUETAHI << endl;
      ibin = 34;
      hcuts->SetBinContent(ibin, MUETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta^{max}(#mu)");
    }


    if (!ok) cout << "==> bmmReader: nothing done with " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}
