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
  cout << "==> bmmSignalReader: Starting analysis..." << endl;
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
void bmmSignalReader::candidateSelection(int mode) {
  int nc0(fCands.size()), nc1(0);

  ((TH1D*)fpHistFile->Get("bnc0"))->Fill(nc0); 
  if (0 == nc0) {
    cout << "no candidate " << TYPE << " found" << endl;
    return; 
  }
  
  fpCand = 0; 
  TAnaCand *pCand; 

  if (0 == nc0) {
    fpCand = fCands[0]; 
  }

  // -- Create list with candidates that fullfil basic cuts on their final state particles
  vector<int> clist; 
  if (nc0 > 1) {
    for (unsigned int iC = 0; iC < fCands.size(); ++iC) {
      pCand = fCands[iC]; 
      
      TAnaTrack *pl1 = fpEvt->getSigTrack(pCand->fSig1); 
      TAnaTrack *pl2 = fpEvt->getSigTrack(pCand->fSig2); 

      if (pl1->fQ*pl2->fQ > 0) continue;

      if (false == fGoodTracks[iC])   continue; 
      if (false == fGoodTracksPT[iC]) continue; 

      if (false == fGoodMuonsID[iC])  continue; 
      if (false == fGoodMuonsPT[iC])  continue; 

      clist.push_back(iC);
    }
  } 

  // -- Select best candidate
  if (clist.size() > 1) {
    double drMin(99.), dr(0.); 
    for (unsigned int i = 0; i < clist.size(); ++i) {
      pCand = fCands[clist[i]]; 
      TAnaTrack *pl1 = fpEvt->getSigTrack(pCand->fSig1); 
      TAnaTrack *pl2 = fpEvt->getSigTrack(pCand->fSig1+1); 

      // -- closest in rphi
      if (0 == mode) {
	dr = pl1->fPlab.DeltaR(pl2->fPlab); 
	if (dr < drMin) {
	  fpCand = pCand; 
	  drMin = dr; 
	}
      }
    }
  } else {
    fpCand = fCands[0];
  }

  ((TH1D*)fpHistFile->Get("bnc1"))->Fill(nc1); 

  if (BLIND && fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX) fpCand = 0; 
  if (fpCand) fillCandidateVariables();
}

// ----------------------------------------------------------------------
void bmmSignalReader::fillHist() {
  bmmReader::fillHist(); 
}

// ---------------------------------------------------------------------- 
void bmmSignalReader::bookHist() {
  bmmReader::bookHist();
  cout << "==> bmmReader: bookHist " << endl;
  TH1D *h1; 
  h1 = new TH1D("h100", "", 100, 0., 100.);
  fTree->Branch("m2",      &fCandM,  "m2/D");
}

// ----------------------------------------------------------------------
void bmmSignalReader::readCuts(TString filename, int dump) {
  bmmReader::readCuts(filename, dump); 

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

  TH1D *hcuts = (TH1D*)gFile->Get("hcuts"); 
  if (0 == hcuts) {
    new TH1D("hcuts", "", 1000, 0., 1000.);
    hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  }
  int ibin; 
  string cstring; 
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!ok) cout << "==> bmmSignalReader: nothing done with  " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}

