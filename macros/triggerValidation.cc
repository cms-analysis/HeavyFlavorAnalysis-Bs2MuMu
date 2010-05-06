#include "triggerValidation.hh"
#include "../interface/HFMasses.hh"

#include "TRandom.h"
#include <cmath>

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/triggerValidation.default.cuts
//           ./runTreeReaders -f test.root 
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
triggerValidation::triggerValidation(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> triggerValidation: constructor..." << endl;
}

// ----------------------------------------------------------------------
triggerValidation::~triggerValidation() {
  cout << "==> triggerValidation: destructor..." << endl;

}

// ----------------------------------------------------------------------
void triggerValidation::startAnalysis() {
  cout << "==> triggerValidation: Starting analysis..." << endl;

  for (int i = 0; i < 256; ++i) {
    fHLTWasRun[i] = fHLTResult[i] = fHLTWasRunResult[i] = 0; 
  }

  for (int i = 0; i < 128; ++i) {
    fL1TWasRun[i] = fL1TResult[i] = fL1TWasRunResult[i] = 0; 
    fLTTWasRun[i] = fLTTResult[i] = fLTTWasRunResult[i] = 0; 
  }
}


// ----------------------------------------------------------------------
void triggerValidation::eventProcessing() {

  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 

  lttValidation();
  l1tValidation();
  hltValidation();
  
}


// ----------------------------------------------------------------------
void triggerValidation::initVariables() {

}

// ----------------------------------------------------------------------
void triggerValidation::lttValidation() {

  int wasrun(0), result(0), word(0), bit(0);
  // cout << "----------------------------------------------------------------------" << endl;
  for (int i = 0; i < 32*NL1TT; ++i) {
    bit  = i%32;
    word = i/32; 
    wasrun = (fpEvt->fL1TTWasRun[word] & 0x1<<bit)?1:0;
    result = (fpEvt->fL1TTWords[word] & 0x1<<bit)?1:0;

    if (fLTTNames[i] != fpEvt->fL1TTNames[i]) {
      cout << "TT: Filling at i = " << i << " the new name " << fpEvt->fL1TTNames[i] << endl;
      fLTTNames[i] = fpEvt->fL1TTNames[i];
    } else {
      fLTTNames[i] = fpEvt->fL1TTNames[i];
    }

    if (wasrun) fLTTWasRun[i] += 1; 
    if (result) fLTTResult[i] += 1; 
    if (wasrun && result) fLTTWasRunResult[i] += 1; 

    // cout << Form("%3i %50s %4i %4i", i, fpEvt->fHLTNames[i].Data(), wasrun, result) << endl;
  }


}


// ----------------------------------------------------------------------
void triggerValidation::l1tValidation() {

  int wasrun(0), result(0), word(0), bit(0);
  // cout << "----------------------------------------------------------------------" << endl;
  for (int i = 0; i < 32*NL1T; ++i) {
    bit  = i%32;
    word = i/32; 
    wasrun = (fpEvt->fL1TWasRun[word] & 0x1<<bit)?1:0;
    result = (fpEvt->fL1TWords[word] & 0x1<<bit)?1:0;

    if (fL1TNames[i] != fpEvt->fL1TNames[i]) {
      cout << "L1: Filling at i = " << i << " the new name " << fpEvt->fL1TNames[i] << endl;
      fL1TNames[i] = fpEvt->fL1TNames[i];
    } else {
      fL1TNames[i] = fpEvt->fL1TNames[i];
    }

    if (wasrun) fL1TWasRun[i] += 1; 
    if (result) fL1TResult[i] += 1; 
    if (wasrun && result) fL1TWasRunResult[i] += 1; 

    // cout << Form("%3i %50s %4i %4i", i, fpEvt->fHLTNames[i].Data(), wasrun, result) << endl;
  }


}


// ----------------------------------------------------------------------
void triggerValidation::hltValidation() {

  int wasrun(0), result(0), word(0), bit(0);
  // cout << "----------------------------------------------------------------------" << endl;
  for (int i = 0; i < 32*NHLT; ++i) {
    if (!strcmp(fpEvt->fHLTNames[i].Data(), "")) break;

    bit  = i%32;
    word = i/32; 

    // cout << bit << " " << word << endl;

    wasrun = (fpEvt->fHLTWasRun[word] & 0x1<<bit)?1:0;
    result = (fpEvt->fHLTWords[word] & 0x1<<bit)?1:0;

    if (fHLTNames[i] != fpEvt->fHLTNames[i]) {
      cout << "HL: Filling at i = " << i << " the new name " << fpEvt->fHLTNames[i] << endl;
      fHLTNames[i] = fpEvt->fHLTNames[i];
    } else {
      fHLTNames[i] = fpEvt->fHLTNames[i];
    }

    if (wasrun) fHLTWasRun[i] += 1; 
    if (result) fHLTResult[i] += 1; 
    if (wasrun && result) fHLTWasRunResult[i] += 1; 

    // cout << Form("%3i %50s %4i %4i", i, fpEvt->fHLTNames[i].Data(), wasrun, result) << endl;
  }
		 
}



// ----------------------------------------------------------------------
void triggerValidation::fillHist() {

  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 

}

// ---------------------------------------------------------------------- 
void triggerValidation::bookHist() {
  cout << "==> triggerValidation: bookHist " << endl;
  
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

}


// ----------------------------------------------------------------------
void triggerValidation::closeHistFile() {
  cout << "==> triggerValidation: closeHist " << endl;
  
  for (int i = 0; i < 256; ++i) {
    if (!strcmp(fHLTNames[i].Data(), "")) break;
    
    cout << Form("%3i %50s %4i %4i %4i", i, fHLTNames[i].Data(), fHLTWasRun[i], fHLTResult[i], fHLTWasRunResult[i]) << endl;
  }

  for (int i = 0; i < 128; ++i) {
    cout << Form("%3i %50s %4i %4i %4i", i, fL1TNames[i].Data(), fL1TWasRun[i], fL1TResult[i], fL1TWasRunResult[i]) << endl;
  }

  for (int i = 0; i < 128; ++i) {
    cout << Form("%3i %50s %4i %4i %4i", i, fLTTNames[i].Data(), fLTTWasRun[i], fLTTResult[i], fLTTWasRunResult[i]) << endl;
  }
  


  fpHistFile->cd();
  fpHistFile->Write();
  fpHistFile->Close();
  delete fpHistFile;

}


// ----------------------------------------------------------------------
void triggerValidation::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> triggerValidation: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> triggerValidation: Cut file  " << fCutFile.Data() << endl;
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

    if (!ok) cout << "==> triggerValidation: ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}
