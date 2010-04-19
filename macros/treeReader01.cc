#include "treeReader01.hh"

#include "TRandom.h"

#include "treeReader01.icc"


// ----------------------------------------------------------------------
// Run with: ./runTreeReader01 -c chains/bg-test -D root 
//           ./runTreeReader01 -f test.root 
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
void treeReader01::startAnalysis() {
  cout << "treeReader01: startAnalysis: ..." << endl;
}


// ----------------------------------------------------------------------
bool treeReader01::goodRun() {
  
  if (fRun < 132440) return true; // assume this is MC

  if (fRun == 132606) return true; 
  if (fRun == 132605) return true; 
  if (fRun == 132602) return true; 
  if (fRun == 132601) return true; 
  if (fRun == 132599) return true; 
  if (fRun == 132598) return true; 
  if (fRun == 132597) return true; 
  if (fRun == 132596) return true; 
  if (fRun == 132572) return true; 
  if (fRun == 132569) return true; 
  if (fRun == 132478) return true; 
  if (fRun == 132477) return true; 
  if (fRun == 132476) return true; 
  if (fRun == 132473) return true; 
  if (fRun == 132471) return true; 
  if (fRun == 132442) return true; 
  if (fRun == 132440) return true; 
  
  return false; 
}



// ----------------------------------------------------------------------
void treeReader01::eventProcessing() {
  TAnaTrack *pTrack;
  TAnaCand *pCand;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "Found " << fpEvt->nGenCands() << " gen cands in event" << endl;
  cout << "Found " << fpEvt->nSigTracks() << " sig tracks in event" << endl;
  cout << "Found " << fpEvt->nRecTracks() << " rec tracks in event" << endl;
  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks()); 
  
  cout << "------------------------------" << endl;
  ((TH1D*)fpHistFile->Get("h20"))->Fill(fpEvt->nRecTracks());
  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {
    pTrack = fpEvt->getRecTrack(it);
    ((TH1D*)fpHistFile->Get("h10"))->Fill(pTrack->fPlab.Perp()); 
    cout << "R: "; pTrack->dump(); 
  }


  cout << "------------------------------" << endl;
  for (int it = 0; it < fpEvt->nCands(); ++it) {
    pCand = fpEvt->getCand(it);
    cout << "C: " << pCand->fType << " "; pCand->dump(); 
    ((TH1D*)fpHistFile->Get("h100"))->Fill(pCand->fMass); 
  }

  fpHistFile->cd();
  fillHist(); 
  fTree->Fill();

}


// ----------------------------------------------------------------------
void treeReader01::fillHist() {


}

// ----------------------------------------------------------------------
void treeReader01::bookHist() {
  TH1 *h;
  cout << "==> treeReader01: bookHist> " << endl;

  h = new TH1D("h1", "nTrk", 40, 0., 40.);
  h = new TH1D("h10", "pT", 40, 0., 20.);
  h = new TH1D("h20", "ntrk", 20, 0, 20.);

  h = new TH1D("h100", "m", 40, 2.8, 3.4);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun ,"run/I");

}


// ----------------------------------------------------------------------
void treeReader01::initVariables() {
  cout << "treeReader01: initVariables: ..." << endl;

  fRun = -1; 


}
