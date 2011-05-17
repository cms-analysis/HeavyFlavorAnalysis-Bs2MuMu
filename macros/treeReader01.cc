#include "treeReader01.hh"

#include "TRandom.h"

using std::string;
using std::cout;
using std::endl;
#include "treeReader01.icc"


// ----------------------------------------------------------------------
// Run with: ./runTreeReader01 -c chains/bg-test -D root 
//           ./runTreeReader01 -f test.root 
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
void treeReader01::startAnalysis() {
  cout << "treeReader01: startAnalysis: ..." << endl;
  fpJSON = new JSON(JSONFILE.c_str()); 
}

// ----------------------------------------------------------------------
void treeReader01::endAnalysis() {
  cout << "treeReader01: endAnalysis: ..." << endl;
}


// ----------------------------------------------------------------------
bool treeReader01::goodRun() {
  return true; 
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
int  treeReader01::numberOfBPixLayers(TAnaTrack *pTrack) {
  bool layer1 = false, layer2 = false, layer3 = false;
  int hits = pTrack->fValidHits;
  //cout<<" muon1 "<<algo1<<" "<<qual1<<" "<<hits1<<hex<<" ";
  if(hits>20) hits=20; // pattern has only 20 locations
  for(int i =0; i<hits; ++i){
    uint32_t pat = pTrack->fHitPattern[i];
    //cout<<pat<<" ";
    if( pat == 0x488 ) layer1 = true;
    else if( pat == 0x490 ) layer2 = true;
    else if( pat == 0x498 ) layer3 = true;
  }
  //cout<<dec<<endl;
  
  int pixHits=0;
  if(layer1) {pixHits++;}
  if(layer2) {pixHits++;}
  if(layer3) {pixHits++;}

  return pixHits;
}
int  treeReader01::numberOfPixLayers(TAnaTrack *pTrack) {
  bool layer1 = false, layer2 = false, layer3 = false, disk1=false, disk2=false;
  int hits = pTrack->fValidHits;
  //cout<<" muon1 "<<algo1<<" "<<qual1<<" "<<hits1<<hex<<" ";
  if(hits>20) hits=20; // pattern has only 20 locations
  for(int i =0; i<hits; ++i){
    uint32_t pat = pTrack->fHitPattern[i];
    //cout<<pat<<" ";
    if( pat == 0x488 ) layer1 = true;
    else if( pat == 0x490 ) layer2 = true;
    else if( pat == 0x498 ) layer3 = true;
    else if( pat == 0x508 ) disk1 = true;
    else if( pat == 0x510 ) disk2 = true;
  }
  //cout<<dec<<endl;
  
  int pixHits=0;
  if(layer1) {pixHits++;}
  if(layer2) {pixHits++;}
  if(layer3) {pixHits++;}
  if(disk1) {pixHits++;}
  if(disk2) {pixHits++;}

  return pixHits;
}
int  treeReader01::numberOfBPixLayer1Hits(TAnaTrack *pTrack) {
  int pixHits=0;
  int hits = pTrack->fValidHits;
  //cout<<" muon1 "<<algo1<<" "<<qual1<<" "<<hits1<<hex<<" ";
  if(hits>20) hits=20; // pattern has only 20 locations
  for(int i =0; i<hits; ++i){
    uint32_t pat = pTrack->fHitPattern[i];
    //cout<<pat<<" ";
    if( pat == 0x488 ) pixHits++;
    //else if( pat == 0x490 ) layer2 = true;
    //else if( pat == 0x498 ) layer3 = true;
  }
  //cout<<dec<<endl;

  return pixHits;
}
