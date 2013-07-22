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
  BLIND = 0; 
}

// ----------------------------------------------------------------------
int  treeReader01::numberOfBPixLayers(TAnaTrack *pTrack) {
  bool layer1 = false, layer2 = false, layer3 = false;
  int hits = pTrack->fValidHits;
  //cout<<" muon1 "<<algo1<<" "<<qual1<<" "<<hits1<<hex<<" ";
  if(hits>20) hits=20; // pattern has only 20 locations
  for(int i =0; i<hits; ++i){
    int pat = pTrack->fHitPattern[i];
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

// ----------------------------------------------------------------------
int  treeReader01::numberOfPixLayers(TAnaTrack *pTrack) {
  bool layer1 = false, layer2 = false, layer3 = false, disk1=false, disk2=false;
  int hits = pTrack->fValidHits;
  //cout<<" muon1 "<<algo1<<" "<<qual1<<" "<<hits1<<hex<<" ";
  if(hits>20) hits=20; // pattern has only 20 locations
  for(int i =0; i<hits; ++i){
    unsigned int pat = pTrack->fHitPattern[i];
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


// ----------------------------------------------------------------------
int treeReader01::numberOfPixelHits(TAnaTrack *pTrack) {
  int hits = pTrack->fValidHits;
  int pixhits(0); 
  if(hits>20) hits=20; // pattern has only 20 locations
  for(int i =0; i<hits; ++i){
    unsigned int pat = pTrack->fHitPattern[i];
    //cout<<pat<<" ";
    if( pat == 0x488 ) ++pixhits;
    else if( pat == 0x490 )  ++pixhits;
    else if( pat == 0x498 )  ++pixhits;
    else if( pat == 0x508 )  ++pixhits;
    else if( pat == 0x510 )  ++pixhits;
  }
  //cout<<dec<<endl;
  
  return pixhits;
}

// ----------------------------------------------------------------------
int  treeReader01::numberOfBPixLayer1Hits(TAnaTrack *pTrack) {
  int pixHits=0;
  int hits = pTrack->fValidHits;
  //cout<<" muon1 "<<algo1<<" "<<qual1<<" "<<hits1<<hex<<" ";
  if(hits>20) hits=20; // pattern has only 20 locations
  for(int i =0; i<hits; ++i){
    unsigned int pat = pTrack->fHitPattern[i];
    //cout<<pat<<" ";
    if( pat == 0x488 ) pixHits++;
    //else if( pat == 0x490 ) layer2 = true;
    //else if( pat == 0x498 ) layer3 = true;
  }
  //cout<<dec<<endl;

  return pixHits;
}


// ----------------------------------------------------------------------
// Hit pattern is the summary information of the hits associated to track in
// AOD.  When RecHits are no longer available, the compact hit pattern should
// allow basic track selection based on the hits in various subdetectors.  The
// hits of a track are saved in unit32_t hitPattern_[28], initialized as
// 0x00000000, ..., 0x00000000.  Set one hit with 10 bits
//
//      +-----+-----+-----+-----+-----+-----+-----+-----+----------------+-----+-----+
//      |tk/mu|  sub-structure  |   sub-sub-structure   |     stereo     |  hit type |
//      +-----+-----+-----+-----+-----+-----+-----+-----+----------------+-----+-----+
//  ... | 10  |   9    8     7  |   6    5     4     3  |        2       |  1     0  | bit
//  ... |  2  |   1    0     3  |   2    1     0     3  |        2       |  1     0  | bit
//
//      |tk = 1      PXB = 1            layer = 1-3                       hit type = 0-3
//      |tk = 1      PXF = 2            disk  = 1-2                       hit type = 0-3
//      |tk = 1      TIB = 3            layer = 1-4      0=rphi,1=stereo  hit type = 0-3
//      |tk = 1      TID = 4            wheel = 1-3      0=rphi,1=stereo  hit type = 0-3
//      |tk = 1      TOB = 5            layer = 1-6      0=rphi,1=stereo  hit type = 0-3
//      |tk = 1      TEC = 6            wheel = 1-9      0=rphi,1=stereo  hit type = 0-3
//      |mu = 0      DT  = 1            4*(stat-1)+superlayer             hit type = 0-3
//      |mu = 0      CSC = 2            4*(stat-1)+(ring-1)               hit type = 0-3
//      |mu = 0      RPC = 3            4*(stat-1)+2*layer+region         hit type = 0-3
//
//      hit type, see DataFormats/TrackingRecHit/interface/TrackingRecHit.h
//      valid    = valid hit                                     = 0
//      missing  = detector is good, but no rec hit found        = 1
//      inactive = detector is off, so there was no hope         = 2
//      bad      = there were many bad strips within the ellipse = 3
// ----------------------------------------------------------------------
int  treeReader01::numberOfTrackerLayers(TAnaTrack *pTrack) {
  bool pixl[3], tibl[4], tobl[6]; 
  bool pixd[2], tidw[3], tecw[9];

  for (int i = 0; i < 3; ++i) pixl[i] = false; 
  for (int i = 0; i < 4; ++i) tibl[i] = false; 
  for (int i = 0; i < 6; ++i) tobl[i] = false; 

  for (int i = 0; i < 2; ++i) pixd[i] = false; 
  for (int i = 0; i < 3; ++i) tidw[i] = false; 
  for (int i = 0; i < 9; ++i) tecw[i] = false; 

  int hits = pTrack->fValidHits;
  if(hits>20) hits=20; // pattern has only 20 locations
  //  cout << "----------------------------------------------------------------------" << endl;
  int hit(0), hitmask(3); 
  int det(0), detpos(7), detmask(0); 
  int lay(0), layerpos(3), layermask(0); 
  detmask = 0x7 << detpos;
  layermask = 0xf << layerpos;
  //  cout << "detmask = " << std::hex << detmask << " laymask = " << layermask << std::dec << endl;
  for(int i =0; i<hits; ++i){
    unsigned int pat = pTrack->fHitPattern[i];
    
    hit = (pat & hitmask);
    det = 0; 
    det = (pat & detmask)>>detpos; 
    lay = 0; 
    lay = (pat & layermask)>>layerpos; 

    //    lay = lay - 1; // FIXME this line is necessary to be correct. But you should use TAnaTrack::fNumberOfValidTrkHits anyway!
    if ((1 == det) && (0 == hit)) pixl[lay] = true; 
    if ((2 == det) && (0 == hit)) pixd[lay] = true; 
	 	       	      
    if ((3 == det) && (0 == hit)) tibl[lay] = true; 
    if ((4 == det) && (0 == hit)) tidw[lay] = true; 
	 	       	      
    if ((5 == det) && (0 == hit)) tobl[lay] = true; 
    if ((6 == det) && (0 == hit)) tecw[lay] = true; 

  }
  
  int trkHits(0); 
  for (int i = 0; i < 3; ++i) {
    if (pixl[i]) {
      ++trkHits;
    }
  }

  for (int i = 0; i < 4; ++i) {
    if (tibl[i]) {
      ++trkHits; 
    }
  }

  for (int i = 0; i < 6; ++i) {
    if (tobl[i]) {
      ++trkHits; 
    }
  }

  for (int i = 0; i < 2; ++i) {
    if (pixd[i]) {
      ++trkHits; 
    }
  }

  for (int i = 0; i < 3; ++i) {
    if (tidw[i]) {
      ++trkHits; 
    }
  }

  for (int i = 0; i < 9; ++i) {
    if (tecw[i]) {
      ++trkHits; 
    }
  }
  return trkHits;
}
