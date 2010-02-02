#include <set>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFTruthCandidate.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>


#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace edm;
using namespace std;

// ----------------------------------------------------------------------
HFTruthCandidate::HFTruthCandidate(const edm::ParameterSet& iConfig):
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", string("goodTracks"))), 
  fMotherID(iConfig.getUntrackedParameter("motherID", 0)), 
  fType(iConfig.getUntrackedParameter("type", 67)) {

  vector<int> defaultIDs;
  defaultIDs.push_back(0);
  fDaughtersID = iConfig.getUntrackedParameter<vector<int> >("daughtersID", defaultIDs);

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFTruthCandidate constructor" << endl;
  cout << "--- motherID:              " << fMotherID << endl;
  fDaughtersSet.clear(); 
  for (unsigned int i = 0; i < fDaughtersID.size(); ++i) {
    cout << "---   daughterID:              " << fDaughtersID[i] << endl;
    fDaughtersSet.insert(TMath::Abs(fDaughtersID[i])); 
    fDaughtersGammaSet.insert(TMath::Abs(fDaughtersID[i])); 
  }    
  fDaughtersGammaSet.insert(22); 
  cout << "----------------------------------------------------------------------" << endl;
}

// ----------------------------------------------------------------------
HFTruthCandidate::~HFTruthCandidate() {  
}


// ----------------------------------------------------------------------
void HFTruthCandidate::beginJob() {
}

// ----------------------------------------------------------------------
void HFTruthCandidate::endJob() {
}


// ----------------------------------------------------------------------
void HFTruthCandidate::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- In generator block, find mother with declared decay channel
  multiset<int> genDaughters; 
  multiset<int> genIndices; 
  multiset<pair<int, int> > genMap; 
  TGenCand *pGen, *pDau;
  int matchedDecay(0); 
  for (int ig = 0; ig < gHFEvent->nGenCands(); ++ig) {
    pGen = gHFEvent->getGenCand(ig);
    if (TMath::Abs(pGen->fID) == fMotherID) {
      //      cout << "mother ";
      //      pGen->dump(); 
      genDaughters.clear(); 
      genIndices.clear();
      for (int id = pGen->fDau1; id <= pGen->fDau2; ++id) {
	pDau = gHFEvent->getGenCand(id);
	//	cout << "  daug: ";
	//	pDau->dump(); 
	genDaughters.insert(TMath::Abs(pDau->fID)); 
	genIndices.insert(id); 
	genMap.insert(make_pair(id, TMath::Abs(pDau->fID))); 
      }
      // -- now check whether this is the decay channel in question
      if (fDaughtersSet == genDaughters) {
	matchedDecay = 1; 
	//	cout << "matched decay" << endl;
      }
      if (fDaughtersGammaSet == genDaughters) {
	matchedDecay = 1; 
	//	cout << "matched decay with bremsstrahlung photon" << endl;
      }
      
      //       if (matchedDecay > 0) {
      // 	for (multiset<int>::iterator i = genDaughters.begin(); i != genDaughters.end(); ++i) {
      // 	  //	  cout << " genDaughter: " << *i << endl;
      // 	}
      //       } else {
      // 	//	cout << "Did not match decay" << endl;
      // 	for (multiset<int>::iterator i = genDaughters.begin(); i != genDaughters.end(); ++i) {
      // 	  //	  cout << " unmatched genDaughter: " << *i << endl;
      // 	}
      //       }
    }
  }

//   // -- Remove decayed particles from daughter list
//   for (multiset<pair<int, int> >::iterator i = genMap.begin(); i != genMap.end(); ++i) {
//     int stable = 0; 
//     if (*i.first == 11)   stable = 1;  
//     if (*i.first == 13)   stable = 1;  
//     if (*i.first == 22)   stable = 1;  
//     if (*i.first == 211)  stable = 1;  
//     if (*i.first == 321)  stable = 1;  
//     if (*i.first == 2212) stable = 1;  
//     if (0 == stable) {
//       pGen = gHFEvent->getGenCand(i.second); 
//       for (int id = pGen->fDau1; id <= pGen->fDau2; ++id) {
// 	pDau = gHFEvent->getGenCand(id);
// 	cout << "  daug: ";
// 	pDau->dump(); 
// 	genDaughters.insert(TMath::Abs(pDau->fID)); 
// 	genIndices.insert(id); 
// 	genMap.insert(make_pair(id, TMath::Abs(pDau->fID))); 
//       }
      
//       cout << " genDaughter: " << *i << endl;
//     }
//   }
  


  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFBu2JpsiKp> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  Vertex dummy; 
  HFKalmanVertexFit  aKal(0, dummy, 1, fType); 
  vector<Track> trackList; 
  vector<int> trackIndices;
  vector<double> trackMasses;

  TAnaTrack *pTrack; 
  if (matchedDecay > 0) {
    for (int it = 0; it < gHFEvent->nRecTracks(); ++it) {
      pTrack = gHFEvent->getRecTrack(it); 
      if (genIndices.find(pTrack->fGenIndex) != genIndices.end()) {
	//	cout << "Found rec track: " << it; 
	//	pTrack->dump(); 
	
	TrackBaseRef TrackView(hTracks, it);
	Track track(*TrackView);
	trackList.push_back(track); 
	trackIndices.push_back(it); 
	double mass = MMUON; 
	if (321  == TMath::Abs(pTrack->fMCID)) mass = MKAON;
	if (211  == TMath::Abs(pTrack->fMCID)) mass = MPION;
	if (13   == TMath::Abs(pTrack->fMCID)) mass = MMUON;
	if (2212 == TMath::Abs(pTrack->fMCID)) mass = MPROTON;
	trackMasses.push_back(mass); 
      }
    }

    if (trackList.size() == fDaughtersSet.size()) {
      aKal.doNotFit(trackList, trackIndices, trackMasses, fType); 
    }
    
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFTruthCandidate);
