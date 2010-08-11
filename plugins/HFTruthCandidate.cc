#include <algorithm>

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

using reco::Track;
using reco::Vertex;
using reco::TrackBaseRef;

// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace edm;
using namespace std;

// ----------------------------------------------------------------------
HFTruthCandidate::HFTruthCandidate(const edm::ParameterSet& iConfig):
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", string("goodTracks"))), 
  fPartialDecayMatching(iConfig.getUntrackedParameter<bool>("partialDecayMatching", false)), 
  fMotherID(iConfig.getUntrackedParameter("motherID", 0)), 
  fType(iConfig.getUntrackedParameter("type", 67)),
  fGenType(iConfig.getUntrackedParameter("GenType", -67)),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)) {

  vector<int> defaultIDs;
  defaultIDs.push_back(0);
  fDaughtersID = iConfig.getUntrackedParameter<vector<int> >("daughtersID", defaultIDs);
  
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFTruthCandidate constructor" << endl;
  cout << "--- verbose:               " << fVerbose << endl;
  cout << "--- tracksLabel:           " << fTracksLabel << endl;
  cout << "--- motherID:              " << fMotherID << endl;
  cout << "--- type:                  " << fType << endl;
  cout << "--- GenType:               " << fGenType << endl;
  fDaughtersSet.clear(); 
  fStableDaughters = 0; 
  for (unsigned int i = 0; i < fDaughtersID.size(); ++i) {
    cout << "---   daughterID:              " << fDaughtersID[i] << endl;
    if (TMath::Abs(fDaughtersID[i]) == 11)   ++fStableDaughters; 
    if (TMath::Abs(fDaughtersID[i]) == 13)   ++fStableDaughters; 
    if (TMath::Abs(fDaughtersID[i]) == 211)  ++fStableDaughters; 
    if (TMath::Abs(fDaughtersID[i]) == 321)  ++fStableDaughters; 
    if (TMath::Abs(fDaughtersID[i]) == 2212) ++fStableDaughters; 
    fDaughtersSet.insert(TMath::Abs(fDaughtersID[i])); 
    fDaughtersGammaSet.insert(TMath::Abs(fDaughtersID[i])); 
    fDaughtersGamma2Set.insert(TMath::Abs(fDaughtersID[i])); 
  }    
  cout << "---    total stable particles: " << fStableDaughters << endl;
  fDaughtersGammaSet.insert(22); 
  fDaughtersGamma2Set.insert(22); 
  fDaughtersGamma2Set.insert(22); 
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
  TGenCand *pGen, *pDau, *pTmp;
  int matchedDecay(0);
  int iMom(-1), motherIndex(-1); 

  vector<int> bla(100); 
  vector<int>::iterator blaIt; 
  
  //   cout << "----------------------------------------------------------------------" << endl;
  for (int ig = 0; ig < gHFEvent->nGenCands(); ++ig) {
    pGen = gHFEvent->getGenCand(ig);
    if (TMath::Abs(pGen->fID) == fMotherID) {
      motherIndex = ig; 
      if (fVerbose > 1) {
	cout << "mother ";
	pGen->dump(); 
      }
      genDaughters.clear(); 
      genIndices.clear();
      genMap.clear();
      // -- version with direct daughters
      //       for (int id = pGen->fDau1; id <= pGen->fDau2; ++id) {
      // 	pDau = gHFEvent->getGenCand(id);
      // 	//	cout << "  daug: ";
      // 	//	pDau->dump(); 
      // 	genDaughters.insert(TMath::Abs(pDau->fID)); 
      // 	genIndices.insert(id); 
      // 	genMap.insert(make_pair(id, TMath::Abs(pDau->fID))); 
      //       }

      // -- version with descendants
      for (int id = ig+1; id < gHFEvent->nGenCands(); ++id) {
	pDau = gHFEvent->getGenCand(id);
	iMom = pDau->fMom1;
	while (iMom > ig) {
	  pTmp = gHFEvent->getGenCand(iMom);
	  iMom = pTmp->fMom1;
	}
	if (iMom == ig) {
	  if (fVerbose > 1) {
	    cout << "  daug: ";
	    pDau->dump(); 
	  }
	  genDaughters.insert(TMath::Abs(pDau->fID)); 
	  genIndices.insert(id); 
	  genMap.insert(make_pair(id, TMath::Abs(pDau->fID))); 
	}
      }

      // -- now check whether this is PARTIALLY the decay channel in question
      if (fPartialDecayMatching) {
	blaIt = set_intersection(genDaughters.begin(), genDaughters.end(), fDaughtersSet.begin(), fDaughtersSet.end(), bla.begin()); 
	if (static_cast<unsigned int>(blaIt - bla.begin()) == fDaughtersSet.size()) {
	  matchedDecay = 1; 
	  if (fVerbose > 0) {
	    cout << "matched partial decay: ";
	    for (vector<int>::iterator it = bla.begin(); it != blaIt; ++it) cout << *it << " "; 
	    cout << endl;
	  }
	  break;
	}
      }

      // -- now check whether this is the decay channel in question
      if (fDaughtersSet == genDaughters) {
	matchedDecay = 1; 
	if (fVerbose > 0) cout << "matched decay" << endl;
	break;
      }
      if (fDaughtersGammaSet == genDaughters) {
	matchedDecay = 1; 
	if (fVerbose > 0) cout << "matched decay with bremsstrahlung photon" << endl;
	break;
      }
      if (fDaughtersGamma2Set == genDaughters) {
	matchedDecay = 1; 
	if (fVerbose > 0) cout << "matched decay with 2 bremsstrahlung photons" << endl;
	break;
      }
    }
  }


  // -- Dump generator candidate made from stable charged particles
  if (matchedDecay > 0) {
    int id(-1), idx(-1); 
    TLorentzVector comp; 
    for (multiset<pair<int, int> >::iterator i = genMap.begin(); i != genMap.end(); ++i) {
      idx = i->first; 
      id  = i->second; 
      if (id == 11 || id == 13 || id == 211 || id == 321 || id ==2212)  {
	comp += gHFEvent->getGenCand(idx)->fP ; 
      }
    }

    TAnaCand *pCand = gHFEvent->addCand();
    pCand->fPlab = comp.Vect();
    pCand->fMass = comp.M();
    pCand->fType = fGenType;
    pCand->fIndex= motherIndex;
    if (fVerbose > 1) {
      char line[200];
      sprintf(line, "p=%8.3f(%+9.3f,%+9.3f,%+9.3f), mass = %f", 
	      pCand->fPlab.Mag(), 
	      pCand->fPlab.X(), pCand->fPlab.Y(), pCand->fPlab.Z(), 
	      pCand->fMass);
      cout << line << endl;
    }

  }


  if (fVerbose > 2 && 0 == matchedDecay)  {
    cout << "Did not match decay" << endl;
    for (multiset<int>::iterator i = genDaughters.begin(); i != genDaughters.end(); ++i) {
      cout << " unmatched genDaughter: " << *i << endl;
    }
  }




  // -- Construct and dump reconstructed candidates matched to generator particles
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

    if (static_cast<int>(trackList.size()) == fStableDaughters) {
      aKal.doNotFit(trackList, trackIndices, trackMasses, fType); 
    }
    
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFTruthCandidate);
