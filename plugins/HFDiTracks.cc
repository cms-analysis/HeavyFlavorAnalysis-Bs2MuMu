#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDiTracks.h"

#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKinematicVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatorics.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace reco;
using namespace edm;

// ----------------------------------------------------------------------
HFDiTracks::HFDiTracks(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fVertexing(iConfig.getUntrackedParameter<int>("vertexing", 1)),
  fTracksLabel(iConfig.getUntrackedParameter<string>("tracksLabel", string("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<string>("PrimaryVertexLabel", string("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 3.0)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 1.0)), 
  fTrackMass(iConfig.getUntrackedParameter<double>("trackMass", 0.1396)), 
  fMassLow(iConfig.getUntrackedParameter<double>("massLow", 0.0)), 
  fMassHigh(iConfig.getUntrackedParameter<double>("massHigh", 12.0)), 
  fType(iConfig.getUntrackedParameter<int>("type", 1300)) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDiTracks constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel.c_str() << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  vertexing                 " << fVertexing << endl;
  cout << "---  muonPt:                   " << fMuonPt << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  trackMass:                " << fTrackMass << endl;
  cout << "---  Type:                     " << fType << endl;
  cout << "---  massLow:                  " << fMassLow << endl;
  cout << "---  massHigh:                 " << fMassHigh << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDiTracks::~HFDiTracks() {
  
}


// ----------------------------------------------------------------------
void HFDiTracks::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel.c_str(), recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFDiTracks> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFDiTracks> No primary vertex found, skipping" << endl;
    return;
  }
  fPV = vertices[gHFEvent->fEventTag]; 
  if (fVerbose > 0) {
    cout << "HFDimuons: Taking vertex " << gHFEvent->fEventTag << " with ntracks = " << fPV.tracksSize() << endl;
  }
  
  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel.c_str(), hTracks);
  if (!hTracks.isValid()) {
    cout << "==>HFDiTracks> No valid TrackCollection with label "<<fTracksLabel.c_str() <<" found, skipping" << endl;
    return;
  }
  
  // -- Transient tracks for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFDiTracks: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  // -- set up vertex fitter 
  HFKalmanVertexFit a(fTTB.product(), fPV, fType, 0); 
  vector<Track> trackList; 
  vector<int> trackIndices;
  vector<double> trackMasses;

  TLorentzVector ditrack, t1, t2;
  for (unsigned int it1 = 0; it1 < hTracks->size()-1; ++it1) {
    TrackBaseRef t1TrackView(hTracks, it1);
    Track tTrack1(*t1TrackView);
    t1.SetPtEtaPhiM(tTrack1.pt(), tTrack1.eta(), tTrack1.phi(), fTrackMass); 
    if (tTrack1.pt() < fMuonPt)  continue;

    for (unsigned int it2 = it1 + 1; it2 < hTracks->size(); ++it2) {
      TrackBaseRef t2TrackView(hTracks, it2);
      Track tTrack2(*t2TrackView);
      if (tTrack2.pt() < fMuonPt)  continue;

      t2.SetPtEtaPhiM(tTrack2.pt(), tTrack2.eta(), tTrack2.phi(), fTrackMass); 
      ditrack = t1 + t2; 

      if (ditrack.M() < fMassLow || ditrack.M() > fMassHigh) {
	if (fVerbose > 0) {
	  cout << "==>HFDiTracks> ditrack mass = " << ditrack.M() << ", skipping" << endl;
	}
	continue; 
      }

      // -- Vertexing, new style
      trackList.clear();
      trackIndices.clear(); 
      trackMasses.clear(); 
      
      trackList.push_back(tTrack1); 
      trackIndices.push_back(it1); 
      trackMasses.push_back(fTrackMass);
      
      trackList.push_back(tTrack2); 
      trackIndices.push_back(it2); 
      trackMasses.push_back(fTrackMass);
      
      if (fVertexing > 0) {
	a.doFit(trackList, trackIndices, trackMasses); 	
      } else {
	a.doNotFit(trackList, trackIndices, trackMasses); 	
      }
      
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void  HFDiTracks::beginJob(const EventSetup& setup) {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDiTracks::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDiTracks);
