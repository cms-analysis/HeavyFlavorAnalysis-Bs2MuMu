#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFMuonAndTrack.h"

#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

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

#define MMUON 0.10566

// ----------------------------------------------------------------------
HFMuonAndTrack::HFMuonAndTrack(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fVertexing(iConfig.getUntrackedParameter<int>("vertexing", 1)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 3.0)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 1.0)), 
  fMassLow(iConfig.getUntrackedParameter<double>("massLow", 8.7)), 
  fMassHigh(iConfig.getUntrackedParameter<double>("massHigh", 10.2)), 
  fType(iConfig.getUntrackedParameter<int>("type", 1300)) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFMuonAndTrack constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  vertexing                 " << fVertexing << endl;
  cout << "---  muonPt:                   " << fMuonPt << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  massLow:                  " << fMassLow << endl;
  cout << "---  massHigh:                 " << fMassHigh << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFMuonAndTrack::~HFMuonAndTrack() {
  
}


// ----------------------------------------------------------------------
void HFMuonAndTrack::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // -- get the primary vertex
  edm::Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFMuonAndTrack> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFMuonAndTrack> No primary vertex found, skipping" << endl;
    return;
  }
  fPV = vertices[gHFEvent->fEventTag]; 
  if (fVerbose > 0) {
    cout << "HFDimuons: Taking vertex " << gHFEvent->fEventTag << " with ntracks = " << fPV.tracksSize() << endl;
  }
  
  // -- get the collection of muons
  Handle<MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  if (!hMuons.isValid()) {
  cout << "==>HFMuonAndTrack> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  edm::Handle<edm::View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFMuonAndTrack> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  // -- get the collection of muons and store their corresponding track indices
  vector<int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index();
    if (fVerbose > 2) cout << "muon->track().index() = "<< muon->track().index()<< endl; 
    if (im >= 0) muonIndices.push_back(im);
  }
  if (fVerbose > 0) {
    cout << "==>HFMuonAndTrack> nMuons = " << hMuons->size() << endl;
    cout << "==>HFMuonAndTrack> nMuonIndices = " << muonIndices.size() << endl;
  }
  if (muonIndices.size() < 1) return;
  
  // -- Transient tracks for vertexing
  try {
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  }  catch (Exception event) {
    // cout << "something bad happened" << endl;
  }

  if (!fTTB.isValid()) {
    cout << " -->HFMuonAndTrack: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  // -- set up vertex fitter 
  HFKalmanVertexFit a(fTTB.product(), fPV, 1300, 0); 
  vector<Track> trackList; 
  vector<int> trackIndices;
  vector<double> trackMasses;

  TLorentzVector dimuon, m1, m2;
  int iMuon1(-1); 
  for (unsigned int imuon1 = 0; imuon1 < muonIndices.size(); ++imuon1) {
    TrackBaseRef mu1TrackView(hTracks, muonIndices[imuon1]);
    Track tMuon1(*mu1TrackView);
    if (tMuon1.pt() < fMuonPt)  continue;

    iMuon1 = muonIndices[imuon1];
    m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON); 

    for (unsigned int itrack2 = 0; itrack2 < hTracks->size(); ++itrack2){    
      if (static_cast<int>(itrack2) == muonIndices[imuon1]) continue; 

      TrackBaseRef rTrackView(hTracks, itrack2);
      Track tTrack2(*rTrackView);
      if (tTrack2.pt() < fTrackPt)  continue;

      m2.SetPtEtaPhiM(tTrack2.pt(), tTrack2.eta(), tTrack2.phi(), MMUON); 
      dimuon = m1 + m2;

      if (dimuon.M() < fMassLow || dimuon.M() > fMassHigh) {
	if (fVerbose > 0) {
	  cout << "==>HFMuonAndTrack> dimuon mass = " << dimuon.M() << ", skipping" << endl;
	}
	continue; 
      }
      
      if (fVerbose > 0) {
	cout << "==>HFMuonAndTrack> dimuon mass = " << dimuon.M() << ", vertexing" << endl;
	cout << "==>HFMuonAndTrack>  tMuon1.charge() = " << tMuon1.charge() << ", tTrack2.charge() " << tTrack2.charge() << endl;
      }
      
      // -- vertexing
      trackList.clear();
      trackIndices.clear(); 
      trackMasses.clear(); 
      
      trackList.push_back(tMuon1); 
      trackIndices.push_back(iMuon1); 
      trackMasses.push_back(MMUON);
      
      trackList.push_back(tTrack2); 
      trackIndices.push_back(itrack2); 
      trackMasses.push_back(MMUON);
      
      if (fVertexing > 0) {
	a.doFit(trackList, trackIndices, trackMasses, fType); 	
      } else {
	a.doNotFit(trackList, trackIndices, trackMasses, fType); 	
      }
      
    } 
  }
}


// ------------ method called once each job just before starting event loop  ------------
void  HFMuonAndTrack::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFMuonAndTrack::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFMuonAndTrack);
