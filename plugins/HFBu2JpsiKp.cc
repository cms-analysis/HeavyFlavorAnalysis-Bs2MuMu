#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBu2JpsiKp.h"

#include <iostream>
#include <utility>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKinematicVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatorics.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDecayTree.h"

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
HFBu2JpsiKp::HFBu2JpsiKp(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)), 
  fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)), 
  fBuWindow(iConfig.getUntrackedParameter<double>("BuWindow", 0.8)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 1.5)), 
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)), 
  fMaxD0(iConfig.getUntrackedParameter<double>("maxD0", 999.)),
  fMaxDz(iConfig.getUntrackedParameter<double>("maxDz", 999.)),
  fVertexing(iConfig.getUntrackedParameter<int>("vertexing", 1)),
  fType(iConfig.getUntrackedParameter<int>("type", 521))  {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBu2JpsiKp constructor" << endl;
  cout << "---  verbose:                  " << fVerbose << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  PrimaryVertexLabel:       " << fPrimaryVertexLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  muonPt :                  " << fMuonPt << endl;
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiWindow:                " << fPsiWindow << endl;
  cout << "---  BuWindow:                 " << fBuWindow << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "---  maxDoca:                  " << fMaxDoca << endl;
  cout << "---  maxD0:                    " << fMaxD0 << endl;
  cout << "---  maxDz:                    " << fMaxDz << endl;
  cout << "---  vertexing:                " << fVertexing << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFBu2JpsiKp::~HFBu2JpsiKp() {
  
}


// ----------------------------------------------------------------------
void HFBu2JpsiKp::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFBu2JpsiKp> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFBu2JpsiKp> No primary vertex found, skipping" << endl;
    return;
  }
  fPV = vertices[gHFEvent->fBestPV]; 
  if (fVerbose > 0) {
    cout << "HFDimuons: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << endl;
  }
  
  // -- get the collection of muons
  Handle<MuonCollection> hMuons;
  iEvent.getByLabel(fMuonsLabel, hMuons);
  if (!hMuons.isValid()) {
  cout << "==>HFBu2JpsiKp> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
    return;
  }

  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFBu2JpsiKp> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  // -- Transient track builder for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFBu2JpsiKp: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  // -- get the collection of muons and store their corresponding track indices
  vector<unsigned int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index(); 
    if (im >= 0) muonIndices.push_back(im);
  }
  if (fVerbose > 0) {
    cout << "==>HFBu2JpsiKp> nMuons = " << hMuons->size() << endl;
    cout << "==>HFBu2JpsiKp> nMuonIndices = " << muonIndices.size() << endl;
  }
  if (muonIndices.size() < static_cast<unsigned int>(fPsiMuons)) return;


  // -- Build muon lists
  TLorentzVector tlv; 
  vector<pair<int, TLorentzVector> > tlist1, tlist2; 
  if (2 == fPsiMuons) {
    tlist1.reserve(10); 
    tlist2.reserve(10); 
  } else {
    tlist1.reserve(100); 
    tlist2.reserve(100); 
  }    
  int isMuon(0); 
  for (unsigned int itrack = 0; itrack < hTracks->size(); ++itrack){    
    TrackBaseRef rTrackView(hTracks, itrack);
    Track tTrack(*rTrackView);
    if (tTrack.d0() > fMaxD0) continue;
    if (tTrack.dz() > fMaxDz) continue;
    tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MMUON); 
    if (2 == fPsiMuons) {
      for (unsigned int im = 0; im < muonIndices.size(); ++im) {
		  if (muonIndices[im] == itrack) {
			  ++isMuon; 
			  tlist1.push_back(make_pair(itrack, tlv)); 
			  tlist2.push_back(make_pair(itrack, tlv)); 
		  }
      } 
    } else if (1 == fPsiMuons) {
      for (unsigned int im = 0; im < muonIndices.size(); ++im) {
		  if (muonIndices[im] == itrack) {
			  ++isMuon;
			  tlist1.push_back(make_pair(itrack, tlv)); 
		  }
      }
      tlist2.push_back(make_pair(itrack, tlv)); 
    } else {
      tlist1.push_back(make_pair(itrack, tlv)); 
      tlist2.push_back(make_pair(itrack, tlv)); 
    }
  }

  if (isMuon < fPsiMuons) {
    if (fVerbose > 0) cout << "==>HFBu2JpsiKp> Not enough muons found for J/psi candidates combinatorics: isMuon =" 
			   << isMuon << ", required: " << fPsiMuons 
			   << endl;
    return;
  }

  HFTwoParticleCombinatorics a(fVerbose); 
  vector<pair<int, int> > psiList; 
  a.combine(psiList, tlist1, tlist2, 2.8, 3.4, 1); 
  if (fVerbose > 0) cout << "==>HFBu2JpsiKp> J/psi list size: " << psiList.size() << endl;
  
  HFKalmanVertexFit     aKal(fTTB.product(), fPV, 100521, fVerbose);  aKal.fMaxDoca     = fMaxDoca; 
  HFKinematicVertexFit  aKin(fTTB.product(), fPV, 300521, fVerbose);  
  HFSequentialVertexFit aSeq(hTracks, fTTB.product(), fPV, 0 /*fVerbose*/);
  vector<Track> trackList; 
  vector<int> trackIndices;
  vector<double> trackMasses;

  // -- Build J/psi + track
  TLorentzVector psi, cpsi, m1, m2, ka, bu;
  for (unsigned int i = 0; i < psiList.size(); ++i) {
    unsigned int iMuon1 = psiList[i].first; 
    unsigned int iMuon2 = psiList[i].second; 
    
    TrackBaseRef mu1TrackView(hTracks, iMuon1);
    Track tMuon1(*mu1TrackView);
    if (tMuon1.pt() < fMuonPt)  continue;
    m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON); 

    TrackBaseRef mu2TrackView(hTracks, iMuon2);
    Track tMuon2(*mu2TrackView);
    if (tMuon2.pt() < fMuonPt)  continue;
    m2.SetPtEtaPhiM(tMuon2.pt(), tMuon2.eta(), tMuon2.phi(), MMUON); 

    psi = m1 + m2; 
    if ((TMath::Abs(psi.M() - MJPSI) > fPsiWindow)) continue;
    
    for (unsigned int iTrack = 0; iTrack < hTracks->size(); ++iTrack){    
      if (iTrack == iMuon1 || iTrack == iMuon2) continue; 
      TrackBaseRef rTrackView(hTracks, iTrack);
      Track tKaon(*rTrackView);
      if (tKaon.d0() > fMaxD0) continue;
      if (tKaon.dz() > fMaxDz) continue;
      if (tKaon.pt() < fTrackPt) continue;
      ka.SetXYZM(tKaon.px(), tKaon.py(), tKaon.pz(), MKAON); 
      if (psi.DeltaR(ka) > fDeltaR) continue; 
      
      bu = ka + psi; 
      if (TMath::Abs(bu.M() - MBPLUS) > fBuWindow) continue; 

      // -- KVF: muon muon kaon
      trackList.clear();
      trackIndices.clear(); 
      trackMasses.clear(); 
      
      trackList.push_back(tMuon1); 
      trackIndices.push_back(iMuon1); 
      trackMasses.push_back(MMUON);
      
      trackList.push_back(tMuon2); 
      trackIndices.push_back(iMuon2); 
      trackMasses.push_back(MMUON);

      trackList.push_back(tKaon); 
      trackIndices.push_back(iTrack); 
      trackMasses.push_back(MKAON);
      
      if (0 == fVertexing) {
		  aKal.doNotFit(trackList, trackIndices, trackMasses, -100521); 	
		  continue; 
      }

      if (fVerbose > 1) {
	cout << "psiList/track: " << i << "/" << iTrack << endl;
      }
      aKal.doFit(trackList, trackIndices, trackMasses, 100521); 	
      aKal.doFit(trackList, trackIndices, trackMasses, 200521, 2); 	

      // -- kinematic fit: J/psi kaon
      trackList.clear();
      trackIndices.clear(); 
      trackMasses.clear(); 
      
      trackList.push_back(tMuon1); 
      trackIndices.push_back(iMuon1); 
      trackMasses.push_back(MMUON);
      
      trackList.push_back(tMuon2); 
      trackIndices.push_back(iMuon2); 
      trackMasses.push_back(MMUON);
      
      trackList.push_back(tKaon); 
      trackIndices.push_back(iTrack); 
      trackMasses.push_back(MKAON);
      
      aKin.doJpsiFit(trackList, trackIndices, trackMasses, 300521); 	
      
      // Sequentialfit without mass constraint
      if (fVerbose > 5) cout << "==>HFBu2JpsiKp> going to sequential fit" << endl;
      HFDecayTree theTree(500521);
      theTree.addTrack(iTrack,321);
    
      HFDecayTreeIterator iterator = theTree.addDecayTree(500443);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);

      if (fVerbose > 5) cout << "==>HFBu2JpsiKp> sequential fit without mass constraint" << endl;
      aSeq.doFit(&theTree);

      // Sequentialfit with mass constraint
      theTree.clear();
      theTree.particleID = 700521;
      theTree.addTrack(iTrack,321);
      iterator = theTree.addDecayTree(700443,MJPSI);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      if (fVerbose > 5) cout << "==>HFBu2JpsiKp> sequential fit with mass constraint" << endl;

      aSeq.doFit(&theTree);
      if (fVerbose > 5) cout << "==>HFBu2JpsiKp> done with fitting for track " << iTrack << endl;
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void  HFBu2JpsiKp::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFBu2JpsiKp::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBu2JpsiKp);
