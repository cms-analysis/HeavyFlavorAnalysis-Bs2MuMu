#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFCharm.h"

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
HFCharm::HFCharm(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fMaxTracks(iConfig.getUntrackedParameter<int>("maxTracks", 150)), 
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fUseMuon(iConfig.getUntrackedParameter<int>("useMuon", 0)), 
  fPhiWindow(iConfig.getUntrackedParameter<double>("phiWindow", 0.3)), 
  fDWindow(iConfig.getUntrackedParameter<double>("DWindow", 0.3)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)), 
  fKaonPt(iConfig.getUntrackedParameter<double>("kaonPt", 1.0)), 
  fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 1.0)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 1.5)),
  fType(iConfig.getUntrackedParameter<int>("type", 1)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFCharm constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  PrimaryVertexLabel:       " << fPrimaryVertexLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  muonPt:                   " << fMuonPt << endl;
  cout << "---  useMuon:                  " << fUseMuon << endl;
  cout << "---  phiWindow:                " << fPhiWindow << endl;
  cout << "---  DWindow:                  " << fDWindow << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  kaonPt:                   " << fKaonPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFCharm::~HFCharm() {
  
}


// ----------------------------------------------------------------------
void HFCharm::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFCharm> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFCharm> No primary vertex found, skipping" << endl;
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
  cout << "==>HFCharm> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFCharm> No valid TrackCollection with label " << fTracksLabel << " found, skipping" << endl;
    return;
  }

  if (hTracks->size() > fMaxTracks) {
    cout << "==>HFCharm> Too many tracks " << hTracks->size() << ", skipping" << endl;
    return;
  }

  // -- Transient tracks for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFCharm: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  // -- get the collection of muons and store their corresponding track indices
  vector<unsigned int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index(); 
    if (im > 0) muonIndices.push_back(im);
  }
  if (fVerbose > 0) {
    cout << "==>HFCharm> nMuons = " << hMuons->size() << endl;
    cout << "==>HFCharm> nMuonIndices = " << muonIndices.size() << endl;
  }

  // -- Build lists
  TLorentzVector muon, dzero, dplus, kaon, pion, track, kaon1, kaon2;
  TLorentzVector tlv; 
  vector<pair<int, TLorentzVector> > kalist, pilist, mulist; 
  mulist.reserve(20); 
  pilist.reserve(200); 
  kalist.reserve(200); 

  int    muMaxIdx(-1), muMaxQ(0); 
  double muMaxPt(-99.), muPt(0.); 
  for (unsigned int itrack = 0; itrack < hTracks->size(); ++itrack){    
    TrackBaseRef rTrackView(hTracks, itrack);
    Track tTrack(*rTrackView);
    if (tTrack.pt() < fPionPt)  continue;
    tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MPION); 
    pilist.push_back(make_pair(itrack, tlv));
    if (tTrack.pt() < fKaonPt)  continue;
    tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MKAON); 
    kalist.push_back(make_pair(itrack, tlv));
    for (unsigned int im = 0; im < muonIndices.size(); ++im) {
      if (muonIndices[im] == itrack) {
	muPt = tTrack.pt();
	if ((muPt > muMaxPt) && (muPt > fMuonPt)) {
	  muMaxPt  = muPt;
	  muMaxQ   = tTrack.charge();
	  muMaxIdx = itrack; 
	  muon.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MMUON); 
	  break;
	}
      }
    }
  }

  // -- skip rest if we want muons, but none satisfy the muon PT cut
  if (fUseMuon && muMaxIdx < 0) return; 

  HFTwoParticleCombinatorics a(fVerbose); 
  vector<pair<int, int> > kapiList; 
  a.combine(kapiList, kalist, pilist, 0.6, 2.1, 0); 
  if (fVerbose > 0) cout << "==>HFCharm> K-pi list size: " << kapiList.size() << endl;
  
  HFKalmanVertexFit  aKal(fTTB.product(), fPV, 0, fVerbose); 
  vector<Track> trackList; 
  vector<int> trackIndices;
  vector<double> trackMasses;

  // -------------------------
  // -- Build (K pi) (+ stuff)
  // -------------------------
  for (unsigned int i = 0; i < kapiList.size(); ++i) {
    unsigned int iKaon = kapiList[i].first; 
    unsigned int iPion = kapiList[i].second; 
    
    TrackBaseRef kaTrackView(hTracks, iKaon);
    Track tKaon(*kaTrackView);
    kaon.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MKAON); 

    TrackBaseRef piTrackView(hTracks, iPion);
    Track tPion(*piTrackView);
    pion.SetPtEtaPhiM(tPion.pt(), tPion.eta(), tPion.phi(), MPION); 

    if (fUseMuon) {
      if (tKaon.charge()*muMaxQ < 0) continue; 
      if (muon.DeltaR(kaon) > fDeltaR) continue; 
      if (muon.DeltaR(pion) > fDeltaR) continue; 
    }
    
    // ----------
    // -- KVF: D0
    // ----------
    dzero = kaon + pion; 
    if ((TMath::Abs(dzero.M() - MD0) < fDWindow)) {
      trackList.clear();
      trackIndices.clear(); 
      trackMasses.clear(); 
      
      trackList.push_back(tKaon); 
      trackIndices.push_back(iKaon); 
      trackMasses.push_back(MKAON);
      
      trackList.push_back(tPion); 
      trackIndices.push_back(iPion); 
      trackMasses.push_back(MPION);
      
      if (fUseMuon) {
	TrackBaseRef muTrackView(hTracks, muMaxIdx);
	Track tMuon(*muTrackView);
	trackList.push_back(tMuon); 
	trackIndices.push_back(muMaxIdx); 
	trackMasses.push_back(MMUON);
      }

      aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+10, 2); 	

    }

    // ---------------
    // -- KVF: K pi pi
    // ---------------
    for (unsigned int iTrack = 0; iTrack < hTracks->size(); ++iTrack){    
      if (iTrack == iKaon) continue; 
      if (iTrack == iPion) continue; 
      TrackBaseRef rTrackView(hTracks, iTrack);
      Track tTrack(*rTrackView);
      if (tTrack.charge() == tKaon.charge()) continue; 

      track.SetPtEtaPhiM(tTrack.pt(), tTrack.eta(), tTrack.phi(), MPION); 
      if (fUseMuon) {
	if (muon.DeltaR(track) > fDeltaR) continue; 
      }

      trackList.clear();
      trackIndices.clear(); 
      trackMasses.clear(); 
      
      trackList.push_back(tKaon); 
      trackIndices.push_back(iKaon); 
      trackMasses.push_back(MKAON);
      
      trackList.push_back(tPion); 
      trackIndices.push_back(iPion); 
      trackMasses.push_back(MPION);
      
      trackList.push_back(tTrack); 
      trackIndices.push_back(iTrack); 
      trackMasses.push_back(MPION);

      if (fUseMuon) {
	TrackBaseRef muTrackView(hTracks, muMaxIdx);
	Track tMuon(*muTrackView);
	trackList.push_back(tMuon); 
	trackIndices.push_back(muMaxIdx); 
	trackMasses.push_back(MMUON);
      }

      // -- D*, with (K pi) in D0 mass window, fitting only D0
      if ((TMath::Abs(dzero.M() - MD0) < 0.1)) {
	aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+20, 2); 	
      }
      
      // -- D+, with fitting of all three tracks
      if (tTrack.pt() < fTrackPt)  continue;
      dplus = kaon + pion + track; 
      if ((TMath::Abs(dplus.M() - MDPLUS) < fDWindow)) {
	aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+30, 3);
      }
    }
  }


  // ---------------------
  // -- Build (K K) + pion
  // ---------------------
  vector<pair<int, int> > phiList; 
  a.combine(phiList, kalist, kalist, 0.9, 1.15, 1); 
  if (fVerbose > 0) cout << "==>HFCharm> KK list size: " << phiList.size() << endl;
  
  for (unsigned int i = 0; i < phiList.size(); ++i) {
    unsigned int iKaon1 = phiList[i].first; 
    unsigned int iKaon2 = phiList[i].second; 
    
    TrackBaseRef ka1TrackView(hTracks, iKaon1);
    Track tKaon1(*ka1TrackView);
    kaon1.SetPtEtaPhiM(tKaon1.pt(), tKaon1.eta(), tKaon1.phi(), MKAON); 

    TrackBaseRef ka2TrackView(hTracks, iKaon2);
    Track tKaon2(*ka2TrackView);
    kaon2.SetPtEtaPhiM(tKaon2.pt(), tKaon2.eta(), tKaon2.phi(), MKAON); 

    if (fUseMuon && muMaxIdx > -1) {
      if (muon.DeltaR(kaon1) > fDeltaR) continue; 
      if (muon.DeltaR(kaon2) > fDeltaR) continue; 
    }


    for (unsigned int iTrack = 0; iTrack < hTracks->size(); ++iTrack){    
      if (iTrack == iKaon1) continue; 
      if (iTrack == iKaon2) continue; 
      TrackBaseRef rTrackView(hTracks, iTrack);
      Track tTrack(*rTrackView);
      if (tTrack.pt() < fTrackPt) continue; 

      track.SetPtEtaPhiM(tTrack.pt(), tTrack.eta(), tTrack.phi(), MPION); 
      if (fUseMuon) {
	if (muon.DeltaR(track) > fDeltaR) continue; 
      }

      trackList.clear();
      trackIndices.clear(); 
      trackMasses.clear(); 
      
      trackList.push_back(tKaon1); 
      trackIndices.push_back(iKaon1); 
      trackMasses.push_back(MKAON);
      
      trackList.push_back(tKaon2); 
      trackIndices.push_back(iKaon2); 
      trackMasses.push_back(MKAON);
      
      trackList.push_back(tTrack); 
      trackIndices.push_back(iTrack); 
      trackMasses.push_back(MPION);

      if (fUseMuon) {
	TrackBaseRef muTrackView(hTracks, muMaxIdx);
	Track tMuon(*muTrackView);
	trackList.push_back(tMuon); 
	trackIndices.push_back(muMaxIdx); 
	trackMasses.push_back(MMUON);
      }

      dplus = kaon1 + kaon2 + track; 

      if (dplus.M() > 1.5 && dplus.M() < 2.5) {
	aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+40, 3);
      }
    }
  }
  
}


// ------------ method called once each job just before starting event loop  ------------
void  HFCharm::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFCharm::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFCharm);
