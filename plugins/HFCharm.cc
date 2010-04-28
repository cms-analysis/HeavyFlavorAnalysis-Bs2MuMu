#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFCharm.h"

#include <iostream>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKinematicVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTwoParticleCombinatorics.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFThreeParticleCombinatorics.hh"
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
  fMaxTracks(iConfig.getUntrackedParameter<int>("maxTracks", 1000)), 
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fUseMuon(iConfig.getUntrackedParameter<int>("useMuon", 0)), 
  fPhiWindow(iConfig.getUntrackedParameter<double>("phiWindow", 0.3)), 
  fDWindow(iConfig.getUntrackedParameter<double>("DWindow", 0.3)), 
  fLcWindow(iConfig.getUntrackedParameter<double>("LcWindow", 0.4)), 
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fProtonPt(iConfig.getUntrackedParameter<double>("protonPt", 1.0)), 
  fKaonPt(iConfig.getUntrackedParameter<double>("kaonPt", 1.0)), 
  fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 1.0)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 1.5)),
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 2.0)),
  fType(iConfig.getUntrackedParameter<int>("type", 1)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFCharm constructor" << endl;
  cout << "---  verbose:                  " << fVerbose << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  PrimaryVertexLabel:       " << fPrimaryVertexLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  useMuon:                  " << fUseMuon << endl;
  cout << "---  phiWindow:                " << fPhiWindow << endl;
  cout << "---  DWindow:                  " << fDWindow << endl;
  cout << "---  muonPt:                   " << fMuonPt << endl;
  cout << "---  protonPt:                 " << fProtonPt << endl;
  cout << "---  kaonPt:                   " << fKaonPt << endl;
  cout << "---  pionPt:                   " << fPionPt << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "---  maxDoca:                  " << fMaxDoca << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFCharm::~HFCharm() {
  
}


// ----------------------------------------------------------------------
void HFCharm::analyze(const Event& iEvent, const EventSetup& iSetup) {

  pid_t pid = getpid();
  char line[100]; 
  sprintf(line, "ps -F %i", pid); 
  if (fVerbose > 0) {
    cout << "==>HFCharm: beginning of analyze():" << endl;
    system(line); 
  }

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
  fPV = vertices[gHFEvent->fBestPV]; 
  if (fVerbose > 0) {
    cout << "HFDimuons: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << endl;
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

  if (hTracks->size() > static_cast<unsigned int>(fMaxTracks)) {
    cout << "==>HFCharm> Too many tracks " << hTracks->size() << ", skipping" << endl;
    return;
  }
  if (fVerbose > 0) {
    cout << "==>HFCharm> ntracks = " << hTracks->size() << endl;
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
  TLorentzVector muon, lambdac, dzero, dplus, kaon, pion, track, kaon1, kaon2, pion1, pion2, pion3;
  TLorentzVector tlv; 
  vector<pair<int, TLorentzVector> > prlist, kalist, pilist, mulist; 
  mulist.reserve(200); 
  pilist.reserve(2000); 
  kalist.reserve(2000); 
  prlist.reserve(2000); 

  int    muMaxIdx(-1), muMaxQ(0); 
  double muMaxPt(-99.), muPt(0.); 
  for (unsigned int itrack = 0; itrack < hTracks->size(); ++itrack){    
    TrackBaseRef rTrackView(hTracks, itrack);
    Track tTrack(*rTrackView);
    if (tTrack.pt() > fPionPt)  {
      tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MPION); 
      pilist.push_back(make_pair(itrack, tlv));
    }

    if (tTrack.pt() > fKaonPt) {
      tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MKAON); 
      kalist.push_back(make_pair(itrack, tlv));
    }

    if (tTrack.pt() > fProtonPt) {
      tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MPROTON); 
      prlist.push_back(make_pair(itrack, tlv));
    }

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
  kapiList.reserve(100000); 
  a.combine(kapiList, kalist, pilist, 0.5, 2.5, 0); 

  vector<pair<int, int> > phiList; 
  a.combine(phiList, kalist, kalist, 0.9, 2.0, 1); 

  if (fVerbose > 0) cout << "==>HFCharm> K-pi list size: " << kapiList.size() << endl;
  if (fVerbose > 0) cout << "==>HFCharm> KK list size: " << phiList.size() << endl;
  
  HFKalmanVertexFit  aKal(fTTB.product(), fPV, 0, fVerbose); 
  aKal.fMaxDoca = fMaxDoca; 
  vector<Track> trackList; 
  vector<int> trackIndices;
  vector<double> trackMasses;

  // ----------------------------
  // -- D0 -> K- Pi+   (with mu-)
  // ----------------------------
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
	trackMasses.push_back(-MMUON); // Note negative mass: fVar1 contains D0 mass WITHOUT muon 
      }

      aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+10, 2); 
    }	
  }


  // ----------------------------------------------
  // -- LambdaC+ -> proton+ kaon- pion+  (with mu-)
  // ----------------------------------------------
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
    
    
    for (unsigned int iTrack = 0; iTrack < hTracks->size(); ++iTrack){    
      if (iTrack == iKaon) continue; 
      if (iTrack == iPion) continue; 
      
      TrackBaseRef rTrackView(hTracks, iTrack);
      Track tTrack(*rTrackView);
      if (tTrack.charge() == tKaon.charge()) continue; 
      
      track.SetPtEtaPhiM(tTrack.pt(), tTrack.eta(), tTrack.phi(), MPROTON); 
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
      trackMasses.push_back(MPROTON);

      lambdac = track + kaon + pion; 
      if (fUseMuon && (TMath::Abs(lambdac.M() - MLAMBDA_C) < fLcWindow)) {
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
	trackMasses.push_back(MPROTON);
	
	TrackBaseRef muTrackView(hTracks, muMaxIdx);
	Track tMuon(*muTrackView);
	trackList.push_back(tMuon); 
	trackIndices.push_back(muMaxIdx); 
	trackMasses.push_back(-MMUON);
      }

      if ((TMath::Abs(lambdac.M() - MLAMBDA_C) < fLcWindow)) {
	aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+50, 3); 	
      }    
    }
  }

  // ---------------------------------
  // -- D*+ -> K- pi+ pi+   (with mu-)
  // ---------------------------------
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
	trackMasses.push_back(-MMUON);
      }

      // -- D*, with (K pi) in D0 mass window, fitting only D0
      dzero = kaon + pion; 
      if ((TMath::Abs(dzero.M() - MD0) < 0.2)) {
	aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+20, 2); 	
      }
      
      //       // -- D+, with fitting of all three tracks: This is done WITHOUT duplicates by Hadi's triplets
      //       dplus = kaon + pion + track; 
      //       if ((TMath::Abs(dplus.M() - MDPLUS) < fDWindow)) {
      // 	//	aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+60, 3);
      //       }


    }
  }


  // -------------------
  // -- KVF: D+ (K pi pi)
  // -------------------

  vector<triplet> kapipiList;

  HFThreeParticleCombinatorics b(fVerbose);
  b.combine(kapipiList, kalist, pilist,1.0, 2.3);
  if (fVerbose > 0) cout << "==>HFCharm> K-pi-pi list size: " << kapipiList.size() << endl;

  for (vector<triplet>::iterator it=kapipiList.begin(); it!=kapipiList.end(); ++it) {
    
    TrackBaseRef kaTrackView(hTracks, it->ka());
    Track tKaon(*kaTrackView);
    kaon.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MKAON); 

    TrackBaseRef pi1TrackView(hTracks, it->pi1());
    Track tPion1(*pi1TrackView);
    pion1.SetPtEtaPhiM(tPion1.pt(), tPion1.eta(), tPion1.phi(), MPION); 

    TrackBaseRef pi2TrackView(hTracks, it->pi2());
    Track tPion2(*pi2TrackView);
    pion2.SetPtEtaPhiM(tPion2.pt(), tPion2.eta(), tPion2.phi(), MPION); 

    if (fUseMuon) {
      if (tKaon.charge()*muMaxQ < 0) continue;                 // muon has same charge as Kaon
      if (muon.DeltaR(kaon) > fDeltaR) continue; 
      if (muon.DeltaR(pion1) > fDeltaR) continue; 
      if (muon.DeltaR(pion2) > fDeltaR) continue; 
    }
    

    if (tPion1.charge() == tKaon.charge()) continue;          // pions have opposite charge from Kaon
    if (tPion2.charge() == tKaon.charge()) continue; 
      
    trackList.clear();
    trackIndices.clear(); 
    trackMasses.clear(); 
      
    trackList.push_back(tKaon); 
    trackIndices.push_back(it->ka()); 
    trackMasses.push_back(MKAON);
      
    trackList.push_back(tPion1); 
    trackIndices.push_back(it->pi1()); 
    trackMasses.push_back(MPION);
      
    trackList.push_back(tPion2); 
    trackIndices.push_back(it->pi2()); 
    trackMasses.push_back(MPION);

    if (fUseMuon) {
	TrackBaseRef muTrackView(hTracks, muMaxIdx);
	Track tMuon(*muTrackView);
	trackList.push_back(tMuon); 
	trackIndices.push_back(muMaxIdx); 
	trackMasses.push_back(-MMUON);
    }

    // -- D+, with fitting of all three tracks
    dplus = kaon + pion1 + pion2; 
    if ((TMath::Abs(dplus.M() - MDPLUS) < fDWindow)) {
      aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+30, 3);
    }
  }
  kapipiList.clear();


  // -------------------------------
  // -- KVF: D+ (Kshort pi --> 3 pi)
  // -------------------------------

  vector<triplet> KshortpiList; 
  b.combine(KshortpiList, pilist, 0.8, 2.3, MKSHORT-0.25, MKSHORT+0.25); 
  if (fVerbose > 0) cout << "==>HFCharm> Kshortpi list size: " << KshortpiList.size() << endl;

  for (vector<triplet>::iterator it=KshortpiList.begin(); it!=KshortpiList.end(); ++it) {
    
    TrackBaseRef pi1TrackView(hTracks, it->pi1());
    Track tPion1(*pi1TrackView);
    pion1.SetPtEtaPhiM(tPion1.pt(), tPion1.eta(), tPion1.phi(), MPION); 

    TrackBaseRef pi2TrackView(hTracks, it->pi2());
    Track tPion2(*pi2TrackView);
    pion2.SetPtEtaPhiM(tPion2.pt(), tPion2.eta(), tPion2.phi(), MPION); 

    TrackBaseRef pi3TrackView(hTracks, it->pi3());
    Track tPion3(*pi3TrackView);
    pion3.SetPtEtaPhiM(tPion3.pt(), tPion3.eta(), tPion3.phi(), MPION); 

    if (fUseMuon) {
      if (tPion3.charge()*muMaxQ > 0) continue;                     // pion coming from D+ vertex has opposite sign from muon
      //      if (muon.DeltaR(pion1) > fDeltaR) continue; 
      //      if (muon.DeltaR(pion2) > fDeltaR) continue; 
      if (muon.DeltaR(pion3) > fDeltaR) continue; 
    }
    
    if (tPion1.charge() == tPion2.charge()) continue;               // pions from Kshort must have opposite charge
      
    trackList.clear();
    trackIndices.clear(); 
    trackMasses.clear(); 
      
    trackList.push_back(tPion1); 
    trackIndices.push_back(it->pi1()); 
    trackMasses.push_back(MPION);
      
    trackList.push_back(tPion2); 
    trackIndices.push_back(it->pi2()); 
    trackMasses.push_back(MPION);
      
    trackList.push_back(tPion3); 
    trackIndices.push_back(it->pi3()); 
    trackMasses.push_back(MPION);

    if (fUseMuon) {
	TrackBaseRef muTrackView(hTracks, muMaxIdx);
	Track tMuon(*muTrackView);
	trackList.push_back(tMuon); 
	trackIndices.push_back(muMaxIdx); 
	trackMasses.push_back(-MMUON);
    }

    // -- Kshort, with fitting of two tracks only
    TLorentzVector kshort = pion1 + pion2; 
    if ((TMath::Abs(kshort.M() - MKSHORT) < fDWindow)) {
	aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+60, 2);
    }
  }
  KshortpiList.clear();


  // ----------------------------------
  // -- D(s)+ -> K+ K- pi+   (with mu-)
  // ----------------------------------

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
	if (tTrack.charge()*muMaxQ > 0) continue; 
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
	trackMasses.push_back(-MMUON);
      }

      dplus = kaon1 + kaon2 + track; 

      // unified very wide window for both D+ and Ds+
      if (dplus.M() > 1.5 && dplus.M() < 2.5) {
	aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+40, 3);
      }
    }
  }

  kapiList.clear();
  kapipiList.clear(); 
  phiList.clear();

  trackList.clear();
  trackIndices.clear();
  trackMasses.clear();
  
  if (fVerbose > 0) {
    cout << "==>HFCharm: end of analyze():" << endl;
    system(line); 
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
