#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFLambdas.h"

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
HFLambdas::HFLambdas(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fMaxTracks(iConfig.getUntrackedParameter<int>("maxTracks", 1000)), 
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fUseMuon(iConfig.getUntrackedParameter<int>("useMuon", 0)), 
  fL0Window(iConfig.getUntrackedParameter<double>("L0Window", 0.4)), 
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fProtonPt(iConfig.getUntrackedParameter<double>("protonPt", 1.0)), 
  fPionPt(iConfig.getUntrackedParameter<double>("pionPt", 1.0)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 1.5)),
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 9999.0)),
  fMaxVtxChi2(iConfig.getUntrackedParameter<double>("maxVtxChi2", 9999.0)),
  fMinVtxSigXY(iConfig.getUntrackedParameter<double>("minVtxSigXY", -1.)),
  fMinVtxSig3d(iConfig.getUntrackedParameter<double>("minVtxSig3d", -1.)),
  fMinCosAngle(iConfig.getUntrackedParameter<double>("minCosAngle", -1.)),
  fMinPtCand(iConfig.getUntrackedParameter<double>("minPtCand", -99.)),
  fType(iConfig.getUntrackedParameter<int>("type", 1)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFLambdas constructor" << endl;
  cout << "---  verbose:                  " << fVerbose << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  PrimaryVertexLabel:       " << fPrimaryVertexLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  useMuon:                  " << fUseMuon << endl;
  cout << "---  phiWindow:                " << fPhiWindow << endl;
  cout << "---  L0Window:                 " << fL0Window << endl;
  cout << "---  muonPt:                   " << fMuonPt << endl;
  cout << "---  protonPt:                 " << fProtonPt << endl;
  cout << "---  pionPt:                   " << fPionPt << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;

  cout << "---  maxDoca:                  " << fMaxDoca << endl;
  cout << "---  maxVtxChi2:               " << fMaxVtxChi2 << endl;
  cout << "---  minVtxSigXY:              " << fMinVtxSigXY << endl;
  cout << "---  minVtxSig3d:              " << fMinVtxSig3d << endl;
  cout << "---  minCosAngle:              " << fMinCosAngle << endl;
  cout << "---  minPtCand:                " << fMinPtCand << endl;

  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFLambdas::~HFLambdas() {
  
}


// ----------------------------------------------------------------------
void HFLambdas::analyze(const Event& iEvent, const EventSetup& iSetup) {

  pid_t pid = getpid();
  char line[100]; 
  sprintf(line, "ps -F %i", pid); 
  if (fVerbose > 0) {
    cout << "==>HFLambdas: beginning of analyze():" << endl;
    system(line); 
  }

  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFLambdas> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFLambdas> No primary vertex found, skipping" << endl;
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
  cout << "==>HFLambdas> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFLambdas> No valid TrackCollection with label " << fTracksLabel << " found, skipping" << endl;
    return;
  }

  if (hTracks->size() > static_cast<unsigned int>(fMaxTracks)) {
    cout << "==>HFLambdas> Too many tracks " << hTracks->size() << ", skipping" << endl;
    return;
  }
  if (fVerbose > 0) {
    cout << "==>HFLambdas> ntracks = " << hTracks->size() << endl;
  }
  
  // -- Transient tracks for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFLambdas: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  // -- get the collection of muons and store their corresponding track indices
  vector<unsigned int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index(); 
    if (im >= 0) muonIndices.push_back(im);
  }
  if (fVerbose > 0) {
    cout << "==>HFLambdas> nMuons = " << hMuons->size() << endl;
    cout << "==>HFLambdas> nMuonIndices = " << muonIndices.size() << endl;
  }

  // -- Build lists
  TLorentzVector muon, lambda0, proton, pion, track;
  TLorentzVector tlv; 
  vector<pair<int, TLorentzVector> > prlist, pilist, mulist; 
  mulist.reserve(200); 
  pilist.reserve(2000); 
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

  /*
  HFTwoParticleCombinatorics a(fVerbose); 
  vector<pair<int, int> > piList; 
  kapiList.reserve(100000); 
  a.combine(kapiList, kalist, pilist, 0.5, 2.5, 0); 

  vector<pair<int, int> > phiList; 
  a.combine(phiList, kalist, kalist, 0.9, 2.0, 1); 

  if (fVerbose > 0) cout << "==>HFLambdas> K-pi list size: " << kapiList.size() << endl;
  if (fVerbose > 0) cout << "==>HFLambdas> KK list size: " << phiList.size() << endl;
  */
  
  HFKalmanVertexFit  aKal(fTTB.product(), fPV, 0, fVerbose); 
  vector<Track> trackList; 
  vector<int> trackIndices;
  vector<double> trackMasses;

  // -------------------
  // -- KVF: L0 -> p pi-
  // -------------------

  typedef pair<int,int> duplet;
  vector<duplet> prpiList;

  HFTwoParticleCombinatorics b(fVerbose);
  b.combine(prpiList, prlist, pilist, MLAMBDA_0-fL0Window, MLAMBDA_0+fL0Window);
  if (fVerbose > 0) cout << "==>HFLambdas> pr-pi list size: " << prpiList.size() << endl;

  aKal.setNoCuts();
  aKal.fMaxDoca     = fMaxDoca; 
  aKal.fVtxChi2     = fMaxVtxChi2; 
  aKal.fVtxSigXY    = fMinVtxSigXY; 
  aKal.fVtxSig3d    = fMinVtxSig3d; 
  aKal.fCosAngle    = fMinCosAngle;
  aKal.fPtCand      = fMinPtCand;

  for (vector<duplet>::iterator it=prpiList.begin(); it!=prpiList.end(); ++it) {
    
    TrackBaseRef prTrackView(hTracks, it->first);
    Track tProton(*prTrackView);
    proton.SetPtEtaPhiM(tProton.pt(), tProton.eta(), tProton.phi(), MPROTON); 

    TrackBaseRef piTrackView(hTracks, it->second);
    Track tPion(*piTrackView);
    pion.SetPtEtaPhiM(tPion.pt(), tPion.eta(), tPion.phi(), MPION); 

    if (tPion.charge() == tProton.charge()) continue;          // pions have opposite charge from proton

    /*
    if (fUseMuon) {
      if (static_cast<unsigned int>(muMaxIdx) == it->first) continue; 
      if (static_cast<unsigned int>(muMaxIdx) == it->second) continue; 
      //if (tProton.charge()*muMaxQ < 0) continue;                 // muon has same charge as Kaon
      if (muon.DeltaR(proton) > fDeltaR) continue; 
      if (muon.DeltaR(pion) > fDeltaR) continue; 
    } */
    

      
    trackList.clear();
    trackIndices.clear(); 
    trackMasses.clear(); 
      
    trackList.push_back(tProton); 
    trackIndices.push_back(it->first); 
    trackMasses.push_back(MPROTON);
      
    trackList.push_back(tPion); 
    trackIndices.push_back(it->second); 
    trackMasses.push_back(MPION);
      
    /*
    if (fUseMuon) {
	TrackBaseRef muTrackView(hTracks, muMaxIdx);
	Track tMuon(*muTrackView);
	trackList.push_back(tMuon); 
	trackIndices.push_back(muMaxIdx); 
	trackMasses.push_back(-MMUON);
    } */

    // -- D+, with fitting of all three tracks
    lambda0 = proton + pion; 
    cout << "lambda0.M(): " << lambda0.M() << endl;
    if ((TMath::Abs(lambda0.M() - MLAMBDA_0) < fL0Window )) {
      cout << "if true" << endl;
      aKal.doFit(trackList, trackIndices, trackMasses, fType*10000+3122, 2);
    }
  }
  prpiList.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void  HFLambdas::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFLambdas::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFLambdas);
