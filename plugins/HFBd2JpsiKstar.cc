#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBd2JpsiKstar.h"

#include <iostream>

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
HFBd2JpsiKstar::HFBd2JpsiKstar(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)), 
  fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)), 
  fKstarWindow(iConfig.getUntrackedParameter<double>("kstarWindow", 0.4)), 
  fBdWindow(iConfig.getUntrackedParameter<double>("BsWindow", 0.8)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 99.)),
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)),
  fMaxD0(iConfig.getUntrackedParameter<double>("maxD0", 999.)),
  fMaxDz(iConfig.getUntrackedParameter<double>("maxDz", 999.)),
  fVertexing(iConfig.getUntrackedParameter<int>("vertexing", 1)),
  fType(iConfig.getUntrackedParameter<int>("type", 511))  {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBd2JpsiKstar constructor" << endl;
  cout << "---  verbose:                  " << fVerbose << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  muonPt :                  " << fMuonPt << endl;
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiWindow:                " << fPsiWindow << endl;
  cout << "---  kstarWindow:              " << fKstarWindow << endl;
  cout << "---  BdWindow:                 " << fBdWindow << endl;
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
HFBd2JpsiKstar::~HFBd2JpsiKstar() {
  
}


// ----------------------------------------------------------------------
void HFBd2JpsiKstar::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFBd2JpsiKstar> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFBd2JpsiKstar> No primary vertex found, skipping" << endl;
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
  cout << "==>HFBd2JpsiKstar> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFBd2JpsiKstar> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  // -- Transient tracks for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFBd2JpsiKstar: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  // -- get the collection of muons and store their corresponding track indices
  vector<unsigned int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index(); 
    if (im >= 0) muonIndices.push_back(im);
  }
  if (fVerbose > 0) {
    cout << "==>HFBd2JpsiKstar> nMuons = " << hMuons->size() << endl;
    cout << "==>HFBd2JpsiKstar> nMuonIndices = " << muonIndices.size() << endl;
  }
  if (muonIndices.size() < static_cast<unsigned int>(fPsiMuons)) return;

  // -- Build muon lists
  TLorentzVector tlv; 
  vector<pair<int, TLorentzVector> > tlist1, tlist2,  klist, pilist; 
  if (2 == fPsiMuons) {
    tlist1.reserve(10); 
    tlist2.reserve(10); 
  } else {
    tlist1.reserve(100); 
    tlist2.reserve(100); 
  }    
  klist.reserve(100); 
  pilist.reserve(100); 
  int isMuon(0); 
  for (unsigned int itrack = 0; itrack < hTracks->size(); ++itrack){    
    TrackBaseRef rTrackView(hTracks, itrack);
    Track tTrack(*rTrackView);
    if (tTrack.pt() < fTrackPt) continue;
    if (tTrack.d0() > fMaxD0) continue;
    if (tTrack.dz() > fMaxDz) continue;
    tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MKAON); 
    klist.push_back(make_pair(itrack, tlv));
    tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MPION); 
    pilist.push_back(make_pair(itrack, tlv));

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
    if (fVerbose > 0) cout << "==>HFBs2JpsiPhi> Not enough muons found for J/psi candidates combinatorics: isMuon =" 
			   << isMuon << ", required: " << fPsiMuons 
			   << endl;
    return;
  }

  HFTwoParticleCombinatorics a(fVerbose); 
  vector<pair<int, int> > psiList; 
  a.combine(psiList, tlist1, tlist2, 2.8, 3.4, 1); 
  if (fVerbose > 0) cout << "==>HFBs2JpsiKp> J/psi list size: " << psiList.size() << endl;
  vector<pair<int, int> > kstarList; 
  a.combine(kstarList, klist, pilist, 0.5, 1.2, 0);
  if (fVerbose > 0) cout << "==>HFBs2JpsiKp> kstar list size: " << kstarList.size() << endl;
  
  HFKalmanVertexFit    aKal(fTTB.product(), fPV, 100511, fVerbose);   aKal.fMaxDoca     = fMaxDoca; 
  HFSequentialVertexFit aSeq(hTracks, fTTB.product(), fPV, fVerbose);
  vector<Track> trackList; 
  vector<int> trackIndices;
  vector<double> trackMasses;

  // -- Build J/psi + kstar
  TLorentzVector psi, kstar, m1, m2, ka1, pi2, bs;
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
    
    for (unsigned int j = 0; j < kstarList.size(); ++j){    
      unsigned int iKaon1 = kstarList[j].first; 
      unsigned int iPion2 = kstarList[j].second; 
      if (iKaon1 == iMuon1 || iKaon1 == iMuon2) continue; 
      if (iPion2 == iMuon1 || iPion2 == iMuon2) continue; 

      TrackBaseRef rTrackView1(hTracks, iKaon1);
      Track tKaon1(*rTrackView1);
      if (tKaon1.pt() < fTrackPt) continue;
      ka1.SetXYZM(tKaon1.px(), tKaon1.py(), tKaon1.pz(), MKAON); 
      if (psi.DeltaR(ka1) > fDeltaR) continue; 

      TrackBaseRef rTrackView2(hTracks, iPion2);
      Track tPion2(*rTrackView2);
      if (tPion2.pt() < fTrackPt) continue;
      pi2.SetXYZM(tPion2.px(), tPion2.py(), tPion2.pz(), MPION); 
      if (psi.DeltaR(pi2) > fDeltaR) continue; 

      kstar = ka1 + pi2; 
      if ((TMath::Abs(kstar.M() - MKSTAR) > fKstarWindow)) continue;
      
      bs = psi + kstar; 
      if (TMath::Abs(bs.M() - MBS) > fBdWindow) continue; 
      
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

      trackList.push_back(tKaon1); 
      trackIndices.push_back(iKaon1); 
      trackMasses.push_back(MKAON);

      trackList.push_back(tPion2); 
      trackIndices.push_back(iPion2); 
      trackMasses.push_back(MPION);
      
      if (0 == fVertexing) {
    	aKal.doNotFit(trackList, trackIndices, trackMasses, -100511); 	
		continue; 
      }

      aKal.doFit(trackList, trackIndices, trackMasses, 100511); 	
      aKal.doFit(trackList, trackIndices, trackMasses, 200511, 2); 	

      // -- sequential fit: J/psi kaon pion
      HFDecayTree theTree(300511);
	  
      HFDecayTreeIterator iterator = theTree.addDecayTree(300443,0); // Don't use kinematic particle for the J/Psi
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      theTree.addTrack(iKaon1,321);
      theTree.addTrack(iPion2,211);
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      aSeq.doFit(&theTree);
      
      // --sequential fit: J/Psi (constraint) kaon pion
      theTree.clear();
      theTree.particleID = 400511;
      
      iterator = theTree.addDecayTree(400443,1,MJPSI);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      theTree.addTrack(iKaon1,321);
      theTree.addTrack(iPion2,211);
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      aSeq.doFit(&theTree);
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void  HFBd2JpsiKstar::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFBd2JpsiKstar::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBd2JpsiKstar);
