#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDstar.h"

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

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

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
HFDstar::HFDstar(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("generalTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 3.0)), 
  fSlowPionPt(iConfig.getUntrackedParameter<double>("slowPionPt", 0.4)), 
  fD0Window(iConfig.getUntrackedParameter<double>("D0Window", 0.1)), 
  fDeltaM(iConfig.getUntrackedParameter<double>("deltaM", 0.03)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 99.0)), 
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)), 
  fMaxD0(iConfig.getUntrackedParameter<double>("maxD0", 999.)),
  fMaxDz(iConfig.getUntrackedParameter<double>("maxDz", 999.)),
  fType(iConfig.getUntrackedParameter<int>("type", 300000))  {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDstar constructor" << endl;
  cout << "---  verbose:                  " << fVerbose << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  PrimaryVertexLabel:       " << fPrimaryVertexLabel << endl;
  cout << "---  D0Window:                 " << fD0Window << endl;
  cout << "---  DeltaM:                   " << fDeltaM << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  slowPionPt:               " << fSlowPionPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "---  maxDoca:                  " << fMaxDoca << endl;
  cout << "---  maxD0:                    " << fMaxD0 << endl;
  cout << "---  maxDz:                    " << fMaxDz << endl;
  cout << "---  type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDstar::~HFDstar() {
  
}


// ----------------------------------------------------------------------
void HFDstar::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  // -- get the magnetic field
  ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  const MagneticField *field = magfield.product();

  // -- Filter on truth candidates
  TAnaCand *pCand(0);
  TAnaTrack *pT(0); 
  int OK(0); 
  for (int iC = 0; iC < gHFEvent->nCands(); ++iC) {
    pCand = gHFEvent->getCand(iC);
    if (54 == pCand->fType) {
      OK = 1; 
      for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
	pT = gHFEvent->getRecTrack(gHFEvent->getSigTrack(i)->fIndex); 
	if (pT->fPlab.Perp() < 0.) {
	  OK = 0; 
	}
      }
      break;
    }
  }
  if (0 == OK) return;
  cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
  for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
    pT = gHFEvent->getRecTrack(gHFEvent->getSigTrack(i)->fIndex); 
    pT->dump();
  }
  cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
  
  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFDstar> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFDstar> No primary vertex found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFDstar> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  // -- Transient track builder for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFDstar: Error: no TransientTrackBuilder found."<<endl;
    return;
  }


  // -- Build lists
  TLorentzVector muon, lambdac, dzero, dplus, kaon, pion, track, kaon1, kaon2, pion1, pion2, pion3;
  TLorentzVector tlv; 
  vector<pair<int, TLorentzVector> > kalist, pilist; 
  pilist.reserve(10000); 
  kalist.reserve(10000); 

  for (unsigned int itrack = 0; itrack < hTracks->size(); ++itrack){    
    TrackBaseRef rTrackView(hTracks, itrack);
    Track tTrack(*rTrackView);

    if (tTrack.d0() > fMaxD0) continue;
    if (tTrack.dz() > fMaxDz) continue;

    if (tTrack.pt() > fTrackPt)  {
      tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MPION); 
      pilist.push_back(make_pair(itrack, tlv));
    }

    if (tTrack.pt() > fTrackPt) {
      tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), MKAON); 
      kalist.push_back(make_pair(itrack, tlv));
    }
  }


  HFTwoParticleCombinatorics a(fVerbose); 
  vector<pair<int, int> > kapiList; 
  kapiList.reserve(100000); 
  a.combine(kapiList, kalist, pilist, MD0-fD0Window, MD0+fD0Window, 0); 
  
  HFSequentialVertexFit aSeq(hTracks, fTTB.product(), recoPrimaryVertexCollection, field, fVerbose);

  // -- Build D0 + track
  TLorentzVector d0, dstar, ka, pi, pis;
  for (unsigned int i = 0; i < kapiList.size(); ++i) {
    unsigned int iKaon = kapiList[i].first; 
    unsigned int iPion = kapiList[i].second; 
    
    TrackBaseRef kaonTrackView(hTracks, iKaon);
    Track tKaon(*kaonTrackView);
    if (tKaon.pt() < fTrackPt)  continue;
    ka.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MKAON); 

    TrackBaseRef pionTrackView(hTracks, iPion);
    Track tPion(*pionTrackView);
    if (tPion.pt() < fTrackPt)  continue;
    pi.SetPtEtaPhiM(tPion.pt(), tPion.eta(), tPion.phi(), MPION); 

    d0 = ka + pi; 
    if ((TMath::Abs(d0.M() - MD0) > fD0Window)) continue;
    
    for (unsigned int iTrack = 0; iTrack < hTracks->size(); ++iTrack){    
      if (iTrack == iKaon || iTrack == iPion) continue; 
      TrackBaseRef rTrackView(hTracks, iTrack);
      Track tSlowPion(*rTrackView);

      if (tSlowPion.d0() > fMaxD0) continue;
      if (tSlowPion.dz() > fMaxDz) continue;
      if (tSlowPion.pt() < fSlowPionPt) continue;

      pis.SetXYZM(tSlowPion.px(), tSlowPion.py(), tSlowPion.pz(), MPION); 
      if (d0.DeltaR(pis) > fDeltaR) continue; 
      
      dstar = d0 + pis; 
      if (TMath::Abs(dstar.M() - MDSTARPLUS) > 0.3) continue; 

      // -- sequential fit: D0 pi_slow
      if (fVerbose > 5) {
	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	cout << "==>HFDstar> going to sequential fit with track indices: " 
	     << iKaon << " " << iPion << " " << iTrack
	     << endl;
      }

      HFDecayTree theTree(300054, true, MDSTARPLUS, false, -1.0, true);
      
      HFDecayTreeIterator iterator = theTree.addDecayTree(300050, true, MD0, false);
      iterator->addTrack(iKaon,321);
      iterator->addTrack(iPion,211);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      theTree.addTrack(iTrack,211);
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      aSeq.doFit(&theTree);
      //       if (fVerbose > 6) {
      // 	cout << " XXXXXX: " << gHFEvent->nCands() << endl;
      // 	cout << " XXXXXX: " << gHFEvent->getCand(gHFEvent->nCands()-2)->fType << endl;
      // 	TAnaCand *pC = gHFEvent->getCand(gHFEvent->nCands()-2); 
      // 	TAnaCand *pD = gHFEvent->getCand(pC->fDau1); 
      // 	cout << " XXXXXX: HFDstarCand: idx = " << pC->fIndex << " type = " << pC->fType
      // 	     << " m* = " << pC->fMass << " m0 = " << pD->fMass << " dm = " << pC->fMass-pD->fMass << endl;
      // 	//  if (fVerbose > -1) cout << "   mdz = " << pD->fMass << endl;
      // 	for (int id = pD->fSig1; id <= pD->fSig2; ++id) {
      // 	  pT = gHFEvent->getRecTrack(gHFEvent->getSigTrack(id)->fIndex);
      // 	  cout << " XXXXXX: " << gHFEvent->getSigTrack(id)->fMCID << " " ; 
      // 	  pT->dump(); 
      // 	}
      // 	// -- slow pion
      // 	pT = gHFEvent->getRecTrack(gHFEvent->getSigTrack(pC->fSig1)->fIndex); 
      // 	cout << " XXXXXX: " << gHFEvent->getSigTrack(pC->fSig1)->fMCID << " " ; 
      // 	pT->dump(); 
      //       }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void  HFDstar::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFDstar::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFDstar);
