#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBd2DstarPi.h"

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
HFBd2DstarPi::HFBd2DstarPi(const ParameterSet& iConfig) :
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
  cout << "--- HFBd2DstarPi constructor" << endl;
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
HFBd2DstarPi::~HFBd2DstarPi() {
  
}


// ----------------------------------------------------------------------
void HFBd2DstarPi::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  // -- get the magnetic field
  ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  const MagneticField *field = magfield.product();

  //   // -- Filter on truth candidates
  //   TAnaCand *pCand(0);
  //   TAnaTrack *pT(0); 
  //   int OK(0); 
  //   for (int iC = 0; iC < gHFEvent->nCands(); ++iC) {
  //     pCand = gHFEvent->getCand(iC);
  //     if (54 == pCand->fType) {
  //       OK = 1; 
  //       for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
  // 	pT = gHFEvent->getRecTrack(gHFEvent->getSigTrack(i)->fIndex); 
  // 	if (pT->fPlab.Perp() < 0.) {
  // 	  OK = 0; 
  // 	}
  //       }
  //       break;
  //     }
  //   }
  //   if (0 == OK) return;
  //   cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
  //   for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
  //     pT = gHFEvent->getRecTrack(gHFEvent->getSigTrack(i)->fIndex); 
  //     pT->dump();
  //   }
  //   cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
  
  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFBd2DstarPi> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFBd2DstarPi> No primary vertex found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFBd2DstarPi> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  // -- Transient track builder for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFBd2DstarPi: Error: no TransientTrackBuilder found."<<endl;
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
  
  vector<pair<int, int> > phiList; 
  a.combine(phiList, pilist, pilist, 0., 6., 1); 
  if (fVerbose > 0) cout << "==>HFBd2DstarPi> phi (haha) list size: " << phiList.size() << endl;
  
  HFSequentialVertexFit aSeq(hTracks, fTTB.product(), recoPrimaryVertexCollection, field, fVerbose);

  // -- Build D0 + track
  TLorentzVector b0, dstar, d0, ka, pi, pis, pif;
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

    for (unsigned int j = 0; j < phiList.size(); ++j){    
      unsigned int iPion1 = phiList[j].first; // slow pion!
      unsigned int iPion2 = phiList[j].second; 
      if (iPion1 == iKaon || iPion2 == iKaon) continue; 
      if (iPion1 == iPion || iPion2 == iPion) continue; 

      TrackBaseRef rTrackView1(hTracks, iPion1);
      Track tPion1(*rTrackView1);
      if (tPion1.pt() < fSlowPionPt) continue;
      pis.SetXYZM(tPion1.px(), tPion1.py(), tPion1.pz(), MPION); 

      TrackBaseRef rTrackView2(hTracks, iPion2);
      Track tPion2(*rTrackView2);
      if (tPion2.pt() < fTrackPt) continue;
      pif.SetXYZM(tPion2.px(), tPion2.py(), tPion2.pz(), MPION); 


      // -- the slow D* pion has the same charge like the fast D0 pion, and other charge correlations
      if (tPion1.charge()*tPion.charge() < 0) continue;       
      if (tPion1.charge()*tPion2.charge() > 0) continue; 
      if (tKaon.charge()*tPion.charge() > 0) continue; 

      if (tPion2.d0() > fMaxD0) continue;
      if (tPion2.dz() > fMaxDz) continue;

      if (d0.DeltaR(pis) > fDeltaR) continue; 
      
      dstar = d0 + pis; 
      double dm = dstar.M() - d0.M();
      if (TMath::Abs(0.145 - dm) > fDeltaM) continue;  // keep 0.135 .. 0.155 
      if (dstar.DeltaR(pif) > fDeltaR) continue; 

      b0 = dstar + pif; 
      if (TMath::Abs(b0.M() - MB_0) > 0.6) continue; 
      if (b0.DeltaR(pif) > fDeltaR) continue; 
      if (b0.DeltaR(d0) > fDeltaR) continue; 
      
      // -- sequential fit: D0 pi_slow pi_fast
      if (fVerbose > 5) {
	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	cout << "==>HFBd2DstarPi> going to sequential fit with track indices: " 
	     << iKaon << " " << iPion << " " << iPion1 << " " << iPion2
	     << endl;
      }

//       HFDecayTree theTree(300030, true, MB_0, false, -1.0, true);
//       theTree.addTrack(iPion2,211); // add the fast pion to the B0 candidate
//       theTree.addTrack(iPion1,211); // add the slow pion to the D*+ candidate
//       theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
//       HFDecayTreeIterator iterator = theTree.addDecayTree(300032, true, MD0, false); // D0
//       iterator->addTrack(iKaon,321);
//       iterator->addTrack(iPion,211);
//       iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
//       aSeq.doFit(&theTree);
      
      HFDecayTree theTree(300030, true, MB_0, false, -1.0, true);
      theTree.addTrack(iPion2,211); // add the fast pion to the B0 candidate
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      HFDecayTreeIterator iterator = theTree.addDecayTree(300031, false, MDSTARPLUS, false); // D*-
      iterator->addTrack(iPion1,211); // add the slow pion to the D*+ candidate
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      iterator = iterator->addDecayTree(300032, true, MD0, false); // D0
      iterator->addTrack(iKaon,321);
      iterator->addTrack(iPion,211);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      
      aSeq.doFit(&theTree);
      
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void  HFBd2DstarPi::beginJob() {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFBd2DstarPi::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBd2DstarPi);
