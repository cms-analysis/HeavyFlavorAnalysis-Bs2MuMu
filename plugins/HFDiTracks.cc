#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDiTracks.h"

#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKalmanVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"
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
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 3.0)), 
  fTrack1Mass(iConfig.getUntrackedParameter<double>("track1Mass", 0.1396)), 
  fTrack2Mass(iConfig.getUntrackedParameter<double>("track2Mass", 0.1396)), 
  fMassLow(iConfig.getUntrackedParameter<double>("massLow", 0.0)), 
  fMassHigh(iConfig.getUntrackedParameter<double>("massHigh", 12.0)), 
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.05)),
  fPvWeight(iConfig.getUntrackedParameter<double>("pvWeight", 0.0)),
  fSameSign(iConfig.getUntrackedParameter<bool>("sameSign", false)), 
  fType(iConfig.getUntrackedParameter<int>("type", 1300)) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDiTracks constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  track1Mass:               " << fTrack1Mass << endl;
  cout << "---  track2Mass:               " << fTrack2Mass << endl;
  cout << "---  massLow:                  " << fMassLow << endl;
  cout << "---  massHigh:                 " << fMassHigh << endl;
  cout << "---  maxDoca:                  " << fMaxDoca << endl;
  cout << "---  pvWeight:                 " << fPvWeight << endl;
  cout << "---  sameSign:                 " << fSameSign << endl;
  cout << "---  Type:                     " << fType << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDiTracks::~HFDiTracks() {
  
}


// ----------------------------------------------------------------------
void HFDiTracks::analyze(const Event& iEvent, const EventSetup& iSetup) {
  // -- get the magnetic field
  ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  const MagneticField *field = magfield.product();
 
  // -- get the primary vertex
  Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFDitracks> No primary vertex collection found, skipping" << endl;
    return;
  }
  const VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFDitracks> No primary vertex found, skipping" << endl;
    return;
  }
  fPV = vertices[gHFEvent->fBestPV]; 
  if (fVerbose > 0) {
    cout << "HFDitracks: Taking vertex " << gHFEvent->fBestPV << " with ntracks = " << fPV.tracksSize() << endl;
  }
  
  // -- get the collection of tracks
  Handle<View<Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFDitracks> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  // -- Build lists
  vector<pair<int, TLorentzVector> > t1list, t2list; 
  TLorentzVector tlv; 
  t1list.reserve(1000); 
  t2list.reserve(1000); 

  for (unsigned int itrack = 0; itrack < hTracks->size(); ++itrack){    
    TrackBaseRef rTrackView(hTracks, itrack);
    Track tTrack(*rTrackView);

    if (tTrack.pt() > fTrackPt)  {
      tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), fTrack1Mass); 
      t1list.push_back(make_pair(itrack, tlv));
    }

    if (tTrack.pt() > fTrackPt) {
      tlv.SetXYZM(tTrack.px(), tTrack.py(), tTrack.pz(), fTrack2Mass); 
      t2list.push_back(make_pair(itrack, tlv));
    }
  }


  HFTwoParticleCombinatorics a(fVerbose); 
  vector<pair<int, int> > ttList; 
  ttList.reserve(100000); 
  if (TMath::Abs(fTrack1Mass - fTrack2Mass) < 0.0001) {
    a.combine(ttList, t1list, t2list, fMassLow, fMassHigh, 1); 
  } else {
    a.combine(ttList, t1list, t2list, fMassLow, fMassHigh, 0); 
  }
  
  // -- Transient tracks for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFDitracks: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  // -- set up vertex fitter 
  HFSequentialVertexFit aSeq(hTracks, fTTB.product(), recoPrimaryVertexCollection, field, fVerbose);

  TLorentzVector ditrack, t1, t2;
  int n1=0, n2=0;
  for (unsigned int i = 0; i < ttList.size(); ++i) {
    unsigned int it1 = ttList[i].first; 
    unsigned int it2 = ttList[i].second; 
    n1++;

    TrackBaseRef t1View(hTracks, it1);
    Track tT1(*t1View);
    if (tT1.pt() < fTrackPt)  continue;
    t1.SetPtEtaPhiM(tT1.pt(), tT1.eta(), tT1.phi(), fTrack1Mass); 

    TrackBaseRef t2View(hTracks, it2);
    Track tT2(*t2View);
    if (tT2.pt() < fTrackPt)  continue;
    t2.SetPtEtaPhiM(tT2.pt(), tT2.eta(), tT2.phi(), fTrack2Mass); 

    // Select opposite charge pairs (or same)
    if(fSameSign) {if( tT1.charge() != tT2.charge() ) continue;}  // select track pairs with same charge
    else          {if( tT1.charge() == tT2.charge() ) continue;}  // select track pairs with opposite charge

    ditrack = t1 + t2; 
    if (ditrack.M() < fMassLow || ditrack.M() > fMassHigh) continue;

    
    // -- Vertexing, with Kinematic Particles
    HFDecayTree theTree(fType, true, 0, false); 
    theTree.addTrack(it1, idFromMass(fTrack1Mass));
    theTree.addTrack(it2, idFromMass(fTrack2Mass));
    //theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
    theTree.setNodeCut(RefCountedHFNodeCut(new HFPvWeightCut(fMaxDoca,fPvWeight)));
    
    aSeq.doFit(&theTree);
    n2++;
  }
  //cout<<n1<<" "<<n2<<endl;
}

// ------------ method called once each job just before starting event loop  ------------
void  HFDiTracks::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDiTracks::endJob() {
}

// ----------------------------------------------------------------------
int HFDiTracks::idFromMass(double mass) {
  if (TMath::Abs(mass - MELECTRON) < 0.0001) return 11;
  if (TMath::Abs(mass - MMUON) < 0.0001) return 13;
  if (TMath::Abs(mass - MPION) < 0.0001) return 211;
  if (TMath::Abs(mass - MKAON) < 0.0001) return 321;
  if (TMath::Abs(mass - MPROTON) < 0.0001) return 2212;
  cout << "%%%> HFDiTracks: mass " << mass << " not associated to any stable particle" << endl;
  return 0; 
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFDiTracks);
