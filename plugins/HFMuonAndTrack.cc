#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFMuonAndTrack.h"

#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"

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
  cout << "---  Type:                     " << fType << endl;
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
  const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());
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
  edm::Handle<edm::View<reco::Track> > hTracks;
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
    //      cout << "%%>   Rec. Track #" << setw(2) << rTrack.index() << " pT: " 
    //           << setprecision(2) << setw(6) << trackView.pt() 
    //           <<  " matched to 0 MC Tracks" << endl;
  }

  if (!fTTB.isValid()) {
    cout << " -->HFMuonAndTrack: Error: no TransientTrackBuilder found."<<endl;
    return;
  }

  std::vector<reco::Track> fitTracks;
  TLorentzVector dimuon, m1, m2;
  int iMuon1(-1); 
  bool found_candidate = false;
  for (unsigned int imuon1 = 0; imuon1 < muonIndices.size(); ++imuon1) {
    TrackBaseRef mu1TrackView(hTracks, muonIndices[imuon1]);
    Track tMuon1(*mu1TrackView);
    iMuon1 = muonIndices[imuon1];
    
    if (tMuon1.pt() < fMuonPt)  continue;
   
    m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON); 
    

    for (unsigned int itrack2 = 0; itrack2 < hTracks->size(); ++itrack2){    
      if (static_cast<int>(itrack2) == muonIndices[imuon1]) continue; 

      TrackBaseRef rTrackView(hTracks, itrack2);
      Track tTrack2(*rTrackView);
      found_candidate = false;

      if (tTrack2.pt() < fTrackPt)  continue;

      //      if ( tMuon1.charge()*tTrack2.charge() < 0  ){
      m2.SetPtEtaPhiM(tTrack2.pt(), tTrack2.eta(), tTrack2.phi(), MMUON); 
      dimuon = m1 + m2;
      found_candidate = true; 
      //      } 
      
      if (found_candidate) {
      	if (dimuon.M() < fMassLow || dimuon.M() > fMassHigh) {
	  if (fVerbose > 0) {
	    cout << "==>HFMuonAndTrack> dimuon mass = " << dimuon.M() << ", skipping" << endl;
	  }
	  found_candidate = false;
	  continue; 
      	}
	
	
	if (fVerbose > 0) {
	  cout << "==>HFMuonAndTrack> dimuon mass = " << dimuon.M() << ", vertexing" << endl;
	  cout << "==>HFMuonAndTrack>  tMuon1.charge() = " << tMuon1.charge() << ", tTrack2.charge() " << tTrack2.charge() << endl;
	}
	
	// -- Vertex the two muons only
	//TAnaCand *pCand = gHFEvent->addCand();
	fitTracks.clear();
	fitTracks.push_back(tMuon1); 
	fitTracks.push_back(tTrack2);
	if ( fVertexing > 0 ) {
	  doVertexFit(fitTracks, iMuon1, itrack2);
	} else if ( fVertexing == 0 ) {
	  fillCandAndSignal(fitTracks, iMuon1, itrack2);
	}
      } 
    }
  }
}

// ----------------------------------------------------------------------
void HFMuonAndTrack::fillCandAndSignal(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2) {

  Track tMuon1 = Tracks[0]; 
  Track tMuon2 = Tracks[1]; 
  
  // -- Build composite
  TLorentzVector comp, M1, M2;
  M1.SetXYZM(Tracks[0].px(), Tracks[0].py(), Tracks[0].pz(), MMUON); 
  M2.SetXYZM(Tracks[1].px(), Tracks[1].py(), Tracks[1].pz(), MMUON); 
  comp = M1 + M2;
  
  // -- fill candidate
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fPlab = comp.Vect();
  pCand->fMass = comp.M();    
  pCand->fType = fType;
  pCand->fDau1 = -1;
  pCand->fDau2 = -1;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + 1;
    
  // -- fill refitted sig tracks
  TAnaTrack *pTrack = gHFEvent->addSigTrack();
  pTrack->fMCID     = tMuon1.charge()*13; 
  pTrack->fGenIndex = gHFEvent->getRecTrack(iMuon1)->fGenIndex; 
  pTrack->fQ        = tMuon1.charge();
  pTrack->fPlab.SetXYZ(Tracks[0].px(),
		       Tracks[0].py(),
		       Tracks[0].pz()
		       ); 
  pTrack->fIndex  = iMuon1;
  
  pTrack            = gHFEvent->addSigTrack();
  pTrack->fMCID     = gHFEvent->getRecTrack(iMuon2)->fMCID;
  pTrack->fMuID     = gHFEvent->getRecTrack(iMuon2)->fMuID; 
  pTrack->fGenIndex = gHFEvent->getRecTrack(iMuon2)->fGenIndex; 
  pTrack->fQ        = tMuon2.charge();
  pTrack->fPlab.SetXYZ(Tracks[1].px(),
		       Tracks[1].py(),
		       Tracks[1].pz()
		       ); 
  pTrack->fIndex  = iMuon2;
}

// ----------------------------------------------------------------------
void HFMuonAndTrack::doVertexFit(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2){

  Track tMuon1 = Tracks[0]; 
  Track tMuon2 = Tracks[1]; 
  
  std::vector<reco::TransientTrack> RecoTransientTrack;
  RecoTransientTrack.clear();
  RecoTransientTrack.push_back(fTTB->build(Tracks[0]));
  RecoTransientTrack.push_back(fTTB->build(Tracks[1]));
  
  // -- Do the vertexing
  KalmanVertexFitter theFitter(true);
  TransientVertex TransSecVtx = theFitter.vertex(RecoTransientTrack); 
  if (TransSecVtx.isValid() ) {
    if (isnan(TransSecVtx.position().x()) 
	|| isnan(TransSecVtx.position().y()) 
	|| isnan(TransSecVtx.position().z()) ) {
      if (fVerbose > 0) cout << "==>HFMuonAndTrack> Something went wrong! SecVtx nan - continue ... " << endl;
      //pCand->fType = -1;
      return; 
    }
  } else {
    if (fVerbose > 0) cout << "==>HFMuonAndTrack> KVF failed! continue ..." << endl;
    //pCand->fType = -1;
    return; 
  }
  
  // -- Get refitted tracks
  std::vector<reco::TransientTrack> refTT = TransSecVtx.refittedTracks();
  std::vector<reco::Track> refT; refT.clear(); 
  for(vector<reco::TransientTrack>::const_iterator i = refTT.begin(); i != refTT.end(); i++) {
    const Track &ftt = i->track();
    refT.push_back(ftt);
  }
  
  // -- Build composite
  TLorentzVector comp, M1, M2;
  M1.SetXYZM(refT[0].px(), refT[0].py(), refT[0].pz(), MMUON); 
  M2.SetXYZM(refT[1].px(), refT[1].py(), refT[1].pz(), MMUON); 
  comp = M1 + M2;
  
  
  // -- Build vertex for ntuple
  TAnaVertex anaVtx;
  ChiSquared chi(TransSecVtx.totalChiSquared(), TransSecVtx.degreesOfFreedom());
  anaVtx.setInfo(chi.value(), chi.degreesOfFreedom(), chi.probability(), 0, 0);
  anaVtx.fPoint.SetXYZ(TransSecVtx.position().x(), 
		       TransSecVtx.position().y(), 
		       TransSecVtx.position().z());
  
  anaVtx.addTrack(iMuon1);
  anaVtx.addTrack(iMuon2);
  
  VertexDistanceXY axy;
  anaVtx.fDxy     = axy.distance(fPV, TransSecVtx).value();
  anaVtx.fDxyE    = axy.distance(fPV, TransSecVtx).error();
  VertexDistance3D a3d;
  anaVtx.fD3d     = a3d.distance(fPV, TransSecVtx).value();
  anaVtx.fD3dE    = a3d.distance(fPV, TransSecVtx).error();
        
  // -- fill candidate
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fPlab = comp.Vect();
  pCand->fMass = comp.M();
  pCand->fVtx  = anaVtx;    
  pCand->fType = fType;
  pCand->fDau1 = -1;
  pCand->fDau2 = -1;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + 1;
  
  // -- fill refitted sig tracks
  TAnaTrack *pTrack = gHFEvent->addSigTrack();
  pTrack->fMCID     = tMuon1.charge()*13; 
  pTrack->fGenIndex = gHFEvent->getRecTrack(iMuon1)->fGenIndex; 
  pTrack->fQ        = tMuon1.charge();
  pTrack->fPlab.SetXYZ(refT[0].px(),
		       refT[0].py(),
		       refT[0].pz()
		       ); 
  pTrack->fIndex  = iMuon1;
  
  pTrack            = gHFEvent->addSigTrack();
  pTrack->fMCID     = gHFEvent->getRecTrack(iMuon2)->fMCID;
  pTrack->fMuID     = gHFEvent->getRecTrack(iMuon2)->fMuID; 
  pTrack->fGenIndex = gHFEvent->getRecTrack(iMuon2)->fGenIndex; 
  pTrack->fQ        = tMuon2.charge();
  pTrack->fPlab.SetXYZ(refT[1].px(),
		       refT[1].py(),
		       refT[1].pz()
		       ); 
  pTrack->fIndex  = iMuon2;
  
}


// ------------ method called once each job just before starting event loop  ------------
void  HFMuonAndTrack::beginJob(const edm::EventSetup& setup) {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFMuonAndTrack::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFMuonAndTrack);
