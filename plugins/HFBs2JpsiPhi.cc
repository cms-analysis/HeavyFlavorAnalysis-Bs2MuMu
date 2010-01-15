#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBs2JpsiPhi.h"

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
#define MKAON 0.49368
#define MJPSI 3.096916
#define MPHI  1.01946
#define MBS0  5.27953

// ----------------------------------------------------------------------
HFBs2JpsiPhi::HFBs2JpsiPhi(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("goodTracks"))), 
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel")),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)), 
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 0.4)), 
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 99.)),
  fType(iConfig.getUntrackedParameter<int>("type", 531)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBs2JpsiPhi constructor" << endl;
  cout << "---  tracksLabel:              " << fTracksLabel << endl;
  cout << "---  muonsLabel:               " << fMuonsLabel << endl;
  cout << "---  trackPt:                  " << fTrackPt << endl;
  cout << "---  muonPt:                   " << fMuonPt << endl;
  cout << "---  deltaR:                   " << fDeltaR << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFBs2JpsiPhi::~HFBs2JpsiPhi() {
  
}


// ----------------------------------------------------------------------
void HFBs2JpsiPhi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // -- get the primary vertex
  edm::Handle<VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(fPrimaryVertexLabel, recoPrimaryVertexCollection);
  if(!recoPrimaryVertexCollection.isValid()) {
    cout << "==>HFBs2JpsiPhi> No primary vertex collection found, skipping" << endl;
    return;
  }
  const reco::VertexCollection vertices = *(recoPrimaryVertexCollection.product());
  if (vertices.size() == 0) {
    cout << "==>HFBs2JpsiPhi> No primary vertex found, skipping" << endl;
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
  cout << "==>HFBs2JpsiPhi> No valid MuonCollection with label "<< fMuonsLabel <<" found, skipping" << endl;
    return;
  }
  
  // -- get the collection of tracks
  edm::Handle<edm::View<reco::Track> > hTracks;
  iEvent.getByLabel(fTracksLabel, hTracks);
  if(!hTracks.isValid()) {
    cout << "==>HFBs2JpsiPhi> No valid TrackCollection with label "<<fTracksLabel <<" found, skipping" << endl;
    return;
  }

  // -- get the collection of muons and store their corresponding track indices
  vector<int> muonIndices;
  for (MuonCollection::const_iterator muon = hMuons->begin(); muon != hMuons->end(); ++muon) {
    int im = muon->track().index(); 
    if (im > 0) muonIndices.push_back(im);
  }
  if (fVerbose > 0) {
    cout << "==>HFBs2JpsiPhi> nMuons = " << hMuons->size() << endl;
    cout << "==>HFBs2JpsiPhi> nMuonIndices = " << muonIndices.size() << endl;
  }
  if (muonIndices.size() < 2) return;

  // -- reconstruct J/Psi, based on 'best' dimuon mass
  TLorentzVector psi, cpsi, m1, m2;
  double dmass(1000.0);
  Track tMuon1, tMuon2; 
  int iMuon1(-1), iMuon2(-1); 
  for (unsigned int imuon1 = 0; imuon1 < muonIndices.size()-1; ++imuon1) {
    TrackBaseRef mu1TrackView(hTracks, muonIndices[imuon1]);
    Track muon1(*mu1TrackView);
    if (muon1.pt() < fMuonPt)  continue;

    for (unsigned int imuon2 = imuon1 + 1; imuon2 < muonIndices.size(); ++imuon2) {
      TrackBaseRef mu2TrackView(hTracks, muonIndices[imuon2]);
      Track muon2(*mu2TrackView);
      if (muon2.pt() < fMuonPt)  continue;

      m1.SetPtEtaPhiM(muon1.pt(), muon1.eta(), muon1.phi(), MMUON); 
      m2.SetPtEtaPhiM(muon2.pt(), muon2.eta(), muon2.phi(), MMUON); 
      cpsi = m1 + m2; 
      if ((TMath::Abs(psi.M()-MJPSI) < dmass) 
	  //&& (muon1.charge()*muon2.charge() < 0) 
	  ) {
	psi    = cpsi;
	dmass  = TMath::Abs(psi.M()-MJPSI);
	tMuon1 = muon1;
	iMuon1 = muonIndices[imuon1]; 
	tMuon2 = muon2;
	iMuon2 = muonIndices[imuon2]; 
      }
    }
  }

  if (iMuon1 < 0 || iMuon2 < 0) {
    //if (fVerbose > 0) cout << " --> HFBs2JpsiPhi: No unlike charge dimuons found, dmass = " << dmass << endl;
    return;
  }
  if (dmass >  0.8) {
    if (fVerbose > 0) cout << " --> HFBs2JpsiPhi: No J/Psi candidate found: " << dmass << endl;
    return;
  }

  // -- Transient tracks for vertexing
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) {
    cout << " -->HFBs2JpsiPhi: Error: no TransientTrackBuilder found."<<endl;
    return;
  }
  std::vector<reco::Track> fitTracks;
  
  TAnaVertex BdVtx;
  TAnaTrack *pTrack; 

  TLorentzVector ka1, ka2, phi, bs0; 
  Track tKaon1, tKaon2; 
  int   iKaon1, iKaon2; 

  // -- Loop over all phi daughters
  for (unsigned int iTrack1 = 0; iTrack1 < hTracks->size(); ++iTrack1){    
    if (static_cast<int>(iTrack1) == iMuon1 || static_cast<int>(iTrack1) == iMuon2) continue; 
    TrackBaseRef rTrackView(hTracks, iTrack1);
    Track track1(*rTrackView);
    if (track1.pt() < fTrackPt)  continue;
    
    for (unsigned int iTrack2 = 0; iTrack2 < hTracks->size(); ++iTrack2){    
      if (iTrack2 == iTrack1) continue; 
      if (static_cast<int>(iTrack2) == iMuon1 || static_cast<int>(iTrack2) == iMuon2) continue; 
      TrackBaseRef sTrackView(hTracks, iTrack2);
      Track track2(*sTrackView);
      
      if (track2.pt() < fTrackPt)  continue;
      if (track1.charge()*track2.charge() > 0) continue;
	
      tKaon1 = track1; iKaon1 = iTrack1; 
      tKaon2 = track2; iKaon2 = iTrack2; 

      ka1.SetXYZM(track1.px(), track1.py(), track1.pz(), MKAON); 
      ka2.SetXYZM(track2.px(), track2.py(), track2.pz(), MKAON); 
      phi = ka1 + ka2;
      if (TMath::Abs(phi.M() - MPHI) > 0.3) continue; 
      bs0 = phi + psi; 
      if (TMath::Abs(bs0.M() - MBS0) > 1.3) continue; 

      fitTracks.clear();
      fitTracks.push_back(tMuon1); 
      fitTracks.push_back(tMuon2); 
      fitTracks.push_back(tKaon1); 
      fitTracks.push_back(tKaon2); 

      RefCountedKinematicTree bsTree; 
      int OK = doVertexFit(fitTracks, bsTree); 
      if (OK == 0 || !bsTree->isValid()) {
	if (fVerbose > 0) {
	  cout << "----------------------------------------" << endl;
	  cout << "==> HFBs2JpsiPhi: Invalid kinematic fit tree !! " << endl;
	  cout << tMuon1.pt() << ", "  << tMuon1.d0() << ", "  << tMuon1.dz() << endl;
	  cout << tMuon2.pt() << ", "  << tMuon2.d0() << ", "  << tMuon2.dz() << endl;
	  cout << tKaon1.pt() << ", "  << tKaon1.d0() << ", "  << tKaon1.dz() << endl;
	  cout << tKaon2.pt() << ", "  << tKaon2.d0() << ", "  << tKaon2.dz() << endl;
	  cout << "----------------------------------------" << endl;
	}
	continue; 
      }
      
      bsTree->movePointerToTheTop();
      RefCountedKinematicVertex bsVertex = bsTree->currentDecayVertex();
      RefCountedKinematicParticle bs = bsTree->currentParticle();

      TVector3 bsPlab = (bs->currentState().globalMomentum().x(), 
			 bs->currentState().globalMomentum().y(), 
			 bs->currentState().globalMomentum().z()
			 );
      AlgebraicVector7 bsParameters = bs->currentState().kinematicParameters().vector();
      double bsMass = bsParameters[6];

      if (bsVertex->vertexIsValid() && bsMass > 4.5 && bsMass < 6.5) {

	if (fVerbose > 0) {
	  cout << "----------------------------------------" << endl;
	  cout << "==> HFBs2JpsiPhi: Filling B0S candidate with mass = " << bsMass << endl;
	  cout << tMuon1.pt() << ", "  << tMuon1.d0() << ", "  << tMuon1.dz() << endl;
	  cout << tMuon2.pt() << ", "  << tMuon2.d0() << ", "  << tMuon2.dz() << endl;
	  cout << tKaon1.pt() << ", "  << tKaon1.d0() << ", "  << tKaon1.dz() << endl;
	  cout << tKaon2.pt() << ", "  << tKaon2.d0() << ", "  << tKaon2.dz() << endl;
	  cout << "----------------------------------------" << endl;
	}

	TAnaVertex anaVt;

	ChiSquared chi(bs->chiSquared(), bs->degreesOfFreedom());
	anaVt.setInfo(bs->chiSquared(), bs->degreesOfFreedom(), chi.probability(), 0, 0);
	
	anaVt.fPoint.SetXYZ(bsVertex->position().x(), 
			    bsVertex->position().y(), 
			    bsVertex->position().z()
			    );
	
	// -- Distance to primary vertex
	VertexDistanceXY axy;
	double dXY      = axy.distance(fPV, bsVertex->vertexState()).value();
	double dXYE     = axy.distance(fPV, bsVertex->vertexState()).error();
	
	VertexDistance3D a3d;
	double d3d      = a3d.distance(fPV, bsVertex->vertexState()).value();
	double d3dE     = a3d.distance(fPV, bsVertex->vertexState()).error();
	
	anaVt.fDxy  = dXY; 
	anaVt.fDxyE = dXYE; 
	
	anaVt.fD3d  = d3d; 
	anaVt.fD3dE = d3dE; 
	anaVt.addTrack(iMuon1); 
	anaVt.addTrack(iMuon2); 
	anaVt.addTrack(iKaon1); 
	anaVt.addTrack(iKaon2); 


	// -- fill candidates
	TAnaCand  *pCand = gHFEvent->addCand();
	
	pCand->fPlab = bsPlab;
	pCand->fMass = bsMass;
	pCand->fVtx  = anaVt;    
	pCand->fType = fType;
	pCand->fDau1 = gHFEvent->nCands();
	pCand->fDau2 = pCand->fDau1+1;
	pCand->fSig1 = -1;
	pCand->fSig2 = -1;

	pCand = gHFEvent->addCand();
	pCand->fPlab = psi.Vect();
	pCand->fMass = psi.M();
	pCand->fType = fType*2000 + 443;
	pCand->fDau1 = -1;
	pCand->fDau2 = -1;
	pCand->fSig1 = iMuon1;
	pCand->fSig2 = iMuon2;

	pCand = gHFEvent->addCand();
	pCand->fPlab = phi.Vect();
	pCand->fMass = phi.M();
	pCand->fType = fType*2000 + 333;
	pCand->fDau1 = -1;
	pCand->fDau2 = -1;
	pCand->fSig1 = iKaon1;
	pCand->fSig2 = iKaon2;


	// -- fill (not yet refitted) sig tracks
	pTrack            = gHFEvent->addSigTrack();
	pTrack->fMCID     = tMuon1.charge()*13; 
	pTrack->fGenIndex = -1; 
	pTrack->fQ        = tMuon1.charge();
	pTrack->fPlab.SetPtEtaPhi(tMuon1.pt(),
				  tMuon1.eta(),
				  tMuon1.phi()
				  ); 
	pTrack->fIndex  = iMuon1;

	pTrack            = gHFEvent->addSigTrack();
	pTrack->fMCID     = tMuon2.charge()*13; 
	pTrack->fGenIndex = -1; 
	pTrack->fQ        = tMuon2.charge();
	pTrack->fPlab.SetPtEtaPhi(tMuon2.pt(),
				  tMuon2.eta(),
				  tMuon2.phi()
				  ); 
	pTrack->fIndex  = iMuon2;

	pTrack            = gHFEvent->addSigTrack();
	pTrack->fMCID     = tKaon1.charge()*321; 
	pTrack->fGenIndex = -1; 
	pTrack->fQ        = tKaon1.charge();
	pTrack->fPlab.SetPtEtaPhi(tKaon1.pt(),
				  tKaon1.eta(),
				  tKaon1.phi()
				  ); 
	pTrack->fIndex  = iKaon1;
	
	pTrack            = gHFEvent->addSigTrack();
	pTrack->fMCID     = tKaon2.charge()*321; 
	pTrack->fGenIndex = -1; 
	pTrack->fQ        = tKaon2.charge();
	pTrack->fPlab.SetPtEtaPhi(tKaon2.pt(),
				  tKaon2.eta(),
				  tKaon2.phi()
				  ); 
	pTrack->fIndex  = iKaon2;


	// -- Vertex the two muons only
	pCand = gHFEvent->addCand(); 
	fitTracks.clear();
	fitTracks.push_back(tMuon1); 
	fitTracks.push_back(tMuon2); 
	doJpsiVertexFit(fitTracks, iMuon1, iMuon2, pCand); 

	
      }
    }

  }
  
}


// ----------------------------------------------------------------------
void HFBs2JpsiPhi::doJpsiVertexFit(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2, TAnaCand *pCand){

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
      cout << "==>HFDimuons> Something went wrong! SecVtx nan - continue ... " << endl;
      pCand->fType = -1;
      return; 
    }
  } else {
    cout << "==>HFDimuons> KVF failed! continue ..." << endl;
    pCand->fType = -1;
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
  pCand->fPlab = comp.Vect();
  pCand->fMass = comp.M();
  pCand->fVtx  = anaVtx;    
  pCand->fType = fType*1000 + 443;
  pCand->fDau1 = -1;
  pCand->fDau2 = -1;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + 1;
  
  // -- fill refitted sig tracks
  TAnaTrack *pTrack = gHFEvent->addSigTrack();
  pTrack->fMCID     = tMuon1.charge()*13; 
  pTrack->fGenIndex = -1; 
  pTrack->fQ        = tMuon1.charge();
  pTrack->fPlab.SetXYZ(refT[0].px(),
		       refT[0].py(),
		       refT[0].pz()
		       ); 
  pTrack->fIndex  = iMuon1;
  
  pTrack            = gHFEvent->addSigTrack();
  pTrack->fMCID     = tMuon2.charge()*13; 
  pTrack->fGenIndex = -1; 
  pTrack->fQ        = tMuon2.charge();
  pTrack->fPlab.SetXYZ(refT[1].px(),
		       refT[1].py(),
		       refT[1].pz()
		       ); 
  pTrack->fIndex  = iMuon2;
}


// ----------------------------------------------------------------------
int HFBs2JpsiPhi::doVertexFit(std::vector<reco::Track> &Tracks, RefCountedKinematicTree &bsTree){

  KinematicParticleFactoryFromTransientTrack pFactory;
  
  ParticleMass muon_mass = MMUON;
  ParticleMass kaon_mass = MKAON;
  float muon_sigma = 0.0000001;
  float kaon_sigma = 0.000016;

  float chi = 0.;
  float ndf = 0.;

  TransientTrack ttMuon1 = fTTB->build(Tracks[0]);
  TransientTrack ttMuon2 = fTTB->build(Tracks[1]);
  TransientTrack ttKaon1 = fTTB->build(Tracks[2]);
  TransientTrack ttKaon2 = fTTB->build(Tracks[3]);
  
  vector<RefCountedKinematicParticle> muonParticles;
  vector<RefCountedKinematicParticle> allParticles;
  muonParticles.push_back(pFactory.particle(ttMuon1, muon_mass, chi, ndf, muon_sigma));
  muonParticles.push_back(pFactory.particle(ttMuon2, muon_mass, chi, ndf, muon_sigma));
  allParticles.push_back(pFactory.particle(ttMuon1,  muon_mass, chi, ndf, muon_sigma));
  allParticles.push_back(pFactory.particle(ttMuon2,  muon_mass, chi, ndf, muon_sigma));
  
  allParticles.push_back(pFactory.particle(ttKaon1, kaon_mass, chi, ndf, kaon_sigma));
  allParticles.push_back(pFactory.particle(ttKaon2, kaon_mass, chi, ndf, kaon_sigma));
  
  //creating the vertex fitter
  KinematicParticleVertexFitter kpvFitter;
  RefCountedKinematicTree jpTree;
  try {
    jpTree = kpvFitter.fit(muonParticles);
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "HFBs2JpsiPhi: exception caught from kin fit 1" << endl;
  }
  if (jpTree->isEmpty()) {
    return 0;
  } 

  KinematicParticleFitter csFitter;
  ParticleMass jpsi = MJPSI;
  float jp_m_sigma = 0.00004;
  KinematicConstraint * jpsi_c2 = new MassKinematicConstraint(jpsi, jp_m_sigma);
  try {
    jpTree = csFitter.fit(jpsi_c2, jpTree);
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "HFBu2JpsiKp: exception caught from kin fit 2" << endl;
    delete jpsi_c2; 
    return 0; 
  }
  if (jpTree->isEmpty()) return 0; 

  try {
    jpTree->movePointerToTheTop();
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "HFBu2JpsiKp: exception caught from kin fit " << endl;
    delete jpsi_c2; 
    return 0; 
  }

  RefCountedKinematicParticle jpsi_part = jpTree->currentParticle();
  allParticles.push_back(jpsi_part);
  bsTree = kpvFitter.fit(allParticles);
  delete jpsi_c2; 
  return 1;
}


// ------------ method called once each job just before starting event loop  ------------
void  HFBs2JpsiPhi::beginJob(const edm::EventSetup& setup) {
}


// ------------ method called once each job just after ending the event loop  ------------
void  HFBs2JpsiPhi::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBs2JpsiPhi);
