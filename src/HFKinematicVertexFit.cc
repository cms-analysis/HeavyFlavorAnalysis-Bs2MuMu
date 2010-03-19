#include <iostream>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFKinematicVertexFit.hh"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"


// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFKinematicVertexFit::HFKinematicVertexFit(const TransientTrackBuilder *TTB, 
					   Vertex  &PV,
					   int type,
					   int verbose) {
  
  fpTTB     = TTB;
  fPV       = PV;
  fType     = type; 
  fVerbose  = verbose;
}

// ----------------------------------------------------------------------
HFKinematicVertexFit::~HFKinematicVertexFit() {
  if (fVerbose > 2) cout << "This is the end" << endl;
}


// ----------------------------------------------------------------------
void HFKinematicVertexFit::doJpsiFit(vector<Track>  &trackList, 
				     vector<int>    &trackIndices,
				     vector<double> &trackMasses, 
				     int type
				     ) {
  
  KinematicParticleFactoryFromTransientTrack pFactory;
  
  float muon_sigma = 0.0000001;
  float kaon_sigma = 0.000016;

  float chi = 0.;
  float ndf = 0.;
  
  TransientTrack ttMuon1 = fpTTB->build(trackList[0]);
  TransientTrack ttMuon2 = fpTTB->build(trackList[1]);
  vector<TransientTrack> otherTracks; 
  for (unsigned int i = 2; i < trackList.size(); ++i) {
    TransientTrack tTrack  = fpTTB->build(trackList[i]);
    otherTracks.push_back(tTrack); 
  }
  
  vector<RefCountedKinematicParticle> muonParticles;
  muonParticles.push_back(pFactory.particle(ttMuon1, trackMasses[0], chi, ndf, muon_sigma));
  muonParticles.push_back(pFactory.particle(ttMuon2, trackMasses[1], chi, ndf, muon_sigma));
  
  vector<RefCountedKinematicParticle> allParticles;
  for (unsigned int i = 0; i < otherTracks.size(); ++i) {
    allParticles.push_back(pFactory.particle(otherTracks[i], trackMasses[i+2], chi, ndf, kaon_sigma));
  }
  
  KinematicParticleVertexFitter kpvFitter;
  RefCountedKinematicTree jpTree; 
  try {
    jpTree = kpvFitter.fit(muonParticles);
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "==>HFKinematicVertexFitter: exception caught from kin fit in J/psi fit" << endl;
  }
  if (jpTree->isEmpty()) {
    if (fVerbose > 0) cout << "==>HFKinematicVertexFitter: empty J/psi tree!! returning!!!" << endl;
    return;
  } 
  
  KinematicParticleFitter csFitter;
  ParticleMass jpsi = MJPSI;
  float jp_m_sigma = 0.00004;
  KinematicConstraint *jpsi_c2 = new MassKinematicConstraint(jpsi, jp_m_sigma);
  try {
    jpTree = csFitter.fit(jpsi_c2, jpTree);
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "==>HFKinematicVertexFitter: exception caught from kin fit in constrained J/psi fit" << endl;
    delete jpsi_c2; 
    return; 
  }
  if (jpTree->isEmpty()) {
    if (fVerbose > 0) cout << "==>HFKinematicVertexFitter: empty constrained J/psi Tree!! returning!!!" << endl;
    return; 
  }
  
  try {
    jpTree->movePointerToTheTop();
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "==>HFKinematicVertexFitter: exception caught from kin fit 3" << endl;
    delete jpsi_c2; 
    return; 
  }

  RefCountedKinematicParticle jpsi_part = jpTree->currentParticle();
  allParticles.push_back(jpsi_part);
  RefCountedKinematicTree buTree = kpvFitter.fit(allParticles);  
  delete jpsi_c2; 

  // -- Fill ntuple
  if (buTree->isEmpty()) {
    if (fVerbose > 0) cout << "==>HFKinematicVertexFitter: empty buTree!! returning!!!" << endl;
    return; 
  }
  try {
    buTree->movePointerToTheTop();
  } catch (cms::Exception &ex) {
    //    cout << ex.explainSelf() << endl;
    if (fVerbose > 0) cout << "==>HFKinematicVertexFitter: exception caught from kin fit 3" << endl;
    return; 
  }
  
  RefCountedKinematicVertex buVertex = buTree->currentDecayVertex();
  RefCountedKinematicParticle Bu = buTree->currentParticle();
  
  TVector3 buPlab(Bu->currentState().globalMomentum().x(), 
		     Bu->currentState().globalMomentum().y(), 
		     Bu->currentState().globalMomentum().z()
		     );
  AlgebraicVector7 buParameters = Bu->currentState().kinematicParameters().vector();
  double buMass = buParameters[6];
  
  if (buVertex->vertexIsValid()) {
    
    if (fVerbose > 0) {
      cout << "----------------------------------------" << endl;
      cout << "==> HFKinematicVertexFit: Filling candidate with mass = " << buMass << endl;
      for (unsigned int i = 0; i < trackList.size(); ++i) {
	cout << trackList[i].px() << ", "  << trackList[i].py() << ", "  << trackList[i].pz() << endl;
      }
      cout << "----------------------------------------" << endl;
    }
    
    TAnaVertex anaVt;
    
    ChiSquared chi(Bu->chiSquared(), Bu->degreesOfFreedom());
    anaVt.setInfo(Bu->chiSquared(), Bu->degreesOfFreedom(), chi.probability(), 0, 0);
    
    anaVt.fPoint.SetXYZ(buVertex->position().x(), 
			buVertex->position().y(), 
			buVertex->position().z()
			);
    
    // -- Distance to primary vertex
    VertexDistanceXY axy;
    double dXY      = axy.distance(fPV, buVertex->vertexState()).value();
    double dXYE     = axy.distance(fPV, buVertex->vertexState()).error();
    
    VertexDistance3D a3d;
    double d3d      = a3d.distance(fPV, buVertex->vertexState()).value();
    double d3dE     = a3d.distance(fPV, buVertex->vertexState()).error();
      
    anaVt.fDxy  = dXY; 
    anaVt.fDxyE = dXYE; 
    
    anaVt.fD3d  = d3d; 
    anaVt.fD3dE = d3dE; 
      
    // -- fill candidate
    TAnaCand  *pCand = gHFEvent->addCand();
    
    pCand->fPlab = buPlab;
    pCand->fMass = buMass;
    pCand->fVtx  = anaVt;    
    pCand->fType = (type != 0? type: fType);
    pCand->fDau1 = -1;
    pCand->fDau2 = -1;
    pCand->fSig1 = gHFEvent->nSigTracks();
    pCand->fSig2 = pCand->fSig1 + trackList.size() - 1;
    
    // -- fill (not yet refitted) sig tracks
    TAnaTrack *pTrack; 
    pTrack            = gHFEvent->addSigTrack();
    pTrack->fMCID     = trackList[0].charge()*13; 
    pTrack->fGenIndex = -1; 
    pTrack->fQ        = trackList[0].charge();
    pTrack->fPlab.SetPtEtaPhi(trackList[0].pt(),
			      trackList[0].eta(),
			      trackList[0].phi()
			      ); 
    pTrack->fIndex  = trackIndices[0];

    pTrack            = gHFEvent->addSigTrack();
    pTrack->fMCID     = trackList[1].charge()*13; 
    pTrack->fGenIndex = -1; 
    pTrack->fQ        = trackList[1].charge();
    pTrack->fPlab.SetPtEtaPhi(trackList[1].pt(),
			      trackList[1].eta(),
			      trackList[1].phi()
			      ); 
    pTrack->fIndex  = trackIndices[1];

    for (unsigned int i = 2; i < otherTracks.size(); ++i) {
      pTrack            = gHFEvent->addSigTrack();
      pTrack->fMCID     = trackList[i].charge()*321; 
      pTrack->fGenIndex = -1; 
      pTrack->fQ        = otherTracks[i].charge();
      pTrack->fPlab.SetPtEtaPhi(trackList[i].pt(),
				trackList[i].eta(),
				trackList[i].phi()
				); 
      pTrack->fIndex  = trackIndices[2+i];
    }

  }

  return;

}






