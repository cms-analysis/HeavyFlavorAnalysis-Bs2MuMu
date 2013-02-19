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
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

using reco::Vertex;
using reco::Track;
using reco::TransientTrack;
using std::vector;
using std::endl;
using std::cout;

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
void HFKinematicVertexFit::doJpsiFit(vector<Track> &trackList,
				     vector<int> &trackIndices,
				     vector<double> &trackMasses,
				     int type) {

  if (fVerbose > 5) cout << "==>HFKinematicVertexFit> doJpsiFit()" << endl;

  KinematicParticleFactoryFromTransientTrack pFactory;
  KinematicParticleVertexFitter kpvFitter;
  KinematicParticleFitter csFitter;
  KinematicConstraint *jpsi_cons = NULL;
  RefCountedKinematicTree jpTree;
  RefCountedKinematicTree buTree;
  RefCountedKinematicVertex buVertex;
  RefCountedKinematicParticle buParticle;
  ParticleMass jpsi = MJPSI;
  
  TVector3 buPlab;
  double buMass;
  
  // fit for J/Psi (first two entries)
  float jp_m_sigma = 0.00004;
  float muon_sigma = 0.0000001;
  float kaon_sigma = 0.000016;
  
  float chi = 0.0;
  float ndf = 0.0;
	
  vector<TransientTrack> muonTracks; // for J/Psi
  vector<TransientTrack> otherTracks; // for total fit
  vector<RefCountedKinematicParticle> muonParticles; // for J/Psi
  vector<RefCountedKinematicParticle> otherParticles; // for total fit
	
  unsigned int j;

  // create the transient tracks
  muonTracks.push_back(fpTTB->build(trackList[0]));
  muonTracks.push_back(fpTTB->build(trackList[1]));
	
  for (j = 2; j < trackList.size(); j++)
    otherTracks.push_back(fpTTB->build(trackList[j]));
	
  // create the kinematic particles from the transient tracks
  for (j = 0; j < muonTracks.size(); j++)
    muonParticles.push_back(pFactory.particle(muonTracks[j],trackMasses[j],chi,ndf,muon_sigma)); // FIXME: replace muon-sigma!!
  for (j = 0; j < otherTracks.size(); j++)
    otherParticles.push_back(pFactory.particle(otherTracks[j],trackMasses[j+2],chi,ndf,kaon_sigma)); // FIXME: replace kaon-sigma!!
	
  try {
    jpTree = kpvFitter.fit(muonParticles);
    if (jpTree->isEmpty()) {
      if (fVerbose > 0) cout << "==> HFKinematicVertexFitter: empty J/Psi tree!! returning!!!" << endl;
      goto bail;
    }
		
    // setup the mass constrain for the J/Psi
    jpsi_cons = new MassKinematicConstraint(jpsi,jp_m_sigma);
    jpTree = csFitter.fit(jpsi_cons, jpTree);
    if (jpTree->isEmpty()) {
      if(fVerbose > 0) cout << "==> HFKinematicVertexFitter: empty constrained J/Psi tree!! returning!!!" << endl;
      goto bail;
    }
		
    jpTree->movePointerToTheTop();
    otherParticles.push_back(jpTree->currentParticle());
		
    buTree = kpvFitter.fit(otherParticles);
    if (buTree->isEmpty()) {
      if (fVerbose > 0) cout << "==> HFKinematicVertexFitter: empty buTree!! returning!!!" << endl;
      goto bail;
    }
		
    buTree->movePointerToTheTop();
  }
  catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==> HFKinematicVertexFitter: cms exception caught: " << ex.what() << endl;
    goto bail;
  }
  catch (VertexException &ex) {
    if (fVerbose > 0) cout << "==> HFKinematicVertexFitter: vertex exception caught: " << ex.what() << endl;
    goto bail;
  }
		
  // -- Fill the ntuple	
  buVertex = buTree->currentDecayVertex();
  buParticle = buTree->currentParticle();
  buPlab = TVector3(buParticle->currentState().globalMomentum().x(),
		    buParticle->currentState().globalMomentum().y(),
		    buParticle->currentState().globalMomentum().z());
  buMass = buParticle->currentState().mass();
	
  if (buVertex->vertexIsValid() && buParticle->chiSquared() > 0) {

    TAnaVertex anaVt;
    TAnaCand  *pCand;
    TAnaTrack *pTrack;
    VertexDistanceXY axy;
    VertexDistance3D a3d;
		
    ChiSquared chi(buParticle->chiSquared(), buParticle->degreesOfFreedom());

    // dump some information if in verbose mode
    if (fVerbose > 0) {
      cout << "----------------------------------------" << endl;
      cout << "==> HFKinematicVertexFitter: Filling candidate with mass = " << buMass << endl;
      for (j = 0; j < trackList.size(); j++)
	cout << trackList[j].px() << ", "  << trackList[j].py() << ", "  << trackList[j].pz() << endl;
      cout << "----------------------------------------" << endl;
    }
		
    anaVt.setInfo(buParticle->chiSquared(), buParticle->degreesOfFreedom(), chi.probability(), 0, -1);
    anaVt.fPoint.SetXYZ(buVertex->position().x(), buVertex->position().y(), buVertex->position().z() );
		
    // -- Distance to primary vertex
    anaVt.fDxy = axy.distance(fPV, buVertex->vertexState()).value();
    anaVt.fDxyE = axy.distance(fPV, buVertex->vertexState()).error();
		
    anaVt.fD3d = a3d.distance(fPV, buVertex->vertexState()).value();
    anaVt.fD3dE = a3d.distance(fPV, buVertex->vertexState()).error();
		
    // -- fill candidate
    pCand = gHFEvent->addCand();
    pCand->fPlab = buPlab;
    pCand->fMass = buMass;
    pCand->fVtx  = anaVt;
    pCand->fType = (type != 0 ? type: fType);
    pCand->fDau1 = -1;
    pCand->fDau2 = -1;
		
    pCand->fSig1 = gHFEvent->nSigTracks();
    pCand->fSig2 = pCand->fSig1 + trackList.size() - 1;
		
    // -- fill (not yet refitted) sig tracks (FIXME)
    pTrack = gHFEvent->addSigTrack();
    pTrack->fMCID		= trackList[0].charge()*13; // FIXME: not always a muon!!
    pTrack->fGenIndex	= gHFEvent->getRecTrack(trackIndices[0])->fGenIndex;
    pTrack->fQ			= trackList[0].charge();
    pTrack->fPlab.SetPtEtaPhi(trackList[0].pt(), trackList[0].eta(), trackList[0].phi() ); 
    pTrack->fIndex  = trackIndices[0];
		
    pTrack = gHFEvent->addSigTrack();
    pTrack->fMCID		= trackList[1].charge()*13; // FIXME: not always a muon!!
    pTrack->fGenIndex	= gHFEvent->getRecTrack(trackIndices[1])->fGenIndex;
    pTrack->fQ			= trackList[1].charge();
    pTrack->fPlab.SetPtEtaPhi(trackList[1].pt(),trackList[1].eta(),trackList[1].phi() );
    pTrack->fIndex  = trackIndices[1];
		
    for (j = 0; j < otherTracks.size(); j++) {
      pTrack = gHFEvent->addSigTrack();
      pTrack->fMCID		= trackList[j+2].charge()*321; // FIXME: not always a kaon!!
      pTrack->fGenIndex	= gHFEvent->getRecTrack(trackIndices[j+2])->fGenIndex;
      pTrack->fQ			= otherTracks[j].charge();
      pTrack->fPlab.SetPtEtaPhi(trackList[j+2].pt(),trackList[j+2].eta(),trackList[j+2].phi() );
      pTrack->fIndex  = trackIndices[j+2];
    }
  }
	
 bail:
  // clean up!!
  if (fVerbose > 5) cout << "==>HFKinematicVertexFit>    trying to clean up" << endl;
  if (jpsi_cons)
    delete jpsi_cons;
  if (fVerbose > 5) cout << "==>HFKinematicVertexFit>    returning" << endl;
}
