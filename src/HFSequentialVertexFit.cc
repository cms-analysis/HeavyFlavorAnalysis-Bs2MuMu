/*
 *  HFSequentialVertexFit.cc
 *
 *  Created by Christoph NŠgeli on 29.4.10.
 */

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "CommonTools/Statistics/interface/ChiSquared.h"

extern TAna01Event *gHFEvent;

struct EmptyTreeError {
  EmptyTreeError() {}
};
struct MassNotFoundException {
  MassNotFoundException() {}
};

using namespace std;
using namespace edm;
using namespace reco;

HFSequentialVertexFit::HFSequentialVertexFit(Handle<View<Track> > hTracks, const TransientTrackBuilder *TTB, Vertex &PV, int verbose) :
  fVerbose(verbose),fPV(PV),fpTTB(TTB),fhTracks(hTracks)
{ } // HFSequentialVertexFit()

HFSequentialVertexFit::~HFSequentialVertexFit()
{} // ~HFSequentialVertexFit()

RefCountedKinematicTree HFSequentialVertexFit::fitTree(HFDecayTree *tree)
{
	KinematicParticleFactoryFromTransientTrack pFactory;
	KinematicParticleVertexFitter kpvFitter;
	RefCountedKinematicTree kinTree;
	RefCountedKinematicTree subTree;
	vector<RefCountedKinematicParticle> kinParticles;
	HFDecayTreeTrackIterator trackIt;
	HFDecayTreeIterator treeIt;
	map<HFDecayTreeIterator, RefCountedKinematicTree> candidateMap;
	map<HFDecayTreeIterator, RefCountedKinematicTree>::iterator mapIt;
	double mass;
	
	// add the particles from the tracks
	for (trackIt = tree->getTrackBeginIterator(); trackIt != tree->getTrackEndIterator(); trackIt++) {
	  mass = getParticleMass(trackIt->second);     
	  float sigma = 0.00001*mass;
	  TrackBaseRef baseRef(fhTracks,(*trackIt).first);
	  kinParticles.push_back(pFactory.particle(fpTTB->build(*baseRef), mass, 0.0f, 0.0f, sigma)); // FIXME: das sigma noch besser anpassen (!!)
	}
	
	// add the particles from the sub_tree.
	for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); treeIt++) {
		subTree = fitTree(&(*treeIt));
		if (subTree->isEmpty()) throw EmptyTreeError();

		subTree->movePointerToTheTop();
		kinParticles.push_back(subTree->currentParticle());
		
		// create the map!
		candidateMap[treeIt] = subTree;
	}
	
	// now, add the candidates using the candidate_map.
	// this ensures that the daughter candidates are in a connected range.
	for (mapIt = candidateMap.begin(); mapIt!=candidateMap.end(); ++mapIt) {
		
		treeIt = mapIt->first;
		subTree = mapIt->second;
		
		// add a candidate to the Event if the vertex is valid
		if(treeIt->particleID && subTree->currentDecayVertex()->vertexIsValid()) {
			
			TAnaVertex anaVt;
			TAnaCand *pCand;
			TAnaTrack *pTrack;
			VertexDistanceXY axy;
			VertexDistance3D a3d;
			TVector3 plab;
			
			RefCountedKinematicParticle particle;
			RefCountedKinematicVertex kinVertex;
			
			vector<int> allTracks = treeIt->getAllTracks();
			double mass;
			unsigned j;
			
			particle = subTree->currentParticle();
			kinVertex = subTree->currentDecayVertex();
			
			plab = TVector3(particle->currentState().globalMomentum().x(),
							particle->currentState().globalMomentum().y(),
							particle->currentState().globalMomentum().z());
			mass = particle->currentState().mass();
			
			ChiSquared chi(particle->chiSquared(),particle->degreesOfFreedom());
			// dump some information if in verbose mode.
			if (fVerbose > 0) {
				cout << "-----------------------------------------" << endl;
				cout << "==> HFSequentialVertexFit: Filling candidate with mass = " << mass << endl;
				cout << "-----------------------------------------" << endl;
			}
			
			anaVt.setInfo(particle->chiSquared(),particle->degreesOfFreedom(),chi.probability(), 0, 0);
			anaVt.fPoint.SetXYZ(kinVertex->position().x(),kinVertex->position().y(),kinVertex->position().z());
			
			// -- Distance to primary vertex
			anaVt.fDxy = axy.distance(fPV, kinVertex->vertexState()).value();
			anaVt.fDxyE = axy.distance(fPV, kinVertex->vertexState()).error();
			
			anaVt.fD3d = a3d.distance(fPV, kinVertex->vertexState()).value();
			anaVt.fD3dE = a3d.distance(fPV, kinVertex->vertexState()).error();
			
			// -- fill the candidate
			pCand = gHFEvent->addCand();
			pCand->fPlab = plab;
			pCand->fMass = mass;
			pCand->fVtx = anaVt;
			pCand->fType = treeIt->particleID;
			pCand->fDau1 = -1; // FIXME: Set here the daughter range!!
			pCand->fDau2 = -1;
			
			pCand->fSig1 = gHFEvent->nSigTracks();
			pCand->fSig2 = pCand->fSig1 + allTracks.size() - 1;
			
			// FIXME: there are still some variables missing in TAnaCand!!!
			
			// -- fill sig tracks. FIXME: fill refitted sig tracks
			for (j = 0; j < allTracks.size(); j++) {
				TAnaTrack *recTrack = gHFEvent->getRecTrack(allTracks[j]);
				
				pTrack = gHFEvent->addSigTrack();
				// FIXME: adjust fMCID to the particle we think it is...
				pTrack->fMCID = recTrack->fMCID;
				pTrack->fGenIndex = recTrack->fGenIndex;
				pTrack->fQ = recTrack->fQ;
				pTrack->fPlab = recTrack->fPlab;
				pTrack->fIndex = allTracks[j];
				// FIXME: fill all the other parameters
			}
		}
	}
	
	kinTree = kpvFitter.fit(kinParticles);
	
	// apply the kinematic constraint.
	if (!kinTree->isEmpty() && tree->massConstraint >= 0) {
		KinematicParticleFitter csFitter;
		auto_ptr<KinematicConstraint> con(new MassKinematicConstraint(tree->massConstraint,tree->massConstraintSigma));
		kinTree = csFitter.fit(&(*con),kinTree);
	}
	
	return kinTree;
} // fitTree()

void HFSequentialVertexFit::doFit(HFDecayTree *tree)
{
  RefCountedKinematicTree kinTree;
  RefCountedKinematicVertex kinVertex;
  RefCountedKinematicParticle kinParticle;
  TVector3 plab;
  double mass;
  
  try {
    kinTree = fitTree(tree);
    kinTree->movePointerToTheTop();
    
    if (tree->particleID && !kinTree->isEmpty() && kinTree->currentDecayVertex()->vertexIsValid()) {
      
      // Fill the n-tuple
      TAnaVertex anaVt;
      TAnaCand *pCand;
      TAnaTrack *pTrack;
      VertexDistanceXY axy;
      VertexDistance3D a3d;
      vector<int> allTracks = tree->getAllTracks();
      unsigned j;
      
      kinParticle = kinTree->currentParticle();
      kinVertex = kinTree->currentDecayVertex();
      plab = TVector3(kinParticle->currentState().globalMomentum().x(),
		      kinParticle->currentState().globalMomentum().y(),
		      kinParticle->currentState().globalMomentum().z());
      mass = kinParticle->currentState().mass();
      
      ChiSquared chi(kinParticle->chiSquared(),kinParticle->degreesOfFreedom());
      
      // dump some information if in verbose mode.
      if (fVerbose > 0) {
		cout << "-----------------------------------------" << endl;
		cout << "==> HFSequentialVertexFit: Filling candidate with mass = " << mass << endl;
		cout << "-----------------------------------------" << endl;
      }
      
      anaVt.setInfo(kinParticle->chiSquared(), kinParticle->degreesOfFreedom(), chi.probability(), 0, 0);
      anaVt.fPoint.SetXYZ(kinVertex->position().x(),kinVertex->position().y(),kinVertex->position().z());
      
      // -- Distance to primary vertex
      anaVt.fDxy = axy.distance(fPV, kinVertex->vertexState()).value();
      anaVt.fDxyE = axy.distance(fPV, kinVertex->vertexState()).error();
      
      anaVt.fD3d = a3d.distance(fPV, kinVertex->vertexState()).value();
      anaVt.fD3dE = a3d.distance(fPV, kinVertex->vertexState()).error();
      
      // -- fill candidate
      pCand = gHFEvent->addCand();
      pCand->fPlab = plab;
      pCand->fMass = mass;
      pCand->fVtx = anaVt;
      pCand->fType = tree->particleID;
      pCand->fDau1 = -1; // FIXME: This could be set if the iterative fit would add candidates, too!
      pCand->fDau2 = -1;
      
      pCand->fSig1 = gHFEvent->nSigTracks();
      pCand->fSig2 = pCand->fSig1 + allTracks.size() - 1;
      
      // FIXME: there are still some variable missing in TAnaCand!!!
      
      // -- fill sig tracks. FIXME: fill refitted sig tracks
      // FIXME: anderen parameter?
      for (j = 0; j < allTracks.size(); j++) {
		TAnaTrack *recTrack = gHFEvent->getRecTrack(allTracks[j]);

		pTrack = gHFEvent->addSigTrack();
		pTrack->fMCID = recTrack->fMCID;
		pTrack->fGenIndex = recTrack->fGenIndex;
		pTrack->fQ = recTrack->fQ;
		pTrack->fPlab = recTrack->fPlab;
		pTrack->fIndex = allTracks[j];
      }
      
    } else if (fVerbose > 0) cout << "==> HFSequentialVertexFit: empty tree or invalid vertex." << endl;
    
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==> HFSequentialVertexFit: cms exception caught: " << ex.what() << endl;
  } catch (VertexException &ex) {
    if (fVerbose > 0) cout << "==> HFSequentialVertexFit: vertex exception caught: " << ex.what() << endl;
  } catch (EmptyTreeError& ex) {
    if (fVerbose > 0) cout << "==> HFSequentialVertexFit: empty tree." << endl;
  }
} // doFit()

double HFSequentialVertexFit::getParticleMass(int particleID)
{
  double mass;
  particleID = abs(particleID);

  switch(particleID) {
  case 11: // electron
    mass = MELECTRON;
    break;
  case 13: // muon
    mass = MMUON;
    break;
  case 211: // pion
    mass = MPION;
    break;
  case 321: // kaon
    mass = MKAON;
    break;
  case 2212: // proton
    mass = MPROTON;
    break;
  default:
    throw MassNotFoundException();
    break;
  }

  return mass;
} // getParticleMass()
