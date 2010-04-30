/*
 *  HFSequentialVertexFit.cc
 *
 *  Created by Christoph Nägeli on 29.4.10.
 */

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"

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

using namespace std;
using namespace edm;
using namespace reco;

HFSequentialVertexFit::HFSequentialVertexFit(Handle<View<Track> > hTracks, const TransientTrackBuilder *TTB, Vertex &PV, int type, int verbose) :
  fType(type),fVerbose(verbose),fPV(PV),fpTTB(TTB),fhTracks(hTracks)
{ } // HFSequentialVertexFit()

HFSequentialVertexFit::~HFSequentialVertexFit()
{} // ~HFSequentialVertexFit()

RefCountedKinematicTree HFSequentialVertexFit::fitTree(HFDecayTree *tree)
{
	KinematicParticleFactoryFromTransientTrack pFactory;
	KinematicParticleVertexFitter kpvFitter;
	RefCountedKinematicTree kinTree;
	vector<RefCountedKinematicParticle> kinParticles;
	HFDecayTreeTrackIterator it;
	HFDecayTreeIterator treeIt;
	
	// add the particles from the tracks
	for (it = tree->getTrackBeginIterator(); it != tree->getTrackEndIterator(); it++) {
	        float sigma = 0.00001*(*it).second;
		TrackBaseRef baseRef(fhTracks,(*it).first);
		kinParticles.push_back(pFactory.particle(fpTTB->build(*baseRef),(*it).second, 0.0f, 0.0f, sigma)); // FIXME: das sigma noch besser anpassen (!!)
	}
	
	// add the particles from the sub_tree
	for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); treeIt++) {
		RefCountedKinematicTree subParticleTree = fitTree(&(*treeIt));
		if (subParticleTree->isEmpty())
			throw EmptyTreeError();
		
		subParticleTree->movePointerToTheTop();
		kinParticles.push_back(subParticleTree->currentParticle());
	}

	kinTree = kpvFitter.fit(kinParticles);
	
	// apply the kinematic constaint.
	if (!kinTree->isEmpty() && tree->massConstraint >= 0) {
		KinematicParticleFitter csFitter;
		auto_ptr<KinematicConstraint> con(new MassKinematicConstraint(tree->massConstraint,tree->massConstraintSigma));
		kinTree = csFitter.fit(&(*con),kinTree);
	}
	
	return kinTree;
} // fitTree()

void HFSequentialVertexFit::doFit(HFDecayTree *tree, int type)
{
	RefCountedKinematicTree kinTree;
	RefCountedKinematicVertex kinVertex;
	RefCountedKinematicParticle kinParticle;
	TVector3 plab;
	double mass;
	
	try {
		kinTree = fitTree(tree);		
	} catch (cms::Exception &ex) {
		if (fVerbose > 0) cout << "==> HFSequentialVertexFit: cms exception caught: " << ex.what() << endl;
	} catch (VertexException &ex) {
		if (fVerbose > 0) cout << "==> HFSequentialVertexFit: vertex exception caught: " << ex.what() << endl;
	} catch (EmptyTreeError& ex) {
		if (fVerbose > 0) cout << "==> HFSequentialVertexFit: empty tree." << endl;
	}
	
	// FIXME: Check that this is really empty in case fitTree did raise an exception!!
	if (!kinTree->isEmpty() && kinTree->currentDecayVertex()->vertexIsValid()) {
		
		// Fill the n-tuple
		TAnaVertex anaVt;
		TAnaCand *pCand;
		TAnaTrack *pTrack;
		VertexDistanceXY axy;
		VertexDistance3D a3d;
		vector<int> allTracks = tree->getAllTracks();
		unsigned j;
		
		kinTree->movePointerToTheTop();
		
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
		pCand->fType = (type != 0) ? type : fType;
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
	
} // doFit()
