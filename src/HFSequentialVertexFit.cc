/*
 *  HFSequentialVertexFit.cc
 *
 *  Created by Christoph NŠgeli <christoph.naegeli@psi.ch> on 29.4.10.
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

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

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

// if this node does not survive the nodeCut, then it returns false and the fitting sequence stops
bool HFSequentialVertexFit::fitTree(HFDecayTree *tree)
{
	KinematicParticleFactoryFromTransientTrack pFactory;
	KinematicParticleVertexFitter kpvFitter;
	vector<RefCountedKinematicParticle> kinParticles;
	vector<pair<int,int> >::const_iterator trackIt;
	vector<pair<int,int> > allTreeTracks;
	RefCountedKinematicTree kinTree;
	HFDecayTreeIterator treeIt;
	RefCountedHFNodeCut nodeCut;
	map<int,int> *kinParticleMap;
	double mass;
	
	// set up the kinParticleMap for 'tree'
	kinParticleMap = tree->getKinParticleMap();
	
	// add the particles from the sub-tree's with vertexing...
	for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
		
		if(!fitTree(&(*treeIt))) return false; // abort if there was some problem
		
		kinTree = *(treeIt->getKinematicTree());
		if (kinTree->isEmpty()) throw EmptyTreeError();
		
		if (treeIt->vertexing) {
			kinTree->movePointerToTheTop();
			kinParticles.push_back(kinTree->currentParticle());
		}
	}
	
	// add the particles from the tracks
	allTreeTracks = tree->getAllTracks(1);
	for (trackIt = allTreeTracks.begin(); trackIt != allTreeTracks.end(); ++trackIt) {
		mass = getParticleMass(trackIt->second);
		float sigma = 0.00001*mass;
		TrackBaseRef baseRef(fhTracks,trackIt->first);
		(*kinParticleMap)[trackIt->first] = kinParticles.size();
		kinParticles.push_back(pFactory.particle(fpTTB->build(*baseRef),mass,0.0f,0.0f,sigma)); // FIXME: das sigma noch besser anpassen
	}
	
	// do the actual fit of this vertex
	kinTree = kpvFitter.fit(kinParticles);
	if (!kinTree->isEmpty() && tree->massConstraint >= 0) {
		KinematicParticleFitter csFitter;
		auto_ptr<KinematicConstraint> con(new MassKinematicConstraint(tree->massConstraint,tree->massConstraintSigma));
		kinTree = csFitter.fit(&(*con),kinTree);
	}
	
	tree->setKinematicTree(kinTree);
	
	// set the node cut variables
	nodeCut = tree->getNodeCut();
	{
		// initialize the node variables
		double maxDoca;
		double vtxChi2;
		TVector3 vtxPos;
		TVector3 ptCand;
		
		RefCountedKinematicVertex kinVertex;
		RefCountedKinematicParticle kinPart;
		
		kinTree->movePointerToTheTop();
		kinPart = kinTree->currentParticle();
		kinVertex = kinTree->currentDecayVertex();
		
		maxDoca = getMaxDoca(kinParticles);
		vtxChi2 = kinPart->chiSquared();
		vtxPos.SetXYZ(kinVertex->position().x(),kinVertex->position().y(),kinVertex->position().z());
		ptCand.SetXYZ(kinPart->currentState().globalMomentum().x(),
					  kinPart->currentState().globalMomentum().y(),
					  kinPart->currentState().globalMomentum().z());
		
		nodeCut->setFields(maxDoca, vtxChi2, vtxPos, ptCand);
	}
	
	return (*nodeCut)();
} // fitTree()

void HFSequentialVertexFit::saveTree(HFDecayTree *tree)
{
  int dau1 = -1,dau2 = -1;
  TAnaCand *pCand,*pMomCand;
  HFDecayTreeIterator treeIt;
  RefCountedKinematicTree subTree;
  VertexState vState; // vertex state of 'tree'

  // create the Ana Candidate of the node if requested and not yet existing
  if(tree->particleID && !tree->getAnaCand())
    tree->setAnaCand(addCandidate(tree,fPV)); // top candidate w.r.t. primary vertex

  // get the current vertex state
  subTree = *(tree->getKinematicTree());
  subTree->movePointerToTheTop();
  vState = subTree->currentDecayVertex()->vertexState();
  
  // create all the requested candidates of the daughters
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
	  if (treeIt->particleID && !treeIt->getAnaCand())
		  treeIt->setAnaCand(addCandidate(&(*treeIt),vState));
  }
  
  // link the candidates
  pMomCand = tree->getAnaCand();
  if(pMomCand) {
    // now link the daughters
    for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
      pCand = treeIt->getAnaCand();
      if(pCand) {
    	// set mother
		pCand->fMom = pMomCand->fIndex;

     	if (dau1 == -1) dau1 = pCand->fIndex;
		else            dau1 = (pCand->fIndex < dau1) ? pCand->fIndex : dau1;
	
		if (dau2 == -1) dau2 = pCand->fIndex;
		else            dau2 = (pCand->fIndex > dau2) ? pCand->fIndex : dau2;
      }
    }

    pMomCand->fDau1 = dau1;
    pMomCand->fDau2 = dau2;
  }

  // recursively continue
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt)
    saveTree(&(*treeIt));

} // saveTree()

template<class T>
TAnaCand *HFSequentialVertexFit::addCand(HFDecayTree *tree, T &toVertex)
{
  TAnaCand *pCand = NULL;
  TAnaVertex anaVtx;
  TAnaTrack *pTrack;
  VertexDistanceXY axy;
  VertexDistance3D a3d;
  vector<pair<int,int> > allTreeTracks = tree->getAllTracks(1);
  map<int,int> *kinParticleMap;
  RefCountedKinematicTree kinTree = *(tree->getKinematicTree());
  RefCountedKinematicParticle kinParticle;
  RefCountedKinematicVertex kinVertex;
  vector<RefCountedKinematicParticle> daughterParticles;
  TVector3 plab;
  double cov[9];
  double mass;
  unsigned j;

  if (tree->particleID == 0) return pCand; // i.e. null
  if (kinTree->isEmpty()) return pCand;
  
  kinTree->movePointerToTheTop();
  kinParticle = kinTree->currentParticle();
  kinVertex = kinTree->currentDecayVertex();
  daughterParticles = kinTree->daughterParticles();

  kinParticleMap = tree->getKinParticleMap();

  if (!kinVertex->vertexIsValid()) return pCand;
  
  plab = TVector3(kinParticle->currentState().globalMomentum().x(),
		  kinParticle->currentState().globalMomentum().y(),
		  kinParticle->currentState().globalMomentum().z());
  mass = kinParticle->currentState().mass();
  
  if (kinParticle->chiSquared() < 0)
	  return pCand;

  ChiSquared chi(kinParticle->chiSquared(),kinParticle->degreesOfFreedom());

  // dump some information if in verbose mode...
  if (fVerbose > 0) {
    cout << "-----------------------------------------" << endl;
    cout << "==> HFSequentialVertexFit: Filling candidate with mass = " << mass << endl;
    cout << "-----------------------------------------" << endl;
  }
  
  anaVtx.setInfo(kinParticle->chiSquared(),kinParticle->degreesOfFreedom(),chi.probability(),0,0);
  anaVtx.fPoint.SetXYZ(kinVertex->position().x(),kinVertex->position().y(),kinVertex->position().z());

  // -- Distance to mother vertex (or primary vertex)
  anaVtx.fDxy = axy.distance(toVertex, kinVertex->vertexState()).value();
  anaVtx.fDxyE = axy.distance(toVertex, kinVertex->vertexState()).error();
  
  anaVtx.fD3d = a3d.distance(toVertex, kinVertex->vertexState()).value();
  anaVtx.fD3dE = a3d.distance(toVertex, kinVertex->vertexState()).error();
  
  // -- set covariance matrix
  cov[0] = kinVertex->error().cxx();
  cov[1] = kinVertex->error().cyx();
  cov[2] = kinVertex->error().czx();
  cov[3] = kinVertex->error().cyx();
  cov[4] = kinVertex->error().cyy();
  cov[5] = kinVertex->error().czy();
  cov[6] = kinVertex->error().czx();
  cov[7] = kinVertex->error().czy();
  cov[8] = kinVertex->error().czz();
  anaVtx.setCovXX(cov);
  
  // -- fill candidate
  pCand = gHFEvent->addCand();
  pCand->fPlab = plab;
  pCand->fMass = mass;
  pCand->fVtx = anaVtx;
  pCand->fType = tree->particleID;

  pCand->fMom = -1; // Mom gets linked later.
  pCand->fDau1 = -1; // Daughters get linked later
  pCand->fDau2 = -1;
  
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + allTreeTracks.size() - 1;

  // FIXME: max / min Doca calculation should NOT use the refitted tracks!!
  pCand->fMaxDoca = getMaxDoca(daughterParticles);
  pCand->fMinDoca = getMinDoca(daughterParticles);
  
  for (j = 0; j < allTreeTracks.size(); j++) {
    
    TransientTrack fitTrack = daughterParticles[(*kinParticleMap)[allTreeTracks[j].first]]->refittedTransientTrack();
    
    pTrack = gHFEvent->addSigTrack();
    pTrack->fIndex = allTreeTracks[j].first;
    pTrack->fMCID = allTreeTracks[j].second;
    pTrack->fPlab = TVector3(fitTrack.track().px(),fitTrack.track().py(),fitTrack.track().pz());
    pTrack->fDof = fitTrack.ndof();
    pTrack->fValidHits = fitTrack.numberOfValidHits();
    pTrack->fChi2 = fitTrack.chi2();
    pTrack->fQ = fitTrack.charge();
  }
  
  return pCand;
} // addCandidate()

void HFSequentialVertexFit::doFit(HFDecayTree *tree)
{
  if (fVerbose > 5) cout << "==>HFSequentialVertexFit> doFit()" << endl;
  
  try {
    tree->resetKinematicTree(1);
	if(fitTree(tree))
		saveTree(tree);
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

double HFSequentialVertexFit::getMaxDoca(vector<RefCountedKinematicParticle> &kinParticles)
{
	double maxDoca = -1.0;
	TwoTrackMinimumDistance md;
	vector<RefCountedKinematicParticle>::iterator in_it, out_it;
	
	for (out_it = kinParticles.begin(); out_it != kinParticles.end(); ++out_it) {
		for (in_it = out_it + 1; in_it != kinParticles.end(); ++in_it) {
			md.calculate((*out_it)->currentState().freeTrajectoryState(),(*in_it)->currentState().freeTrajectoryState());
			if (md.distance() > maxDoca)
				maxDoca = md.distance();
		}
	}
	
	return maxDoca;
} // getMaxDoca()

double HFSequentialVertexFit::getMinDoca(vector<RefCountedKinematicParticle> &kinParticles)
{
  double minDoca = 99999.9;
  TwoTrackMinimumDistance md;
  unsigned j,k,n;

  n = kinParticles.size();
  for (j = 0; j < n; j++) {
    for (k = j+1; k < n; k++) {
      md.calculate(kinParticles[j]->currentState().freeTrajectoryState(),kinParticles[k]->currentState().freeTrajectoryState());
      if (md.distance() < minDoca)
	minDoca = md.distance();
    }
  }

  return minDoca;
} // getMinDoca()
