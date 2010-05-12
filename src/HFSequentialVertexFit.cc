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

void HFSequentialVertexFit::fitTree(HFDecayTree *tree)
{
  KinematicParticleFactoryFromTransientTrack pFactory;
  KinematicParticleVertexFitter kpvFitter;
  vector<RefCountedKinematicParticle> kinParticles;
  RefCountedKinematicTree kinTree;
  HFDecayTreeTrackIterator trackIt;
  HFDecayTreeIterator treeIt;
  double mass;

  // add the particles from the tracks
  for (trackIt = tree->getTrackBeginIterator(); trackIt != tree->getTrackEndIterator(); trackIt++) {
    mass = getParticleMass(trackIt->second);     
    float sigma = 0.00001*mass;
    TrackBaseRef baseRef(fhTracks,(*trackIt).first);
    kinParticles.push_back(pFactory.particle(fpTTB->build(*baseRef), mass, 0.0f, 0.0f, sigma)); // FIXME: das sigma noch besser anpassen (!!)
  }

  // add the particles from the sub tree
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); treeIt++) {
    fitTree(&(*treeIt));
    kinTree = *(treeIt->getKinematicTree());
    if (kinTree->isEmpty()) throw EmptyTreeError();
    
    kinTree->movePointerToTheTop();
    kinParticles.push_back(kinTree->currentParticle());
  }

  // do the actual fit
  kinTree = kpvFitter.fit(kinParticles);
  if(!kinTree->isEmpty() && tree->massConstraint >= 0) {
    KinematicParticleFitter csFitter;
    auto_ptr<KinematicConstraint> con(new MassKinematicConstraint(tree->massConstraint,tree->massConstraintSigma));
    kinTree = csFitter.fit(&(*con),kinTree);
  }

  tree->setKinematicTree(kinTree);
} // fitTree()

void HFSequentialVertexFit::saveTree(HFDecayTree *tree)
{
  int dau1 = -1,dau2 = -1;
  TAnaCand *pCand,*pMomCand;
  HFDecayTreeIterator treeIt;
  RefCountedKinematicTree subTree;
  VertexState vState;

  // create the Ana Candidate of the node if requested and not yet existing
  if(tree->particleID && !tree->getAnaCand())
    tree->setAnaCand(addCandidate(tree,fPV)); // top candidate w.r.t. primary vertex

  // iterate the daughters
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
    if(treeIt->particleID && !treeIt->getAnaCand()) {
      subTree = (*(tree->getKinematicTree()));
      subTree->movePointerToTheTop();
      vState = subTree->currentDecayVertex()->vertexState();
      treeIt->setAnaCand(addCandidate(&(*treeIt),vState));
    }
  }

  pMomCand = tree->getAnaCand();
  if(pMomCand) {
    // now link the daughters
    for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
      pCand = treeIt->getAnaCand();
      if(pCand) {
	// set mother
	pCand->fMom = pMomCand->fVar2;

	if (dau1 == -1) dau1 = pCand->fVar2;
	else            dau1 = (pCand->fVar2 < dau1) ? pCand->fVar2 : dau1;
	
	if (dau2 == -1) dau2 = pCand->fVar2;
	else            dau2 = (pCand->fVar2 > dau2) ? pCand->fVar2 : dau2;
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
  TAnaTrack *pTrack,*recTrack;
  VertexDistanceXY axy;
  VertexDistance3D a3d;
  vector<pair<int,int> > allTracks = tree->getAllTracks();
  RefCountedKinematicTree kinTree = *(tree->getKinematicTree());
  RefCountedKinematicParticle kinParticle;
  RefCountedKinematicVertex kinVertex;
  TVector3 plab;
  double mass;
  unsigned j;

  if (tree->particleID == 0) return pCand; // i.e. null
  if (kinTree->isEmpty()) return pCand;
  
  kinTree->movePointerToTheTop();
  kinParticle = kinTree->currentParticle();
  kinVertex = kinTree->currentDecayVertex();

  if (!kinVertex->vertexIsValid()) return pCand;
  
  plab = TVector3(kinParticle->currentState().globalMomentum().x(),
		  kinParticle->currentState().globalMomentum().y(),
		  kinParticle->currentState().globalMomentum().z());
  mass = kinParticle->currentState().mass();
  
  if (kinParticle->chiSquared() < 0) return pCand; 
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
  
  // -- fill candidate
  pCand = gHFEvent->addCand();
  pCand->fVar2 = gHFEvent->nCands()-1; // pointer to index. required for linking of daughters.
  pCand->fPlab = plab;
  pCand->fMass = mass;
  pCand->fVtx = anaVtx;
  pCand->fType = tree->particleID;

  pCand->fMom = -1; // Mom gets linked later.
  pCand->fDau1 = -1; // Daughters get linked later
  pCand->fDau2 = -1;
  
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + allTracks.size() - 1;

  // FIXME: calculate doca's
  // FIXME: take refitted tracks!!
  
  for (j = 0; j < allTracks.size(); j++) {
    recTrack = gHFEvent->getRecTrack(allTracks[j].first);
    
    pTrack = gHFEvent->addSigTrack();
    pTrack->fMCID = allTracks[j].second;
    pTrack->fGenIndex = recTrack->fGenIndex;
    pTrack->fQ = recTrack->fQ;
    pTrack->fPlab = recTrack->fPlab;
    pTrack->fIndex = allTracks[j].first;
  }
  
  return pCand;
} // addCandidate()

void HFSequentialVertexFit::doFit(HFDecayTree *tree)
{
  RefCountedKinematicTree kinTree;
  RefCountedKinematicVertex kinVertex;
  RefCountedKinematicParticle kinParticle;
  TVector3 plab;
  
  try {
    tree->resetKinematicTree(1);
    fitTree(tree);
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
