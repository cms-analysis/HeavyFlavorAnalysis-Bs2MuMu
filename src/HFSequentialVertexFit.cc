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
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include "CommonTools/Statistics/interface/ChiSquared.h"

const static unsigned kNBR_CLOSED_TRACKS = 20;

extern TAna01Event *gHFEvent;

struct EmptyTreeError {
  EmptyTreeError() {}
};
struct MassNotFoundException {
  MassNotFoundException() {}
};

struct ImpactParameters {
	ImpactParameters() {
		lip = Measurement1D();
		tip = Measurement1D();
	}
	Measurement1D lip;
	Measurement1D tip;
};

using namespace std;
using namespace edm;
using namespace reco;

HFSequentialVertexFit::HFSequentialVertexFit(Handle<View<Track> > hTracks, const TransientTrackBuilder *TTB, Handle<VertexCollection> pvCollection, const MagneticField *field, int verbose) :

	fVerbose(verbose),
	fpTTB(TTB),
	fhTracks(hTracks),
	fPVCollection(pvCollection),
	magneticField(field)
{} // HFSequentialVertexFit()

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
		float sigma;
		float mass = getParticleMass(trackIt->second,&sigma);
		TrackBaseRef baseRef(fhTracks,trackIt->first);
		(*kinParticleMap)[trackIt->first] = kinParticles.size();
		kinParticles.push_back(pFactory.particle(fpTTB->build(*baseRef),mass,0.0f,0.0f,sigma));
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
		
		maxDoca = tree->maxDoca = getMaxDoca(kinParticles);
		tree->minDoca = getMinDoca(kinParticles);
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
    tree->setAnaCand(addCandidate(tree)); // top candidate w.r.t. primary vertex

  // get the current vertex state
  subTree = *(tree->getKinematicTree());
  subTree->movePointerToTheTop();
  vState = subTree->currentDecayVertex()->vertexState();
  
  // create all the requested candidates of the daughters
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
	  if (treeIt->particleID && !treeIt->getAnaCand())
		  treeIt->setAnaCand(addCandidate(&(*treeIt),&vState));
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

// Utility routine to sort the mindoca array of the candidate...
static bool doca_less(pair<int,pair<double,double> > x,pair<int,pair<double,double> > y)
{ return x.second.first < y.second.first; } // doca_less()


TAnaCand *HFSequentialVertexFit::addCandidate(HFDecayTree *tree, VertexState *wrtVertexState)
{
  TAnaCand *pCand = NULL;
  TAnaVertex anaVtx;
  TAnaTrack *pTrack;
  VertexDistanceXY axy;
  VertexDistance3D a3d;
  vector<pair<int,int> > allTreeTracks = tree->getAllTracks(1);
  set<int> allUsedTrackIndices = tree->getAllTracksIndices();
  map<int,int> *kinParticleMap;
  RefCountedKinematicTree kinTree = *(tree->getKinematicTree());
  RefCountedKinematicParticle kinParticle;
  RefCountedKinematicVertex kinVertex;
  vector<RefCountedKinematicParticle> daughterParticles;
  TVector3 plab;
  double cov[9];
  double mass;
  unsigned j;
  int pvIx = -1; // PV index of this candidate
  ImpactParameters pvImpParams;
  AnalyticalImpactPointExtrapolator extrapolator(magneticField);
  TransverseImpactPointExtrapolator transverseExtrapolator(magneticField);
  TrajectoryStateOnSurface tsos;

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

  if (!wrtVertexState) {
	  
	  VertexCollection::const_iterator vertexIt;
	  // calculate the impact parameters for all primary vertices
	  if (fVerbose > 0)
		  cout << "==> HFSequentialVertexFit: Number of PV vertices to compare is " << fPVCollection->size() << endl;
	  
	  // iterate through all PVs and estimate the impact parameters of this particle
	  for (vertexIt = fPVCollection->begin(), j = 0; vertexIt != fPVCollection->end(); ++vertexIt,++j) {
		  
		  std::pair<bool,Measurement1D> currentIp;
		  
		  // extrapolate to PCA
		  tsos = extrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),RecoVertex::convertPos(vertexIt->position()));
		  
		  // compute with iptools
		  currentIp = IPTools::signedDecayLength3D(tsos,GlobalVector(0,0,1),*vertexIt);
		  if (!currentIp.first) {
			  if (fVerbose > 0) cout << "==>HFSequentialVertexFit: Unable to compute lip to vertex at index " << j << endl;
			  continue;
		  }
		  
		  // store?
		  if (pvIx >= 0 && fabs(currentIp.second.value()) >= fabs(pvImpParams.lip.value())) continue;
		  
		  pvIx = j;
		  pvImpParams.lip = currentIp.second;
	  }
	  
	  // now, compute the tip w.r.t. PV
	  tsos = transverseExtrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),RecoVertex::convertPos((*fPVCollection)[pvIx].position()));
	  pvImpParams.tip = axy.distance(VertexState(tsos.globalPosition(),tsos.cartesianError().position()),VertexState(RecoVertex::convertPos((*fPVCollection)[pvIx].position()),RecoVertex::convertError((*fPVCollection)[pvIx].error())));
  }
  
  if (wrtVertexState) {
	  // -- Distance to mother vertex
	  anaVtx.fDxy = axy.distance(*wrtVertexState, kinVertex->vertexState()).value();
	  anaVtx.fDxyE = axy.distance(*wrtVertexState, kinVertex->vertexState()).error();
	  
	  anaVtx.fD3d = a3d.distance(*wrtVertexState, kinVertex->vertexState()).value();
	  anaVtx.fD3dE = a3d.distance(*wrtVertexState, kinVertex->vertexState()).error();
  } else if (pvIx >= 0) {
	  // -- Distance w.r.t primary vertex
	  Vertex currentPV = (*fPVCollection)[pvIx];
	  anaVtx.fDxy = axy.distance(currentPV,kinVertex->vertexState()).value();
	  anaVtx.fDxyE = axy.distance(currentPV,kinVertex->vertexState()).error();
	  
	  anaVtx.fD3d = a3d.distance(currentPV,kinVertex->vertexState()).value();
	  anaVtx.fD3dE = a3d.distance(currentPV,kinVertex->vertexState()).error();
	  
  } else if (fVerbose > 0)
	  cout << "==> HFSequentialVertexFit: No idea what distance to compute in TAnaVertex.fDxy and TAnaVertex.fD3d" << endl;
  
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
  
  pCand->fMaxDoca = tree->maxDoca;
  pCand->fMinDoca = tree->minDoca;
  
  pCand->fPvIdx = pvIx;
  pCand->fPvLip = pvImpParams.lip.value();
  pCand->fPvLipE = pvImpParams.lip.error();
  pCand->fPvTip = pvImpParams.tip.value();
  pCand->fPvTipE = pvImpParams.tip.error();
  
  for (j = 0; j < allTreeTracks.size(); j++) {
    
    TransientTrack fitTrack = daughterParticles[(*kinParticleMap)[allTreeTracks[j].first]]->refittedTransientTrack();
    
    pTrack = gHFEvent->addSigTrack();
    pTrack->fIndex = allTreeTracks[j].first;
    pTrack->fMCID = allTreeTracks[j].second; // Here, we use the MCID of the sigTrack to store the assumed particle ID for the mass hypothesis
    pTrack->fPlab = TVector3(fitTrack.track().px(),fitTrack.track().py(),fitTrack.track().pz());
    pTrack->fDof = fitTrack.ndof();
    pTrack->fValidHits = fitTrack.numberOfValidHits();
    pTrack->fChi2 = fitTrack.chi2();
    pTrack->fQ = fitTrack.charge();
  }
  
  // fill the closest approaching tracks -- only if this is supposed to be a SV
  if (!wrtVertexState) {
	  for (j = 0; j < fhTracks->size(); j++) {
		  
		  if (allUsedTrackIndices.count(j)>0) continue; // this tracks belongs to the candidate
		  
		  TrackBaseRef baseRef(fhTracks,j);
		  TransientTrack transTrack = fpTTB->build(*baseRef);
		  
		  tsos = extrapolator.extrapolate(transTrack.initialFreeState(),kinVertex->position());
		  
		  // measure the distance...
		  Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(),tsos.cartesianError().position()),kinVertex->vertexState());
		  
		  // add it to the candidate...
		  pCand->fNstTracks.push_back(make_pair(j,make_pair(doca.value(),doca.error())));
	  }
	  
	  // sort the vector & keep only the first ten
	  sort(pCand->fNstTracks.begin(),pCand->fNstTracks.end(),doca_less);
	  if (pCand->fNstTracks.size() > kNBR_CLOSED_TRACKS) // erase the elements if bigger than kNBR_CLOSED_TRACKS
		  pCand->fNstTracks.erase(pCand->fNstTracks.begin() + kNBR_CLOSED_TRACKS,pCand->fNstTracks.end());
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

float HFSequentialVertexFit::getParticleMass(int particleID, float *mass_sigma)
{
  float mass;
  float sigma = 0.0;
  particleID = abs(particleID);
  
  // sigma corresponds to standard uncertainty as can be found in the PDG
  switch(particleID) {
  case 11: // electron
    mass = MELECTRON;
	sigma = 0.013E-9f;
    break;
  case 13: // muon
    mass = MMUON;
	sigma = 4E-9f;
    break;
  case 211: // pion
    mass = MPION;
	sigma = 3.5E-7f;
    break;
  case 321: // kaon
    mass = MKAON;
	sigma = 1.6E-5f;
    break;
  case 2212: // proton
    mass = MPROTON;
	sigma = 8E-8f;
    break;
  default:
    throw MassNotFoundException();
    break;
  }
  
  if (mass_sigma) *mass_sigma = sigma;
  
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
