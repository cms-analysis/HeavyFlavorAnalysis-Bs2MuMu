/*
 *  HFSequentialVertexFit.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 29.4.10.
 */

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFSequentialVertexFit.h"
#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFMasses.hh"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include "CommonTools/Statistics/interface/ChiSquared.h"

#include <algorithm>

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
	ImpactParameters(Measurement1D plip, Measurement1D ptip) {
	    lip = plip;
	    tip = ptip;
	}
	Measurement1D lip;
	Measurement1D tip;
};

using namespace std;
using namespace edm;
using namespace reco;

HFSequentialVertexFit::HFSequentialVertexFit(Handle<View<Track> > hTracks, const TransientTrackBuilder *TTB, Handle<VertexCollection> pvCollection, const MagneticField *field, int verbose, bool removeCandTracksFromVtx) :

	fVerbose(verbose),
	fpTTB(TTB),
	fhTracks(hTracks),
	fPVCollection(pvCollection),
	magneticField(field),
	removeCandTracksFromVtx_(removeCandTracksFromVtx)
{} // HFSequentialVertexFit()

HFSequentialVertexFit::~HFSequentialVertexFit()
{} // ~HFSequentialVertexFit()

// if this node does not survive the nodeCut, then it returns false and the fitting sequence stops
bool HFSequentialVertexFit::fitTree(HFDecayTree *tree)
{
	KinematicParticleFactoryFromTransientTrack pFactory;
	vector<RefCountedKinematicParticle> kinParticles;
	vector<track_entry_t>::const_iterator trackIt;
	vector<track_entry_t> allTreeTracks;
	RefCountedKinematicTree kinTree;
	HFDecayTreeIterator treeIt;
	//RefCountedHFNodeCut nodeCut;
	map<int,int> *kinParticleMap;
	int mass_constrained_tracks = 0;
	
	// set up the kinParticleMap for 'tree'
	kinParticleMap = tree->getKinParticleMap();
	
	// add the particles from the tracks, we have to add them first for the KinematicConstrainedVertexFitter
	allTreeTracks = tree->getAllTracks(1);
	sort(allTreeTracks.begin(), allTreeTracks.end()); // sort such that all mass fit tracks are first
	for (trackIt = allTreeTracks.begin(); trackIt != allTreeTracks.end(); ++trackIt) {
		float sigma;
		float mass = getParticleMass(trackIt->particleID,&sigma);
		TrackBaseRef baseRef(fhTracks,trackIt->trackIx);
		(*kinParticleMap)[trackIt->trackIx] = kinParticles.size();
		kinParticles.push_back(pFactory.particle(fpTTB->build(*baseRef),mass,0.0f,0.0f,sigma));
		if (trackIt->massFit) mass_constrained_tracks++;
	}
	
	// add the particles from the sub-tree's with vertexing...
	for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
		
		if(!fitTree(&(*treeIt))) return false; // abort if there was some problem
		
		kinTree = *(treeIt->getKinematicTree());
		if (kinTree->isEmpty()) throw EmptyTreeError();
		
		if (treeIt->vertexing()) {
			kinTree->movePointerToTheTop();
			kinParticles.push_back(kinTree->currentParticle());
		}
	}
	
	// do the actual fit of this vertex
	if (tree->mass_tracks() > 0 && mass_constrained_tracks > 0) {
		KinematicConstrainedVertexFitter kcvFitter;
		auto_ptr<MultiTrackKinematicConstraint> tr_c(new MultiTrackMassKinematicConstraint(tree->mass_tracks(),mass_constrained_tracks));
		kinTree = kcvFitter.fit(kinParticles,&(*tr_c));
	} else {
		KinematicParticleVertexFitter kpvFitter;
		kinTree = kpvFitter.fit(kinParticles);
	}
	
	if (!kinTree->isEmpty() && tree->massConstraint()) {
		KinematicParticleFitter csFitter;
		auto_ptr<KinematicConstraint> con(new MassKinematicConstraint(tree->mass(),tree->massSigma()));
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
		tree->set_maxDoca(maxDoca);
		tree->set_minDoca(getMinDoca(kinParticles));
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
  if(tree->particleID() && !tree->getAnaCand())
    tree->setAnaCand(addCandidate(tree)); // top candidate w.r.t. primary vertex

  // get the current vertex state
  subTree = *(tree->getKinematicTree());
  subTree->movePointerToTheTop();
  vState = subTree->currentDecayVertex()->vertexState();
  
  // create all the requested candidates of the daughters
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
	  if (treeIt->particleID() && !treeIt->getAnaCand())
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

  // override the dxy of the daughters if requested
  if (tree->daughtersToPV()) computeDaughterDistance(tree);

  // recursively continue
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt)
    saveTree(&(*treeIt));

} // saveTree()

void HFSequentialVertexFit::computeDaughterDistance(HFDecayTree *tree)
{
  TAnaCand *mom, *dau;
  HFDecayTreeIterator it;
  VertexDistanceXY axy;
  VertexDistance3D a3d;
  RefCountedKinematicTree dauTree;
  RefCountedKinematicVertex dauVertex;

  // Load the Primary Vertex Collection

  mom = tree->getAnaCand();
  for(it = tree->getVerticesBeginIterator(); it != tree->getVerticesEndIterator(); ++it) {
    dau = it->getAnaCand();
    if (!dau) continue;

    dauTree = *(it->getKinematicTree());
    dauTree->movePointerToTheTop();
    dauVertex = dauTree->currentDecayVertex();

    // Vertex Distanz neu berechnen zum PV: xy
    dau->fVtx.fDxy  = axy.distance( (*fPVCollection)[mom->fPvIdx], dauVertex->vertexState() ).value();
    dau->fVtx.fDxyE = axy.distance( (*fPVCollection)[mom->fPvIdx], dauVertex->vertexState() ).error();

    // Vertex Distanz neu berechnen zum PV: 3d
    dau->fVtx.fD3d  = a3d.distance( (*fPVCollection)[mom->fPvIdx], dauVertex->vertexState() ).value();
    dau->fVtx.fD3dE = a3d.distance( (*fPVCollection)[mom->fPvIdx], dauVertex->vertexState() ).error();
  }
}

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
  vector<track_entry_t> allTreeTracks = tree->getAllTracks(1);
  set<int> allUsedTrackIndices = tree->getAllTracksIndices();
  map<int,int> *kinParticleMap;
  RefCountedKinematicTree kinTree = *(tree->getKinematicTree());
  RefCountedKinematicParticle kinParticle;
  RefCountedKinematicVertex kinVertex;
  vector<RefCountedKinematicParticle> daughterParticles;
  TVector3 plab;
  double cov[9];
  double mass;
  unsigned int j;
  int pvIx = -1, pvIx2 = -1; // PV index of this candidate
  ImpactParameters pvImpParams;
  ImpactParameters pvImpParams2nd(Measurement1D(9999.,9999.),Measurement1D(9999.,9999.)); // stores just tip of second best PV to detect pile-up problems
  AnalyticalImpactPointExtrapolator extrapolator(magneticField);
  TransverseImpactPointExtrapolator transverseExtrapolator(magneticField);
  TrajectoryStateOnSurface tsos;

  if (tree->particleID() == 0) return pCand; // i.e. null
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

  if (!wrtVertexState) {   // do for the candidate (no VertexState)
	  
	  VertexCollection::const_iterator vertexIt;
	  // calculate the impact parameters for all primary vertices
	  if (fVerbose > 0)
		  cout << "==> HFSequentialVertexFit: Number of PV vertices to compare is " << fPVCollection->size() << endl;
	  double pvWeightCut = nodeCut->getPvWeightCut();
	  
	  // iterate through all PVs and estimate the impact parameters of this particle
	  unsigned int nGoodVtx;
	  for (vertexIt = fPVCollection->begin(), j = 0, nGoodVtx = 0; vertexIt != fPVCollection->end(); ++vertexIt,++j) {
		  
		  std::pair<bool,Measurement1D> currentIp;
		  
		  // extrapolate to PCA
		  tsos = extrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),RecoVertex::convertPos(vertexIt->position()));
		  
		  // compute with iptools
		  currentIp = IPTools::signedDecayLength3D(tsos,GlobalVector(0,0,1),*vertexIt);

		  // Compute the PV weight
		  double weight = (vertexIt->ndof()+2.)/(2.*vertexIt->tracksSize());
// 		  cout<<j<<" "<<nGoodVtx<<" "<<vertexIt->position()<<" "<<vertexIt->chi2()<<" "
// 		   <<vertexIt->ndof()<<" "<<vertexIt->tracksSize()<<" "<<vertexIt->nTracks()<<" "
// 		   <<vertexIt->isFake()<<" "<<vertexIt->isValid()<<" "<<vertexIt->normalizedChi2()<<" "<<weight<<" "
// 		   <<currentIp.second.value()<<" "<<pvWeightCut<<endl;

		  if( weight < pvWeightCut) { // Check the PV weight and skip it if lower than the cut 
		    if (fVerbose > 2) cout<<"==>HFSequentialVertexFit: PV "<<j<<" rejected because of too low weight "<<weight<<endl; 
		    continue;
		  }

		  if (!currentIp.first) {
			  if (fVerbose > 0) cout << "==>HFSequentialVertexFit: Unable to compute lip to vertex at index " << j << endl;
			  continue;
		  }
		  
		  // store?
		  if (nGoodVtx == 0 ) // the first PV, this is currently the best one ;)
		  {
		      pvIx = j;
		      pvImpParams.lip = currentIp.second;
		  }
		  else if (nGoodVtx == 1) // now the second PV
		  {
		      if (fabs(currentIp.second.value()) >= fabs(pvImpParams.lip.value())) // not the best but the second best
			  pvImpParams2nd.lip = currentIp.second;
		      else // the best, the previous one is the current 2nd best
		      {
			  pvIx2 = pvIx;
			  pvIx = j;
			  pvImpParams2nd.lip = pvImpParams.lip;
			  pvImpParams.lip = currentIp.second;
		      }
		  }
		  else // we have more than 2 PV
		  {
		      if (fabs(currentIp.second.value()) >= fabs(pvImpParams.lip.value())) // not the best
		      {
			  if (fabs(currentIp.second.value()) < fabs(pvImpParams2nd.lip.value())) // but the second best
			    {pvImpParams2nd.lip = currentIp.second; pvIx2 = j; }
		      }
		      else // this is currently the best one, keep it and put the old best one to 2nd best
		      {
			  pvIx2 = pvIx;
			  pvIx = j;
			  pvImpParams2nd.lip = pvImpParams.lip;
			  pvImpParams.lip = currentIp.second;
		      }
		  }
		  nGoodVtx++; // Count the no. of good vertices
	  }
	  
	  // now, compute the tip w.r.t. PV
	  tsos = transverseExtrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),RecoVertex::convertPos((*fPVCollection)[pvIx].position()));
	  pvImpParams.tip = axy.distance(VertexState(tsos.globalPosition(),tsos.cartesianError().position()),VertexState(RecoVertex::convertPos((*fPVCollection)[pvIx].position()),RecoVertex::convertError((*fPVCollection)[pvIx].error())));
	  if (fVerbose > 2) cout<<"==> HFSequentialVertexFit: Selected best PVs "<<pvIx<<" "<<pvIx2<<endl;
  }
  
  double vtxDistanceCosAlphaPlab(0); // will be used later while calculating the lifetime but easier to determine here.
  cov99_t vtxDistanceCov;
  jac9_t vtxDistanceJac3d, vtxDistanceJac2d;
  if (wrtVertexState) {
	  // -- Distance to mother vertex
	  anaVtx.fDxy = axy.distance(*wrtVertexState, kinVertex->vertexState()).value();
	  anaVtx.fDxyE = axy.distance(*wrtVertexState, kinVertex->vertexState()).error();
	  
	  anaVtx.fD3d = a3d.distance(*wrtVertexState, kinVertex->vertexState()).value();
	  anaVtx.fD3dE = a3d.distance(*wrtVertexState, kinVertex->vertexState()).error();
	  // -- get covariance matrix for error propagation in lifetime calculation
	  vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(wrtVertexState->error()),
		  kinParticle->currentState().kinematicParametersError().matrix());
	  vtxDistanceJac3d = makeJacobianVector3d(wrtVertexState->position(), kinVertex->vertexState().position(), plab);
	  vtxDistanceJac2d = makeJacobianVector2d(wrtVertexState->position(), kinVertex->vertexState().position(), plab);
          // -- get sign of distance
	  const GlobalVector diff = kinVertex->vertexState().position() - wrtVertexState->position() ;
	  const TVector3 tv3diff = TVector3(diff.x(),diff.y(),diff.z());
	  vtxDistanceCosAlphaPlab = plab.Dot(tv3diff) / (plab.Mag() * tv3diff.Mag());
  } else if (pvIx >= 0) {
	  // -- Distance w.r.t primary vertex
	  Vertex currentPV = (*fPVCollection)[pvIx];

	  // refit the vertex without the tracks of the candidate
	  if (removeCandTracksFromVtx_)
	  {
	      vector<TransientTrack> vrtxRefit;
	      vector<track_entry_t> completeTrackList = tree->getAllTracks(0);
	      bool removedTracks(false); // to check if we really need to perform the refit
	      for (std::vector<TrackBaseRef>::const_iterator itTBR = currentPV.tracks_begin();  itTBR != currentPV.tracks_end(); itTBR++)
	      {   // loop over all tracks in the PV to check if tracks from the candidate were used for fittingthe PV
		  TrackRef tref = itTBR->castTo<TrackRef>();
		  bool trkFound(false);
		  for (vector<track_entry_t>::const_iterator trackIt = completeTrackList.begin(); trackIt != completeTrackList.end(); ++trackIt) 
		  {
		      TrackBaseRef curTr(fhTracks,trackIt->trackIx);
		      TrackRef curTref = curTr.castTo<TrackRef>();
		      if (tref == curTref) trkFound = true;
		  }
		  if (!trkFound)
		  {
		      TransientTrack tTrk = fpTTB->build(*(*itTBR));
		      vrtxRefit.push_back(tTrk);
		  }
		  else
		  { // track contained in PV, don't add to list
		      removedTracks = true;
		  }
	      }
	      if (removedTracks)
	      { // so we need to fit a new PV
		  AdaptiveVertexFitter avf;
		  TransientVertex newVtx = avf.vertex(vrtxRefit);
		  if (newVtx.isValid()) currentPV = reco::Vertex(newVtx);
	      }
	  }

	  anaVtx.fDxy = axy.distance(currentPV,kinVertex->vertexState()).value();
	  anaVtx.fDxyE = axy.distance(currentPV,kinVertex->vertexState()).error();
	  
	  anaVtx.fD3d = a3d.distance(currentPV,kinVertex->vertexState()).value();
	  anaVtx.fD3dE = a3d.distance(currentPV,kinVertex->vertexState()).error();

	  // -- get covariance matrix for error propagation in lifetime calculation
	  vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(currentPV.error()),
		  kinParticle->currentState().kinematicParametersError().matrix());
	  vtxDistanceJac3d = makeJacobianVector3d(currentPV.position(), kinVertex->vertexState().position(), plab);
	  vtxDistanceJac2d = makeJacobianVector2d(currentPV.position(), kinVertex->vertexState().position(), plab);
          // -- get sign of distance
	  const TVector3 p1(currentPV.position().x(), currentPV.position().y(), currentPV.position().z());
	  const TVector3 p2(kinVertex->vertexState().position().x(), kinVertex->vertexState().position().y(), kinVertex->vertexState().position().z());
	  const TVector3 pDiff = p2-p1;
	  vtxDistanceCosAlphaPlab = plab.Dot(pDiff) / (plab.Mag() * pDiff.Mag());
	  
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
  pCand->fType = tree->particleID();

  pCand->fMom = -1; // Mom gets linked later.
  pCand->fDau1 = -1; // Daughters get linked later
  pCand->fDau2 = -1;
  
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + allTreeTracks.size() - 1;
  
  pCand->fMaxDoca = tree->maxDoca();
  pCand->fMinDoca = tree->minDoca();
  
  pCand->fPvIdx = pvIx;
  pCand->fPvIdx2 = pvIx2;
  pCand->fPvLip = pvImpParams.lip.value();
  pCand->fPvLipE = pvImpParams.lip.error();
  pCand->fPvTip = pvImpParams.tip.value();
  pCand->fPvTipE = pvImpParams.tip.error();
  pCand->fPvLip2 = pvImpParams2nd.lip.value();
  pCand->fPvLipE2 = pvImpParams2nd.lip.error();
  pCand->fPvTip2 = pvImpParams2nd.tip.value();
  pCand->fPvTipE2 = pvImpParams2nd.tip.error();
  
  // -- calculate lifetime
  {
    // TMath::Ccgs() is to convert from cm to s (speed of light in cgs system, CMS uses cm)
    if (pCand->fPlab.Mag() > 0)
    {
	const double massOverC = tree->mass() / TMath::Ccgs();
	// from 3d vertexing
	pCand->fTau3d = anaVtx.fD3d / pCand->fPlab.Mag() * vtxDistanceCosAlphaPlab * massOverC;
	pCand->fTau3dE = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac3d)) * massOverC;
	//const double tEsimple = anaVtx.fD3dE / pCand->fPlab.Mag() * vtxDistanceCosAlphaPlab * massOverC;
	//cout << tree->particleID() << " t: " << pCand->fTau3d << " tE: " << pCand->fTau3dE << " tE/t: " << pCand->fTau3dE/pCand->fTau3d << 
	//    " tEsimple: " << tEsimple << " tEsimple/tE: " << tEsimple/pCand->fTau3dE << " m/c: " << massOverC << endl;
	// from 2d vertexing
	const double sinTheta = TMath::Sin(pCand->fPlab.Theta());
	const double flightlength2d = sinTheta != 0 ? anaVtx.fDxy / sinTheta : 0;
	pCand->fTauxy = flightlength2d / pCand->fPlab.Mag() * vtxDistanceCosAlphaPlab * massOverC;
	pCand->fTauxyE = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac2d)) * massOverC;
    }
    else
    {
	pCand->fTau3d = pCand->fTauxy = -99.;
    }
    // just for convenience extract the error from cov into const variables (will be optimized away by the compiler)
    // covariance matrix has order (x,y,z,p_x,p_y,p_z,m)
    // see RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h
    const double pxE = kinParticle->currentState().kinematicParametersError().matrix()(3,3);
    const double pyE = kinParticle->currentState().kinematicParametersError().matrix()(4,4);
    const double pzE = kinParticle->currentState().kinematicParametersError().matrix()(5,5);
    pCand->fTau3dE = TMath::Sqrt(anaVtx.fD3dE*anaVtx.fD3dE + pxE*pxE + pyE*pyE + pzE*pzE ) / TMath::Ccgs();
  }

  for (j = 0; j < allTreeTracks.size(); j++) {
    
    TransientTrack fitTrack = daughterParticles[(*kinParticleMap)[allTreeTracks[j].trackIx]]->refittedTransientTrack();
    
    pTrack = gHFEvent->addSigTrack();
    pTrack->fIndex = allTreeTracks[j].trackIx;
    pTrack->fMCID = allTreeTracks[j].particleID; // Here, we use the MCID of the sigTrack to store the assumed particle ID for the mass hypothesis
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

// -------------------------------------------------------------
HFSequentialVertexFit::cov33_t HFSequentialVertexFit::GlobalError2SMatrix_33(GlobalError m_in)
{
    cov33_t m_out;
    for(int i=0; i!=3; i++)
	for(int j=i; j!=3; j++)
	{
	    m_out(i,j) = m_in.matrix()(i+1,j+1);
	}
    return m_out;
}

HFSequentialVertexFit::cov99_t
    HFSequentialVertexFit::makeCovarianceMatrix(const HFSequentialVertexFit::cov33_t cov_vtx1,
						const HFSequentialVertexFit::cov77_t cov_vtx2)
{
    cov99_t cov;
    cov.Place_at(cov_vtx1,0,0);
    cov.Place_at(cov_vtx2.Sub<cov66_t>(0,0),3,3);
    return cov;
}

// -------------------------------------------------------------
HFSequentialVertexFit::jac9_t HFSequentialVertexFit::makeJacobianVector3d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2, const AlgebraicVector3 &momentum)
{
    HFSequentialVertexFit::jac9_t jac;
    const AlgebraicVector3 dist = vtx2 - vtx1;
    const double factor2 = 1. / ROOT::Math::Mag2(momentum);
    const double lifetime = ROOT::Math::Dot(dist, momentum) * factor2;
    jac.Place_at(-momentum*factor2,0);
    jac.Place_at( momentum*factor2,3);
    jac.Place_at( factor2*(dist-2*lifetime*momentum*factor2),6);
    return jac;
}

HFSequentialVertexFit::jac9_t HFSequentialVertexFit::makeJacobianVector3d(const GlobalPoint &vtx1, const GlobalPoint &vtx2, const TVector3 &tv3momentum)
{   // TODO: Update 2d calculation to projected version as in 3d
    return makeJacobianVector3d(AlgebraicVector3(vtx1.x(),vtx1.y(),vtx1.z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}

HFSequentialVertexFit::jac9_t HFSequentialVertexFit::makeJacobianVector3d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
	const GlobalPoint &vtx2, const TVector3 &tv3momentum)
{
    return makeJacobianVector3d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}

// -------------------------------------------------------------
HFSequentialVertexFit::jac9_t HFSequentialVertexFit::makeJacobianVector2d(const GlobalPoint &vtx1, const GlobalPoint &vtx2, const TVector3 &tv3momentum)
{
    return makeJacobianVector2d(AlgebraicVector3(vtx1.x(),vtx1.y(),vtx1.z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}

HFSequentialVertexFit::jac9_t HFSequentialVertexFit::makeJacobianVector2d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
	const GlobalPoint &vtx2, const TVector3 &tv3momentum)
{
    return makeJacobianVector2d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}


HFSequentialVertexFit::jac9_t HFSequentialVertexFit::makeJacobianVector2d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2, const AlgebraicVector3 &momentum)
{
    HFSequentialVertexFit::jac9_t jac;
    const double momentumMag = ROOT::Math::Mag(momentum);
    const AlgebraicVector3 dist = vtx2 - vtx1;
    const double distMag = ROOT::Math::Mag(dist);
    const double factorPositionComponent = 1./(distMag*momentumMag);
    const double factorMomentumComponent = 1./pow(momentumMag,3);
    jac(0)=-dist(0)*factorPositionComponent;
    jac(1)=-dist(1)*factorPositionComponent;
    jac(3)= dist(0)*factorPositionComponent;
    jac(4)= dist(1)*factorPositionComponent;
    jac(6)= momentum(0)*factorMomentumComponent;
    jac(7)= momentum(1)*factorMomentumComponent;
    return jac;
}

