/*
 *  HFDecayTree.cp
 *  HFDecayTree
 *
 *  Created by Christoph on 28.4.10.
 *
 */
#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDecayTree.h"

using namespace std;

// constants used to decide which particle depending on mass.
const static double PROTON_BOUND	=	0.750;
const static double KAON_BOUND		=	0.250;
const static double PION_BOUND		=	0.120;
const static double MUON_BOUND		=	0.050;

HFNodeCut::HFNodeCut() : fMaxDoca(0.0), fVtxChi2(0.0), fVtxPos(), fPtCand()
{} // HFNodeCut()

void HFNodeCut::setFields(double maxDoca, double vtxChi2, TVector3 vtxPos, TVector3 ptCand)
{
	fMaxDoca = maxDoca;
	fVtxChi2 = vtxChi2;
	fVtxPos = vtxPos;
	fPtCand = ptCand;
} // setFields()

bool HFNodeCut::operator()() {return true;}

HFDecayTree::HFDecayTree(int pID, int doVertexing, double constraint, double constraintSigma) :
  vertexing(doVertexing),particleID(pID),massConstraint(constraint),massConstraintSigma(constraintSigma),maxDoca(0),minDoca(0),kinTree(0)
{
  if(massConstraintSigma <= 0.0 && massConstraint > 0.0)
    massConstraintSigma = 0.0001 * massConstraint;

  nodeCut = RefCountedHFNodeCut(new HFNodeCut);
} // HFDecayTree()

void HFDecayTree::addTrack(int trackIx, double trackMass)
{
	int type;
	
	// try to figure out the particle ID by it's mass.
	// NOTE: cannot distinguish between particle and antiparticle (=> deprecated)
	if		(trackMass >= PROTON_BOUND)	type = 2212;	// p+
	else if	(trackMass >= KAON_BOUND)	type = 321;		// K+
	else if (trackMass >= PION_BOUND)	type = 211;		// Pi+
	else if (trackMass >= MUON_BOUND)	type = 13;		// mu-
	else								type = 11;		// e-
	
	addTrack(trackIx, type);
} // addTrack()

void HFDecayTree::addTrack(int trackIx, int trackID)
{
	trackIndices[trackIx] = trackID;
} // addTrack()

void HFDecayTree::appendDecayTree(HFDecayTree subTree)
{
  subVertices.push_back(subTree);
} // appendDecayTree()

HFDecayTreeIterator HFDecayTree::addDecayTree(int pID, int doVertexing, double mass, double mass_sigma)
{
  return subVertices.insert(subVertices.end(),HFDecayTree(pID,doVertexing,mass,mass_sigma));
} // addDecayTree()

void HFDecayTree::clear()
{
  vertexing = 1;
  particleID = -1.0;
  massConstraint = -1.0;
  massConstraintSigma = -1.0;
  maxDoca = -1.0;
  minDoca = -1.0;
  
  // clear the containers
  trackIndices.clear();
  kinParticleMap.clear();
  subVertices.clear();

  // clear the kinematic tree
  delete kinTree;
  kinTree = NULL;

  anaCand = NULL;

  nodeCut = RefCountedHFNodeCut(new HFNodeCut);
} // clear()

HFDecayTreeTrackIterator HFDecayTree::getTrackBeginIterator()
{
	return trackIndices.begin();
} // getTrackBeginIterator()

HFDecayTreeTrackIterator HFDecayTree::getTrackEndIterator()
{
	return trackIndices.end();
} // getTrackBeginIterator()

HFDecayTreeIterator HFDecayTree::getVerticesBeginIterator()
{
	return subVertices.begin();
} // getVerticesBeginIterator()

HFDecayTreeIterator HFDecayTree::getVerticesEndIterator()
{
	return subVertices.end();
} // getVerticesEndIterator()

void HFDecayTree::getAllTracks(vector<pair<int,int> > *out_vector, int onlyThisVertex)
{
	HFDecayTreeTrackIterator trackIt;
	HFDecayTreeIterator treeIt;
	
	for (trackIt = trackIndices.begin(); trackIt!=trackIndices.end(); ++trackIt)
		out_vector->push_back(*trackIt);
	
	for (treeIt = subVertices.begin(); treeIt!=subVertices.end(); ++treeIt) {
		if (!treeIt->vertexing || !onlyThisVertex)
			treeIt->getAllTracks(out_vector,onlyThisVertex);
	}
} // getAllTracks()

vector<pair<int,int> > HFDecayTree::getAllTracks(int onlyThisVertex)
{
  vector<pair<int,int> > tracks;
  getAllTracks(&tracks,onlyThisVertex);
  return tracks;
} // getAllTracks()

set<int> HFDecayTree::getAllTracksIndices(int onlyThisVertex)
{
	vector<pair<int,int> > tracks;
	set<int> result;
	getAllTracks(&tracks,onlyThisVertex);
	
	for(vector<pair<int,int> >::const_iterator it = tracks.begin(); it != tracks.end();++it)
		result.insert(it->first);
	
	return result;
} // getAllTracksIndices()

map<int,int> *HFDecayTree::getKinParticleMap()
{
	return &kinParticleMap;
} // getKinParticleMap()

void HFDecayTree::setKinParticleMap(map<int,int> newMap)
{
	kinParticleMap = newMap;
} // setKinParticleMap()

RefCountedKinematicTree* HFDecayTree::getKinematicTree()
{
  return kinTree;
} // getKinematicTree()

void HFDecayTree::setKinematicTree(RefCountedKinematicTree newTree)
{
  if(!kinTree) kinTree = new RefCountedKinematicTree;
  
  *kinTree = newTree; // make a copy from the reference counting pointer
} // setKinematicTree()

void HFDecayTree::resetKinematicTree(int recursive)
{
  HFDecayTreeIterator treeIt;

  if(recursive) {
    for(treeIt = getVerticesBeginIterator(); treeIt!=getVerticesEndIterator(); ++treeIt)
      treeIt->resetKinematicTree(recursive);
  }
  
  kinParticleMap.clear();
  delete kinTree;
  kinTree = NULL;
  anaCand = NULL;
} // resetKinematicTree()

TAnaCand *HFDecayTree::getAnaCand()
{
  return anaCand;
} // getAnaCand()

void HFDecayTree::setAnaCand(TAnaCand *cand)
{
  anaCand = cand;
} // setAnaCand()

RefCountedHFNodeCut HFDecayTree::getNodeCut()
{ return nodeCut; }

void HFDecayTree::setNodeCut(RefCountedHFNodeCut newNodeCut)
{ nodeCut = newNodeCut; }

void HFDecayTree::dump(unsigned indent)
{
	HFDecayTreeIterator treeIt;
	HFDecayTreeTrackIterator trackIt;
	
	dumpTabs(indent);
	cout << "HFDecayTree (particleID = " << particleID << ", vertexing = " << vertexing
		<< ", mass = " << massConstraint << ", massSigma = " << massConstraintSigma << ") {" << endl;
	
	dumpTabs(indent+1);
	for (trackIt = trackIndices.begin(); trackIt!=trackIndices.end(); ++trackIt)
		cout << '(' << "trackIx = " << trackIt->first << ", trackParticleID = " << trackIt->second << ")\t";
	cout << endl;
	
	for (treeIt = subVertices.begin(); treeIt != subVertices.end(); ++treeIt)
		treeIt->dump(indent+1);
	
	dumpTabs(indent);
	cout << '}' << endl;
} // dump()


void HFDecayTree::dumpTabs(unsigned indent)
{
	for (unsigned j = 0; j < indent; j++) cout << '\t';
} // dumpTabs()
