/*
 *  HFDecayTree.cp
 *  HFDecayTree
 *
 *  Created by Christoph on 28.4.10.
 *  Modified by Frank on 31.3.11: added flag for massConstraint, sign of mass no longer determines behaviour.
 *
 */
#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDecayTree.h"

using namespace std;

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

HFDecayTree::HFDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma, bool daughtersToPV) :
  particleID_(pID), vertexing_(doVertexing), mass_(mass), massConstraint_(massConstraint), massSigma_(massSigma), maxDoca_(0), minDoca_(0), daughtersToPV_(daughtersToPV), kinTree_(0)
{
  if(massConstraint && massSigma <= 0.0) massSigma_ = 0.0001 * mass;

  nodeCut_ = RefCountedHFNodeCut(new HFNodeCut);
} // HFDecayTree()

// the following constructor is DEPRECATED, attribute set in header will issue a warning
// initializers adjusted for new data member structure
HFDecayTree::HFDecayTree(int pID, int doVertexing, double constraint, double constraintSigma) :
  particleID_(pID), vertexing_(doVertexing>0?true:false),
  mass_(constraint>=0?constraint:0.0), massConstraint_(constraint<0?false:true), massSigma_(constraintSigma),
  maxDoca_(0), minDoca_(0), kinTree_(0)
{
  if(massConstraint_ && massSigma_ <= 0.0) massSigma_ = 0.0001 * mass_;

  nodeCut_ = RefCountedHFNodeCut(new HFNodeCut);
} // HFDecayTree()

void HFDecayTree::addTrack(int trackIx, int trackID)
{
	trackIndices_[trackIx] = trackID;
} // addTrack()

void HFDecayTree::appendDecayTree(HFDecayTree subTree)
{
  subVertices_.push_back(subTree);
} // appendDecayTree()

HFDecayTreeIterator HFDecayTree::addDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma, bool daughtersToPV)
{
  return subVertices_.insert(subVertices_.end(), HFDecayTree(pID, doVertexing, mass, massConstraint, massSigma, daughtersToPV));
} // addDecayTree()

// this version is DEPRECATED
HFDecayTreeIterator HFDecayTree::addDecayTree(int pID, int doVertexing, double mass, double mass_sigma)
{
  return subVertices_.insert(subVertices_.end(), HFDecayTree(pID, doVertexing>0?true:false, mass, mass<0?false:true, mass_sigma));
} // addDecayTree()

void HFDecayTree::clear()
{
  clear(-1, true, 0.0, false, -1.0);
  // particleID_ = -1.0;
  // vertexing_ = true;
  // mass_ = 0.0;
  // massConstraint_ = false;
  // massSigma_ = -1.0;
}

void HFDecayTree::clear(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma, bool daughtersToPV)
{
  particleID_ = pID;
  vertexing_ = doVertexing;
  mass_ = mass;
  massConstraint_ = massConstraint;
  massSigma_ = massSigma;
  maxDoca_ = -1.0;
  minDoca_ = -1.0;
  daughtersToPV_ = daughtersToPV;

  
  // clear the containers
  trackIndices_.clear();
  kinParticleMap_.clear();
  subVertices_.clear();

  // clear the kinematic tree
  delete kinTree_;
  kinTree_ = NULL;

  anaCand_ = NULL;

  nodeCut_ = RefCountedHFNodeCut(new HFNodeCut);
} // clear()

HFDecayTreeTrackIterator HFDecayTree::getTrackBeginIterator()
{
	return trackIndices_.begin();
} // getTrackBeginIterator()

HFDecayTreeTrackIterator HFDecayTree::getTrackEndIterator()
{
	return trackIndices_.end();
} // getTrackBeginIterator()

HFDecayTreeIterator HFDecayTree::getVerticesBeginIterator()
{
	return subVertices_.begin();
} // getVerticesBeginIterator()

HFDecayTreeIterator HFDecayTree::getVerticesEndIterator()
{
	return subVertices_.end();
} // getVerticesEndIterator()

void HFDecayTree::getAllTracks(vector<pair<int,int> > *out_vector, int onlyThisVertex)
{
	HFDecayTreeTrackIterator trackIt;
	HFDecayTreeIterator treeIt;
	
	for (trackIt = trackIndices_.begin(); trackIt!=trackIndices_.end(); ++trackIt)
		out_vector->push_back(*trackIt);
	
	for (treeIt = subVertices_.begin(); treeIt!=subVertices_.end(); ++treeIt) {
		if (!treeIt->vertexing_ || !onlyThisVertex)
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
	return &kinParticleMap_;
} // getKinParticleMap()

void HFDecayTree::setKinParticleMap(map<int,int> newMap)
{
	kinParticleMap_ = newMap;
} // setKinParticleMap()

RefCountedKinematicTree* HFDecayTree::getKinematicTree()
{
  return kinTree_;
} // getKinematicTree()

void HFDecayTree::setKinematicTree(RefCountedKinematicTree newTree)
{
  if(!kinTree_) kinTree_ = new RefCountedKinematicTree;
  
  *kinTree_ = newTree; // make a copy from the reference counting pointer
} // setKinematicTree()

void HFDecayTree::resetKinematicTree(int recursive)
{
  HFDecayTreeIterator treeIt;

  if(recursive) {
    for(treeIt = getVerticesBeginIterator(); treeIt!=getVerticesEndIterator(); ++treeIt)
      treeIt->resetKinematicTree(recursive);
  }
  
  kinParticleMap_.clear();
  delete kinTree_;
  kinTree_ = NULL;
  anaCand_ = NULL;
} // resetKinematicTree()

TAnaCand *HFDecayTree::getAnaCand()
{
  return anaCand_;
} // getAnaCand()

void HFDecayTree::setAnaCand(TAnaCand *cand)
{
  anaCand_ = cand;
} // setAnaCand()

RefCountedHFNodeCut HFDecayTree::getNodeCut()
{ return nodeCut_; }

void HFDecayTree::setNodeCut(RefCountedHFNodeCut newNodeCut)
{ nodeCut_ = newNodeCut; }

void HFDecayTree::dump(unsigned indent)
{
	HFDecayTreeIterator treeIt;
	HFDecayTreeTrackIterator trackIt;
	
	dumpTabs(indent);
	cout << "HFDecayTree (particleID = " << particleID_ << ", vertexing = " << vertexing_
		<< ", massConstraint = " << massConstraint_ << ", mass = " << mass_ << ", massSigma = " << massSigma_ << ") {" << endl;
	
	dumpTabs(indent+1);
	for (trackIt = trackIndices_.begin(); trackIt!=trackIndices_.end(); ++trackIt)
		cout << '(' << "trackIx = " << trackIt->first << ", trackParticleID = " << trackIt->second << ")\t";
	cout << endl;
	
	for (treeIt = subVertices_.begin(); treeIt != subVertices_.end(); ++treeIt)
		treeIt->dump(indent+1);
	
	dumpTabs(indent);
	cout << '}' << endl;
} // dump()


void HFDecayTree::dumpTabs(unsigned indent)
{
	for (unsigned j = 0; j < indent; j++) cout << '\t';
} // dumpTabs()
