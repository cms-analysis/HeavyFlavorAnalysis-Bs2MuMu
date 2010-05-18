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

HFDecayTree::HFDecayTree(int pID, double constraint, double constraintSigma) : particleID(pID), massConstraint(constraint), massConstraintSigma(constraintSigma), kinTree(0)
{
  if(massConstraintSigma <= 0.0 && massConstraint > 0.0)
    massConstraintSigma = 0.0001 * massConstraint;
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

HFDecayTreeIterator HFDecayTree::addDecayTree(int pID, double mass, double mass_sigma)
{
  return subVertices.insert(subVertices.end(),HFDecayTree(pID,mass,mass_sigma));
} // addDecayTree()

void HFDecayTree::clear()
{
  particleID = -1.0;
  massConstraint = -1.0;
  massConstraintSigma = -1.0;
  
  // clear the containers
  trackIndices.clear();
  subVertices.clear();

  // clear the kinematic tree
  delete kinTree;
  kinTree = NULL;

  anaCand = NULL;
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

void HFDecayTree::getAllTracks(vector<pair<int,int> > *out_vector)
{
	HFDecayTreeTrackIterator trackIt;
	HFDecayTreeIterator treeIt;
	
	for (trackIt = trackIndices.begin(); trackIt!=trackIndices.end(); ++trackIt)
		out_vector->push_back(*trackIt);
	
	for (treeIt = subVertices.begin(); treeIt!=subVertices.end(); ++treeIt)
		treeIt->getAllTracks(out_vector);
} // getAllTracks()

vector<pair<int,int> > HFDecayTree::getAllTracks()
{
  vector<pair<int,int> > tracks;
  getAllTracks(&tracks);
  return tracks;
} // getAllTracks()

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

void HFDecayTree::dump(unsigned indent)
{
	HFDecayTreeIterator treeIt;
	HFDecayTreeTrackIterator trackIt;
	
	dumpTabs(indent);
	cout << "HFDecayTree (particleID = " << particleID << ", mass = "
		<< massConstraint << ", massSigma = " << massConstraintSigma << ") {" << endl;
	
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
