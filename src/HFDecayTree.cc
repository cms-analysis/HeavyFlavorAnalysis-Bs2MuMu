/*
 *  HFDecayTree.cp
 *  HFDecayTree
 *
 *  Created by Christoph on 28.4.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDecayTree.h"

using namespace std;

void HFDecayTree::addTrack(int trackIx, double trackMass)
{
	if (trackMass <= 0)
		cerr << "HFDecayTree: Adding Track with non-positive mass!!!" << endl;
	
	if(!trackIndices.insert(make_pair(trackIx, trackMass)).second)
		cerr << "HFDecayTree: Adding the same Trackindex again. Ignoring!!!" << endl;
} // addTrack()

void HFDecayTree::appendDecayTree(HFDecayTree subTree)
{
	subVertices.push_back(subTree);
} // appendDecayTree()

HFDecayTreeIterator HFDecayTree::addDecayTree()
{
	return subVertices.insert(subVertices.end(),HFDecayTree());
} // addDecayTree()

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


std::vector<int> HFDecayTree::getAllTracks()
{
	vector<int> tracks;
	HFDecayTreeTrackIterator trackIt;
	HFDecayTreeIterator treeIt;
	
	for (trackIt = trackIndices.begin(); trackIt!=trackIndices.end(); ++trackIt)
		tracks.push_back(trackIt->first);
	
	for (treeIt = subVertices.begin(); treeIt!=subVertices.end(); ++treeIt) {
		vector<int> tmp = treeIt->getAllTracks();
		tracks.insert(tracks.end(),tmp.begin(),tmp.end());
	}
	
	return tracks;
} // getAllTracks()

// FIXME: this has to be completed for the new variables introduced.
void HFDecayTree::dump(unsigned indent)
{
	set<pair<int,double> >::const_iterator it;
	HFDecayTreeIterator ptr;
	unsigned j;
	
	for (j = 0; j < indent; j++) cout << '\t';
	cout << "Dumping HFDecayTree (" << massConstraint << ") {" << endl;
	
	// one more then above!!
	for (j = 0; j <= indent; j++) cout << '\t';
	for (it = trackIndices.begin(); it != trackIndices.end(); ++it)
		cout << '(' << it->first << ',' << it->second << ")\t";
	cout << endl;
	
	for (ptr = subVertices.begin(); ptr!=subVertices.end(); ++ptr)
		ptr->dump(indent + 1);
	
	for (j = 0; j < indent; j++) cout << '\t';
	cout << '}' << endl;
} // dump()
