/*
 *  HFDecayTree.h
 *  HFDecayTree
 *
 *  Created by Christoph on 28.4.10.
 */

#ifndef HFDECAYTREE_H
#define HFDECAYTREE_H

// STL
#include <map>
#include <vector>
#include <utility>

class HFDecayTree;
typedef std::vector<HFDecayTree>::iterator HFDecayTreeIterator;
typedef std::map<int,double>::iterator HFDecayTreeTrackIterator;

class HFDecayTree
{
	public:
		HFDecayTree(double constraint = -1.0, double constraintSigma = -1.0);
		virtual ~HFDecayTree() {}
		
		// Constructing the tree structure
		void addTrack(int trackIx, double trackMass); // add a track with a given mass to this node. NOTE: does not maintain the order of insertion!!
		void appendDecayTree(HFDecayTree subTree); // to append an already constructed decay tree
		HFDecayTreeIterator addDecayTree(double mass = -1.0, double mass_sigma = -1.0); // to get a reference to a new tree to be constructed
		
		void clear();

		// Accessing the track data
		HFDecayTreeTrackIterator getTrackBeginIterator();
		HFDecayTreeTrackIterator getTrackEndIterator();
		
		HFDecayTreeIterator getVerticesBeginIterator();
		HFDecayTreeIterator getVerticesEndIterator();
		
		std::vector<int> getAllTracks();
		
		// Debugging!
		void dump(unsigned indent = 0);
		
		double massConstraint; // if < 0, then no massconstraint at this vertex
		double massConstraintSigma;
	private:
		std::map<int,double> trackIndices;
		std::vector<HFDecayTree> subVertices;
};

#endif
