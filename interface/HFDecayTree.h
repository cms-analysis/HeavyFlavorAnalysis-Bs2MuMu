/*
 *  HFDecayTree.h
 *  HFDecayTree
 *
 *  Created by Christoph on 28.4.10.
 */

#ifndef HFDECAYTREE_H
#define HFDECAYTREE_H

// STL
#include <vector>
#include <set>
#include <utility>

struct lt_track
{
	bool operator()(const std::pair<int,double> a, const std::pair<int,double> b) const {
		return a.first < b.first;
	}
};

class HFDecayTree;
typedef std::vector<HFDecayTree>::iterator HFDecayTreeIterator;
typedef std::set<std::pair<int,double>,lt_track>::iterator HFDecayTreeTrackIterator;

class HFDecayTree
{
	public:
                HFDecayTree() : massConstraint(-1), massConstraintSigma(-1) {}
		HFDecayTree(double constraint, double constraintSigma) : massConstraint(constraint),massConstraintSigma(constraintSigma) {}
		virtual ~HFDecayTree() {}
		
		// Constructing the tree structure
		void addTrack(int trackIx, double trackMass); // add a track with a given mass to this node. NOTE: does not maintain the order of insertion!!
		void appendDecayTree(HFDecayTree subTree); // to append an already constructed decay tree
		HFDecayTreeIterator addDecayTree(); // to get a reference to a new tree to be constructed
		
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
		std::set<std::pair<int,double>,lt_track> trackIndices;
		std::vector<HFDecayTree> subVertices;
};

#endif
