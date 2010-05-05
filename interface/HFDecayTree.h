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
typedef std::map<int,int>::iterator HFDecayTreeTrackIterator;

class HFDecayTree
{
	public:
		HFDecayTree(int pID = -1.0, double constraint = -1.0, double constraintSigma = -1.0);
		virtual ~HFDecayTree() {}
		
		// Constructing the tree structure
		void addTrack(int trackIx, double trackMass) __attribute__((deprecated)); // DEPRECATED, use addTrack(trackIx,trackID) instead!
		void addTrack(int trackIx, int trackID); // Add a track with a given type.
		
		void appendDecayTree(HFDecayTree subTree); // to append an already constructed decay tree
		HFDecayTreeIterator addDecayTree(int pID = 0, double mass = -1.0, double mass_sigma = -1.0); // to get a reference to a new tree to be constructed
		
		void clear();

		// Accessing the track data
		HFDecayTreeTrackIterator getTrackBeginIterator();
		HFDecayTreeTrackIterator getTrackEndIterator();
		
		HFDecayTreeIterator getVerticesBeginIterator();
		HFDecayTreeIterator getVerticesEndIterator();
		
		void getAllTracks(std::vector<int> *out_vector);
		std::vector<int> getAllTracks();
		
		// Debugging!
		void dump(unsigned indent = 0);
		
		double particleID; // if < 0, then no TAnaCandidate should be created.
		double massConstraint; // if < 0, then no massconstraint at this vertex
		double massConstraintSigma;
	private:
		void dumpTabs(unsigned indent); // used by dump()
		
		std::map<int,int> trackIndices;
		std::vector<HFDecayTree> subVertices;
};

#endif
