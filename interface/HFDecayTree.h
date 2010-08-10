/*
 *  HFDecayTree.h
 *  HFDecayTree
 *
 *  Created by Christoph on 28.4.10.
 */

#ifndef HFDECAYTREE_H
#define HFDECAYTREE_H

// CMSSW
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/GeometrySurface/interface/ReferenceCounted.h"

// STL
#include <map>
#include <vector>
#include <utility>

// ROOT
#include <TVector3.h>

/* The Functor Cut Object base class */
class HFNodeCut : public ReferenceCounted {
	public:
		HFNodeCut();
		
		void setFields(double maxDoca, double vtxChi2, TVector3 vtxPos, TVector3 ptCand);
		
		/* This function has to be overwritten. True if the particle passes the test */
		virtual bool operator()();
		
        public:
		double fMaxDoca;
		double fVtxChi2;
		TVector3 fVtxPos;
		TVector3 fPtCand;
};
typedef ReferenceCountingPointer<HFNodeCut> RefCountedHFNodeCut;

class HFMaxDocaCut : public HFNodeCut {
	public:
		HFMaxDocaCut(double docaCut) : fDocaCut(docaCut) {}
		virtual bool operator()() { return (fMaxDoca < fDocaCut);}
	protected:
		double fDocaCut;
};

class HFDecayTree;
typedef std::vector<HFDecayTree>::iterator HFDecayTreeIterator;
typedef std::map<int,int>::iterator HFDecayTreeTrackIterator;

class HFDecayTree
{
	public:
		HFDecayTree(int pID = 0.0, int doVertexing = 1, double constraint = -1.0, double constraintSigma = -1.0);
		virtual ~HFDecayTree() { delete kinTree; }
		
		// Constructing the tree structure
		void addTrack(int trackIx, double trackMass) __attribute__((deprecated)); // DEPRECATED, use addTrack(trackIx,trackID) instead!
		void addTrack(int trackIx, int trackID); // Add a track with a given type.
		
		
		void appendDecayTree(HFDecayTree subTree); // to append an already constructed decay tree
		HFDecayTreeIterator addDecayTree(int pID = 0, int doVertexing = 1, double mass = -1.0, double mass_sigma = -1.0); // to get a reference to the subvertex
		
		void clear();
		
		// Accessing the track data
		HFDecayTreeTrackIterator getTrackBeginIterator();
		HFDecayTreeTrackIterator getTrackEndIterator();
		
		HFDecayTreeIterator getVerticesBeginIterator();
		HFDecayTreeIterator getVerticesEndIterator();
		
		void getAllTracks(std::vector<std::pair<int,int> > *out_vector, int onlyThisVertex = 0);
		std::vector<std::pair<int,int> > getAllTracks(int onlyThisVertex = 0);
		
		RefCountedKinematicTree *getKinematicTree();
		void setKinematicTree(RefCountedKinematicTree newTree);
		void resetKinematicTree(int recursive = 0);
		
		// Reconstruction
		TAnaCand *getAnaCand();
		void setAnaCand(TAnaCand *cand);
		
		RefCountedHFNodeCut getNodeCut();
		void setNodeCut(RefCountedHFNodeCut newNodeCut);

		// Debugging!
		void dump(unsigned indent = 0);
		
	public:
		// Tree Variables...
		int vertexing; // do a vertexing at this node
		double particleID; // if == 0, then no TAnaCandidate should be created.
		double massConstraint; // if <= 0, then no massconstraint at this vertex
		double massConstraintSigma;
	private:
		void dumpTabs(unsigned indent); // used by dump()
		
		std::map<int,int> trackIndices; // map: trackIx -> particleTyp
		std::vector<HFDecayTree> subVertices;
		RefCountedKinematicTree *kinTree;
		TAnaCand *anaCand;
		RefCountedHFNodeCut nodeCut;
};

#endif
