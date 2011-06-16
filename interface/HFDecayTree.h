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
		HFDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false);
		HFDecayTree(int pID, int doVertexing, double constraint, double constraintSigma = -1.0) __attribute__((deprecated)); // DEPRECATED use new constructor instead
		virtual ~HFDecayTree() { delete kinTree_; }
		
		// Constructing the tree structure
		void addTrack(int trackIx, int trackID); // Add a track with a given type.
		
		void appendDecayTree(HFDecayTree subTree); // to append an already constructed decay tree
		HFDecayTreeIterator addDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false); // to get a reference to the subvertex
		//HFDecayTreeIterator addDecayTree(int pID = 0, int doVertexing = 1, double mass = -1.0, double mass_sigma = -1.0) __attribute__ ((deprecated)); // to get a reference to the subvertex
		HFDecayTreeIterator addDecayTree(int pID, int doVertexing, double mass, double mass_sigma = -1.0) __attribute__ ((deprecated)); // to get a reference to the subvertex
		
		void clear();
		void clear(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false); // variant to clear and initialize the tree with same signature as constructor
		
		// Accessing the track data
		HFDecayTreeTrackIterator getTrackBeginIterator();
		HFDecayTreeTrackIterator getTrackEndIterator();
		
		HFDecayTreeIterator getVerticesBeginIterator();
		HFDecayTreeIterator getVerticesEndIterator();
		
		void getAllTracks(std::vector<std::pair<int,int> > *out_vector, int onlyThisVertex = 0);
		std::vector<std::pair<int,int> > getAllTracks(int onlyThisVertex = 0);
		std::set<int> getAllTracksIndices(int onlyThisVertex = 0);
		
		// Kinematic Tree associated stuff
		std::map<int,int> *getKinParticleMap();
		void setKinParticleMap(std::map<int,int> newMap);
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

		// Accessors to previously public data members
		bool vertexing() { return vertexing_; };
		double particleID() { return particleID_; };
		bool massConstraint() { return massConstraint_; };
		double mass() { return mass_; };
		double massSigma() { return massSigma_; };
		double maxDoca() { return maxDoca_; };
		double minDoca() { return minDoca_; };
		bool daughtersToPV() { return daughtersToPV_; }
		
		void set_vertexing(bool vertexing) { vertexing_ = vertexing; };
		void set_particleID(double particleID) { particleID_ = particleID; };
		void set_massConstraint(bool massConstraint) { massConstraint_ = massConstraint; };
		void set_mass(double mass) { mass_ = mass; };
		void set_massSigma(double massSigma) { massSigma_ = massSigma; };
		void set_maxDoca(double maxDoca) { maxDoca_ = maxDoca; };
		void set_minDoca(double minDoca) { minDoca_ = minDoca; };
		void set_daughtersToPV(bool daughtersToPV) { daughtersToPV_ = daughtersToPV; }

	private:
		// Tree Variables...
		double particleID_; // if == 0, then no TAnaCandidate should be created.
		bool vertexing_; // do a vertexing at this node
		double mass_;
		bool massConstraint_; // false: no massconstraint at this vertex
		double massSigma_;
		double maxDoca_;
		double minDoca_;
		bool daughtersToPV_;

		void dumpTabs(unsigned indent); // used by dump()
		
		std::map<int,int> trackIndices_; // map: trackIx -> particleTyp
		std::map<int,int> kinParticleMap_; // map: trackIx -> entry in the daughter kinematic particles...
		std::vector<HFDecayTree> subVertices_;
		RefCountedKinematicTree *kinTree_;
		TAnaCand *anaCand_;
		RefCountedHFNodeCut nodeCut_;
};

#endif
