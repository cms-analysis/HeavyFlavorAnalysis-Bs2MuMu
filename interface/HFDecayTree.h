/*
 *  HFDecayTree.h
 *  HFDecayTree
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 28.4.10.
 */

#ifndef HFDECAYTREE_H
#define HFDECAYTREE_H

// CMSSW
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/GeometrySurface/interface/ReferenceCounted.h"

// STL
#include <set>
#include <map>
#include <vector>
#include <utility>

// ROOT
#include <TVector3.h>

struct track_entry_t {
	track_entry_t(int ix, int pid, bool mfit) : trackIx(ix), particleID(pid), massFit(mfit) {}
	int trackIx;
	int particleID;
	bool massFit;
};
bool operator<(const track_entry_t &t1, const track_entry_t &t2);

/* The Functor Cut Object base class */
class HFNodeCut : public ReferenceCounted {
	public:
		HFNodeCut();
		
		void setFields(double maxDoca, double vtxChi2, TVector3 vtxPos, TVector3 ptCand);
		
		/* This function has to be overwritten. True if the particle passes the test */
		virtual bool operator()();
		virtual double getPvWeightCut(void) {return 0.0;}
		
	public:
		double fMaxDoca;
		double fVtxChi2;
		//double fPvWeight;
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
class HFPvWeightCut : public HFNodeCut {
	public:
                HFPvWeightCut(double docaCut, double pvWeightCut) : fDocaCut(docaCut), fPvWeightCut(pvWeightCut) {}
		//virtual bool operator()() { return ( (fMaxDoca < fDocaCut) || (fPvWeight < fPvWeightCut) );}
		virtual bool operator()() { return ( fMaxDoca < fDocaCut );}
		virtual double getPvWeightCut(void) {return fPvWeightCut;}
	protected:
		double fDocaCut, fPvWeightCut;
};

class HFDecayTree;
typedef std::vector<HFDecayTree>::iterator HFDecayTreeIterator;
typedef std::set<track_entry_t>::iterator HFDecayTreeTrackIterator;

class HFDecayTree
{
	public:
		HFDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false);
		virtual ~HFDecayTree() { delete kinTree_;}
		
		// Constructing the tree structure
		void addTrack(int trackIx, int trackID, bool massFit = true); // Add a track with a given type and massFit
		
		void appendDecayTree(HFDecayTree subTree); // to append an already constructed decay tree
		HFDecayTreeIterator addDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false); // to get a reference to the subvertex
		
		void clear();
		void clear(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false); // variant to clear and initialize the tree with same signature as constructor
		
		// Accessing the track data
		HFDecayTreeTrackIterator getTrackBeginIterator();
		HFDecayTreeTrackIterator getTrackEndIterator();
		
		HFDecayTreeIterator getVerticesBeginIterator();
		HFDecayTreeIterator getVerticesEndIterator();
		
		void getAllTracks(std::vector<track_entry_t> *out_vector, int onlyThisVertex = 0);
		std::vector<track_entry_t> getAllTracks(int onlyThisVertex = 0);
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
		double mass_tracks() { return mass_tracks_; }
		double massSigma() { return massSigma_; };
		double maxDoca() { return maxDoca_; };
		double minDoca() { return minDoca_; };
		bool daughtersToPV() { return daughtersToPV_; }
		
		void set_vertexing(bool vertexing) { vertexing_ = vertexing; };
		void set_particleID(double particleID) { particleID_ = particleID; };
		void set_massConstraint(bool massConstraint) { massConstraint_ = massConstraint; };
		void set_mass(double mass) { mass_ = mass; };
		void set_mass_tracks(double mass_tracks) { mass_tracks_ = mass_tracks; }
		void set_massSigma(double massSigma) { massSigma_ = massSigma; };
		void set_maxDoca(double maxDoca) { maxDoca_ = maxDoca; };
		void set_minDoca(double minDoca) { minDoca_ = minDoca; };
		void set_daughtersToPV(bool daughtersToPV) { daughtersToPV_ = daughtersToPV; }

	private:
		// Tree Variables...
		double particleID_; // if == 0, then no TAnaCandidate should be created.
		bool vertexing_; // do a vertexing at this node
		double mass_;
		double mass_tracks_; // apply a mass constraint to the tracks
		bool massConstraint_; // false: no massconstraint at this vertex
		double massSigma_;
		double maxDoca_;
		double minDoca_;
		bool daughtersToPV_;

		void dumpTabs(unsigned indent); // used by dump()
		
		std::set<track_entry_t> trackIndices_; // added tracks
		std::map<int,int> kinParticleMap_; // map: trackIx -> entry in the daughter kinematic particles...
		std::vector<HFDecayTree> subVertices_;
		RefCountedKinematicTree *kinTree_;
		TAnaCand *anaCand_;
		RefCountedHFNodeCut nodeCut_;
};

#endif
