#ifndef MASS_READER_H
#define MASS_READER_H

#include "treeReader01.hh"
#include <set>

class massReader : public treeReader01 {
	
	public:
		massReader(TChain *tree, TString evtClassName);
		~massReader();

		virtual void bookHist();
		virtual void eventProcessing();

	protected:
		// For subclasses
		TTree *reduced_tree;
		multiset<int> trueDecay;
		
		virtual int loadCandidateVariables(TAnaCand *pCand);
		
		// loading variables
		virtual int checkTruth(TAnaCand *cand); // check if all are originating from the same particle
		virtual int countMuons(TAnaCand *cand); // count the number of identified muons

		// creates the decay of the TGenCand
		void buildDecay(TGenCand *gen, multiset<int> *particles);
	
	protected:
		const char *fTreeName;
		
	protected:
		int fCandidate;
		TVector3 fMomentum;
		TVector3 *fMomentumPtr;
		float fMass;
		int fTruth;		// is this background or a true candidate?
                float fPt; // pt of the top particle
		float fNbrMuons;  // number of muons in the muon list.
		float fD3;
		float fD3E;
		float fDxy;	// distance to originating vertex (2d)
		float fDxyE;	// error (2d)
		float fAlpha; // angle between momentum and dist(vertex, motherVertex)
		float fAlphaXY; // Angle between momentum and dist(vertex, motherVertex) in xy-plane
		float fChi2; // chi2 of the vertex
		float fNdof; // number of degrees of freedom of vertex
		float fMaxDoca; // max doca
};

#endif
