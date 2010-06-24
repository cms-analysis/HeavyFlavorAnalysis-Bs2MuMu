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

	private:
		// Private variables
		int fCandidate;
		TVector3 fMomentum;
		TVector3 *fMomentumPtr;
		double fMass;
		int fTruth;		// is this background or a true candidate?
		int fNbrMuons;  // number of muons in the muon list.
		double fD3;
		double fD3E;
		double fDxy;	// distance to originating vertex (2d)
		double fDxyE;	// error (2d)
		double fAlpha; // angle between momentum and dist(vertex, motherVertex)
		double fChi2; // chi2 of the vertex
		double fNdof; // number of degrees of freedom of vertex
};

#endif
