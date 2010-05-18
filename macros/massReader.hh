#ifndef MASS_READER_H
#define MASS_READER_H

#include "treeReader01.hh"

class massReader : public treeReader01 {
	
	public:
		massReader(TChain *tree, TString evtClassName);
		~massReader();

		virtual void bookHist();
		virtual void eventProcessing();

	protected:
		// For subclasses
		TTree *reduced_tree;
		
		virtual int loadCandidateVariables(TAnaCand *pCand);

	private:
		// Private variables
		int fCandidate;
		TVector3 fMomentum;
		TVector3 *fMomentumPtr;
		double fMass;
		int fTruth;		// is this background or a true candidate?
		int fTwoMuon; // are both muons fromthe muon list?
		double fDxy;	// distance to originating vertex
		double fDxyE;	// error
		double fAlpha; // angle between momentum and dist(vertex, motherVertex)
		double fChi2; // chi2 of the vertex
		double fNdof; // number of degrees of freedom of vertex
				
		// Private functions
		int checkTruth(TAnaCand *cand, int truth_type);
		int checkMuons(TAnaCand *cand);
};

#endif
