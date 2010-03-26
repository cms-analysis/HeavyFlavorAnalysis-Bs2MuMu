#ifndef MASS_READER_H
#define MASS_READER_H

#include "treeReader01.hh"

using namespace std;

class massReader : public treeReader01 {
	
public:
	massReader(TChain *tree, TString evtClassName);
	~massReader();

	void bookHist();
	void eventProcessing();

private:
	// Private variables
	TTree *reduced_tree;
	int fCandidate;
	TVector3 fMomentum;
	TVector3 *fMomentumPtr;
	double fMass;
	int fTruth;
	
	// Private functions
	int checkTruth(TAnaCand *cand, int truth_type);
};

#endif
