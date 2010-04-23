#ifndef DZERO_READER_H
#define DZERO_READER_H

#include "treeReader01.hh"

using namespace std;

class d0Reader : public treeReader01 {
	
public:
	d0Reader(TChain *tree, TString evtClassName);
	~d0Reader();

	void bookHist();
	void eventProcessing();
        void startAnalysis();
        void readCuts(TString filename, int dump);

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
        int CheckGoodRun(int run);
        std::vector<int> goodRuns;

  double DCAMAX, SIGVTX, TYPE;

  int goodrun, nPV;
  int runNr;
  double dppt, dpm, dpdxy, dpdxyE, dpd3d, dpd3dE;
  double bpt, bm, bdxy, bdxyE, bd3d, bd3dE;
  double fdistPV, fdistDp, fdeltaR, fdeltaRMax, fdeltaRMu;
  double dcamin, dcamax;
  double ptmin, etamax, ptmu, etamu;
  double chi2vtx, dalpha;
};

#endif
