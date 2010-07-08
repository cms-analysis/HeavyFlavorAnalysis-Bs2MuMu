#ifndef LAMBDA_READER_H
#define LAMBDA_READER_H

#include "treeReader01.hh"

class lambdaReader : public treeReader01 {

public:
    lambdaReader(TChain *tree, TString evtClassName);
    ~lambdaReader();

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

    // Private functions
    int checkTruth(TAnaCand *cand, int truth_type);
    int CheckGoodRun(int run);

    // list of good runs
    std::vector<int> goodRuns;

    double DCAMAX, SIGVTX, TYPE;

    // data members for TTree's
    int goodrun, nPV;
    int runNr;
    // J/Psi
    int jpsi_truth;
    double jpsi_m, jpsi_pt, jpsi_dxy, jpsi_dxyE, jpsi_d3d, jpsi_d3dE;
    double jpsi_chi2vtx, jpsi_deltaR;
    double jpsi_dca, jpsi_mu1_pt, jpsi_mu2_pt, jpsi_mu1_eta, jpsi_mu2_eta;
    double fdistPV, fdistDp, jpsi_mu1_dR, jpsi_mu2_dR, fdeltaRMu;
};

#endif
