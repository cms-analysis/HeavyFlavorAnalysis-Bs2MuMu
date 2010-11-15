#ifndef PHIISOLATIONREADER_H
#define PHIISOLATIONREADER_H

#include "phiReader.hh"

class phiIsolationReader : public phiReader {
	
	public:
		phiIsolationReader(TChain *tree, TString evtClassName);
		virtual ~phiIsolationReader();
		
		virtual void bookHist();
	
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
	
	protected:
		TH2D *pt_ix;
		TH2D *ip_ix;
		TH2D *pt_ip;
		TH2D *ptrel_ip;
};

#endif
