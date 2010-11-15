#ifndef KPISOLATIONREADER_H
#define KPISOLATIONREADER_H

#include "kpReader.hh"

class kpIsolationReader : public kpReader {
	
	public:
		kpIsolationReader(TChain *tree, TString evtClassName);
		virtual ~kpIsolationReader();
		
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
