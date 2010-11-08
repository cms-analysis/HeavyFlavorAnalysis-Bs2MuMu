#ifndef IMPACT_READER_H
#define IMPACT_READER_H

#include <set>
#include <map>

#include "treeReader01.hh"

class impactReader : public treeReader01 {
	
	public:
		impactReader(TChain *tree, TString evtClassName);
		
		virtual void eventProcessing();
		virtual void bookHist();
		virtual void closeHistFile();
	
	private:
		TTree *reduced_tree;
		float fLip_CMSSW;
		float fTip_CMSSW;
		float fLip_geom;
		float fTip_geom;
		int fIx_CMSSW;
		int fIx_geom;
		int fNbrPV;
	private:
		int calculatePVIx(TAnaCand *pCand);
		TVector3 calculatePVDist(TAnaCand *pCand, int pvIx);
};

#endif
