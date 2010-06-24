#ifndef KPREADER_H
#define KPREADER_H

#include "massReader.hh"

class kpReader : public massReader {
	
	public:
		kpReader(TChain *tree, TString evtClassName);
		
		virtual void bookHist();
	
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand);
	
	private:
		double fMassJPsi;
		
		TVector3 fPlabMu1;
		TVector3 fPlabMu2;
		TVector3 fPlabKp;
		
		TVector3 *fPlabMu1Ptr;
		TVector3 *fPlabMu2Ptr;
		TVector3 *fPlabKpPtr;
};

#endif
