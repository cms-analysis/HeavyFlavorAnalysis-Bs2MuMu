#ifndef PHIREADER_H
#define PHIREADER_H

#include "massReader.hh"

class phiReader : public massReader {
	
	public:
		phiReader(TChain* tree, TString evtClassName);
		
		virtual void bookHist();
		
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand);
	
	private:
		double fMassJPsi;
		double fMassPhi;
		
		TVector3 fPlabMu1;
		TVector3 fPlabMu2;
		TVector3 fPlabKp1;
		TVector3 fPlabKp2;
		
		TVector3 *fPlabMu1Ptr;
		TVector3 *fPlabMu2Ptr;
		TVector3 *fPlabKp1Ptr;
		TVector3 *fPlabKp2Ptr;
};

#endif
