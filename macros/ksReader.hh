#ifndef KSREADER_H
#define KSREADER_H

#include "massReader.hh"

class ksReader : public massReader {
	
	public:
		ksReader(TChain *tree, TString evtClassName);
		
		virtual void bookHist();

	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand);
	
	private:
		double fMassKs;
		double fAlphaKs;
		double fDxyKs;
		double fDxyeKs;
		double fChi2Ks;
		
		double fDeltaKs;
		double fDeltaKsTrue;
		
		double fMassJPsi;
		
		TVector3 fPlabMu1;
		TVector3 fPlabMu2;
		TVector3 fPlabPi1;
		TVector3 fPlabPi2;
		TVector3 fPlabKs;
		TVector3 fPlabKsTrue;
		
		TVector3 *fPlabMu1Ptr;
		TVector3 *fPlabMu2Ptr;
		TVector3 *fPlabPi1Ptr;
		TVector3 *fPlabPi2Ptr;
		TVector3 *fPlabKsPtr;
		TVector3 *fPlabKsTruePtr;
	private:
		int getGenIndex(TAnaCand *pCand, int candID);
};

#endif
