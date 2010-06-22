#ifndef KSREADER_H
#define KSREADER_H

#include "massReader.hh"

#include <set>

class ksReader : public massReader {
	
	public:
		ksReader(TChain *tree, TString evtClassName);
		~ksReader();
		
		virtual void bookHist();
		virtual void eventProcessing();

	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand, int truth_type);
	
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
		
		multiset<int> trueDecay;
	
	private:
		unsigned unrecoverableDecays;
		unsigned recoverableCounter;
	private:
		unsigned buildDecay(int genIx, multiset<int> *particles, unsigned *nbrMuons = NULL);
		int getGenIndex(TAnaCand *pCand, int candID);
};

#endif
