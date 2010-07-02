#ifndef KSREADER_H
#define KSREADER_H

#include "massReader.hh"

class ksReader : public massReader {
	
	public:
		ksReader(TChain *tree, TString evtClassName);
		~ksReader();
		
		virtual void bookHist();
		virtual void eventProcessing();

	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand);
	
	private:
		float fMassKs;
		float fAlphaKs;
		float fDxyKs;
		float fDxyeKs;
		float fChi2Ks;
		
		float fDeltaKs;
		float fDeltaKsTrue;
		
		float fMaxDocaKs;
		
		float fMassJPsi;
		
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
	
	private:
		set<int> decay_indices;
		unsigned long long total_counter;
		unsigned long long reco_counter;
};

#endif
