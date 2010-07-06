#ifndef KPREADER_H
#define KPREADER_H

#include "massReader.hh"

class kpReader : public massReader {
	
	public:
		kpReader(TChain *tree, TString evtClassName);
		~kpReader();
		
		virtual void bookHist();
		virtual void eventProcessing();
	
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand);
	
	private: // reduced Tree variables
		float fMassJPsi;
		float fDeltaR; // deltaR of J/Psi and Kp

		float fPtMu1;
		float fPtMu2;
		float fPtKp;

		TVector3 fPlabMu1;
		TVector3 fPlabMu2;
		TVector3 fPlabKp;
		
		TVector3 *fPlabMu1Ptr;
		TVector3 *fPlabMu2Ptr;
		TVector3 *fPlabKpPtr;
	
	private:
		map<int,int> decay_indices; // (genIx, ident_muons)
		unsigned long long total_counter;
		unsigned long long reco_counter;
		unsigned long long reco_single;
		unsigned long long reco_double;
};

#endif
