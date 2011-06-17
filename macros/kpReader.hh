#ifndef KPREADER_H
#define KPREADER_H

#include "massReader.hh"

class kpReader : public massReader {
	
	public:
		kpReader(TChain *tree, TString evtClassName);
		virtual ~kpReader();
		
		virtual void bookHist();
	
	protected:
		virtual void clearVariables();
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand);
		
		virtual bool parseCut(char *cutName, float cutLow, float cutHigh, int dump = 1);
		virtual bool applyCut();
	
	private: // reduced Tree variables
		float fMassJPsi;
		float fMassJPsiRec;
		
		float fChi2Jpsi; // chi^2 of the j/psi daughter
		float fPtKp;
		float fPtKp_Gen; // pt of generator kaon
		float fEtaKp_Gen; // eta of generator kaon
		float fEtaKp;
		
		int fTrackQual_kp;
		
		int fQ_kp;
		
		float fD3_BpJpsi;
		float fD3e_BpJpsi;
	
	private:
		// Additional cut variables
		int fCutTrackQual_kp;
		double fCutMass_JPsiLow;
		double fCutMass_JPsiHigh;
		double fCutPt_Kaon;
};

#endif
