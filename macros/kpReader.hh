#ifndef KPREADER_H
#define KPREADER_H

#include "massReader.hh"

class kpReader : public massReader {
	
	public:
		kpReader(TChain *tree, TString evtClassName);
		virtual ~kpReader();
		
		virtual void bookHist();
		virtual void eventProcessing();
	
	protected:
		virtual void clearVariables();
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand);
		
		virtual bool parseCut(char *cutName, float cutLow, float cutHigh, int dump = 1);
		virtual bool applyCut();
	
	private: // reduced Tree variables
		float fMassJPsi;
		float fMassJPsiRec;
		float fDeltaR; // deltaR of J/Psi and Kp
		
		float fChi2Jpsi; // chi^2 of the j/psi daughter

		float fPtMu1;
		float fPtMu2;
		float fPtKp;
		
		float fPtKp_Gen; // pt of generator kaon
		float fEtaKp_Gen; // eta of generator kaon
		
		int fMuID1,fMuID2;
		float fEtaMu1,fEtaMu2;
		float fEtaKp;
		
		int fTrackQual_mu1;
		int fTrackQual_mu2;
		int fTrackQual_kp;
		
		int fQ_mu1;
		int fQ_mu2;
		int fQ_kp;
		
		float fD3_BpJpsi;
		float fD3e_BpJpsi;
	
	private:
		// Additional cut variables
		int fCutTrackQual_mu1;
		int fCutTrackQual_mu2;
		int fCutTrackQual_kp;
		int fCutMuID_mask;
		bool fCutMuID_reqall;
		bool fCutOppSign_mu;
		double fCutMass_JPsiLow;
		double fCutMass_JPsiHigh;
		double fCutPt_Kaon;
	
	private:
		std::map<int,int> decay_indices; // (genIx, ident_muons)
		unsigned long long total_counter;
		unsigned long long reco_counter;
		unsigned long long reco_single;
		unsigned long long reco_double;
};

#endif
