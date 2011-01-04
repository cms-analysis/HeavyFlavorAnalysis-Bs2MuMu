#ifndef PHIREADER_H
#define PHIREADER_H

#include "massReader.hh"
#include <map>

class phiReader : public massReader {
	
	public:
		phiReader(TChain* tree, TString evtClassName);
		virtual ~phiReader();
		
		virtual void bookHist();
		virtual void eventProcessing();
		
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *cand);
		
		virtual bool parseCut(char *cutName, float cutLow, float cutHigh, int dump = 1);
		virtual bool applyCut();
	
	private: // additional tree variables
		float fMassJPsi;
		float fMassJPsiRec;
		float fMassPhi;
		float fDeltaR; // deltaR of J/Psi Phi
		float fDeltaR_Kaons; // delta R of the two Kaons
		
		float fPtMu1;
		float fPtMu2;
		float fPtKp1;
		float fPtKp2;
		
		int fMuID1,fMuID2;
		float fEtaMu1,fEtaMu2;
		
		int fTrackQual_mu1;
		int fTrackQual_mu2;
		int fTrackQual_kp1;
		int fTrackQual_kp2;
		
		int fQ_mu1;
		int fQ_mu2;
		int fQ_kp1;
		int fQ_kp2;
		
		float fD3_BsJpsi;
		float fD3e_BsJpsi;
	
	private:
		// Additional cut variables
		int fCutTrackQual_mu1;
		int fCutTrackQual_mu2;
		int fCutTrackQual_kp1;
		int fCutTrackQual_kp2;
		bool fCutOppSign_mu;
		bool fCutOppSign_kp;
		double fCutMass_JPsiLow;
		double fCutMass_JPsiHigh;
		double fCutMass_PhiLow;
		double fCutMass_PhiHigh;
		double fCutPt_Kp2;
		double fCutDeltaR_Kaons;
	
	private:
		std::map<int,int> decay_indices; // (genIx, ident_muons)
		unsigned long long total_counter;
		unsigned long long reco_single;
		unsigned long long reco_double;
};

#endif
