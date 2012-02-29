#ifndef PHIREADER_H
#define PHIREADER_H

#include "massReader.hh"
#include <map>

class phiReader : public massReader {
	
	public:
		phiReader(TChain* tree, TString evtClassName);
		virtual ~phiReader();
		
		virtual void bookHist();
		
	protected:
		virtual void clearVariables();
		virtual int loadCandidateVariables(TAnaCand *pCand);
		
		virtual bool parseCut(char *cutName, float cutLow, float cutHigh, int dump = 1);
		virtual bool applyCut();
	
	private: // additional tree variables
		float fMassJPsi;
		float fMassJPsiRec;
		float fMassPhi;
		float fDeltaR_Kaons; // delta R of the two Kaons
		
		float fPtKp1;
		float fPtKp2;
		
		int fTrackQual_kp1;
		int fTrackQual_kp2;
		
		int fQ_kp1;
		int fQ_kp2;
		
		float fD3_BsJpsi;
		float fD3e_BsJpsi;
	
	private:
		// Additional cut variables
		int fCutTrackQual_kp1;
		int fCutTrackQual_kp2;
		bool fCutOppSign_kp;
		double fCutMass_JPsiLow;
		double fCutMass_JPsiHigh;
		double fCutMass_PhiLow;
		double fCutMass_PhiHigh;
		double fCutPt_Kp2;
		double fCutDeltaR_Kaons;
};

#endif
