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
	
	private: // additional tree variables
		float fMassJPsi;
		float fMassPhi;
		float fDeltaR; // deltaR of J/Psi Phi
		
		float fPtMu1;
		float fPtMu2;
		float fPtKp1;
		float fPtKp2;
		
		TVector3 fPlabMu1;
		TVector3 fPlabMu2;
		TVector3 fPlabKp1;
		TVector3 fPlabKp2;
		
		TVector3 *fPlabMu1Ptr;
		TVector3 *fPlabMu2Ptr;
		TVector3 *fPlabKp1Ptr;
		TVector3 *fPlabKp2Ptr;
		
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
		std::map<int,int> decay_indices; // (genIx, ident_muons)
		unsigned long long total_counter;
		unsigned long long reco_single;
		unsigned long long reco_double;
};

#endif
