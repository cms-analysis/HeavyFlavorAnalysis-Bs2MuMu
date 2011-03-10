/*
 *  mumuReader.h
 *  macros
 *
 *  Created by Christoph on 09.01.11.
 *  Copyright 2011 PSI. All rights reserved.
 *
 */

#ifndef MUMUREADER_H
#define MUMUREADER_H

#include "massReader.hh"

class mumuReader : public massReader {
	
	public:
		mumuReader(TChain *tree, TString evtClassName);
		
		virtual void bookHist();
	
	protected:
		virtual void clearVariables();
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *pCand);

		virtual bool parseCut(char *cutName, float cutLow, float cutHigh, int dump = 1);
		virtual bool applyCut();
	
	private:
		float fPtMu1,fPtMu2;
		int fMuID1,fMuID2;
		float fDeltaR;
		float fEtaMu1,fEtaMu2;
		int fTrackQual_mu1;
		int fTrackQual_mu2;
		int fQ_mu1;
		int fQ_mu2;	
	private:
		int fCutTrackQual_mu1;
		int fCutTrackQual_mu2;
		int fCutMuID_mask;
		bool fCutMuID_reqall;
		bool fCutOppSign_mu;
};

#endif
