/*
 *  otherDecays.h
 *  macros
 *
 *  Created by Christoph on 06.01.11.
 *  Copyright 2011 PSI. All rights reserved.
 *
 */

#include "phiReader.hh"

class otherDecays : public phiReader {
	public:
		otherDecays(TChain* tree, TString evtClassName);
	
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
	
};
