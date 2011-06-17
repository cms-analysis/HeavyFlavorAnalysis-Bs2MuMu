/*
 *  mumuReader.h
 *  macros
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 09.01.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef MUMUREADER_H
#define MUMUREADER_H

#include "massReader.hh"

class mumuReader : public massReader {
	
	public:
		mumuReader(TChain *tree, TString evtClassName);
	
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
		virtual int checkTruth(TAnaCand *pCand);

};

#endif
