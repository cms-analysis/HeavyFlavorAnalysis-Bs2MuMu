/*
 *  mumuPsiReader.h
 *  macros
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 04.07.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#include "mumuReader.hh"

class mumuPsiReader : public mumuReader {
	
	public:
		mumuPsiReader(TChain *tree, TString evtClassName);
};