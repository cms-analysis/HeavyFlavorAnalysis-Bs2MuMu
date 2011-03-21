/*
 *  mumuBdReader.h
 *  macros
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */
#ifndef MUMUBDREADER_H
#define MUMUBDREADER_H

#include "mumuReader.hh"

/*
 * The only difference to the mumuReader is that we set the trueDecay a fTruthType variable
 * such that we consider the decay Bd -> mumu.
 * This reader is only used to compute efficiencies from signal MC.
 */

class mumuBdReader : public mumuReader {
	
	public:
		mumuBdReader(TChain *tree, TString evtClassName);
};

#endif
