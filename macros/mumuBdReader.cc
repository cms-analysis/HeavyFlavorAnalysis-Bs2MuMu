/*
 *  mumuBdReader.cc
 *  macros
 *
 *  Created by Christoph on 17.03.11.
 *  Copyright 2011 PSI. All rights reserved.
 *
 */

#include "mumuBdReader.hh"

mumuBdReader::mumuBdReader(TChain *tree, TString evtClassName) :
	mumuReader(tree, evtClassName)
{
	// override the behaviour
	std::cout << "Setting mumuReader to Bd -> MuMu!" << std::endl;
	fTruthType = 511;
	
	trueDecay.clear();
	trueDecay.insert(13); // mu
	trueDecay.insert(13); // mu
	trueDecay.insert(fTruthType); // Bd
} // mumuBdReader()
