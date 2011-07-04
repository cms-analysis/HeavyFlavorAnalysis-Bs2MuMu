/*
 *  mumuPsiReader.cpp
 *  macros
 *
 *  Created by Christoph on 04.07.11.
 *  Copyright 2011 PSI. All rights reserved.
 *
 */

#include "mumuPsiReader.hh"

mumuPsiReader::mumuPsiReader(TChain *tree, TString evtClassName) :
	mumuReader(tree,evtClassName)
{
	// override the std behaviour
	std::cout << "Setting mumuReader to J/Psi -> MuMu!" << std::endl;
	fTruthType = 443;
	
	trueDecay.clear();
	trueDecay.insert(13); // mu
	trueDecay.insert(13); // mu
	trueDecay.insert(fTruthType); // J/Psi
} // mumuPsiReader()
