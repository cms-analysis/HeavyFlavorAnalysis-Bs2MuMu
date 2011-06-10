/*
 *  bplus_estimator.h
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 09.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef BPLUS_ESTIMATOR
#define BPLUS_ESTIMATOR

#include "external_constants.h"

// STL
#include <map>
#include <string>

// ROOT headers
#include <TCut.h>
#include <TTree.h>

// Copied from massReader.h
enum
{
	kGeneratorCand	= 1 << 0,
	kAcceptance		= 1 << 1,
	kEffMuon		= 1 << 2
	//	kEffTrig is not needed
	//	kEffCand is not needed
	//	kEffAna is not needed
};

void estimate_bplus(std::map<bmm_param,measurement_t> *bplus, TTree *dataTree, TTree *mcTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut);


#endif
