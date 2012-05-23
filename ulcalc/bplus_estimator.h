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

void estimate_bplus(std::map<bmm_param,measurement_t> *bplus, TTree *dataTree, TTree *mcTree, TTree *accTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut, std::map<systematics_t,double> *systematics_table);

#endif
