/*
 *  bmm_estimator.h
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef BMM_ESTIMATOR
#define BMM_ESTIMATOR

#include "external_constants.h"
#include "bplus_estimator.h"

// STL
#include <map>
#include <string>

// ROOT headers
#include <TCut.h>
#include <TTree.h>

// routines for B->mumu
void estimate_bmm(std::map<bmm_param,measurement_t> *bmm, TTree *dataTree, TTree *mcTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut, std::pair<double,double> bd_window, std::pair<double,double> bs_window, bool is_bstomumu, double eff_filter = 1, bool enable_systematics = true);

#endif
