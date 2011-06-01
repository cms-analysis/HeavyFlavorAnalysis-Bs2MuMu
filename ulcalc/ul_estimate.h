/*
 *  ul_estimate.h
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef UL_ESTIMATE_H
#define UL_ESIIMATE_H

#include "external_constants.h"

// ROOT headers
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooStats/ConfInterval.h>

void add_channels(std::map<bmm_param,measurement_t> *bmm, std::set<int> *channels);
RooWorkspace *build_model_nchannel(std::map<bmm_param,measurement_t> *bsmm, std::map<bmm_param,measurement_t> *bdmm, bool no_errors, int verbosity);
void estimate_start_values(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity);

RooStats::ConfInterval *est_ul_fc(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, uint32_t nbins = 20, std::pair<double,double> *rg = NULL, double err = 0.0, double *ulLimit = NULL, double *cpuUsed = NULL);
RooStats::ConfInterval *est_ul_bc(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, double *ulLimit = NULL, double *cpuUsed = NULL);
RooStats::ConfInterval *est_ul_cls(RooWorkspace *wspace, RooDataSet *data, double cLevel, double *ulLimit = NULL, bool splitModel = true, double *cpuUsed = NULL);

void compute_vars(std::map<bmm_param,measurement_t> *bmm, bool bstomumu);
measurement_t compute_efftot_bplus(std::map<bmm_param,measurement_t> *bmm, int channel);
measurement_t compute_efftot_bmm(std::map<bmm_param,measurement_t> *bmm, int channel);

#endif
