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
RooWorkspace *build_model_nchannel(std::map<bmm_param,measurement_t> *bsmm, std::map<bmm_param,measurement_t> *bdmm, bool no_errors, int verbosity, bool compute_bd_ul);
RooWorkspace *build_model_light(std::map<bmm_param,measurement_t> *bsmm, int verbosity);
void measure_params(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity);

void est_ul_clb(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity, double err, double *pvalue);
RooStats::ConfInterval *est_ul_cls(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, double err = 0.0, double *ulLimit = NULL, double *cpuUsed = NULL);
RooStats::ConfInterval *est_ul_fc(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, uint32_t nbins = 20, std::pair<double,double> *rg = NULL, double *ulLimit = NULL, double *loLimit = NULL, double *cpuUsed = NULL);
RooStats::ConfInterval *est_ul_bc(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, double *ulLimit = NULL, double *cpuUsed = NULL);
RooStats::ConfInterval *est_ul_bc_light(RooWorkspace *wspace, RooDataSet *data, double cLevel, int verbosity, double *upperLimit, double *cpuUsed);

void compute_vars(std::map<bmm_param,measurement_t> *bmm, bool bstomumu);

#endif
