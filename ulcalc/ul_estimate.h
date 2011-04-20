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

void add_channels(std::map<bmm_param,double> *bmm, std::set<int> *channels);

RooWorkspace *build_model_light(std::map<bmm_param_tag,double> *bsmmVars, bool silent);
RooWorkspace *build_model(std::map<bmm_param_tag,double> *bsmmVars, std::map<bmm_param_tag,double> *bdmmVars, bool silent);
RooWorkspace *build_model_split(std::map<bmm_param_tag,double> *bsmmBar, std::map<bmm_param_tag,double> *bsmmEnd, std::map<bmm_param_tag,double> *bdmmBar, std::map<bmm_param_tag,double> *bdmmEnd, bool silent);
RooWorkspace *build_model_nchannel(std::map<bmm_param,double> *bsmm, std::map<bmm_param,double> *bdmm, bool silent);

RooDataSet *build_data(RooWorkspace *wspace, double nsObs, double ndObs, double nbObs);
RooDataSet *build_data_split(RooWorkspace *wspace,double nsObsB, double nsObsE, double ndObsB, double ndObsE, double nbObsB, double nbObsE);
void estimate_start_values(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels);

RooStats::ConfInterval *est_ul_fc(RooWorkspace *wspace, RooDataSet *data, double cLevel, double *ulLimit = NULL, bool splitModel = true, double *cpuUsed = NULL);
RooStats::ConfInterval *est_ul_bc(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, double *ulLimit = NULL, double *cpuUsed = NULL);
RooStats::ConfInterval *est_ul_mc(RooWorkspace *wspace, RooDataSet *data, double cLevel, double *ulLimit = NULL, bool splitModel = true, double *cpuUsed = NULL);
RooStats::ConfInterval *est_ul_cls(RooWorkspace *wspace, RooDataSet *data, double cLevel, double *ulLimit = NULL, bool splitModel = true, double *cpuUsed = NULL);

void compute_vars(std::map<bmm_param,double> *bmm, bool bstomumu);
double compute_efftot_bplus(std::map<bmm_param,double> *bmm, int channel);
double compute_efftot_bmm(std::map<bmm_param,double> *bmm, int channel);

#endif
