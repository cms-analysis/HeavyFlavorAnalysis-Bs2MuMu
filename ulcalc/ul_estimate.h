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
#include <RooStats/HypoTestResult.h>

void add_channels(std::map<bmm_param,measurement_t> *bmm, std::set<int> *channels);
RooWorkspace *build_model_nchannel(std::map<bmm_param,measurement_t> *bsmm, std::map<bmm_param,measurement_t> *bdmm, bool no_errors, int verbosity, bool compute_bd_ul, bool fixed_bkg, bool floatPoissonians, bool smCrossFeed, bool bothPoi);
void measure_params(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity);

void est_ul_clb(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity, double *pvalue, uint32_t nbrProof, int nToys);
void est_ul_clb_hybrid(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity, double *pvalue, uint32_t nbrProof, int nToys, bool bdmm, bool fixedBkg, bool smExpectation, bool smCrossFeed);
void est_int_hybrid(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, double *ulLimit, std::pair<double,double> *rg, uint32_t* inBins, double *cpuUsed, uint32_t nbrProof, int nToys, bool bdmm, bool fixedBkg, bool smExpectation, bool smCrossFeed);
void est_ul_hybrid(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, double *ulLimit, std::pair<double,double> *rg, uint32_t* inBins, double *cpuUsed, uint32_t nbrProof, int nToys, bool bdmm, bool fixedBkg, bool smExpectation, bool smCrossFeed);
void est_ul_cls(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, double *ulLimit, std::pair<double,double> *rg, uint32_t *npts, double *cpuUsed, uint32_t nbrProof, int nToys);
void est_ul_fc(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, double *ulLimit, double *loLimit, std::pair<double,double> *rg, uint32_t *inBins, double *cpuUsed, uint32_t nbrProof, int nToys);
void est_ul_bc(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, double *ulLimit, double *cpuUsed);
void est_ul_zbi(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, bool bdlimit, double *ul);
void est_sign_bkg(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity, double *pvalue, uint32_t nbrProof, int nToys, bool fixedBkg);

void compute_vars(std::map<bmm_param,measurement_t> *bmm, bool bstomumu);

#endif
