/*
 *  NCRootUtils.h
 *  NCRootUtils
 *
 *  Created by Christoph on 4.3.10.
 *  Copyright 2010 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef NCRootUtils_
#define NCRootUtils_

#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TTree.h>
#include <TCut.h>
#include <TString.h>
#include <TF1.h>

/* The classes below are exported */
#pragma GCC visibility push(default)

/* Utility routines for Histograms */
void normalize_to_bins(TH1 *h);
void cut_negative_bins(TH1D *h, double min = 0.0);

void dump_branch_int(TTree *tree, const char *name, unsigned max = 100u);
void dump_branch_float(TTree *tree, const char *name, unsigned max = 100u);
void dump_branch_double(TTree *tree, const char *name, unsigned max = 100u);
void dump_tree_formula(TTree *tree, const char *formula, unsigned max = 100u);

/* Canvas utility functions */
void draw_to_pad(TH1* histo, TCanvas *canvas = new TCanvas, int pad_nbr = 1, const char *option = "");
void draw_to_pad(TGraph *graph, TCanvas *canvas = new TCanvas, int pad_nbr = 1, const char *option = "ap");

/* Graph functions */
void set_graph_appearance(TGraph* graph, int mstyle, const char *title = "", const char *xname = "", const char *yname = "");

/* Routines for fitting Histograms: ROOT fitting */
typedef void (*init_fit_fct)(TH1D*,TF1*,std::map<int,double> *, double, double);

struct _fit_t {
	init_fit_fct init_fit;
	TF1 *fit_fct;
	int mean_param;
	int sigma_param;
	int height_param;
};
typedef struct _fit_t fit_t;

void adjust_parameter_gaussian(TH1D* h, TF1* fct, std::map<int,double>* fixed = NULL, double a = 0.0, double b = 0.0);
void adjust_parameter_gauss_expo(TH1D* h, TF1* fct, std::map<int,double>* fixed = NULL, double a = 0.0, double b = 0.0);
void adjust_parameter_gauss_linear(TH1D* h, TF1* fct, std::map<int,double>* fixed = NULL, double a = 0.0, double b = 0.0);
void adjust_parameter_double_gauss_linear(TH1D* h, TF1* fct, std::map<int,double>* fixed = NULL, double a = 0.0, double b = 0.0);
void adjust_parameter_gauss_p2(TH1D* h, TF1 *fct, std::map<int,double> *fixed = NULL, double a = 0.0, double b = 0.0);
void adjust_parameter_crystall_ball(TH1D* h, TF1* fct, std::map<int,double>* fixed = NULL, double a = 0.0, double b = 0.0);

void init_gauss(fit_t *ft, double xmin, double xmax);
void init_gauss_expo(fit_t *ft, double xmin, double xmax);
void init_gauss_linear(fit_t *ft, double xmin, double xmax);
void init_gauss_p2(fit_t *ft, double xmin, double xmax);
void init_double_gauss_linear(fit_t *ft, double xmin, double xmax);
void init_crystal_ball(fit_t *ft, double xmin, double xmax);

/* Calculate the Number of events */
double signal_events(TF1 *signalFct, double mu, double sigma, double bin_width);
double signal_events_gauss_linear(TF1 *fct, double bin_width);
double signal_events_gauss_linear_histo(TH1D *h, double xmin, double xmax);
double signal_events_gauss_p2(TF1 *fct, double bin_width, double *signal_err = NULL);
double signal_events_gauss_p2_histo(TH1D* h, double xmin, double xmax);
double signal_events_gauss_expo(TF1 *fct, double bin_width);
double signal_events_gauss_expo_histo(TH1D *h, double xmin, double xmax);
double background_events_gauss_linear(TF1 *fct, double bin_width);
double background_events_gauss_p2(TF1 *fct, double bin_width);
double background_events_gauss_expo(TF1 *fct, double bin_width);

/* Calculate the Significance */
double significance(TF1 *signalFct, TF1 *backFct, double mu, double sigma, double bin_width);
double significance_gauss_linear(TF1 *fct,double bin_width);
double significance_gauss_p2(TF1 *fct,double bin_width);
double significance_gauss_expo(TF1 *fct,double bin_width);

/* Variable plotting */

struct _function_range_t {
	TString name;
	std::vector<double> vals;
	std::string options;
};
typedef struct _function_range_t function_range_t;
typedef void (*process_fct)(function_range_t *, void*);

void parse_variable_config(const char *config_file, process_fct f, void *param = NULL);

void plot_vars(TTree *tree, const char *var_file, TCut preCut = TCut(), const char *ending = "eps");
void plot_tree_with_config(TTree *tree, const char *config_file, const char *histoname, const char *dir = NULL);

#pragma GCC visibility pop

#endif
