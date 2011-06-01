/*
 *  NCRootUtilsLinkDef.h
 *  NCRootUtils
 *
 *  Created by Christoph on 5.3.10.
 *  Copyright 2010 Christoph Nägeli. All rights reserved.
 *
 */

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class fit_t;

#pragma link C++ function normalize_to_bins(TH1*);
#pragma link C++ function cut_negative_bins(TH1D*);
#pragma link C++ function cut_negative_bins(TH1D*,double);

#pragma link C++ function dump_branch_int(TTree*,const char*);
#pragma link C++ function dump_branch_int(TTree*,const char*,unsigned);
#pragma link C++ function dump_branch_float(TTree*,const char*);
#pragma link C++ function dump_branch_float(TTree*,const char*,unsigned);
#pragma link C++ function dump_branch_double(TTree*,const char*);
#pragma link C++ function dump_branch_double(TTree*,const char*,unsigned);
#pragma link C++ function dump_tree_formula(TTree*,const char*);
#pragma link C++ function dump_tree_formula(TTree*,const char*,unsigned);

#pragma link C++ function draw_to_pad(TH1*);
#pragma link C++ function draw_to_pad(TH1*,TCanvas*);
#pragma link C++ function draw_to_pad(TH1*,TCanvas*,int);
#pragma link C++ function draw_to_pad(TH1*,TCanvas*,int,const char*);

#pragma link C++ function draw_to_pad(TGraph*);
#pragma link C++ function draw_to_pad(TGraph*,TCanvas*);
#pragma link C++ function draw_to_pad(TGraph*,TCanvas*,int);
#pragma link C++ function draw_to_pad(TGraph*,TCanvas*,int,const char*);

#pragma link C++ function set_graph_appearance(TGraph*,int,const char *,const char*,const char*);

#pragma link C++ function adjust_parameter_gaussian(TH1D*, TF1*);
#pragma link C++ function adjust_parameter_gaussian(TH1D*, TF1*, map<int,double>*);
#pragma link C++ function adjust_parameter_gaussian(TH1D*, TF1*, map<int,double>*, double);
#pragma link C++ function adjust_parameter_gaussian(TH1D*, TF1*, map<int,double>*, double, double);
#pragma link C++ function adjust_parameter_gauss_expo(TH1D*, TF1*);
#pragma link C++ function adjust_parameter_gauss_expo(TH1D*, TF1*, map<int,double>*);
#pragma link C++ function adjust_parameter_gauss_expo(TH1D*, TF1*, map<int,double>*, double);
#pragma link C++ function adjust_parameter_gauss_expo(TH1D*, TF1*, map<int,double>*, double, double);
#pragma link C++ function adjust_parameter_gauss_linear(TH1D*, TF1*);
#pragma link C++ function adjust_parameter_gauss_linear(TH1D*, TF1*, map<int,double>*);
#pragma link C++ function adjust_parameter_gauss_linear(TH1D*, TF1*, map<int,double>*, double);
#pragma link C++ function adjust_parameter_gauss_linear(TH1D*, TF1*, map<int,double>*, double, double);
#pragma link C++ function adjust_parameter_double_gauss_linear(TH1D*, TF1*);
#pragma link C++ function adjust_parameter_double_gauss_linear(TH1D*, TF1*, map<int,double>*);
#pragma link C++ function adjust_parameter_double_gauss_linear(TH1D*, TF1*, map<int,double>*, double);
#pragma link C++ function adjust_parameter_double_gauss_linear(TH1D*, TF1*, map<int,double>*, double, double);
#pragma link C++ function adjust_parameter_gauss_p2(TH1D*, TF1 *);
#pragma link C++ function adjust_parameter_gauss_p2(TH1D*, TF1 *, map<int,double> *);
#pragma link C++ function adjust_parameter_gauss_p2(TH1D*, TF1 *, map<int,double> *, double);
#pragma link C++ function adjust_parameter_gauss_p2(TH1D*, TF1 *, map<int,double> *, double, double);
#pragma link C++ function adjust_parameter_crystall_ball(TH1D*, TF1*);
#pragma link C++ function adjust_parameter_crystall_ball(TH1D*, TF1*, map<int,double>*);
#pragma link C++ function adjust_parameter_crystall_ball(TH1D*, TF1*, map<int,double>*, double);
#pragma link C++ function adjust_parameter_crystall_ball(TH1D*, TF1*, map<int,double>*, double, double);

#pragma link C++ function init_gauss(fit_t*,double,double);
#pragma link C++ function init_gauss_expo(fit_t*,double,double);
#pragma link C++ function init_gauss_linear(fit_t*,double,double);
#pragma link C++ function init_gauss_p2(fit_t*,double,double);
#pragma link C++ function init_double_gauss_linear(fit_t*,double,double);
#pragma link C++ function init_crystal_ball(fit_t*,double,double);

#pragma link C++ function signal_events(TF1*,double,double,double);
#pragma link C++ function signal_events_gauss_linear(TF1*,double);
#pragma link C++ function signal_events_gauss_linear_histo(TH1D*,double,double);
#pragma link C++ function signal_events_gauss_p2(TF1*,double);
#pragma link C++ function signal_events_gauss_p2(TF1*,double,double*);
#pragma link C++ function signal_events_gauss_p2_histo(TH1D*,double,double);
#pragma link C++ function signal_events_gauss_expo(TF1*,double);
#pragma link C++ function signal_events_gauss_expo_histo(TH1D*,double,double);
#pragma link C++ function background_events_gauss_linear(TF1*,double);
#pragma link C++ function background_events_gauss_p2(TF1 *, double);
#pragma link C++ function background_events_gauss_expo(TF1*,double);
#pragma link C++ function significance(TF1*,TF1*,double,double,double);
#pragma link C++ function significance_gauss_linear(TF1*,double);
#pragma link C++ function significance_gauss_p2(TF1 *,double);
#pragma link C++ function significance_gauss_expo(TF1*,double);

#pragma link C++ function parse_variable_config(const char *, process_fct);
#pargma link C++ function parse_variable_config(const char *, process_fct, void*);

#pragma link C++ function plot_vars(TTree*,const char*);
#pragma link C++ function plot_vars(TTree*,const char*,TCut);
#pragma link C++ function plot_vars(TTree*,const char*,TCut,const char*);

#pragma link C++ function plot_tree_with_config(TTree*,const char*,const char*);
#pragma link C++ function plot_tree_with_config(TTree*,const char*,const char*,const char*);

#endif
