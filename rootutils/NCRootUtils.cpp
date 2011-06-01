/*
 *  NCRootUtils.cp
 *  NCRootUtils
 *
 *  Created by Christoph on 4.3.10.
 *  Copyright 2010 Christoph NŠgeli. All rights reserved.
 *
 */

// C++ headers
#include <iostream>
#include <cmath>
#include <set>
#include <map>

// ROOT headers
#include <TCollection.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TTreeFormula.h>

// boost headers
#include <boost/scoped_ptr.hpp>

// my headers
#include "NCRootUtils.h"
#include "NCFunctions.h"

using std::set;
using std::map;
using boost::scoped_ptr;

static const double signal_factor = 2.39258;

void normalize_to_bins(TH1 *h)
{
	TAxis *ax = h->GetXaxis();
	int32_t j,nbins;
	
	nbins = ax->GetNbins();
	for (j = 1; j <= nbins; j++)
		h->SetBinContent(j, h->GetBinContent(j)/(ax->GetBinUpEdge(j)-ax->GetBinLowEdge(j)));
} // normalize_to_bins()

void cut_negative_bins(TH1D *h, double min)
{
	Int_t j;
	double val;
	
	for (j = 1; j <= h->GetNbinsX(); j++) {
		val = h->GetBinContent(j);
		val = (val < min) ? min : val;
		h->SetBinContent(j, val);
	}
} // cut_negative_bins()

template<class T>
static void dump_branch(TTree *tree, const char *name, unsigned max, T val)
{
	using namespace std;
	TBranch *branch;
	set<T> vals;
	typename set<T>::iterator it;
	unsigned j,nentries;
	
	branch = tree->FindBranch(name);
	if(!branch) {
		cerr << "No Branch name " << name << " in the tree." << endl;
		goto bail;
	}
	branch->SetAddress(&val);
	
	nentries = tree->GetEntries();
	for (j = 0; j < nentries && vals.size() < max; j++) {
		branch->GetEvent(j);
		vals.insert(val);
	}
	
	// dump the set
	cout << "The Tree has the following values in the branch " << name << endl;
	for (it = vals.begin(); it != vals.end(); ++it)
		cout << '\t' << *it;
	cout << endl;
bail:
	return;
} // dump_values()

void dump_branch_int(TTree *tree, const char *name, unsigned max)
{
	dump_branch(tree,name,max,int());
} // dump_values_int()

void dump_branch_float(TTree *tree, const char *name, unsigned max)
{
	dump_branch(tree,name,max,float());
} // dump_branch_float()

void dump_branch_double(TTree *tree, const char *name, unsigned max)
{
	dump_branch(tree,name,max,double());
} // dump_values_double()

void dump_tree_formula(TTree *tree, const char *formula, unsigned max)
{
	set<Double_t> vals;
	set<Double_t>::const_iterator it;
	Long64_t j,nentries;
	Double_t var;
	TTreeFormula *treeFormula = new TTreeFormula("",formula,tree);
	
	
	nentries = tree->GetEntries();
	for (j = 0; j < nentries && vals.size() < max; j++) {
		
		tree->GetEntry(j);
		var = treeFormula->EvalInstance();
		vals.insert(var);
	}
	
	std::cout << "The formula evaluates to the following values:" << std::endl;
	for (it = vals.begin(); it != vals.end(); ++it)
		std::cout << '\t' << *it;
	std::cout << std::endl;
	
	delete treeFormula;
} // dump_formula()

void draw_to_pad(TH1* histo, TCanvas *canvas, int pad_nbr, const char *option)
{
	canvas->cd(pad_nbr);
	histo->Draw(option);
} // draw_to_pad()

void draw_to_pad(TGraph *graph, TCanvas *canvas, int pad_nbr, const char *option)
{
	canvas->cd(pad_nbr);
	graph->Draw(option);
} // draw_to_pad()

void set_graph_appearance(TGraph* graph, int mstyle, const char *title, const char *xname, const char *yname)
{
	graph->SetTitle(title);
	graph->SetMarkerStyle(mstyle);
	graph->GetXaxis()->SetTitle(xname);
	graph->GetYaxis()->SetTitle(yname);
} // set_graph_appearance()


static void fix_parameters(TF1 *fct, map<int,double> *fixed)
{
	if(!fixed) return;
	for(map<int,double>::const_iterator it = fixed->begin(); it != fixed->end(); ++it) {
		fct->SetParameter(it->first, it->second);
		fct->SetParLimits(it->first, it->second, it->second);
	}
} // fix_parameters()

// Initializes a fit sturecture used to fit with a crystal ball function
void init_crystal_ball(fit_t *ft, double xmin, double xmax)
{
	static unsigned uuid = 1;
	ft->init_fit = adjust_parameter_crystall_ball;
	ft->fit_fct = new TF1(Form("crystall_ball_fct_uuid%u",uuid++),f_cb,xmin,xmax,5);
	ft->mean_param = 0;
	ft->sigma_param = 1;
	ft->height_param = 4;
}

// Initializes a fit sturecture used to fit with a gaussian
void init_gauss(fit_t *ft, double xmin, double xmax)
{
	static unsigned uuid = 1; // make the name unique
	ft->init_fit = adjust_parameter_gaussian;
	ft->fit_fct = new TF1(Form("gauss_fct_uuid%u",uuid++),f_gauss,xmin,xmax,3);
	ft->mean_param = 1;
	ft->sigma_param = 2;
	ft->height_param = 0;
} // init_gauss()

void init_gauss_linear(fit_t *ft, double xmin, double xmax)
{
	static unsigned uuid = 1; // make the name unique
	ft->init_fit = adjust_parameter_gauss_linear;
	ft->fit_fct = new TF1(Form("gauss_linear_fct_uuid%u",uuid++),f_gauss_linear,xmin,xmax,5);
	ft->mean_param = 1;
	ft->sigma_param = 2;
	ft->height_param = 0;
	
	// set the names of the variables
	ft->fit_fct->SetParNames("c", "mu", "sigma", "a", "m");
} // init_gauss_linear()

void init_gauss_p2(fit_t *ft, double xmin, double xmax)
{
	ft->init_fit = adjust_parameter_gauss_p2;
	ft->fit_fct = new TF1("",f_p2ag,xmin,xmax,6);
	ft->height_param = 0;
	ft->mean_param = 1;
	ft->sigma_param = 2;
	
	ft->fit_fct->SetParNames("c","mu","sigma","a0","a1","a2");
} // init_gauss_p2()

void init_double_gauss_linear(fit_t *ft, double xmin, double xmax)
{
	static unsigned uuid = 1; // make then name unique
	ft->init_fit = adjust_parameter_double_gauss_linear;
	ft->fit_fct = new TF1(Form("double_gauss_linear_fct_uuid%u",uuid++),f_double_gauss_linear,xmin,xmax,8);
	ft->mean_param = 1;
	ft->sigma_param = 2;
	ft->height_param = 0;
} // init_double_gauss_linear()

// Initializes a fit sturecture used to fit with a gaussian+exponential
void init_gauss_expo(fit_t *ft, double xmin, double xmax)
{
	static unsigned uuid = 1; // make the name unique
	ft->init_fit = adjust_parameter_gauss_expo;
	ft->fit_fct = new TF1(Form("gauss_expo_fct_uuid%u",uuid++),f_gauss_expo,xmin,xmax,5);
	ft->mean_param = 1;
	ft->sigma_param = 2;
	ft->height_param = 0;
	
	// set the names of the variables
	ft->fit_fct->SetParNames("c","mu","sigma","a","lambda");
} // init_gauss_expo()

// adjust the parameters given the function f is a f_Gauss()+f_expo()
void adjust_parameter_gauss_expo(TH1D* h, TF1* fct, map<int,double>* fixed, double a, double b)
{
	double x0,x1,f0,f1;
	double left,right,mu,c,sigma;
	double middle,width;
	int b0,b1,nbins,max_bin;
	TAxis *ax;
	
	fct->GetRange(x0, x1);
	b0 = h->FindBin(x0);
	b1 = h->FindBin(x1);
	nbins = h->GetXaxis()->GetNbins();
	
	// f0 und f1 bestimmen durch Mittelung Ÿber 3 bins
	if(--b0<1) b0 = 1;
	if(++b1>nbins) b1 = nbins;
	
	f0 = h->Integral(b0, b0+2) / 3.0;
	f1 = h->Integral(b1-2, b1) / 3.0;
	if (f0 == 0.0 || f1 == 0.0) {
		left = right = 0.0;
	} else {
		right = log(f0/f1)/(x1-x0);
		left = f0*exp(right*x0);
	}
	
	// adjust the gaussian parameters
	ax = h->GetXaxis();
	middle = (ax->GetBinUpEdge(ax->GetNbins()) + ax->GetBinLowEdge(1)) / 2.0;
	width = (ax->GetBinUpEdge(ax->GetNbins()) - ax->GetBinLowEdge(1)) / 2.0;
	width /= 6.0; // wir suchen das maximum im mittleren drittel
	ax->SetRangeUser(middle - width, middle + width);
	max_bin = h->GetMaximumBin();
	ax->UnZoom();
	
	mu = h->GetBinCenter(max_bin);
	c = h->GetBinContent(max_bin) - left*exp(-right*mu);
	
	sigma = h->Integral(b0,b1,"width");
	if (right)
		sigma -= ( left/right*(exp(-right*x0) - exp(-right*x1)) );
	sigma *= 1.0/(sqrt(2.0*M_PI)*c);
	
	fct->SetParameters(c,mu,sigma,left,right);
	// sigma and gaussian normalization have to be positive
	fct->SetParLimits(0,0.0,99999.0);
	fct->SetParLimits(2,0.0,99999.0);
	
	fix_parameters(fct,fixed);
} // adjust_parameter_gauss_expo()

void adjust_parameter_crystall_ball(TH1D* h, TF1* fct, std::map<int,double>* fixed, double a, double b)
{
	double mean,sigma,alpha,n,N;
	int max_bin;
	map<int,double>::const_iterator it;
	
	max_bin = h->GetMaximumBin();
	
	N = h->GetBinContent(max_bin);
	mean = h->GetBinCenter(max_bin);
	sigma = h->Integral("width") / (TMath::Sqrt(TMath::TwoPi())*N);
	
	// does this work?
	alpha = 1.0;
	n = 1.0;
	
	fct->SetParameters(mean, sigma, alpha, n, N);
	fct->SetParLimits(2,0.0,99999.0);
	
	fix_parameters(fct, fixed);
} // adjust_parameter_crystall_ball()

// adjust the parameters given the fucntion f is a gaussian (c.f NCFunctions.cp)
void adjust_parameter_gaussian(TH1D* h, TF1* fct, map<int,double>* fixed, double a, double b)
{
	double sigma = h->GetRMS();
	
	fct->SetParameter(0, h->Integral("width")/(TMath::Sqrt(TMath::TwoPi()) * sigma)); // normalization constant
	fct->SetParameter(1, h->GetMean()); // mean
	fct->SetParameter(2, sigma); // standard deviation
	fct->SetParLimits(2, 0.0, 99999.0); // confine sigma to be positive.
	
	fix_parameters(fct, fixed);
} // adjust_parameter_gaussian()

void adjust_parameter_gauss_linear(TH1D* h, TF1* fct, map<int,double>* fixed, double a, double b)
{
	int b0,b1,nbins,max_bin;
	double x0,x1,f0,f1;
	double c,mu,sigma,m,a0; // steigung und offset von polynom
	double width,middle;
	TAxis *ax;
	
	fct->GetRange(x0, x1);
	
	b0 = h->FindBin(x0);
	b1 = h->FindBin(x1);
	nbins = h->GetXaxis()->GetNbins();
	
	if (--b0 < 1)		b0 = 1;
	if (++b1 > nbins)	b1 = nbins;
	
	// mitteln Ÿber 3 bins
	f0 = h->Integral(b0,b0+2) / 3.0;
	f1 = h->Integral(b1-2,b1) / 3.0;
	
	// lineare parameter bestimmen
	m = (f1 - f0)/(x1 - x0);
	a0 = f1 - m*x1;
	
	// gaussche parameter setzen
	ax = h->GetXaxis();
	middle = ( ax->GetBinUpEdge(ax->GetNbins()) + ax->GetBinLowEdge(1) )/ 2.0;
	width = ax->GetBinUpEdge(ax->GetNbins()) - ax->GetBinLowEdge(1);
	width /= 3.0; // one third
	ax->SetRangeUser(middle - width / 2.0, middle + width / 2.0);
	max_bin = h->GetMaximumBin();
	ax->UnZoom();
	
	mu = h->GetBinCenter(max_bin);
	c = h->GetBinContent(max_bin) - (m*mu + a0);
	
	sigma = h->Integral(b0,b1,"width");
	sigma -= 0.5*m*(x1*x1 - x0*x0) + a0*(x1-x0); // subtract integral of background
	sigma *= 1/(sqrt(2.0*M_PI)*c);
	
	fct->SetParameters(c,mu,sigma,a0,m);
	fct->SetParLimits(2,0.0,99999.0);
	
	fix_parameters(fct, fixed);
} // adjust_parameter_gauss_linear()

void adjust_parameter_gauss_p2(TH1D* h, TF1 *fct, std::map<int,double> *fixed, double a, double b)
{
	double x0,x1;
	int j,b0,b1,range,bin_max;
	double c,mu,sigma;
	TH1D *clone_hist =  NULL;
	TF1 *pol2 = NULL;
	TCanvas *canvas = new TCanvas;
	
	fct->GetRange(x0, x1);
	b0 = h->FindBin(x0);
	b1 = h->FindBin(x1);
	
	if (b0 < 1) b0 = 1;
	if (b1 > h->GetNbinsX()) b1 = h->GetNbinsX();
	
	clone_hist = (TH1D*)h->Clone("");
	range = (b1 - b0) / 4;
	for (j = b0+range; j <= b1-range; j++) {
		clone_hist->SetBinContent(j, 0.0);
		clone_hist->SetBinError(j, 0.0);
	}
	
	// check if the histogram is empty
	if (clone_hist->Integral() == 0) {
		pol2 = new TF1("",f_p2,x0,x1,3);
		pol2->SetParameters(0.0, 0.0, 0.0);
	} else {
		// do a pol2 fit with this histogram
		clone_hist->Fit("pol2","Q0");
		pol2 = clone_hist->GetFunction("pol2");
	}
	
	// go and look for the maximum in the middle third
	if (a < b)	h->GetXaxis()->SetRangeUser(a, b);
	else		h->GetXaxis()->SetRangeUser(x0 + (x1 - x0)/3.0, x1 - (x1 - x0)/3.0);
	bin_max = h->GetMaximumBin();
	h->GetXaxis()->UnZoom();
	
	mu = h->GetBinCenter(bin_max);
	if (pol2) {
		c = h->GetBinContent(bin_max) - pol2->Eval(mu);
	} else {
		c = 0.0;
	}
	
	if (c == 0.0) {
		sigma = 1.0;
	} else {
		sigma = 1/( sqrt(2*M_PI) * c ) * (h->Integral(b0,b1,"width") - pol2->Integral(x0, x1));
	}
	if (sigma <= 0) sigma = 0.04; // about resolution of CMS detector in our channels
	if (sigma > 0.1) sigma = 0.04;
	
	
	if (pol2)
		fct->SetParameters(c, mu, sigma, pol2->GetParameter(0), pol2->GetParameter(1), pol2->GetParameter(2));
	else
		fct->SetParameters(c, mu, sigma, 0.0, 0.0, 0.0);
	
	delete pol2;
	delete clone_hist;
	delete canvas;
} // adjust_parameter_gauss_p2()

void adjust_parameter_double_gauss_linear(TH1D* h, TF1* fct, map<int,double>* fixed, double a, double b)
{
	int b0,b1,nbins,max_bin;
	double x0,x1,f0,f1;
	double c,mu,sigma,m,a0; // steigung und offset von polynom
	double width,middle;
	TAxis *ax;
	
	fct->GetRange(x0, x1);
	
	b0 = h->FindBin(x0);
	b1 = h->FindBin(x1);
	nbins = h->GetXaxis()->GetNbins();
	
	if (--b0 < 1)		b0 = 1;
	if (++b1 > nbins)	b1 = nbins;
	
	// mitteln Ÿber 3 bins
	f0 = h->Integral(b0,b0+2) / 3.0;
	f1 = h->Integral(b1-2,b1) / 3.0;
	
	// lineare parameter bestimmen
	m = (f1 - f0)/(x1 - x0);
	a0 = f1 - m*x1;
	
	// gaussche parameter setzen
	ax = h->GetXaxis();
	middle = ( ax->GetBinUpEdge(ax->GetNbins()) + ax->GetBinLowEdge(1) )/ 2.0;
	width = ax->GetBinUpEdge(ax->GetNbins()) - ax->GetBinLowEdge(1);
	width /= 3.0; // one third
	ax->SetRangeUser(middle - width / 2.0, middle + width / 2.0);
	max_bin = h->GetMaximumBin();
	ax->UnZoom();
	
	mu = h->GetBinCenter(max_bin);
	c = h->GetBinContent(max_bin) - (m*mu + a0);
	
	sigma = h->Integral(b0,b1,"width");
	sigma -= 0.5*m*(x1*x1 - x0*x0) + a0*(x1-x0); // subtract integral of background
	sigma *= 1/(sqrt(2.0*M_PI)*c);
	
	fct->SetParameters(c/2.0,mu,sigma,c/2.0,mu,sigma,a0,m);
	fct->SetParLimits(2, 0.0, 1.0E7);
	fct->SetParLimits(5, 0.0, 1.0E7);
	
	fix_parameters(fct, fixed);
} // adjust_parameter_double_gauss_linear()

double signal_events(TF1 *signalFct, double mu, double sigma, double bin_width)
{
	return signalFct->Integral(mu - 2*sigma, mu + 2*sigma) / bin_width;
} // signal_events()

double signal_events_gauss_linear(TF1 *fct, double bin_width)
{
	double params[5];
	double result = 0.0;
	scoped_ptr<TF1> gauss;
	
	fct->GetParameters(params);
	
	gauss.reset(new TF1("",f_gauss,fct->GetXmin(),fct->GetXmax(),3));
	gauss->SetParameters(params);
	
	result = signal_events(&(*gauss), params[1], params[2], bin_width);
	
	return result;
} // signal_events_gauss_linear()

double signal_events_gauss_linear_histo(TH1D *h, double xmin, double xmax)
{
	fit_t fit;
	
	init_gauss_linear(&fit, xmin, xmax);
	adjust_parameter_double_gauss_linear(h, fit.fit_fct);
	h->Fit(fit.fit_fct, "LR");
	
	return signal_events_gauss_linear(fit.fit_fct, h->GetBinWidth(1));
} // signal_events_gauss_linear_histo()

double signal_events_gauss_p2(TF1 *fct, double bin_width, double *signal_err)
{
	double params[6];
	double result;
	double error;
	scoped_ptr<TF1> gauss;
	
	fct->GetParameters(params);
	
	gauss.reset(new TF1("",f_gauss,fct->GetXmin(),fct->GetXmax(),3));
	gauss->SetParameters(params);
	result = signal_events(&(*gauss), params[1], params[2], bin_width);
	
	params[0] = fct->GetParError(0); // error on C
	gauss->SetParameters(params);
	error = signal_events(&(*gauss), params[1], params[2], bin_width);
	
	if (signal_err)
		*signal_err = error;
	
	return result;
} // signal_events_gauss_p2()

double signal_events_gauss_p2_histo(TH1D* h, double xmin, double xmax)
{
	fit_t fit;
	
	init_gauss_p2(&fit, xmin, xmax);
	adjust_parameter_gauss_p2(h, fit.fit_fct);
	h->Fit(fit.fit_fct, "LR");
	
	return signal_events_gauss_p2(fit.fit_fct, h->GetBinWidth(1));
} // signal_events_gauss_p2_histo()

double signal_events_gauss_expo(TF1 *fct, double bin_width)
{
	double params[5];
	double result = 0.0;
	scoped_ptr<TF1> gauss;
	
	fct->GetParameters(params);
	gauss.reset(new TF1("",f_gauss,fct->GetXmin(),fct->GetXmax(),3));
	gauss->SetParameters(params);
	
	result = signal_events(&(*gauss), params[1], params[2], bin_width);
	
	return result;
} // signal_events_gauss_expo()

double signal_events_gauss_expo_histo(TH1D *h, double xmin, double xmax)
{
	fit_t fit;
	
	init_gauss_expo(&fit, xmin, xmax);
	adjust_parameter_gauss_expo(h, fit.fit_fct);
	h->Fit(fit.fit_fct,"LR");
	
	return signal_events_gauss_expo(fit.fit_fct, h->GetBinWidth(1));
} // signal_events_gauss_expo_histo()

double background_events_gauss_linear(TF1 *fct, double bin_width)
{
	double params[5];
	double result = 0.0;
	scoped_ptr<TF1> pol1;
	
	fct->GetParameters(params);
	pol1.reset(new TF1("",f_p1,fct->GetXmin(),fct->GetXmax(),2));
	pol1->SetParameters(params+3);
	
	result = signal_events(&(*pol1), params[1], params[2], bin_width);
	
	return result;
} // background_events_gauss_linear()

double background_events_gauss_p2(TF1 *fct, double bin_width)
{
	double params[6];
	double result = 0.0;
	scoped_ptr<TF1> pol2;
	
	fct->GetParameters(params);
	pol2.reset(new TF1("",f_p2,fct->GetXmin(),fct->GetXmax(),3));
	pol2->SetParameters(params + 3);
	
	result = signal_events(&(*pol2), params[1], params[2], bin_width);
	
	return result;
} // background_events_gauss_p2()

double background_events_gauss_expo(TF1 *fct, double bin_width)
{
	double params[5];
	double result = 0.0;
	scoped_ptr<TF1> expo;
	
	fct->GetParameters(params);	
	expo.reset(new TF1("",f_expo,fct->GetXmin(),fct->GetXmax(),2));
	expo->SetParameters(params+3);
	
	result = signal_events(&(*expo), params[1], params[2], bin_width);
	
	return result;
} // background_events_gauss_expo()

double significance(TF1 *signalFct, TF1 *backFct, double mu, double sigma, double bin_width)
{
	double s = signal_events(signalFct, mu, sigma, bin_width);
	double b = signal_events(backFct, mu, sigma, bin_width);
	
	return s / sqrt(s + b);
} // significance()

double significance_gauss_linear(TF1 *fct, double bin_width)
{
	double params[5];
	double result = 0.0;
	scoped_ptr<TF1> gauss;
	scoped_ptr<TF1> pol1;
	
	fct->GetParameters(params);
	
	gauss.reset(new TF1("",f_gauss,fct->GetXmin(),fct->GetXmax(),3));
	gauss->SetParameters(params);
	
	pol1.reset(new TF1("",f_p1,fct->GetXmin(),fct->GetXmax(),2));
	pol1->SetParameters(params+3);
	
	result = significance( &(*gauss), &(*pol1), params[1], params[2], bin_width);
	
	return result;
} // significance_gauss_linear()

double significance_gauss_p2(TF1 *fct, double bin_width)
{
	double result = 0.0;
	scoped_ptr<TF1> gauss(new TF1("",f_gauss,fct->GetXmin(),fct->GetXmax(),3));
	scoped_ptr<TF1> pol2(new TF1("",f_p2,fct->GetXmin(),fct->GetXmax(),3));
	
	gauss->SetParameters(fct->GetParameters());
	pol2->SetParameters(fct->GetParameters()+3);
	
	result = significance(&(*gauss), &(*pol2), fct->GetParameter(1), fct->GetParameter(2), bin_width);
	
	return result;
} // background_events_gauss_p2()

double significance_gauss_expo(TF1 *fct,double bin_width)
{
	double params[5];
	double result = 0.0;
	scoped_ptr<TF1> gauss;
	scoped_ptr<TF1> expo;
	
	fct->GetParameters(params);
	gauss.reset(new TF1("",f_gauss,fct->GetXmin(),fct->GetXmax(),3));
	gauss->SetParameters(params);
	
	expo.reset(new TF1("",f_expo,fct->GetXmin(),fct->GetXmax(),2));
	expo->SetParameters(params+3);
	
	result = significance(&(*gauss), &(*expo), params[1], params[2], bin_width);
	
	return result;
} // significance_gauss_expo()

void parse_variable_config(const char *config_file, process_fct f, void *param)
{
	// parse the config file and call the function with each entry!!
	FILE *file = fopen(config_file,"r");
	char buffer[1024];
	TString line;
	TString delim("\t");
	TObjArray *tokens;
	function_range_t tok;
	Int_t j;
	
	while (fgets(buffer, sizeof(buffer), file)) {
		
		if (buffer[strlen(buffer)-1] != '\n') {
			std::cerr << "parse_variable_config(): line not ending with newline. bailling out..." << std::endl;
			break;
		} else {
			buffer[strlen(buffer)-1] = '\0'; // cut of the trailing new line
		}
		
		if (buffer[0] == '#') continue; // comment line
		line = TString(buffer);
		tokens = line.Tokenize(delim);
		
		tok.name = ((TObjString*)((*tokens)[0]))->GetString(); // first one is variable name
		tok.vals.clear();
		tok.options = "";
		
		// last one is options string, all other fill as double in the vals vector
		for (j = 1; j < tokens->GetEntries() - 1; j++) {
			
			if (((TObjString*)((*tokens)[j]))->GetString().CompareTo("DBL_MAX") == 0)
				tok.vals.push_back(DBL_MAX);
			else if (((TObjString*)((*tokens)[j]))->GetString().CompareTo("DBL_MIN") == 0)
				tok.vals.push_back(DBL_MIN);
			else
				tok.vals.push_back(((TObjString*)((*tokens)[j]))->GetString().Atof());
		}
		
		// option string
		if (tokens->GetEntries() > 1)
			tok.options = std::string(((TObjString*)tokens->Last())->GetString().Data());
		
		// execute the function with the tokens!
		f(&tok,param);
	}
} // parse_variable_config()

// plot all the variables in the config file...
// Format for plotting is
//	variable	min	max
//	pt	0	60
// where the range is optional
void plot_vars(TTree *tree, const char *var_file, TCut preCut, const char *ending)
{
	// CAN BE REWRITTEN USING THE ABOVE PARSER!!
	using std::cout;
	using std::endl;
	FILE *config_file = NULL;
	TString str;
	TObjArray *tokens;
	TString name;
	char buffer[1024];
	double min,max;
	bool with_range;
	TString delim("\t ");
	TH1 *h = NULL;
	TString outfilename;
	TString drawVar,histName;
	TCanvas *c = new TCanvas;
	
	config_file = fopen(var_file,"r");
	
	// read each line
	while (fgets(buffer,sizeof(buffer),config_file)) {
		
		if (buffer[0] == '#') // comment line
			continue;
		
		str.Append(buffer);
		if (str[str.Length()-1] == '\n') {
			
			// end of line reached
			str.Remove(str.Length()-1);
			
			// tokenize by tabs
			tokens = str.Tokenize(delim);
			
			// First one is the variable name...
			name = ((TObjString*)((*tokens)[0]))->GetString();
			cout << "Processing variable: " << name.Data() << endl;
			
			if (tokens->GetEntries() >= 3) {
				with_range = true;
				min = ((TObjString*)((*tokens)[1]))->GetString().Atof();
				max = ((TObjString*)((*tokens)[2]))->GetString().Atof();
				cout << "\tWith forced range [" << min << ", " << max << "]" << endl;
			} else
				with_range = false;
			
			
			if (with_range) {
				histName = name;
				histName.ReplaceAll("/",".");
				h = new TH1D(histName.Data(),name.Data(),100,min,max);
			}
			
			// draw the
			drawVar = TString(Form("%s%s%s",name.Data(),with_range ? " >> " : "",with_range ? histName.Data() : ""));
			cout << "tree->Draw(" << drawVar.Data() << ", " << preCut.GetTitle() << " );" << endl;
			tree->Draw(drawVar.Data(),preCut);
						
			if (!with_range)
				h = (TH1*)gROOT->FindObject("htemp");
			
			// change drawing options
			h->SetFillColor(kRed);
			h->SetLineColor(kRed);
			h->SetFillStyle(3004);
			h->Draw();
			
			// save
			outfilename = name.ReplaceAll("/",".");
			c->SaveAs(Form("%s.%s",outfilename.Data(),ending));
			
			if (with_range)
				delete h;
						
			// cleanup
			str = TString(); // empty the string!
			tokens->Delete(); // free the memory
			delete tokens; // delete the array
		}
	}
	
	delete c;
	
	fclose(config_file);
} // plot_vars

static void plot_tree_with_token(TTree *tree,TString &token, const char *histoname, const char *tofile, TCanvas *c)
{
	TObjArray *lines = NULL;
	TString cut;
	Int_t j;
	
	// parse each line...
	lines = token.Tokenize(TString("\n\r"));
	for (j = 0; j < lines->GetEntries(); j++) {
		
		TString line = ((TObjString*)(*lines)[j])->GetString();
		
		if (j != 0)
			cut.Append(" && ");
		cut.Append(line);
	}
	
	// now, draw the histogram
	tree->Draw(Form("mass >> %s",histoname), cut.Data());
	c->SaveAs(tofile);
	
	// clean up
	if (lines)
		lines->Delete();
	
	delete lines;
} // plot_tree_with_token()

// plot a number of histograms using a config file
// This is the same plot with different cuts!! (e.g. supplied by TMVA)
void plot_tree_with_config(TTree *tree, const char *config_file, const char *histoname, const char *dir)
{
	FILE *cfile = fopen(config_file,"r");
	char buffer[1024];
	TString line;
	TString token;
	int counter = 0;
	char filename[256];
	TCanvas *c = new TCanvas;
	
	if(!cfile)
		goto bail;
	
	// read the config file and separate by lines beginning with ---
	while (fgets(buffer,sizeof(buffer),cfile)) {
		
		line.Append(buffer);
		
		// append to the 'token' until we reach a line beginning with -
		if (line[line.Length()-1] != '\n') // not yet the end of line reached
			continue;
		
		if (line[0] == '-') {
			// we have read a whole token. process the token and reset
			sprintf(filename, "%s%s%s_%d.pdf", dir ? dir : "", dir ? "/" : "", histoname, counter++);
			plot_tree_with_token(tree,token,histoname,filename,c);
			
			// clean the strings
			token = TString();
			line = TString();
		}
		else { // append the line to the token
			token.Append(line);
			line = TString();
		}
	}
	
	sprintf(filename, "%s_%d.pdf", histoname, counter++);
	plot_tree_with_token(tree,token,histoname,filename,c);
bail:
	if(cfile)
		fclose(cfile);
	
	delete c;
} // plot_tree_with_config()
