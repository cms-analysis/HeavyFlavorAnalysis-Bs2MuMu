/*
 *  ul_main.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

// Standard Headers
#include <iostream>
#include <string>
#include <vector>

// ROOT headers
#include <TMath.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooRandom.h>
#include <RooStats/ConfInterval.h>

#include "external_constants.h"
#include "ul_estimate.h"

using namespace std;

enum algo_t {
	kAlgo_Bayesian = 1,
	kAlgo_FeldmanCousins,
	kAlgo_CLs,
	kAlgo_CLb,
	kAlgo_Hybrid,
	kAlgo_CLb_Hybrid,
	kAlgo_Zbi,
	kAlgo_Interval_Hybrid,
	kAlgo_Bkg_Hybrid,
	kAlgo_None
};

static const char *workspace_path = NULL;
static const char *configfile_path = NULL;
static double gCLLevel = 0.95;
static pair<double,double> gMuSRange(-1,-1);
static uint32_t gFCSteps = 20;
static algo_t gAlgorithm = kAlgo_Bayesian;
static bool gDisableErrors = false;
static int gVerbosity = 1; // 0 = quite, 1 = standard, 2 = verbose
static const char *output_path = NULL;
static bool gBdToMuMu = false;
static bool gSeed = false;
static uint32_t gProofWorkers = 0;
static int gToys = 50000;
static bool gFixBackground = false;
static bool gFloatingData = false;
static bool gSMExpectation = false;
static bool gSMCrossFeed = false;

static const char *algo_name(algo_t a)
{
	const char *result;
	
	switch (a) {
		case kAlgo_Bayesian:
			result = "Bayesian";
			break;
		case kAlgo_FeldmanCousins:
			result = "Feldman-Cousins";
			break;
		case kAlgo_CLs:
			result = "CLs";
			break;
		case kAlgo_CLb:
			result = "CLb";
			break;
		case kAlgo_Hybrid:
			result = "Hybrid";
			break;
		case kAlgo_CLb_Hybrid:
			result = "CLb using Hybrid";
			break;
		case kAlgo_Interval_Hybrid:
			result = "Hybrid Two Sided Interval";
			break;
		case kAlgo_Zbi:
			result = "Zbi";
			break;
		case kAlgo_None:
			result = "None";
			break;
		case kAlgo_Bkg_Hybrid:
			result = "Bkg using Hybrid";
			break;
		default:
			cerr << "Unknown algorithm selected: " << a << endl;
			result = NULL;
			abort();
			break;
	}
	
	return result;
} // algo_name()

static bool validate_bmm(map<bmm_param,measurement_t> *bmm, int ch)
{
	bool valid;
	
	// check normalization decay info
	valid = bmm->count(make_pair(kTot_bplus, ch)) > 0;
	valid = valid || ( bmm->count(make_pair(kAcc_bplus, ch)) > 0 && bmm->count(make_pair(kEff_mu_bplus, ch)) > 0
					&& bmm->count(make_pair(kEff_trig_bplus, ch)) > 0 && bmm->count(make_pair(kEff_cand_bplus, ch)) > 0
					&& bmm->count(make_pair(kEff_ana_bplus, ch)) > 0 );
	if (!valid) {
		cerr << "Insufficient information in config file for channel " << ch << ": Supply either 'TOT_BPLUS' or ('ACC_BPLUS' and 'EFF_MU_BPLUS' and 'EFF_TRIG_BPLUS' and 'EFF_CAND_BPLUS' and 'EFF_ANA_BPLUS')" << endl;
		goto bail;
	}
	
	valid = bmm->count(make_pair(kObsBkg_bmm, ch)) > 0;
	if (!valid) {
		cerr << "Insufficient information in config file. Please specify 'OBS_BKG' for channel " << ch << endl;
		goto bail;
	}
	
bail:
	return valid;
} // validate_bmm()

static bool process_line(string name, vector<double> *values, map<bmm_param,measurement_t> *bsmm, map<bmm_param,measurement_t> *bdmm)
{
	bmm_param p;
	bool bsparam;
	bool success = false;
	
	if (values->size() <= 1)
		goto bail;
	
	p.first = find_bmm_param_by_name(name, &bsparam);
	if (p.first == kUnknownParam) goto bail;
	
	// Channel Index
	p.second = (int)((*values)[0]);
	
	// Value
	if (bsparam)	(*bsmm)[p].setVal((*values)[1]);
	else			(*bdmm)[p].setVal((*values)[1]);
	
	// error on variable
	if (values->size() == 3) {
		if (bsparam) {
			(*bsmm)[p].setErrHi((*values)[2]);
			(*bsmm)[p].setErrLo((*values)[2]);
		} else {
			(*bdmm)[p].setErrHi((*values)[2]);
			(*bdmm)[p].setErrLo((*values)[2]);
		}
	} else if (values->size() >= 4) {
		if (bsparam) {
			(*bsmm)[p].setErrHi((*values)[2]);
			(*bsmm)[p].setErrLo((*values)[3]);
		} else {
			(*bdmm)[p].setErrHi((*values)[2]);
			(*bdmm)[p].setErrLo((*values)[3]);
		}
	}
	
	success = true;
bail:
	return success;
} // process_line()

static void parse_input(const char *path, map<bmm_param,measurement_t> *bsmm, map<bmm_param,measurement_t> *bdmm)
{
	FILE *file = fopen(path, "r");
	char* input_line;
	size_t line_len;
	string line;
	string delim(" \t:");
	vector<double> values;
	string name;
	string::iterator it,next;
	
	while ( (input_line = fgetln(file, &line_len)) ) {
		
		// if comment line, read next!
		if (input_line[0] == '#')
			continue;
		
		// remove endline at end of line (if any)
		if (input_line[line_len-1] == '\n' || input_line[line_len-1] == '\r')
			line_len--;
		
		if (line_len == 0) // empty line
			continue;
		
		// create a c++ string out of input line
		line = string(input_line,line_len);
		
		// tokenize by delim
		it = find_first_of(line.begin(), line.end(), delim.begin(), delim.end());
		name = string(line,0,it-line.begin());
		
		// now, search for all tokens
		values.clear();
		
		do {
			next = find_first_of(it + 1, line.end(), delim.begin(), delim.end());
			if ((next - it) > 1) {
				values.push_back(atof( string(line, it-line.begin(), next - it).c_str() ));
			}
			it = next;
		} while (next != line.end());
		
		if (!process_line(name,&values,bsmm,bdmm)) {
			cerr << "Unable to process line: '" << line << "'." << endl;
			abort();
		}
	}
	
	fclose(file);
} // parse_input()

static void usage()
{
	cerr << "ulcalc [--SM-crossfeed][--SM-exp] [--float-data] [--fixed-bkg] [--toys <NbrMCToys>] [--proof <nbr_workers>] [--seed] [--bdtomumu] [--disable-errors] [[-n <nbr steps>] -r x,y] [-l cl] [-w workspace_outfile.root] [-a <\"bayes\"|\"fc\"|\"cls\"|\"clb\"|\"hybrid\"|\"int_hybrid\"|\"clb_hybrid\"|\"zbi\"|\"bkg\"|\"none\">] [-q] [-v] [-o <outputfile>] <configfile>" << endl;
} // usage()

static bool parse_arguments(const char **first, const char **last)
{
	bool ok = true;
	const char *arg;
	string s;
	string::iterator col;
	
	while (first != last) {
		
		arg = *first++;
		if (arg[0] == '-') {
			if (strcmp(arg, "-l") == 0) {
				if (first == last) {
					usage();
					abort();
				}
				gCLLevel = atof(*first++); // read the confidence level
			} else if (strcmp(arg, "-w") == 0) {
				if (first == last) {
					usage();
					abort();
				}
				workspace_path = *first++; // set the workspace output path
			} else if (strcmp(arg, "-n") == 0) {
				if (first == last) {
					usage();
					abort();
				}
				gFCSteps = (uint32_t)atoi(*first++);
			} else if (strcmp(arg, "-r") == 0) {
				if (first == last) {
					usage();
					abort();
				}
				s = *first++;
				if ( (col = find(s.begin(), s.end(), ',')) == s.end() ) {
					usage();
					abort();
				}
				gMuSRange.first = atof(string("").append(s.begin(),col).c_str());
				gMuSRange.second = atof(string("").append(col+1,s.end()).c_str());
			} else if (strcmp(arg, "-a") == 0) {
				if (first == last) {
					usage();
					abort();
				}
				s = *first++;
				if (s.compare("bayes") == 0)	gAlgorithm = kAlgo_Bayesian;
				else if (s.compare("fc") == 0) {
					gAlgorithm = kAlgo_FeldmanCousins;
				} else if (s.compare("cls") == 0) {
					gAlgorithm = kAlgo_CLs;
				} else if (s.compare("clb") == 0) {
					gAlgorithm = kAlgo_CLb;
				} else if (s.compare("hybrid") == 0) {
					gAlgorithm = kAlgo_Hybrid;
				} else if (s.compare("clb_hybrid") == 0) {
					gAlgorithm = kAlgo_CLb_Hybrid;
				} else if (s.compare("int_hybrid") == 0) {
					gAlgorithm = kAlgo_Interval_Hybrid;
				} else if (s.compare("zbi") == 0) {
					gAlgorithm = kAlgo_Zbi;
				} else if (s.compare("bkg") == 0) {
					gAlgorithm = kAlgo_Bkg_Hybrid;
				} else if (s.compare("none") == 0) {
					gAlgorithm = kAlgo_None;
				} else {
					usage();
					abort();
				}
			} else if (strcmp(arg, "--disable-errors") == 0) {
				gDisableErrors = true;
			} else if (strcmp(arg, "-v") == 0) {
				gVerbosity = 2; // verbose
			} else if (strcmp(arg, "-q") == 0) {
				gVerbosity = 0; // quite mode
			} else if (strcmp(arg, "-o") == 0) {
				if (first == last) {
					usage();
					abort();
				}
				output_path = *first++;
			} else if (strcmp(arg, "--bdtomumu") == 0) {
				gBdToMuMu = true;
			} else if (strcmp(arg, "--seed") == 0) {
				gSeed = true;
			} else if (strcmp(arg, "--proof") == 0) {
				if (first == last) {
					usage();
					abort();
				}
				gProofWorkers = atoi(*first++);
			} else if (strcmp(arg, "--toys") == 0) {
				if (first == last) {
					usage();
					abort();
				}
				gToys = atoi(*first++);
			} else if (strcmp(arg, "--fixed-bkg") == 0) {
				gFixBackground = true;
			} else if (strcmp(arg, "--float-data") == 0) {
				gFloatingData = true;
			} else if (strcmp(arg, "--SM-exp") == 0) {
				gSMExpectation = true;
			} else if (strcmp(arg, "--SM-crossfeed") == 0) {
				gSMCrossFeed = true;
			} else {
				cerr << "Unknown option '" << arg << "'." << endl;
				usage();
				abort();
			}
		} else
			configfile_path = arg;
	}
	
	if (!configfile_path) {
		ok = false;
		goto bail;
	}
	
	// Output current configuration
	cout << "Reading configuration file: " << configfile_path << endl;
	cout << "Confidence Level: " << gCLLevel << endl;
	if (output_path)
		cout << "Output File: " << output_path << endl;
	if (workspace_path)
		cout << "Output Workspace into: " << workspace_path << endl;
	if (gAlgorithm == kAlgo_FeldmanCousins || gAlgorithm == kAlgo_CLs || gAlgorithm == kAlgo_Hybrid || gAlgorithm == kAlgo_Interval_Hybrid) {
		cout << "mu_s range: (" << gMuSRange.first << ", " << gMuSRange.second << ")." << endl;
		cout << "Number of Steps " << gFCSteps << endl;
	}
	cout << "Algorithm is " << algo_name(gAlgorithm) << endl;
	cout << "Errors on Variables are " << (gDisableErrors ? "disabled" : "enabled") << endl;
	cout << "Verbosity level " << gVerbosity << endl;
	cout << "Computing Upper limit for decay " << (gBdToMuMu ? "Bd -> mumu" : "Bs -> mumu") << endl;
	if (gSeed) cout << "Unique Random seed" << endl;
	cout << "Number of MC Toys used: " << gToys << endl;
	if (gProofWorkers > 0) cout << "Number of Proof Workers: " << gProofWorkers << endl;
	if (gFloatingData) cout << "Non Integer Observation allowed" << endl;
	if (gSMCrossFeed) cout << "Cross Feed fixed to SM" << endl;
	cout << "-------------------------------------" << endl;
	if (gFixBackground)	cout << "Toys with fixed sideband window" << endl;
	else				cout << "Toys with floating sideband window" << endl;
	
bail:
	return ok;
} // parse_arguments()

static void dump_params(pair<bmm_param,measurement_t> p)
{
	cout << '\t' << find_bmm_name(p.first.first) << '[' << p.first.second << "]: " << p.second.getVal();
	if (p.second.getErrHi() > 0 || p.second.getErrLo() > 0) {
		cout << "+" << p.second.getErrHi() << "-" << p.second.getErrLo();
	}
	cout << endl;
} // dump_params()

static void recursive_calc(RooWorkspace *wspace, RooArgSet *obs, set<int> *channels, set<int>::iterator a,set<int>::iterator b, map<bmm_param,measurement_t> *bsmm, map<bmm_param,measurement_t> *bdmm, double *upperLimit, double *lowerLimit)
{
	RooDataSet *data = NULL;
	double bkg;
	double obsBs,obsBd,bbb;
	set<int>::const_iterator it;
	RooAbsReal *bs_mean;
	RooAbsReal *bd_mean;
	RooArgSet *saved_vars = NULL;
	int ch;
	
	if (a == b) { // end of recursion over channels
		
		data = new RooDataSet("data","",*obs);
		data->add(*obs);
		if (gVerbosity > 0)
			data->Print("v");
		
		wspace->import(*data);
		
		// save the variable state s.t. it does not get lost for the weight computation in the next round
		saved_vars = dynamic_cast<RooArgSet*>(wspace->allVars().snapshot());
		
		// set the gamma bkg priors if there are...
		for (it = channels->begin(); it != channels->end(); ++it) {
			if (wspace->var(Form("gamma_%d",*it)))
				wspace->var(Form("gamma_%d",*it))->setVal(((RooRealVar&)((*data->get(0))[Form("NbObs_%d",*it)])).getVal()+1);
		}
		
		switch (gAlgorithm) {
			case kAlgo_Bayesian:
				est_ul_bc(wspace, data, channels, gCLLevel, gVerbosity, upperLimit, NULL);
				break;
			case kAlgo_FeldmanCousins:
				est_ul_fc(wspace, data, channels, gCLLevel, gVerbosity, upperLimit, lowerLimit, ((gMuSRange.first >= 0) ? &gMuSRange : NULL), &gFCSteps, NULL, gProofWorkers, gToys);
				break;
			case kAlgo_CLs:
				est_ul_cls(wspace, data, channels, gCLLevel, gVerbosity, upperLimit, ((gMuSRange.first >= 0) ? &gMuSRange : NULL), &gFCSteps, NULL, gProofWorkers, gToys);
				break;
			case kAlgo_CLb:
				// note, here upper limit represents p-value of background model.
				est_ul_clb(wspace, data, channels, gVerbosity, upperLimit, gProofWorkers, gToys);
				break;
			case kAlgo_Hybrid:
				est_ul_hybrid(wspace, data, channels, gCLLevel, gVerbosity, upperLimit, ((gMuSRange.first >= 0) ? &gMuSRange : NULL), &gFCSteps, NULL, gProofWorkers, gToys, gBdToMuMu, gFixBackground, gSMExpectation, gSMCrossFeed);
				break;
			case kAlgo_Interval_Hybrid:
				est_int_hybrid(wspace, data, channels, gCLLevel, gVerbosity, upperLimit, ((gMuSRange.first >= 0) ? &gMuSRange : NULL), &gFCSteps, NULL, gProofWorkers, gToys, gBdToMuMu, gFixBackground, gSMExpectation, gSMCrossFeed);
				break;
			case kAlgo_CLb_Hybrid:
				// note, here upper limit represents p-value of background model.
				est_ul_clb_hybrid(wspace, data, channels, gVerbosity, upperLimit, gProofWorkers, gToys, gBdToMuMu, gFixBackground, gSMExpectation, gSMCrossFeed);
				break;
			case kAlgo_Zbi:
				est_ul_zbi(wspace,data,channels,gCLLevel,gBdToMuMu,upperLimit);
				break;
			case kAlgo_Bkg_Hybrid:
				est_sign_bkg(wspace, data, channels, gVerbosity, upperLimit, gProofWorkers, gToys, gFixBackground);
				break;
			case kAlgo_None:
				measure_params(wspace, data, channels, gVerbosity);
				break;
			default:
				cerr << "Unknown algorithm selected to determine the upper limit..." << endl;
				abort();
				break;
		}
		wspace->allVars() = *saved_vars; // restore the values
		
		delete data;
		goto bail;
	}
	
	// save current channel and advance iterator
	ch = *a++;
	
	// Set current value and do recursive...
	bkg = ((*bsmm)[make_pair(kObsBkg_bmm, ch)]).getVal();
	bbb = bkg - wspace->var(Form("PeakBkgSB_%d",ch))->getVal();
	if (bbb < 0) bbb = 0;
	if (!gFixBackground) ((RooRealVar&)((*obs)[Form("NbObs_%d",ch)])).setVal(bkg);
	wspace->var(Form("NbObs_%d",ch))->setVal(bkg);
	wspace->var(Form("nu_b_%d",ch))->setVal(bbb);
	
	
	if (bsmm->count(make_pair(kObsB_bmm, ch)) > 0) {
		obsBs = ((*bsmm)[make_pair(kObsB_bmm, ch)]).getVal();
		if (!gFloatingData) obsBs = (uint32_t)obsBs;
	} else {
		bs_mean = wspace->function(Form("bs_mean_%d",ch));
		obsBs = bs_mean->getVal();
		if (!gFloatingData) obsBs = round(obsBs);
	}
	
	if (bdmm->count(make_pair(kObsB_bmm, ch)) > 0) {
		obsBd = ((*bdmm)[make_pair(kObsB_bmm, ch)]).getVal();
		if (!gFloatingData) obsBd = (uint32_t)obsBd;
	} else {
		bd_mean = wspace->function(Form("bd_mean_%d",ch));
		obsBd = bd_mean->getVal();
		if (!gFloatingData) obsBd = round(obsBd);
	}
	
	// iterate through these variables
	((RooRealVar&)((*obs)[Form("NsObs_%d",ch)])).setVal(obsBs);
	((RooRealVar&)((*obs)[Form("NdObs_%d",ch)])).setVal(obsBd);
	recursive_calc(wspace,obs,channels,a,b,bsmm,bdmm,upperLimit,lowerLimit);
bail:
	if(saved_vars)
		delete saved_vars;
	
	return;
} // recursive_calc()

int main(int argc, const char *argv [])
{
	RooWorkspace *wspace = NULL;
	RooArgSet observables;
	map<bmm_param,measurement_t> bsmm,bdmm;
	set<int> channels;
	set<int>::const_iterator ch;
	double upperLimit = 0,lowerLimit = 0;
	
	if(!parse_arguments(&argv[1], &argv[argc])) {
		usage();
		abort();
	}
	
	parse_input(configfile_path,&bsmm,&bdmm);
	
	// some variables have to be redundant for convenience
	add_channels(&bsmm, &channels);
	add_channels(&bdmm, &channels);
	for (set<int>::const_iterator c = channels.begin(); c!=channels.end();++c) {
		bdmm[make_pair(kObsBkg_bmm, *c)] = bsmm[make_pair(kObsBkg_bmm, *c)];
		
		if (bsmm.count(make_pair(kTot_bplus, *c)) > 0)
			bdmm[make_pair(kTot_bplus, *c)] = bsmm[make_pair(kTot_bplus, *c)];
		if (bsmm.count(make_pair(kAcc_bplus, *c)) > 0)
			bdmm[make_pair(kAcc_bplus, *c)] = bsmm[make_pair(kAcc_bplus, *c)];
		if (bsmm.count(make_pair(kEff_mu_bplus, *c)) > 0)
			bdmm[make_pair(kEff_mu_bplus, *c)] = bsmm[make_pair(kEff_mu_bplus, *c)];
		if (bsmm.count(make_pair(kEff_trig_bplus, *c)) > 0)
			bdmm[make_pair(kEff_trig_bplus, *c)] = bsmm[make_pair(kEff_trig_bplus, *c)];
		if (bsmm.count(make_pair(kEff_cand_bplus, *c)) > 0)
			bdmm[make_pair(kEff_cand_bplus, *c)] = bsmm[make_pair(kEff_cand_bplus, *c)];
		if (bsmm.count(make_pair(kEff_ana_bplus, *c)) > 0)
			bdmm[make_pair(kEff_ana_bplus, *c)] = bsmm[make_pair(kEff_ana_bplus, *c)];
		if (bsmm.count(make_pair(kObs_bplus, *c)) > 0)
			bdmm[make_pair(kObs_bplus, *c)] = bsmm[make_pair(kObs_bplus, *c)];
		if (bsmm.count(make_pair(kPeakBkgOff_bmm, *c)) > 0)
			bdmm[make_pair(kPeakBkgOff_bmm, *c)] = bsmm[make_pair(kPeakBkgOff_bmm, *c)];
		
		// check given parameters
		if(!validate_bmm(&bsmm,*c) || !validate_bmm(&bdmm,*c)) {
			cerr << "Insufficient parameters in the parameter file!!" << endl;
			abort();
		}
	}
	
	// show the input parsed
	if (gVerbosity > 0) {
		cout << "---------------------------" << endl;
		cout << "Bs -> mumu parameters read:" << endl;
		for_each(bsmm.begin(), bsmm.end(), dump_params);
		cout << "---------------------------" << endl;
		cout << "Bd -> mumu parameters read:" << endl;
		for_each(bdmm.begin(), bdmm.end(), dump_params);
		cout << "---------------------------" << endl;
	}
	
	compute_vars(&bsmm,true);
	compute_vars(&bdmm,false);
	wspace = build_model_nchannel(&bsmm, &bdmm, gDisableErrors, gVerbosity, gBdToMuMu, gFixBackground, gFloatingData, gSMCrossFeed, (gAlgorithm == kAlgo_Bkg_Hybrid));
	
	// initialize random number generator seed
	if (gSeed) RooRandom::randomGenerator()->SetSeed(0);
	
	// set the measured candidates...
	observables.addClone(*wspace->set("obs"));
	recursive_calc(wspace,&observables,&channels,channels.begin(),channels.end(),&bsmm,&bdmm,&upperLimit,&lowerLimit);
	
	if (gAlgorithm == kAlgo_CLb || gAlgorithm == kAlgo_CLb_Hybrid) {
		cout << "p value for background model: " << upperLimit << " corresponding to " << sqrt(2.)*TMath::ErfInverse(1-2*upperLimit) << " sigmas." << endl;
	} else {
		cout << (gBdToMuMu ? "Bd -> mumu" : "Bs -> mumu") << " upper limit for config file '" << configfile_path << "' using algorithm " << algo_name(gAlgorithm) << ": " << upperLimit*(gBdToMuMu ? bdtomumu() : bstomumu()) << "\t(" << upperLimit << ") @ " << (int)(gCLLevel*100.) << " % CL" << endl;
		cout << (gBdToMuMu ? "Bd -> mumu" : "Bs -> mumu") << " lower limit for config file '" << configfile_path << "' using algorithm " << algo_name(gAlgorithm) << ": " << lowerLimit*(gBdToMuMu ? bdtomumu() : bstomumu()) << "\t(" << lowerLimit << ") @ " << (int)(gCLLevel*100.) << " % CL" << endl;
	}
	
	if (workspace_path) wspace->writeToFile(workspace_path,kTRUE);
	
	if (output_path) {
		FILE *outFile = fopen(output_path, "w");
		if (gAlgorithm == kAlgo_CLb || gAlgorithm == kAlgo_CLb_Hybrid) {
			fprintf(outFile, "p value of background model: %f\n corresponding to %f sigmas\n", upperLimit, sqrt(2.)*TMath::ErfInverse(1 - 2*upperLimit));
		} else {
			for (ch = channels.begin(); ch != channels.end(); ++ch) {
				fprintf(outFile, "NbObs_%d=%d, NsObs_%d=%d, NdObs_%d=%d\n", *ch, (int)(bsmm[make_pair(kObsBkg_bmm, *ch)].getVal()), *ch, (int)(bsmm[make_pair(kObsB_bmm, *ch)].getVal()),*ch, (int)(bdmm[make_pair(kObsB_bmm, *ch)].getVal()));
			}
			fprintf(outFile, "%s upper limit with algorithm %s: %.5e (%.3f) @ %d %% CL\n", (gBdToMuMu ? "Bd->mumu" : "Bs->mumu"), algo_name(gAlgorithm), upperLimit*(gBdToMuMu ? bdtomumu() : bstomumu()), upperLimit, (int)(100.*gCLLevel));
			fprintf(outFile, "%s lower limit with algorithm %s: %.5e (%.3f) @ %d %% CL\n", (gBdToMuMu ? "Bd->mumu" : "Bs->mumu"), algo_name(gAlgorithm), lowerLimit*(gBdToMuMu ? bdtomumu() : bstomumu()), lowerLimit, (int)(100.*gCLLevel));
		}
		fclose(outFile);
	}
	delete wspace;
	
	return 0;
} // main()
