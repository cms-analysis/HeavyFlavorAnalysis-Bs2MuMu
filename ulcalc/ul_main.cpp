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
#include <RooWorkspace.h>
#include <RooStats/ConfInterval.h>

#include "external_constants.h"
#include "ul_estimate.h"

using namespace std;

static const char *workspace_path = NULL;
static const char *configfile_path = NULL;
static double gCLLevel = 0.9;

#ifdef __linux__
/* Linux does not have that nice function fgetln
 * Use a very simple approach in that case */
static char *fgetln(FILE *f, size_t *len)
{
  static char buffer[2048];
  char *result;

  result = fgets(buffer,sizeof(buffer),f);
  if (result) *len = strlen(buffer);
  else *len = 0;

  return result;
} // fgetln()
#endif

static bool validate_bmm(map<bmm_param,double> *bmm, int ch)
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
	
	// FIXME: complete this list...
	
bail:
	return valid;
} // validate_bmm()

static bool process_line(string name, vector<double> *values, map<bmm_param,double> *bsmm, map<bmm_param,double> *bdmm)
{
	bmm_param p;
	bool bsparam;
	bool success = false;
	
	if (values->size() <= 1)
		goto bail;
	
	p.first = find_bmm_param_by_name(name, &bsparam);
	if (p.first == kUnknownParam) goto bail;
	
	p.second = (int)(*values)[0]; // channel index
	
	if (bsparam)	(*bsmm)[p] = (*values)[1];
	else			(*bdmm)[p] = (*values)[1];
	
	success = true;
bail:
	return success;
} // process_line()

static void parse_input(const char *path, map<bmm_param,double> *bsmm, map<bmm_param,double> *bdmm)
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
	cerr << "ulcalc [-l cl] [-w workspace_outfile.root] <configfile>" << endl;
} // usage()

static bool parse_arguments(const char **first, const char **last)
{
	bool ok = true;
	const char *arg;
	
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
	if (workspace_path)
		cout << "Output Workspace into: " << workspace_path << endl;
	
bail:
	return ok;
} // parse_arguments()

static void dump_params(pair<bmm_param,double> p)
{
	cout << '\t' << find_bmm_name(p.first.first) << '[' << p.first.second << "]: " << p.second << endl;
} // dump_params()

int main(int argc, const char *argv [])
{
	RooStats::ConfInterval *inter;
	RooWorkspace *wspace = NULL;
	RooDataSet *data;
	map<bmm_param,double> bsmm,bdmm;
	set<int> channels;
	bool allData;
	double upper_limit;
	
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
		
		// check given parameters
		if(!validate_bmm(&bsmm,*c) || !validate_bmm(&bdmm,*c)) {
			cerr << "Insufficient parameters in the parameter file!!" << endl;
			abort();
		}
	}
	
	// show the input parsed
	cout << "---------------------------" << endl;
	cout << "Bs -> mumu parameters read:" << endl;
	for_each(bsmm.begin(), bsmm.end(), dump_params);
	cout << "---------------------------" << endl;
	cout << "Bd -> mumu parameters read:" << endl;
	for_each(bdmm.begin(), bdmm.end(), dump_params);
	cout << "---------------------------" << endl;
	
	compute_vars(&bsmm,true);
	compute_vars(&bdmm,false);
	wspace = build_model_nchannel(&bsmm,&bdmm,false);
	
	// set the measured candidates...
	allData = true;
	for (set<int>::const_iterator it = channels.begin(); allData && it != channels.end(); ++it)
		allData = (bsmm.count(make_pair(kObsB_bmm, *it)) > 0) && (bdmm.count(make_pair(kObsB_bmm, *it)) > 0);
	
	if (allData) {
		cout << "Values found for all Observables. Computing this upper limit..." << endl;
		data = build_data_nchannel(wspace,&bsmm,&bdmm);
		data->Print("v");
		delete (inter = est_ul_bc(wspace, data, channels, gCLLevel, &upper_limit, NULL));
		delete data;
		
		cout << "Found Bayesian upper limit: " << upper_limit*BFSM << endl;
	} else {
		cout << "Computing exptected upper limit" << endl;
		cerr << "Not yet implemented... Aborting" << endl;
		abort();
	}
	
	
	if (workspace_path) wspace->writeToFile(workspace_path,kTRUE);
	delete wspace;
	
	return 0;
} // main()
