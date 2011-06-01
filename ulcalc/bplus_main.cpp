/*
 *  bplus_main.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christph Nägeli. All rights reserved.
 *
 */

// My headers
#include "external_constants.h"
#include "bplus_estimator.h"

// Standard headers
#include <cmath>
#include <iostream>
#include <map>
#include <float.h>

// ROOT headers
#include <TFile.h>

using namespace std;

/* Space for global configuration variables */
static map<double,TCut> g_eta_cuts;
static const char *g_mc_file;
static const char *g_data_file;
static const char *g_outputfile;
static bool g_append = false;

static void usage()
{
	cout << "bplus [-a] -c <cutsfile> -m <mcfile> -d <datafile> <outputfile>" << endl;
	abort();
} // usage()


static void parse_arguments(const char **first, const char **last)
{
	const char *arg;
	string s;
	string::iterator colon;
	
	while (first != last) {
		
		arg = *first++;
		if (arg[0] == '-') {
			
			if (strcmp(arg, "-c") == 0) {
			
				g_eta_cuts.clear(); // renew the cuts
				if (first == last) usage();
				parse_cuts(*first++, &g_eta_cuts);
			} else if (strcmp(arg, "-m") == 0) {
				if (first == last) usage();
				g_mc_file = *first++;
			} else if (strcmp(arg, "-d") == 0) {
				if (first == last) usage();
				g_data_file = *first++;
			} else if (strcmp(arg, "-a") == 0) {
				g_append = true;
			}
		} else {
			g_outputfile = arg;
		}
	}
	
	if (!(g_mc_file && g_data_file && g_outputfile) || g_eta_cuts.size() == 0)
		usage();
	
	// print the sets
	cout << "MC File: " << g_mc_file << endl;
	cout << "Datafile: " << g_data_file << endl;
	cout << "Outputfile: " << g_outputfile << endl;
	cout << "Channels: " << endl;
	for (map<double,TCut>::const_iterator it = g_eta_cuts.begin(); it != g_eta_cuts.end(); ++it)
		cout << '\t' << it->first << ": " << it->second.GetTitle() << endl;
} // parse_arguments()

int main(int argc, const char *argv [])
{
	map<double,TCut>::const_iterator it;
	map<bmm_param,measurement_t> bplus;
	TFile *dataFile = NULL, *mcFile = NULL;
	TTree *dataTree = NULL, *mcTree = NULL;
	FILE *outputFile = NULL;
	uint32_t channelIx = 0;
	measurement_t m;
	double last_eta = -1;
	
	parse_arguments(&argv[1], &argv[argc]);
	
	// prepare the input
	dataFile = new TFile(g_data_file);
	mcFile = new TFile(g_mc_file);
	
	dataTree = (TTree*)dataFile->Get("T");
	mcTree = (TTree*)mcFile->Get("T");
	
	// prepare output file
	outputFile = fopen(g_outputfile, (g_append ? "a": "w") );
	
	for (it = g_eta_cuts.begin(); it != g_eta_cuts.end(); ++it) {
		
		bplus.clear();
		estimate_bplus(&bplus, dataTree, mcTree, last_eta, it->first, channelIx, it->second);
		last_eta = it->first;
		
		// add the numbers to the file
		fprintf(outputFile, "#############################################\n");
		fprintf(outputFile, "# B+ -> J/psi K+ (channel = %u, eta < %.2f) #\n", channelIx, it->first);
		fprintf(outputFile, "#############################################\n");
		
		if (bplus.count(make_pair(kTot_bplus, channelIx)) > 0) {
			m = bplus[make_pair(kTot_bplus, channelIx)];
			fprintf(outputFile, "TOT_BPLUS\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
		} else {
			// acceptance
			m = bplus[make_pair(kAcc_bplus, channelIx)];
			fprintf(outputFile, "ACC_BPLUS\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
			// muon efficiency
			m = bplus[make_pair(kEff_mu_bplus, channelIx)];
			fprintf(outputFile, "EFF_MU_BPLUS\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
			// trigger efficiency
			m = bplus[make_pair(kEff_trig_bplus, channelIx)];
			fprintf(outputFile, "EFF_TRIG_BPLUS\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
			// cand efficiency
			m = bplus[make_pair(kEff_cand_bplus, channelIx)];
			fprintf(outputFile, "EFF_CAND_BPLUS\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
			// analysis efficiency
			m = bplus[make_pair(kEff_ana_bplus, channelIx)];
			fprintf(outputFile, "EFF_ANA_BPLUS\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
			// observed bpluses
			m = bplus[make_pair(kObs_bplus, channelIx)];
			fprintf(outputFile, "OBS_BPLUS\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
		}
		
		channelIx++;
	}
	
	if (outputFile)
		fclose(outputFile);
	
	return 0;
} // main()
