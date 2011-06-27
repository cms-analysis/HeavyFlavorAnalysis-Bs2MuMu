/*
 *  bmm_main.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#include "bplus_estimator.h"
#include "bmm_estimator.h"

// Standard headers
#include <stdlib.h>
#include <iostream>
#include <map>

// ROOT headers
#include <TFile.h>

using namespace std;

/* Space for global configuration variables */
static map<double,TCut> g_eta_cuts;
static const char *g_bs_file;
static const char *g_bd_file;
static const char *g_data_file;
static const char *g_outputfile;

static void usage()
{
	cout << "bmm_est -c <cutsfile> -s <bs_mcfile> -0 <b0_mcfile> -d <datafile> <complement_file>" << endl;
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
			} else if (strcmp(arg, "-s") == 0) {
				if (first == last) usage();
				g_bs_file = *first++;
			} else if (strcmp(arg, "-0") == 0) {
				if (first == last) usage();
				g_bd_file = *first++;
			} else if (strcmp(arg, "-d") == 0) {
				if (first == last) usage();
				g_data_file = *first++;
			} else {
				cerr << "Unknown option: " << arg << endl;
				usage();
			}
		} else {
			g_outputfile = arg;
		}
	}
	
	if (!(g_bd_file && g_bs_file && g_data_file && g_outputfile) || g_eta_cuts.size() == 0)
		usage();
	
	// print the sets
	cout << "Bs MC File: " << g_bs_file << endl;
	cout << "Bd MC File: " << g_bd_file << endl;
	cout << "Datafile: " << g_data_file << endl;
	cout << "Outputfile: " << g_outputfile << endl;
	cout << "Channels: " << endl;
	for (map<double,TCut>::const_iterator it = g_eta_cuts.begin(); it != g_eta_cuts.end(); ++it)
		cout << '\t' << it->first << ": " << it->second.GetTitle() << endl;
} // parse_arguments()

int main(int argc, const char *argv [])
{
	FILE *outputFile = NULL;
	TFile *bsFile = NULL, *bdFile = NULL, *dataFile = NULL;
	TTree *bsTree, *bdTree, *dataTree;
	map<double,TCut>::const_iterator it;
	map<bmm_param,measurement_t> bsmm;
	map<bmm_param,measurement_t> bdmm;
	double last_eta = -1.0;
	pair<double,double> bdWindow(5.2,5.32);
	pair<double,double> bsWindow(5.32,5.45);
	measurement_t m;
	
	static uint32_t channelIx = 0;
	
	parse_arguments(&argv[1], &argv[argc]);
	
	// prepare the input file
	dataFile = new TFile(g_data_file);
	bsFile = new TFile(g_bs_file);
	bdFile = new TFile(g_bd_file);
	
	dataTree = (TTree*)dataFile->Get("T");
	bsTree = (TTree*)bsFile->Get("T");
	bdTree = (TTree*)bdFile->Get("T");
	
	// prepare the output
	outputFile = fopen(g_outputfile, "a");
	
	for (it = g_eta_cuts.begin(); it != g_eta_cuts.end(); ++it) {
		
		bsmm.clear();
		estimate_bmm(&bsmm, dataTree, bsTree, last_eta, it->first, channelIx, it->second, bdWindow, bsWindow, true);
		
		bdmm.clear();
		estimate_bmm(&bdmm, dataTree, bdTree, last_eta, it->first, channelIx, it->second, bdWindow, bsWindow, false);
		
		// add the numbers to the file
		fprintf(outputFile, "#######################################\n");
		fprintf(outputFile, "# B -> mumu (channel = %u, eta < %.2f) #\n", channelIx, it->first);
		fprintf(outputFile, "#######################################\n");
		fprintf(outputFile, "LOW_BD\t%u\t%f\n", channelIx, bdWindow.first);
		fprintf(outputFile, "HIGH_BD\t%u\t%f\n", channelIx, bdWindow.second);
		fprintf(outputFile, "LOW_BS\t%u\t%f\n", channelIx, bsWindow.first);
		fprintf(outputFile, "HIGH_BS\t%u\t%f\n", channelIx, bsWindow.second);
		
		// crossfeed
		m = bsmm[make_pair(kProb_swind_bmm, channelIx)];
		fprintf(outputFile, "PSS\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		m = bdmm[make_pair(kProb_swind_bmm, channelIx)];
		fprintf(outputFile, "PSD\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		m = bsmm[make_pair(kProb_dwind_bmm, channelIx)];
		fprintf(outputFile, "PDS\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		m = bdmm[make_pair(kProb_dwind_bmm, channelIx)];
		fprintf(outputFile, "PDD\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		
		// Bs -> mumu efficiencies
		// acceptance
		m = bsmm[make_pair(kAcc_bmm, channelIx)];
		fprintf(outputFile, "ACC_BSMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// muon efficiency
		m = bsmm[make_pair(kEff_mu_bmm, channelIx)];
		fprintf(outputFile, "EFF_MU_BSMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// trigger efficiency
		m = bsmm[make_pair(kEff_trig_bmm, channelIx)];
		fprintf(outputFile, "EFF_TRIG_BSMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// cand efficiency
		m = bsmm[make_pair(kEff_cand_bmm, channelIx)];
		fprintf(outputFile, "EFF_CAND_BSMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// ana effeciency
		m = bsmm[make_pair(kEff_ana_bmm, channelIx)];
		fprintf(outputFile, "EFF_ANA_BSMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// total efficiency for convenience
		m = compute_efftot_bmm(&bsmm, channelIx);
		fprintf(outputFile, "# EFF_TOT_BSMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		
		// Bd -> mumu efficiencies
		// acceptance
		m = bdmm[make_pair(kAcc_bmm, channelIx)];
		fprintf(outputFile, "ACC_BDMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// muon efficiency
		m = bdmm[make_pair(kEff_mu_bmm, channelIx)];
		fprintf(outputFile, "EFF_MU_BDMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// trigger efficiency
		m = bdmm[make_pair(kEff_trig_bmm, channelIx)];
		fprintf(outputFile, "EFF_TRIG_BDMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// cand efficiency
		m = bdmm[make_pair(kEff_cand_bmm, channelIx)];
		fprintf(outputFile, "EFF_CAND_BDMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// ana efficiency
		m = bdmm[make_pair(kEff_ana_bmm, channelIx)];
		fprintf(outputFile, "EFF_ANA_BDMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		// total efficiency for convenience
		m = compute_efftot_bmm(&bdmm, channelIx);
		fprintf(outputFile, "# EFF_TOT_BDMM\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
		
		// Observations
		fprintf(outputFile, "OBS_BKG\t%u\t%f\n", channelIx, bsmm[make_pair(kObsBkg_bmm, channelIx)].getVal());
		fprintf(outputFile, "OBS_BSMM\t%u\t%f\n", channelIx, bsmm[make_pair(kObsB_bmm, channelIx)].getVal());
		fprintf(outputFile, "OBS_BDMM\t%u\t%f\n", channelIx, bdmm[make_pair(kObsB_bmm, channelIx)].getVal());
		
		last_eta = it->first;
		channelIx++;
	}
	
	delete bdFile;
	delete bsFile;
	delete dataFile;
	
	fclose(outputFile);
	
	return 0;
} // main()
