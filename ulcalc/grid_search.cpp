/*
 *  grid_search.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 23.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

// my headers
#include "bmm_estimator.h"
#include "ul_estimate.h"

// standard headers
#include <iostream>
#include <algorithm>
#include <cmath>

// ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TMath.h>

// RooFit headers
#include <RooRealVar.h>

// FIXME: Change this!!!
#define TOTAL_BPLUS 800000

#ifdef __linux__
#define INT64_MAX 0x7fffffffffffffffLL
#endif

using namespace std;

static const char *grid_tree_filename = NULL;
static const char *bsmumu_mc_filename = NULL;
static const char *bmumu_data_filename = NULL;
static const char *output_filename = NULL;

static double gCLLevel = 0.9;
static int64_t gStart = 0;
static int64_t gEnd = INT64_MAX;

static map<pair<int,int>,double> ulValues; // pair = (nbr_bkg,nbr_sig)
static double Pss0;
static double eff_tot0;

static void usage()
{
	cout << "usage: grid_search [-r start,end] [-l <conf_level>] -s <bsmcfile> -f <datafile> -o <output_tree> <gridtree>" << endl;
	abort();
} // main()

static void parse_arguments(const char **first, const char **last)
{
	while (first != last) {
		
		if (strcmp(*first, "-s") == 0) {
			if (++first == last)
				usage();
			bsmumu_mc_filename = *first;
		} else if (strcmp(*first, "-f") == 0) {
			if (++first == last)
				usage();
			bmumu_data_filename = *first;
		} else if (strcmp(*first, "-l") == 0) {
			if (++first == last)
				usage();
			gCLLevel = atof(*first);
		} else if (strcmp(*first, "-r") == 0) {
			if (++first == last)
				usage();
#ifdef __linux__
			sscanf(*first,"%ld,%ld",&gStart,&gEnd);
#else
			sscanf(*first,"%lld,%lld",&gStart,&gEnd);
#endif
		} else if (strcmp(*first, "-o") == 0) {
			if (++first == last)
				usage();
			output_filename = *first;
		}
		
		else {
			grid_tree_filename = *first;
		}
		
		first++;
	}
	
	// check if everything is configured.
	if (!grid_tree_filename || !bsmumu_mc_filename || !bmumu_data_filename || !output_filename)
		usage();
	
	// show the arguments parsed...
	cout << "--------------------------------------" << endl;
	cout << "Grid Search using the following configuration" << endl;
	cout << "Grid Tree: " << grid_tree_filename << endl;
	cout << "Data File: " << bmumu_data_filename << endl;;
	cout << "Bs MC File: " << bsmumu_mc_filename << endl;
	cout << "Conf. Level: " << gCLLevel << endl;
	cout << "Output File: " << output_filename << endl;
	cout << "--------------------------------------" << endl;
} // parse_arguments()

static double exp_bs(map<bmm_param,measurement_t> *bsmm, double bkg, int channel)
{
	return compute_tau(bsmm, bsmm, channel, true)*bkg + (*bsmm)[make_pair(kProb_swind_bmm, channel)].getVal()*(*bsmm)[make_pair(kExp_bmm, channel)].getVal();
} // exp_bs()

static RooDataSet *build_observation(RooWorkspace *wspace, double nbr_bkg, double nbr_bs)
{
	RooRealVar *ns = wspace->var("NsObs");
	RooRealVar *nb = wspace->var("NbObs");
	RooArgSet vars(*ns,*nb);
	RooDataSet *data = new RooDataSet("data","",vars);
	
	ns->setVal((double)nbr_bs);
	nb->setVal((double)nbr_bkg);
	
	data->add(vars);
	
	return data;
} // build_observation()

static double process_cut(RooWorkspace *wspace, TTree *dataTree, TTree *bsmcTree, TCut anaCut, float_t *nbr_bkg = NULL, float_t *eff_bs = NULL)
{
	map<bmm_param,measurement_t> bsmm;
	RooDataSet *data;
	RooStats::ConfInterval *confI;
	int64_t n_bkg,n_s,j;
	double mean_bs;
	double weight_ns;
	double std_ul;
	
	double avg_weight = 0.0;
	double avg_ul = 0.0;
	
	pair<int,int> key;
	
	// get the efficiencies for this cut (fast version w/o bplus estimate...)
	estimate_bmm(&bsmm, dataTree, bsmcTree, -1.0, 100.0, 0, anaCut, make_pair(5.15,5.32), make_pair(5.32,5.45), true);
	bsmm[make_pair(kLow_signal_window_bmm, 0)] = measurement_t(5.32,0);
	bsmm[make_pair(kHigh_signal_window_bmm, 0)] = measurement_t(5.45,0);
	bsmm[make_pair(kTot_bplus, 0)] = measurement_t(TOTAL_BPLUS);
	compute_vars(&bsmm, true);
	
	// try in signal window
	n_bkg = (int64_t)(bsmm[make_pair(kObsBkg_bmm, 0)].getVal());
	mean_bs = exp_bs(&bsmm, (double)n_bkg, 0);
	j = 0;
	while (abs(j) < 10 && avg_weight < 0.95) {
		
       	        n_s = (int64_t)round(mean_bs) + j;
		if (n_s >= 0) {
			weight_ns = TMath::PoissonI(n_s, mean_bs);
			
			key = make_pair((int)n_bkg, (int)n_s);
			cout << "Computing bkg = " << n_bkg << ", n_s = " << n_s << endl;
			if (ulValues.count(key) == 0) {
				
				data = build_observation(wspace, n_bkg, n_s);
				confI = est_ul_bc_light(wspace, data, gCLLevel, 0, &std_ul, NULL);
				
				ulValues[key] = std_ul;
				
				delete confI;
				delete data;
			}
			std_ul = ulValues[key];
			
			avg_weight += weight_ns;
			avg_ul += weight_ns * std_ul*(Pss0*eff_tot0) / ( (bsmm[make_pair(kProb_swind_bmm, 0)].getVal()) * compute_efftot_bmm(&bsmm, 0).getVal() );
		}
		
		if (j == 0) j++;
		else if (j > 0) j = -j;
		else j = -j + 1;
	}
	
	// renormalize...
	avg_ul /= avg_weight;
	
	
	// OUTPUT upper limit
	cout << "------------------------------------------------" << endl;
	cout << "Expected Upper limit " << avg_ul << " corresponding to " << avg_ul*bstomumu() << " found for cuts" << endl;
	cout << anaCut.GetTitle() << endl;
	cout << "------------------------------------------------" << endl;
	
	if (nbr_bkg) *nbr_bkg = (float_t)n_bkg;
	if (eff_bs) *eff_bs = (float_t)(compute_efftot_bmm(&bsmm, 0).getVal());
	
	return avg_ul;
} // process_cut()


int main(int argc, const char *argv[])
{
	TFile *outFile = NULL;
	TFile *file = NULL;
	TFile *dataFile = NULL;
	TFile *bsmcFile = NULL;
	TTree *gridTree = NULL;
	TTree *outTree = NULL;
	TTree *bsmcTree;
	TTree *dataTree;
	TObjArray *branches;
	TBranch *br;
	int64_t j;
	size_t k;
	vector<pair<string,float_t>*> theCuts;
	pair<string,float_t> *p;
	TCut anaCut("alpha < 0.05  && chi2/Ndof < 2.4 && d3/d3e > 9 &&  iso10_pt9_pv > 0.55 && pt_mu2 > 2.8 && pt > 4.4");
	float_t upper_limit,nbr_bkg;
	float_t eff_bs;
	map<bmm_param,measurement_t> bsmm;
	RooWorkspace *wspace = NULL;
	set<int> channels;
	
	// consider only channel 0
	channels.insert(0);
	
	// parse the arguments
	parse_arguments(&argv[1],&argv[argc]);
	
	cout << "Base cut:" << endl;
	cout << '\t' << anaCut.GetTitle() << endl;
	cout << "------------------------------------------------" << endl;
	
	// prepare the grid tree
	file = new TFile(grid_tree_filename);
	gridTree = (TTree*)file->Get("GT");
	outFile = new TFile(output_filename,"recreate");
	outTree = new TTree("OT","Grid search UL tree");
	outTree->Branch("upper_limit", &upper_limit, "upper_limit/F");
	outTree->Branch("nbr_bkg",&nbr_bkg,"nbr_bkg/F");
	outTree->Branch("eff_bs", &eff_bs, "eff_bs/F");
	
	// setup the addresses
	branches = gridTree->GetListOfBranches();
	for (j = 0; j < branches->GetEntries(); j++) {
		string name;
		string by("by");
		string zzz("_zzz");
		string::iterator pos;
		
		br = dynamic_cast<TBranch*>((*branches)[j]);
		
		// replace 'by' with '/'
		name = string(br->GetName());
		while ( (pos = search(name.begin(), name.end(), by.begin(), by.end())) != name.end() ) {
			// ersetze dieses 'by' durch '/'
			*pos++ = '/';
			pos = copy(pos + by.size()-1, name.end(), pos);
			name.erase(pos, name.end());
		}
		
		// replace '_zzz' with '[0]'
		if ( (pos = search(name.begin(), name.end(), zzz.begin(), zzz.end())) != name.end()) {
			*pos++ = '[';
			*pos++ = '0';
			*pos++ = ']';
			name.erase(pos,name.end());
		}
		
		p = new pair<string,float_t>(name,0.0);
		theCuts.push_back(p);
		br->SetAddress(&p->second);
		
		outTree->Branch(br->GetName(), &p->second, Form("%s/F",br->GetName()));
	}
	
	// loop over the tree
	if (gridTree->GetEntries() < gEnd)
		gEnd = gridTree->GetEntries();
	
	dataFile = new TFile(bmumu_data_filename);
	bsmcFile = new TFile(bsmumu_mc_filename);
	dataTree = (TTree*)dataFile->Get("T");
	bsmcTree = (TTree*)bsmcFile->Get("T");
	bsmm[make_pair(kTot_bplus, 0)] = measurement_t(TOTAL_BPLUS);
	bsmm[make_pair(kLow_signal_window_bmm, 0)] = measurement_t(5.32,0);
	bsmm[make_pair(kHigh_signal_window_bmm, 0)] = measurement_t(5.45,0);
	estimate_bmm(&bsmm, dataTree, bsmcTree, -1.0, 100.0, 0, anaCut, make_pair(5.15, 5.32), make_pair(5.32, 5.45), true);
	compute_vars(&bsmm, true);
	wspace = build_model_light(&bsmm, 1);
	Pss0 = bsmm[make_pair(kProb_swind_bmm, 0)].getVal();
	eff_tot0 = compute_efftot_bmm(&bsmm, 0).getVal();	
	cout << "--------------------------------" << endl;
	cout << "Processing tree range [" << gStart << ", " << gEnd << ")." << endl;
	for (j = gStart; j < gEnd; j++) {
		
		// load the cuts...
		gridTree->GetEntry(j);
		
		// build the cut...
		anaCut = TCut("");
		for (k = 0; k < theCuts.size(); k++) {
			p = theCuts[k];
			anaCut = anaCut && TCut(Form("%s %c %f", p->first.c_str(), (p->second < 0) ? '<' : '>', fabs(p->second) ));
		}
		
		// compute the upper limit parameters...
		upper_limit = (float_t)process_cut(wspace, dataTree, bsmcTree, anaCut, &nbr_bkg, &eff_bs);
		outTree->Fill();
	}
	
	outFile->Write();
	outFile->Close();
	while (theCuts.size() > 0) {
		p = theCuts.back();
		theCuts.pop_back();
		delete p;
	}
	// clean up
	delete wspace;
	delete file;
	delete outFile;
	delete dataFile;
	delete bsmcFile;
	
	return 0;
} // main()
