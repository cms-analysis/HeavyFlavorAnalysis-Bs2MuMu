/*
 *  ulc_ext.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 19.03.12.
 *
 */

#include "external_constants.h"
#include "bplus_estimator.h"
#include "bmm_estimator.h"

// Standard headers
#include <iostream>
#include <map>
#include <stdlib.h>

// ROOT headers
#include <TCut.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

using namespace std;

static map<double,TCut> g_eta_cuts;
static const char *g_acc_file = NULL;
static const char *g_sig_mc_file = NULL;
static const char *g_rare_file = NULL;
static const char *g_data_file = NULL;
static const char *g_outputfile = NULL;
static map<systematics_t,double> g_systematics_table;

static void insert_systematics(map<systematics_t,double> *g_systematics_table)
{
	g_systematics_table->insert(make_pair(g_sys_acc_efftrack, 0.04));
	g_systematics_table->insert(make_pair(g_sys_acc_ppro_barrel, 0.035));
	g_systematics_table->insert(make_pair(g_sys_acc_ppro_endcap, 0.05));
	g_systematics_table->insert(make_pair(g_sys_effana, 0.04));
	g_systematics_table->insert(make_pair(g_sys_massscale, 0.03));
	g_systematics_table->insert(make_pair(g_sys_effmu_barrel, 0.04));
	g_systematics_table->insert(make_pair(g_sys_effmu_endcap, 0.08));
	g_systematics_table->insert(make_pair(g_sys_efftrig_barrel, 0.03));
	g_systematics_table->insert(make_pair(g_sys_efftrig_endcap, 0.06));
	g_systematics_table->insert(make_pair(g_sys_normfit, 0.05));
	g_systematics_table->insert(make_pair(g_sys_shapecombbkg, 0.04));
} // insert_systematics()

static void usage()
{
	cout << "ulc_ext [--enable-systematics] -c <cutsfile> -a <accfile> -r <rare_mcfile> -m <sig_mcfile> -d <datafile> <outputfile>" << endl;
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
			
			// command line option
			if			(strcmp(arg, "--enable-systematics") == 0) {
				insert_systematics(&g_systematics_table);
			} else if	(strcmp(arg, "-c") == 0) {
				g_eta_cuts.clear(); // renew cuts with the file
				if (first == last) usage();
				parse_cuts(*first++, &g_eta_cuts);
			} else if	(strcmp(arg, "-a") == 0) {
				if (first == last) usage();
				g_acc_file = *first++;
			} else if	(strcmp(arg, "-m") == 0) {
				if (first == last) usage();
				g_sig_mc_file = *first++;
			} else if	(strcmp(arg, "-r") == 0) {
				if (first == last) usage();
				g_rare_file = *first++;
			} else if	(strcmp(arg, "-d") == 0) {
				if (first == last) usage();
				g_data_file = *first++;
			}
		} else
			g_outputfile = arg;
	}
	
	// check consistency
	if( !(g_acc_file && g_sig_mc_file && g_rare_file && g_data_file && g_outputfile) || g_eta_cuts.size() == 0 )
		usage();
	
	// print the configuration
	cout << "Acc File: " << g_acc_file << endl;
	cout << "Signal MC File: " << g_sig_mc_file << endl;
	cout << "Rare Bkg File: " << g_rare_file << endl;
	cout << "Data File: " << g_data_file << endl;
	cout << "Output ulc: " << g_outputfile << endl;
	cout << "Channels: " << endl;
	for (map<double,TCut>::const_iterator it = g_eta_cuts.begin(); it != g_eta_cuts.end(); ++it)
		cout << '\t' << it->first << ": " << it->second.GetTitle() << endl;
	cout << "Systematics " << (g_systematics_table.size() > 0 ? "enabled" : "disabled") << endl;	
} // parse_arguments()

int main(int argc, const char *argv [])
{
	map<systematics_t,double>::const_iterator syst_it;
	map<double,TCut>::const_iterator it;
	map<bmm_param,measurement_t> bplus;
	map<bmm_param,measurement_t> bsmm;
	map<bmm_param,measurement_t> bdmm;
	FILE *outputFile = NULL;
	TFile *dataFile = NULL, *sigMcFile = NULL, *rareMcFile = NULL, *accFile = NULL;
	TTree *dataTree, *sigTree, *rareTree, *accTree;
	uint32_t channelIx = 0;
	double last_eta = -1;
	measurement_t m;
	TFile rootfile("ulc.root","recreate");
	auto_ptr<map<int,triplet<measurement_t> > > rare_effs;
	
	// run in batch mode
	gROOT->SetBatch(kTRUE);
	
	// parse the command line options
	parse_arguments(&argv[1], &argv[argc]);
	
	// prepare the input
	dataFile = new TFile(g_data_file, "update");
	sigMcFile = new TFile(g_sig_mc_file, "update");
	rareMcFile = new TFile(g_rare_file, "update");
	accFile = new TFile(g_acc_file, "update");
	
	dataTree	= (TTree*)dataFile->Get("T");
	sigTree		= (TTree*)sigMcFile->Get("T");
	rareTree	= (TTree*)rareMcFile->Get("T");
	accTree		= (TTree*)accFile->Get("T");
	
	// prepare the output file
	outputFile = fopen(g_outputfile,"w");
	
	// process all channels
	for (it = g_eta_cuts.begin(); it != g_eta_cuts.end(); ++it) {
		
		// get the efficiencies first...
		bsmm.clear();
		bdmm.clear();
		
		rootfile.cd(); // save histograms here
		rare_effs.reset(estimate_bmm_eff(&bsmm, accTree, sigTree, rareTree, last_eta, it->first, channelIx, it->second, bd_range, bs_range, true, &g_systematics_table));
		estimate_bmm_eff(&bdmm, accTree, sigTree, rareTree, last_eta, it->first, channelIx, it->second, bd_range, bs_range, false, &g_systematics_table);
		
		// measure bplus related stuff
		bplus.clear();
		estimate_bplus(&bplus, dataTree, sigTree, accTree, last_eta, it->first, channelIx, it->second, &g_systematics_table);
		
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
			
			// for convenience, compute eff_tot and tot_bplus in comment line...
			m = compute_efftot_bplus(&bplus, channelIx);
			fprintf(outputFile, "# EFF_TOT_BPLUS\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
			m = bplus[make_pair(kObs_bplus, channelIx)] / m;
			fprintf(outputFile, "# TOT_BPLUS\t%u\t%f\t%f\n", channelIx, m.getVal(), m.getErr());
			
			// save the total bpluse -> j/psi k+
			bplus[make_pair(kTot_bplus, channelIx)] = m;
		}
		
		m = bplus[make_pair(kTot_bplus, channelIx)] / (bf_Bu2JpsiKp() * bf_PsiToMuMu());
		estimate_bmm_obs(&bsmm, dataTree, last_eta, it->first, channelIx, it->second, bd_range, bs_range, true, &(*rare_effs), m);
		estimate_bmm_obs(&bdmm, dataTree, last_eta, it->first, channelIx, it->second, bd_range, bs_range, false, &(*rare_effs), m);
		
		// compute histogram specific stuff
		bsmm[make_pair(kLow_signal_window_bmm, channelIx)] = measurement_t(bs_range.first);
		bsmm[make_pair(kHigh_signal_window_bmm, channelIx)] = measurement_t(bs_range.second);
		
		bdmm[make_pair(kLow_signal_window_bmm, channelIx)] = measurement_t(bd_range.first);
		bdmm[make_pair(kHigh_signal_window_bmm, channelIx)] = measurement_t(bd_range.second);
		
		m = measurement_t(compute_tau(&bsmm, &bdmm, channelIx, true),0);
		if( (syst_it = g_systematics_table.find(g_sys_shapecombbkg)) != g_systematics_table.end() )
			m.setErr(m.getVal() * syst_it->second);
		bsmm[make_pair(kTau_bmm, channelIx)] = m;
		
		m = measurement_t(compute_tau(&bsmm, &bdmm, channelIx, false),0);
		if ( (syst_it = g_systematics_table.find(g_sys_shapecombbkg)) != g_systematics_table.end() )
			m.setErr(m.getVal() * syst_it->second);
		bdmm[make_pair(kTau_bmm, channelIx)] = m;
		
		
		// add the numbers to the file
		fprintf(outputFile, "#######################################\n");
		fprintf(outputFile, "# B -> mumu (channel = %u, eta < %.2f) #\n", channelIx, it->first);
		fprintf(outputFile, "#######################################\n");
		fprintf(outputFile, "LOW_BD\t%u\t%f\n", channelIx, bd_range.first);
		fprintf(outputFile, "HIGH_BD\t%u\t%f\n", channelIx, bd_range.second);
		fprintf(outputFile, "LOW_BS\t%u\t%f\n", channelIx, bs_range.first);
		fprintf(outputFile, "HIGH_BS\t%u\t%f\n", channelIx, bs_range.second);		
		
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
		
		// Peaking Background
		m = bsmm[make_pair(kPeakBkgOff_bmm, channelIx)];
		fprintf(outputFile, "PEAK_BKG_OFF\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
		m = bsmm[make_pair(kPeakBkgOn_bmm, channelIx)];
		fprintf(outputFile, "PEAK_BKG_BS\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
		m = bdmm[make_pair(kPeakBkgOn_bmm, channelIx)];
		fprintf(outputFile, "PEAK_BKG_BD\t%u\t%f\t%f\n",channelIx,m.getVal(),m.getErr());
		
		// Print tau
		fprintf(outputFile, "TAU_BS\t%u\t%f\t%f\n", channelIx, bsmm[make_pair(kTau_bmm, channelIx)].getVal(),bsmm[make_pair(kTau_bmm, channelIx)].getErr());
		fprintf(outputFile, "TAU_BD\t%u\t%f\t%f\n", channelIx, bdmm[make_pair(kTau_bmm, channelIx)].getVal(),bdmm[make_pair(kTau_bmm, channelIx)].getErr());
		
		// Print extra stuff as comment
		// needed for MVA training
		{
			measurement_t tot_bplus; // actuall B+ -> J/psi K+
			measurement_t exp_bmm;
			measurement_t gen_bmm;
			
			fprintf(outputFile, "# Scaling factor for Background given by TAU_BS and TAU_BD\n");
			
			if (bplus.count(make_pair(kTot_bplus, channelIx)) > 0)
				tot_bplus = bplus[make_pair(kTot_bplus, channelIx)];
			else {
				tot_bplus = bplus[make_pair(kObs_bplus, channelIx)];
				tot_bplus = tot_bplus / compute_efftot_bplus(&bplus, channelIx);
			}
			
			// scaling bsmm
			exp_bmm = f_ratio() * c_s_theory() * tot_bplus;
			gen_bmm = measurement_t((double)sigTree->Draw("", "candidate == -80"),0);
			fprintf(outputFile, "# SCALE_BS = %f\n", exp_bmm.getVal()/gen_bmm.getVal());
			
			// scaling bdmm
			exp_bmm = c_d_theory() * tot_bplus;
			gen_bmm = measurement_t((double)sigTree->Draw("", "candidate == -90"),0);
			fprintf(outputFile, "# SCALE_BD = %f\n", exp_bmm.getVal()/gen_bmm.getVal());
			
			// dump output, all B+
			exp_bmm = tot_bplus / (bf_Bu2JpsiKp() * bf_PsiToMuMu());
			fprintf(outputFile, "# TOTAL PRODUCED B+ (not B+ -> J/psi K+) = %.2f +/- %.2f\n", exp_bmm.getVal(), exp_bmm.getErr());
		}
		
		last_eta = it->first; // for next round
		channelIx++;
	}
	
	// cleanupt
	delete dataFile;
	delete sigMcFile;
	delete rareMcFile;
	delete accFile;
	if (outputFile) fclose(outputFile);
	
	return 0;
} // main()
