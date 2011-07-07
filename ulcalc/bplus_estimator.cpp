/*
 *  bplus_estimator.cpp
 *  final_calculator
 *
 *  Created by Christoph on 09.03.11.
 *  Copyright 2011 PSI. All rights reserved.
 *
 */

#include "bplus_estimator.h"

// Helper Classes
#include "../rootutils/NCRootUtils.h"

// Standard headers
#include <iostream>
#include <cmath>

// ROOT headers
#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>

using namespace std;

const static double system_prod = 0.04; // 4 % uncertainty due to unknown production processes
const static double system_ana = 0.071; // 7.1 % uncertainty on analysis cuts
const static double system_muon = 0.0354; // 5 % on the ratio
const static double system_trig = 0.0141421; // 2 % on the ratio
const static double system_fit = 0.05; // 5 % uncertainty due to unknow shape of background
const static double system_track = 0.039; // 3.9 % hadron tracking uncertainty

void estimate_bplus(map<bmm_param,measurement_t> *bplus, TTree *dataTree, TTree *mcTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut, double eff_filter, bool enable_systematics)
{
	double nbr_gens;
	double nbr_acc;
	double nbr_trig;
	double nbr_cand;
	double nbr_ana;
	double nbr_mu;
	double lambda;
	fit_t fit;
	TCut cut;
	measurement_t v;
	
	TH1D *m = new TH1D("mbplus","",20,5.18,5.36);
	TCanvas *can = new TCanvas;
	
	static int plot_nbr = 1;
	
	// display statistics
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(111);
	
	anaCut = anaCut && TCut(bmmBaseCut) && TCut("TMath::Abs(eta_kp) < 2.4 && q_mu1*delta_phi < 0");
	anaCut = anaCut && TCut(Form("TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f",maxEta,maxEta));
	anaCut = anaCut && TCut(Form("!(TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f)",minEta,minEta));
	
	// compute the efficiencies
	
	// generated (and filter efficiency)
	cut = TCut(Form("candidate == 0 && (eff_flags & %d)",kGeneratorCand));
	nbr_gens = (double)mcTree->Draw("",cut);
	
	cut = cut && Form("(eff_flags & %d)",kAcceptance);
	nbr_acc = (double)mcTree->Draw("",cut);
	
	cut = cut && Form("(eff_flags & %d)", kEffMuon);
	nbr_mu = (double)mcTree->Draw("",cut);
	
	cut = cut && "(triggers & (1 << 7))";
	nbr_trig = (double)mcTree->Draw("",cut);
	
	cut = TCut(Form("%s && %s && candidate == 300521 && truth == 1 && pt_kp_gen > 0.4 && TMath::Abs(eta_kp_gen) < 2.5 && TMath::Abs(eta_kp) < 2.4 && pt_kp > 0.5",bmmGeneratorCuts,bmmBaseCut));
	nbr_cand = (double)mcTree->Draw("",cut);
	
	cut = cut && anaCut;
	nbr_ana = (double)mcTree->Draw("", cut);
	
	// store the efficiencies in the map
	// The error is of statistical nature only.
	
	// acceptance
	lambda = nbr_acc / (nbr_gens / eff_filter);
	(*bplus)[make_pair(kAcc_bplus, channelIx)] = measurement_t(lambda, std_dev_binomail(lambda, nbr_gens));
	if (enable_systematics) {
		v = (*bplus)[make_pair(kAcc_bplus, channelIx)];
		(*bplus)[make_pair(kAcc_bplus, channelIx)] = v + measurement_t(0, v.getVal()*system_track);
		
		v = (*bplus)[make_pair(kAcc_bplus, channelIx)];
		(*bplus)[make_pair(kAcc_bplus, channelIx)] = v + measurement_t(0, v.getVal()*system_prod);
	}
	
	// muon efficiency
	lambda = nbr_mu / nbr_acc;
	(*bplus)[make_pair(kEff_mu_bplus, channelIx)] = measurement_t(lambda, std_dev_binomail(lambda, nbr_acc));
	if (enable_systematics) {
		v = (*bplus)[make_pair(kEff_mu_bplus, channelIx)];
		(*bplus)[make_pair(kEff_mu_bplus, channelIx)] = v + measurement_t(0,v.getVal()*system_muon);
	}
	
	// trigger efficiency
	lambda = nbr_trig / nbr_mu;
	(*bplus)[make_pair(kEff_trig_bplus, channelIx)] = measurement_t(lambda, std_dev_binomail(lambda, nbr_mu));
	if (enable_systematics) {
		v = (*bplus)[make_pair(kEff_trig_bplus, channelIx)];
		(*bplus)[make_pair(kEff_trig_bplus, channelIx)] = v + measurement_t(0,v.getVal()*system_trig);
	}
	
	// cand efficiency
	lambda = nbr_cand / nbr_trig;
	(*bplus)[make_pair(kEff_cand_bplus, channelIx)] = measurement_t(lambda, std_dev_binomail(lambda, nbr_trig));
	
	// ana efficiency
	lambda = nbr_ana / nbr_cand;
	(*bplus)[make_pair(kEff_ana_bplus, channelIx)] = measurement_t(lambda, std_dev_binomail(lambda, nbr_cand));
	if (enable_systematics) {
		v = (*bplus)[make_pair(kEff_ana_bplus, channelIx)];
		(*bplus)[make_pair(kEff_ana_bplus, channelIx)] = v + measurement_t(0, v.getVal()*system_ana);
	}
	
	// get the number of Bp -> J/psi Kp
	dataTree->Draw("mass_c >> mbplus", anaCut);
	init_gauss_linear(&fit, 5.20, 5.35);
	adjust_parameter_gauss_linear(m, fit.fit_fct);
	m->Fit(fit.fit_fct, "LR");
	
	// observed bplus
	// FIXME: Better solution of uncertainty
	lambda = signal_events_gauss_linear(fit.fit_fct, m->GetBinWidth(1));
	fit.fit_fct->SetParameter(0, fit.fit_fct->GetParError(0));
	(*bplus)[make_pair(kObs_bplus, channelIx)] = measurement_t(lambda,signal_events_gauss_linear(fit.fit_fct, m->GetBinWidth(1)));
	if (enable_systematics) {
		v = (*bplus)[make_pair(kObs_bplus, channelIx)];
		(*bplus)[make_pair(kObs_bplus, channelIx)] = v + measurement_t(0, v.getVal()*system_fit);
	}
	
	// save the histogram
	m->Draw("E1");
	can->SaveAs(Form("obs_bplus_%d.eps",plot_nbr++));
	
	delete can;
	delete fit.fit_fct;
	delete m;
} // estimate_bplus()
