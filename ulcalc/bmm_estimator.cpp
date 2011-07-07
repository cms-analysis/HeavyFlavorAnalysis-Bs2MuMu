/*
 *  bmm_estimator.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#include "bmm_estimator.h"

#include <iostream>

// ROOT headers
#include <TFile.h>

using namespace std;

// Systematics
const static double system_prod = 0.04; // 4 % uncertainty due to unknown production processes
const static double system_ana = 0.113; // 11.3 % systematics on ana efficiency
const static double system_muon = 0.0354; // 5 % on the ratio
const static double system_trig = 0.0141421; // 2 % on the ratio

void estimate_bmm(map<bmm_param,measurement_t> *bmm, TTree *dataTree, TTree *mcTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut, pair<double,double> bd_window, pair<double,double> bs_window, bool is_bstomumu, double eff_filter, bool enable_systematics)
{
	TCut cut;
	TCut histo_cut(Form("%f < mass && mass < %f",low_histo_bound,high_histo_bound));
	TCut bd_mass_cut(Form("%f < mass && mass < %f",bd_window.first,bd_window.second));
	TCut bs_mass_cut(Form("%f < mass && mass < %f",bs_window.first,bs_window.second));
	double eff,tot;
	double nbr_gens,nbr_acc,nbr_trig,nbr_cand,nbr_ana,nbr_mu;
	measurement_t m;
	
	/* Extend the anaCut with the base cut */
	anaCut = anaCut && TCut(bmmBaseCut);
	anaCut = anaCut && TCut(Form("TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f", maxEta,maxEta));
	anaCut = anaCut && TCut(Form("!(TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f)", minEta,minEta));
	
	/* Compute the cross feed */
	cut = anaCut && TCut(bmmGeneratorCuts) && TCut("truth == 1");
	tot = (double)mcTree->Draw("", cut);
	
	/* Pdx */
	eff = (double)mcTree->Draw("", cut && TCut(Form("%f < mass && mass < %f",bd_window.first,bd_window.second)));
	eff /= tot;
	(*bmm)[make_pair(kProb_dwind_bmm,channelIx)] = measurement_t(eff, std_dev_binomail(eff, tot));
	
	/* Psx */
	eff = (double)mcTree->Draw("", cut && TCut(Form("%f < mass && mass < %f",bs_window.first,bs_window.second)));
	eff /= tot;
	(*bmm)[make_pair(kProb_swind_bmm,channelIx)] = measurement_t(eff, std_dev_binomail(eff, tot));
	
	// generate (and filter efficiency)
	cut = TCut(Form("candidate == 0 && (eff_flags & %d)",kGeneratorCand));
	nbr_gens = (double)mcTree->Draw("",cut);
	
	cut = cut && TCut(Form("eff_flags & %d",kAcceptance));
	nbr_acc = (double)mcTree->Draw("",cut);
	
	cut = cut && TCut(Form("eff_flags & %d",kEffMuon));
	nbr_mu = (double)mcTree->Draw("",cut);
	
	cut = cut && TCut("triggers & (1 << 5)");
	nbr_trig = (double)mcTree->Draw("",cut);
	
	cut = TCut(Form("%s && %s && candidate == 301313 && truth == 1",bmmGeneratorCuts,bmmBaseCut));
	nbr_cand = (double)mcTree->Draw("",cut);
	
	cut = cut && anaCut;
	nbr_ana = (double)mcTree->Draw("",cut);
	
	// store the efficiencies in the map
	// the error is of statistical nature only
	
	// acceptance
	eff = nbr_acc / (nbr_gens / eff_filter);
	(*bmm)[make_pair(kAcc_bmm, channelIx)] = measurement_t(eff, std_dev_binomail(eff, nbr_gens));
	if (enable_systematics) {
		m = (*bmm)[make_pair(kAcc_bmm, channelIx)];
		(*bmm)[make_pair(kAcc_bmm, channelIx)] = m + measurement_t(0,m.getVal()*system_prod);
	}
	
	// muon eff
	eff = nbr_mu / nbr_acc;
	(*bmm)[make_pair(kEff_mu_bmm, channelIx)] = measurement_t(eff, std_dev_binomail(eff, nbr_acc));
	if (enable_systematics) {
		m = (*bmm)[make_pair(kEff_mu_bmm, channelIx)];
		(*bmm)[make_pair(kEff_mu_bmm, channelIx)] = m + measurement_t(0,m.getVal()*system_muon);
	}
	
	// trigger efficiency
	eff = nbr_trig / nbr_mu;
	(*bmm)[make_pair(kEff_trig_bmm, channelIx)] = measurement_t(eff, std_dev_binomail(eff, nbr_mu));
	if (enable_systematics) {
		m = (*bmm)[make_pair(kEff_trig_bmm, channelIx)];
		(*bmm)[make_pair(kEff_trig_bmm, channelIx)] = m + measurement_t(0,m.getVal()*system_trig);
	}
	
	// cand efficiency
	eff = nbr_cand / nbr_trig;
	(*bmm)[make_pair(kEff_cand_bmm, channelIx)] = measurement_t(eff, std_dev_binomail(eff, nbr_trig));
	
	// ana efficiency
	eff = nbr_ana / nbr_cand;
	(*bmm)[make_pair(kEff_ana_bmm, channelIx)] = measurement_t(eff, std_dev_binomail(eff, nbr_cand));
	if (enable_systematics) {
		m = (*bmm)[make_pair(kEff_ana_bmm, channelIx)];
		(*bmm)[make_pair(kEff_ana_bmm, channelIx)] = m + measurement_t(0, m.getVal()*system_ana);
	}
	
	/* Observations */
	(*bmm)[make_pair(kObsBkg_bmm, channelIx)] = measurement_t((double)dataTree->Draw("", anaCut && histo_cut && !bd_mass_cut && !bs_mass_cut));
	
	if (!is_bstomumu)
		bs_mass_cut = bd_mass_cut;
	(*bmm)[make_pair(kObsB_bmm, channelIx)] = measurement_t((double)dataTree->Draw("", anaCut && bs_mass_cut));
	
} // estimate_bmm()
