/*
 *  external_constants.h
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef EXTERNAL_CONSTANTS_H
#define EXTERNAL_CONSTANTS_H

#include <map>
#include <string>
#include <utility>

/* Different filter efficiencies */
const double eff_filter_bplus = 0.24;
const double eff_filter_bmm = 0.63;
const double c_s_theory = 1.80896e-5;
const double c_d_theory = 4.96761e-7;
const double Ntot_Bplus = 435000; // c.f. Lab Journal

/* Eta value for split analysis:
 *	Barrel is < eta_barrel
 *	while we consider endcap to be >= eta_barrel
 */
const double eta_barrel = 1.4;

/* cut to be applied for analysis */
extern const char *bmmGeneratorCuts;
extern const char *bmmBaseCut;
extern const char *bmmAnaBarrelCut;
extern const char *bmmAnaEndcapCut;

/* Standard Model Branching fraction */
const static double BFSM = 3.86e-9;

/* Parameters estimated, saved in map */
enum bmm_param_tag {
	/* Bp -> J/psi Kp */
	kAcc_bplus = 1,
	kEff_mu_bplus,
	kEff_trig_bplus,
	kEff_cand_bplus,
	kEff_ana_bplus,
	kObs_bplus,
	kTot_bplus,
	/* B -> mumu */
	kAcc_bmm,
	kEff_mu_bmm,
	kEff_trig_bmm,
	kEff_cand_bmm,
	kEff_ana_bmm,
	kEff_total_bmm,	// Used for shortcut in Grid search
	kExp_bmm,
	kProb_swind_bmm,
	kProb_dwind_bmm,
	kLow_signal_window_bmm,
	kHigh_signal_window_bmm,
	kObsBkg_bmm,
	kObsB_bmm,
	/* Unknown param */
	kUnknownParam
};
typedef std::pair<bmm_param_tag,int> bmm_param;


/* Model constants:
 *	Total histogram region [4.7,6.2]
 */
const double low_histo_bound = 4.7;
const double high_histo_bound = 6.2;
double compute_tau(std::map<bmm_param,double> *bsmm, std::map<bmm_param,double> *bdmm, int channel, bool tau_s);


/* conversion routines */
std::string find_bmm_name(bmm_param_tag p);
bmm_param_tag find_bmm_param_by_name(std::string name, bool *bsparam);


#endif
