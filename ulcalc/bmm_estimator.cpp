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
#include <TKey.h>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TEfficiency.h>
#include <TEventList.h>

using namespace std;

// Systematics
const static double system_prod = 0.04; // 4 % uncertainty due to unknown production processes
const static double system_ana = 0.113; // 11.3 % systematics on ana efficiency
const static double system_muon = 0.0354; // 5 % on the ratio
const static double system_trig = 0.0141421; // 2 % on the ratio

const static measurement_t misid_kaon_pion(0.001,0.0002);
const static measurement_t misid_proton(0.0005,0.0001);

// enumeration of known decays for truth matching
// corresponds to HFTruthCandidates_cff.py
enum {
	kDecay_BsToMuMu			= 1,
	kDecay_BsToMuMuGa		= 2,
	kDecay_BsToKK			= 3,
	kDecay_BsToKPi			= 4,
	kDecay_BsToPiPi			= 5,
	kDecay_BsToPiMuNu		= 6,
	kDecay_BsToKMuNu		= 7,
	kDecay_BdToMuMu			= 8,
	kDecay_BdToPiPi			= 9,
	kDecay_BdToKPi			= 10,
	kDecay_BdToKK			= 11,
	kDecay_BdToMuMuPi0		= 12,
	kDecay_BdToPiMuNu		= 13,
	kDecay_BuTo3MuNu		= 14,
	kDecay_LambdaBToPPi		= 15,
	kDecay_LambdaBToPK		= 16,
	kDecay_LambdaBToPMuNu	= 17,
	kDecay_Bs2JpsiPhi		= 18,
	kDecay_Bu2JpsiKp		= 19,
	kDecay_Bd2JpsiKstar		= 20,
	kDecay_Bd2JpsiKs		= 21,
	kDecay_PsiToMuMu		= 22,
	kDecay_Psi2SToMuMu		= 23,
	kDecay_Ups1SToMuMu		= 24,
	kDecay_Ups2SToMuMu		= 25,
	kDecay_Ups3SToMuMu		= 26
};

static triplet<measurement_t> eff_raredecay(TTree *mcTree, TCut accCut, TCut cut, TCut bdCut, TCut bsCut, int trueCand, measurement_t mu_misid1, measurement_t mu_misid2, double eff_trigger, measurement_t filterEff, measurement_t bfEff, int channelIx)
{
	Long64_t nbr_reco,nbr_off,nbr_bs,nbr_bd;
	TEfficiency effCalc("eff","",3,0,1); // we only need one bin
	TCut candCut(Form("candidate == %d", 1000000+trueCand));
	triplet<measurement_t> result;
	measurement_t m;
	string histoName(Form("h_pkg_ch%d_%d",channelIx,trueCand));
	TH1D* histo = new TH1D(histoName.c_str(),"",50,4.9,5.9);
	
	// configure the efficiency calculation
	effCalc.SetConfidenceLevel(0.68);
	effCalc.SetStatisticOption(TEfficiency::kFCP);
	
	nbr_reco = mcTree->Draw("", candCut && accCut);
	nbr_off = mcTree->Draw("", candCut && accCut && cut && !bdCut && !bsCut);
	nbr_bd = mcTree->Draw("", candCut && accCut && cut && bdCut);
	nbr_bs = mcTree->Draw("", candCut && accCut && cut && bsCut);
	
	if (nbr_reco > 0) {
		
		effCalc.SetTotalEvents(0, nbr_reco);
		effCalc.SetPassedEvents(0, nbr_off);
		m = measurement_t(effCalc.GetEfficiency(0), (effCalc.GetEfficiencyErrorLow(0) + effCalc.GetEfficiencyErrorUp(0))/2.0);
		m = m * filterEff * bfEff;
		m = m * mu_misid1 * mu_misid2;
		m = m * measurement_t(eff_trigger);
		result.a = m;
		
		effCalc.SetTotalEvents(1, nbr_reco);
		effCalc.SetPassedEvents(1, nbr_bd);
		m = measurement_t(effCalc.GetEfficiency(1), (effCalc.GetEfficiencyErrorLow(1) + effCalc.GetEfficiencyErrorUp(1))/2.0);
		m = m * filterEff * bfEff;
		m = m * mu_misid1 * mu_misid2;
		m = m * measurement_t(eff_trigger);
		result.b = m;
		
		effCalc.SetTotalEvents(2, nbr_reco);
		effCalc.SetPassedEvents(2, nbr_bs);
		m = measurement_t(effCalc.GetEfficiency(2), (effCalc.GetEfficiencyErrorLow(2) + effCalc.GetEfficiencyErrorUp(2))/2.0);
		m = m * filterEff * bfEff;
		m = m * mu_misid1 * mu_misid2;
		m = m * measurement_t(eff_trigger);
		result.c = m;		
	} else
		cerr << "===>ul_ext: Error in computing rare decay with true candidate " << trueCand << endl;
	
	// save histogram to file
	mcTree->Draw(Form("mass>>%s",histoName.c_str()), candCut && accCut && cut);
	histo->Scale( (result.a + result.b + result.c).getVal() / histo->Integral());
	histo->Write(histo->GetName(), TObject::kOverwrite);
	
	return result;
} // eff_raredecay()

static map<int,triplet<measurement_t> > * compute_peaking_eff(TTree *accTree, TTree *mcTree, TCut cut, TCut bdCut, TCut bsCut, double eff_muons, double eff_trigger, int channelIx)
{
	TEfficiency effCalc("eff","",2,0,1); // one bin only
	Long64_t nbr_gens,nbr_acc;
	TCut accCut("TMath::Abs(eta_mu1_gen) < 2.5 && TMath::Abs(eta_mu2_gen) < 2.5 && pt_mu1_gen > 3.5 && pt_mu2_gen > 3.5 && TMath::Abs(eta_mu1) < 2.4 && TMath::Abs(eta_mu2) < 2.4 && (track_qual_mu1 & 4) && (track_qual_mu2 & 4)");
	measurement_t eff_single_mu(TMath::Sqrt(eff_muons),0); // FIXME: set correct error
	measurement_t filterHadronic,filterSemilep;
	map<int,triplet<measurement_t> > *rares = new map<int,triplet<measurement_t > >;
	
	effCalc.SetConfidenceLevel(0.68);
	effCalc.SetStatisticOption(TEfficiency::kFCP);
	
	// compute hadronic filter efficiency based on bs->mumu
	nbr_gens = mcTree->Draw("","candidate == -80");
	nbr_acc = mcTree->Draw("", TCut("candidate == 80") && accCut);
	effCalc.SetTotalEvents(0, nbr_gens);
	effCalc.SetPassedEvents(0, nbr_acc);
	filterHadronic = measurement_t(effCalc.GetEfficiency(0), (effCalc.GetEfficiencyErrorLow(0) + effCalc.GetEfficiencyErrorUp(0))/2.0);
	
	// compute semileptonic filter efficiency based on Bd -> pi mu nu
	nbr_gens = accTree->Draw("", "candidate == -95");
	nbr_acc = accTree->Draw("", TCut("candidate == 1000095") && accCut);
	effCalc.SetTotalEvents(1, nbr_gens);
	effCalc.SetPassedEvents(1, nbr_acc);
	filterSemilep = measurement_t(effCalc.GetEfficiency(1), (effCalc.GetEfficiencyErrorLow(1) + effCalc.GetEfficiencyErrorUp(1))/2.0);

	(*rares)[kDecay_BsToKK]			= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 82, misid_kaon_pion, misid_kaon_pion, eff_trigger, filterHadronic, bf_BsToKK()*f_ratio(), channelIx);
	
	(*rares)[kDecay_BsToKPi]		= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 83, misid_kaon_pion, misid_kaon_pion, eff_trigger, filterHadronic, bf_BsToKPi()*f_ratio(), channelIx);
	
	(*rares)[kDecay_BsToPiPi]		= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 84, misid_kaon_pion, misid_kaon_pion, eff_trigger, filterHadronic, bf_BsToPiPi()*f_ratio(), channelIx);
	
	(*rares)[kDecay_BsToKMuNu]		= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 86, misid_kaon_pion, eff_single_mu, eff_trigger, filterSemilep, bf_BsToKMuNu()*f_ratio(), channelIx);
	
	(*rares)[kDecay_BdToPiPi]		= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 91, misid_kaon_pion, misid_kaon_pion, eff_trigger, filterHadronic, bf_BdToPiPi(), channelIx);
	
	(*rares)[kDecay_BdToKPi]		= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 92, misid_kaon_pion, misid_kaon_pion, eff_trigger, filterHadronic, bf_BdToKPi(), channelIx);
	
	(*rares)[kDecay_BdToKK]			= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 93, misid_kaon_pion, misid_kaon_pion, eff_trigger, filterHadronic, bf_BdToKK(), channelIx);
	
	(*rares)[kDecay_BdToPiMuNu]		= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 95, misid_kaon_pion, eff_single_mu, eff_trigger, filterSemilep, bf_BdToPiMuNu(), channelIx);
	
	(*rares)[kDecay_LambdaBToPPi]	= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 60, misid_proton, misid_kaon_pion, eff_trigger, filterHadronic, bf_LambdaBToPPi()*f_ratio_lb(), channelIx);
	
	(*rares)[kDecay_LambdaBToPK]	= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 61, misid_proton, misid_kaon_pion, eff_trigger, filterHadronic, bf_LambdaBToPK()*f_ratio_lb(), channelIx);
	
	(*rares)[kDecay_LambdaBToPMuNu]	= eff_raredecay(mcTree, accCut, cut, bdCut, bsCut, 62, misid_proton, eff_single_mu, eff_trigger, filterSemilep, bf_LambdaBToPMuNu()*f_ratio_lb(), channelIx);
	
	return rares;
} // compute_peaking_eff()

static triplet<measurement_t> compute_peaking_exp(map<int,triplet<measurement_t> > rare_effs, measurement_t tot_bu)
{
	triplet<measurement_t> expected;
	int j;
	
	for (j = 0; j < 3; j++) {
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_BsToKK].getIx(j)			* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_BsToKPi].getIx(j)		* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_BsToPiPi].getIx(j)		* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_BsToKMuNu].getIx(j)		* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_BdToPiPi].getIx(j)		* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_BdToKPi].getIx(j)		* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_BdToKK].getIx(j)			* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_BdToPiMuNu].getIx(j)		* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_LambdaBToPPi].getIx(j)	* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_LambdaBToPK].getIx(j)	* tot_bu);
		expected.setIx(j, expected.getIx(j) + rare_effs[kDecay_LambdaBToPMuNu].getIx(j)	* tot_bu);
	}
	
	return expected;
} // compute_peaking_exp()

std::map<int,triplet<measurement_t> > *estimate_bmm_eff(std::map<bmm_param,measurement_t> *bmm, TTree *accTree, TTree *sigTree, TTree *rareTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut, std::pair<double,double> bd_window, std::pair<double,double> bs_window, bool is_bstomumu, std::map<systematics_t,double> *systematics_table)
{
	TEfficiency effCalc("eff","",1,0,1); // we only need one bin
	int cand = is_bstomumu ? 80 : 90;
	map<systematics_t,double>::const_iterator syst_it;
	measurement_t m;
	map<int,triplet<measurement_t> > *result = NULL;
	TEventList *elist;
	TEventList treeList;
	string prefix( (is_bstomumu ? "sig" : "dig") );
	
	// Cuts
	TCut histo_cut(Form("%f < mass && mass < %f",low_histo_bound,high_histo_bound));
	TCut bd_mass_cut(Form("%f < mass && mass < %f",bd_window.first,bd_window.second));
	TCut bs_mass_cut(Form("%f < mass && mass < %f",bs_window.first,bs_window.second));
	TCut channelCut = TCut(Form("TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f",maxEta,maxEta)) && TCut(Form("!(TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f)",minEta,minEta));
	TCut acceptanceCutData("TMath::Abs(eta_mu1) < 2.4 && TMath::Abs(eta_mu2) < 2.4 && (track_qual_mu1 & 4) && (track_qual_mu2 & 4) && (q_mu1*q_mu2 < 0)");
	TCut acceptanceCutGen("TMath::Abs(eta_mu1_gen) < 2.5 && TMath::Abs(eta_mu2_gen) < 2.5 && pt_mu1_gen > 1. && pt_mu2_gen > 1.");
	TCut truthCut(Form("true_decay == %d", (is_bstomumu ? 1 : 8))); // see massReader.hh
	TCut muonCut("tight_mu1 && tight_mu2");
	TCut triggerCut("triggered_bs");
	TCut cut;
	
	// candidate counting
	Long64_t nbr_gens,nbr_reco,nbr_ana,nbr_mu,nbr_trig;
	Long64_t nbr_bs,nbr_bd;
	
	// print progress
	cout << "===> Processing efficiencies for channel " << channelIx << Form(": B%c -> mumu", (is_bstomumu ? 's' : 'd')) << endl;
	
	// configure the efficiency calculation
	effCalc.SetConfidenceLevel(0.68);
	effCalc.SetStatisticOption(TEfficiency::kFCP);
	
	/***********************************
	 * Working with the signal MC Tree *
	 ***********************************/
	sigTree->SetEventList(0);
	sigTree->GetDirectory()->cd();
	
	// generated
	cout << "\tCounting generator candidates..." << flush;
	cut = TCut(Form("candidate == %d", -cand));
	sigTree->Draw(">>elist",cut);
	elist = (TEventList*)gDirectory->Get("elist");
	nbr_gens = elist->GetN();
	cout << '\t' << nbr_gens << endl;
	elist->Write(Form("%sGens",prefix.c_str()), TObject::kOverwrite);
	
	// accepted
	cout << "\tCounting reco candidates..." << flush;
	cut = TCut(Form("candidate == %d", cand)) && acceptanceCutData && acceptanceCutGen && channelCut;
	sigTree->Draw(">>elist",cut);
	elist = (TEventList*)gDirectory->Get("elist");
	nbr_reco = elist->GetN();
	cout << '\t' << nbr_reco << endl;
	elist->Write(Form("%sGensAcc_%d",prefix.c_str(),channelIx), TObject::kOverwrite);
	
	// analysis
	cout << "\tCounting ana candiates..." << flush;
	cut = TCut("candidate == 301313") && acceptanceCutData && acceptanceCutGen && channelCut && anaCut && truthCut;
	sigTree->Draw(">>elist",cut);
	elist = (TEventList*)gDirectory->Get("elist");
	nbr_ana = elist->GetN();
	treeList.Clear();
	treeList.Add(elist);
	sigTree->SetEventList(&treeList);
	cout << '\t' << nbr_ana << endl;
	elist->Write(Form("%sAna_%d",prefix.c_str(),channelIx), TObject::kOverwrite);
	
	// muon
	cout << "\tCounting muon candidates..." << flush;
	sigTree->Draw(">>elist",muonCut);
	elist = (TEventList*)gDirectory->Get("elist");
	nbr_mu = elist->GetN();
	treeList.Clear();
	treeList.Add(elist);
	sigTree->SetEventList(&treeList);
	cout << '\t' << nbr_mu << endl;
	elist->Write(Form("%sMuon_%d",prefix.c_str(),channelIx), TObject::kOverwrite);
	
	// trigger
	cout << "\tCounting trigger candidates..." << flush;
	sigTree->Draw(">>elist",triggerCut);
	elist = (TEventList*)gDirectory->Get("elist");
	nbr_trig = elist->GetN();
	treeList.Clear();
	treeList.Add(elist);
	sigTree->SetEventList(&treeList);
	cout << '\t' << nbr_trig << endl;
	elist->Write(Form("%sTrig_%d",prefix.c_str(),channelIx), TObject::kOverwrite);
	
	// cross feed
	cout << "\tComputing cross feed..." << flush;
	nbr_bs = sigTree->Draw("",bs_mass_cut);
	nbr_bd = sigTree->Draw("",bd_mass_cut);
	cout << "\tdone" << endl;
	
	// reset the event list
	sigTree->SetEventList(NULL);
	
	effCalc.SetTotalEvents(0, nbr_trig);
	effCalc.SetPassedEvents(0, nbr_bd);
	(*bmm)[make_pair(kProb_dwind_bmm,channelIx)] = measurement_t(effCalc.GetEfficiency(0),(effCalc.GetEfficiencyErrorLow(0) + effCalc.GetEfficiencyErrorUp(0))/2.);
	
	effCalc.SetTotalEvents(0, nbr_trig);
	effCalc.SetPassedEvents(0, nbr_bs);
	(*bmm)[make_pair(kProb_swind_bmm,channelIx)] = measurement_t(effCalc.GetEfficiency(0),(effCalc.GetEfficiencyErrorLow(0) + effCalc.GetEfficiencyErrorUp(0))/2.);
	
	if( (syst_it = systematics_table->find(g_sys_massscale)) != systematics_table->end() ) {
		// Systematic on Pdx
		m = (*bmm)[make_pair(kProb_dwind_bmm,channelIx)];
		m = measurement_t(0,m.getVal() * syst_it->second);
		(*bmm)[make_pair(kProb_dwind_bmm,channelIx)] = m + (*bmm)[make_pair(kProb_dwind_bmm,channelIx)];
		
		// Systematic on Psx
		m = (*bmm)[make_pair(kProb_swind_bmm,channelIx)];
		m = measurement_t(0,m.getVal() * syst_it->second);
		(*bmm)[make_pair(kProb_swind_bmm,channelIx)] = m + (*bmm)[make_pair(kProb_swind_bmm,channelIx)];
	}
	
	// acceptance
	effCalc.SetTotalEvents(0, nbr_gens);
	effCalc.SetPassedEvents(0, nbr_reco);
	(*bmm)[make_pair(kAcc_bmm, channelIx)] = measurement_t(effCalc.GetEfficiency(0),(effCalc.GetEfficiencyErrorLow(0) + effCalc.GetEfficiencyErrorUp(0))/2.);
	
	// systematics on acceptance ratio
	switch (channelIx) {
		case 0:	// Barrel
			syst_it = systematics_table->find(g_sys_acc_ppro_barrel);
			break;
		case 1: // Endcap
			syst_it = systematics_table->find(g_sys_acc_ppro_endcap);
			break;
		default:
			syst_it = systematics_table->end();
			break;
	}
	if (syst_it != systematics_table->end()) {
		m = (*bmm)[make_pair(kAcc_bmm, channelIx)];
		m = measurement_t(0, m.getVal() * syst_it->second);
		(*bmm)[make_pair(kAcc_bmm, channelIx)] = m + (*bmm)[make_pair(kAcc_bmm, channelIx)];
	}
	
	// historical reason
	(*bmm)[make_pair(kEff_cand_bmm, channelIx)] = measurement_t(1.,0);
	
	// analysis efficiency
	effCalc.SetTotalEvents(0, nbr_reco);
	effCalc.SetPassedEvents(0, nbr_ana);
	(*bmm)[make_pair(kEff_ana_bmm, channelIx)] = measurement_t(effCalc.GetEfficiency(0),(effCalc.GetEfficiencyErrorLow(0) + effCalc.GetEfficiencyErrorUp(0))/2.);
	if( (syst_it = systematics_table->find(g_sys_effana)) != systematics_table->end() ) {
		m = (*bmm)[make_pair(kEff_ana_bmm, channelIx)];
		m = measurement_t(0, m.getVal() * syst_it->second);
		(*bmm)[make_pair(kEff_ana_bmm, channelIx)] = m + (*bmm)[make_pair(kEff_ana_bmm, channelIx)];
	}
	
	// muon efficiency
	effCalc.SetTotalEvents(0, nbr_ana);
	effCalc.SetPassedEvents(0, nbr_mu);
	(*bmm)[make_pair(kEff_mu_bmm, channelIx)] = measurement_t(effCalc.GetEfficiency(0),(effCalc.GetEfficiencyErrorLow(0) + effCalc.GetEfficiencyErrorUp(0))/2.);
	
	switch (channelIx) {
		case 0: // Barrel
			syst_it = systematics_table->find(g_sys_effmu_barrel);
			break;
		case 1: // Endcap
			syst_it = systematics_table->find(g_sys_effmu_endcap);
			break;
		default:
			syst_it = systematics_table->end();
			break;
	}
	if ( syst_it != systematics_table->end() ) {
		m = (*bmm)[make_pair(kEff_mu_bmm, channelIx)];
		m = measurement_t(0, m.getVal() * syst_it->second);
		(*bmm)[make_pair(kEff_mu_bmm, channelIx)] = m + (*bmm)[make_pair(kEff_mu_bmm, channelIx)];
	}
	
	// trigger efficiency
	effCalc.SetTotalEvents(0, nbr_mu);
	effCalc.SetPassedEvents(0, nbr_trig);
	(*bmm)[make_pair(kEff_trig_bmm, channelIx)] = measurement_t(effCalc.GetEfficiency(0),(effCalc.GetEfficiencyErrorLow(0) + effCalc.GetEfficiencyErrorUp(0))/2.);
	
	switch (channelIx) {
		case 0:
			syst_it = systematics_table->find(g_sys_efftrig_barrel);
			break;
		case 1:
			syst_it = systematics_table->find(g_sys_efftrig_endcap);
			break;
		default:
			syst_it = systematics_table->end();
			break;
	}
	if ( syst_it != systematics_table->end() ) {
		m = (*bmm)[make_pair(kEff_trig_bmm, channelIx)];
		m = measurement_t(0, m.getVal() * syst_it->second);
		(*bmm)[make_pair(kEff_trig_bmm, channelIx)] = m + (*bmm)[make_pair(kEff_trig_bmm, channelIx)];
	}
	
	if (is_bstomumu) {
		cout << "\tComputing efficiency of peaking background..." << flush;
		accTree->SetEventList(0);
		rareTree->SetEventList(0);
		result = compute_peaking_eff(accTree, rareTree, histo_cut && acceptanceCutData && channelCut && anaCut, bd_mass_cut, bs_mass_cut, ((*bmm)[make_pair(kEff_mu_bmm, channelIx)]).getVal(), ((*bmm)[make_pair(kEff_trig_bmm,channelIx)]).getVal(),channelIx);
		cout << "\tdone" << endl;
	}
	
	return result;
} // estimate_bmm_eff()

void estimate_bmm_obs(std::map<bmm_param,measurement_t> *bmm, TTree *dataTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut, std::pair<double,double> bd_window, std::pair<double,double> bs_window, bool is_bstomumu, std::map<int,triplet<measurement_t> > *rare_effs, measurement_t tot_bu)
{
	TCut histo_cut(Form("%f < mass && mass < %f",low_histo_bound,high_histo_bound));
	TCut bd_mass_cut(Form("%f < mass && mass < %f",bd_window.first,bd_window.second));
	TCut bs_mass_cut(Form("%f < mass && mass < %f",bs_window.first,bs_window.second));
	TCut acceptanceCutData("TMath::Abs(eta_mu1) < 2.4 && TMath::Abs(eta_mu2) < 2.4 && (track_qual_mu1 & 4) && (track_qual_mu2 & 4) && (q_mu1*q_mu2 < 0)");
	TCut channelCut = TCut(Form("TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f",maxEta,maxEta)) && TCut(Form("!(TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f)",minEta,minEta));
	TCut muonCut("tight_mu1 && tight_mu2");
	TCut triggerCut("triggered_bs");
	TCut cut;
	triplet<measurement_t> exp_pkg;
	
	// peaking background
	cout << "\tComputing peaking background..." << flush;
	exp_pkg = compute_peaking_exp(*rare_effs, tot_bu);
	
	(*bmm)[make_pair(kPeakBkgOff_bmm, channelIx)] = exp_pkg.a;
	(*bmm)[make_pair(kPeakBkgOn_bmm, channelIx)] = (is_bstomumu ? exp_pkg.c : exp_pkg.b);
	cout << "\tdone" << endl;
	
	// observations
	cout << "\tComputing observation..." << flush;
	cut = TCut("candidate == 301313") && acceptanceCutData && channelCut && anaCut && muonCut && triggerCut && histo_cut && !bd_mass_cut && !bs_mass_cut;
	(*bmm)[make_pair(kObsBkg_bmm, channelIx)] = measurement_t((double)dataTree->Draw("",cut),0.0);
	
	if (!is_bstomumu)
		bs_mass_cut = bd_mass_cut;
	cut = TCut("candidate == 301313") && acceptanceCutData && channelCut && anaCut && muonCut && triggerCut && bs_mass_cut;
	(*bmm)[make_pair(kObsB_bmm, channelIx)] = measurement_t((double)dataTree->Draw("",cut),0.0);
	cout << "\tdone" << endl;
} // estiamte_bmm_obs()
