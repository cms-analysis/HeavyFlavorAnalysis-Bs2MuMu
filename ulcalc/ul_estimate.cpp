/*
 *  ul_estimate.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#include "ul_estimate.h"

// ROOT headers
#include <TStopwatch.h>

// RooFit headers
#include <RooProdPdf.h>
#include <RooRealVar.h>

// RooStats headers
#include <RooStats/ModelConfig.h>
#include <RooStats/FeldmanCousins.h>
#include <RooStats/BayesianCalculator.h>
#include <RooStats/MCMCCalculator.h>
#include <RooStats/HybridCalculator.h>

/* Add all the channels present in bmm to channels */
void add_channels(map<bmm_param,double> *bmm, set<int> *channels)
{
	map<bmm_param,double>::const_iterator it;
	for (it = bmm->begin(); it != bmm->end(); ++it)
		channels->insert(it->first.second);
} // add_channels()

RooWorkspace *build_model_nchannel(map<bmm_param,double> *bsmm, map<bmm_param,double> *bdmm, bool silent)
{
	RooStats::ModelConfig *splusbModel = NULL;
	RooStats::ModelConfig *bModel = NULL;
	RooWorkspace *wspace = new RooWorkspace;
	RooProdPdf *totalPdf = NULL;
	RooArgSet observables;
	RooArgSet nuisanceParams;
	set<int> channels;
	set<int>::const_iterator chan;
	
	add_channels(bsmm,&channels);
	add_channels(bdmm,&channels);
	
	// make sure we cover the entire physical range
	wspace->factory("mu_s[1,0,1000]");	// initialize to standard model
	wspace->factory("mu_d[1,0,1000]");	// initialize to standard model
	
	// Create channel specific variables
	for (chan = channels.begin(); chan != channels.end(); ++chan) {
		
		// Observables
		wspace->factory(Form("NbObs_%d[0,1000]",*chan)); // Observed Background
		wspace->factory(Form("NsObs_%d[0,1000]",*chan)); // Observed Evts in Bs window
		wspace->factory(Form("NdObs_%d[0,1000]",*chan)); // Observed Evts in Bd window
		observables.add(*wspace->var(Form("NbObs_%d",*chan)));
		observables.add(*wspace->var(Form("NsObs_%d",*chan)));
		observables.add(*wspace->var(Form("NdObs_%d",*chan)));
		
		// nuisance parameter
		wspace->factory(Form("nu_b_%d[0,0,1000]",*chan)); // background strength
		
		// build the constants
		wspace->factory(Form("TauS_%d[%f]",*chan,compute_tau(bsmm, bdmm, *chan, true)));
		wspace->factory(Form("TauD_%d[%f]",*chan,compute_tau(bdmm, bdmm, *chan, false)));
		wspace->factory(Form("NuS_%d[%f]", *chan, (*bsmm)[make_pair(kExp_bmm, *chan)]));
		wspace->factory(Form("NuD_%d[%f]", *chan, (*bdmm)[make_pair(kExp_bmm, *chan)]));
		
		wspace->factory(Form("Pss_%d[%f]", *chan, (*bsmm)[make_pair(kProb_swind_bmm, *chan)]));
		wspace->factory(Form("Psd_%d[%f]", *chan, (*bdmm)[make_pair(kProb_swind_bmm, *chan)]));
		wspace->factory(Form("Pds_%d[%f]", *chan, (*bsmm)[make_pair(kProb_dwind_bmm, *chan)]));
		wspace->factory(Form("Pdd_%d[%f]", *chan, (*bdmm)[make_pair(kProb_dwind_bmm, *chan)]));
		
		wspace->factory(Form("Poisson::bkg_window_%d(NbObs_%d,nu_b_%d)",*chan,*chan,*chan));
		wspace->factory(Form("Poisson::bs_window_%d(NsObs_%d,FormulaVar::bs_mean_%d(\"TauS_%d*nu_b_%d + Pss_%d*NuS_%d*mu_s + Psd_%d*NuD_%d*mu_d\",{TauS_%d,nu_b_%d,Pss_%d,NuS_%d,mu_s,Psd_%d,NuD_%d,mu_d}))",*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan));
		wspace->factory(Form("Poisson::bd_window_%d(NdObs_%d,FormulaVar::bd_mean_%d(\"TauD_%d*nu_b_%d + Pds_%d*NuS_%d*mu_s + Pdd_%d*NuD_%d*mu_d\",{TauD_%d,nu_b_%d,Pds_%d,NuS_%d,mu_s,Pdd_%d,NuD_%d,mu_d}))",*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan));
	}
	
	// build the product
	totalPdf = new RooProdPdf("total_pdf","Total Model PDF",RooArgList(wspace->allPdfs()));
	wspace->import(*totalPdf);
	
	// uniform prior in case of bayesian code
	wspace->factory("Uniform::prior_mus(mu_s)");
	
	// define the sets
	wspace->defineSet("obs", observables);
	wspace->defineSet("poi", "mu_s"); // Parameter of interest (at the moment, just mu_s)
	nuisanceParams = wspace->allVars();
	nuisanceParams.remove(*wspace->set("obs"));
	nuisanceParams.remove(*wspace->set("poi"));
	RooStats::RemoveConstantParameters(&nuisanceParams);
	wspace->defineSet("nui", nuisanceParams);
	
	// Model configuration
	splusbModel = new RooStats::ModelConfig("splusbConfig");
	splusbModel->SetWorkspace(*wspace); // set the workspace
	splusbModel->SetPdf(*wspace->pdf("total_pdf"));
	splusbModel->SetParametersOfInterest(*wspace->set("poi"));
	splusbModel->SetObservables(*wspace->set("obs"));
	splusbModel->SetNuisanceParameters(*wspace->set("nui"));
	splusbModel->SetPriorPdf(*wspace->pdf("prior_mus"));
	wspace->import(*splusbModel);
	
	bModel = new RooStats::ModelConfig("bConfig");
	bModel->SetWorkspace(*wspace); // set the workspace
	bModel->SetPdf(*wspace->pdf("total_pdf"));
	bModel->SetParametersOfInterest(*wspace->set("poi"));
	bModel->SetObservables(*wspace->set("obs"));
	bModel->SetNuisanceParameters(*wspace->set("nui"));
	bModel->SetPriorPdf(*wspace->pdf("prior_mus"));
	wspace->import(*bModel);
	
	if (!silent) {
		cout << "-------------------------------------" << endl;
		cout << "Workspace Configuration:" << endl;
		wspace->Print();
		cout << "-------------------------------------" << endl;
		cout << "S+B ModelConfig configuration:" << endl;
		splusbModel->Print();
		cout << "-------------------------------------" << endl;
		cout << "B ModelConfig configuration:" << endl;
		bModel->Print();
		cout << "-------------------------------------" << endl;
		cout << "Variables:" << endl;
		wspace->allVars().Print("v");
		cout << "-------------------------------------" << endl;
	}
	
	delete bModel;
	delete splusbModel;
	delete totalPdf;
	
	return wspace;
} // build_model_nchannel()

RooDataSet *build_data(RooWorkspace *wspace, double nsObs, double ndObs, double nbObs)
{
	RooRealVar *ns = wspace->var("NsObs");
	RooRealVar *nd = wspace->var("NdObs");
	RooRealVar *nb = wspace->var("NbObs");
	RooArgSet vars(*ns,*nd,*nb);
	RooDataSet *data = new RooDataSet("data","",vars);
	
	// set the values...
	ns->setVal(nsObs);
	nd->setVal(ndObs);
	nb->setVal(nbObs);
	
	// add to dataset
	data->add(vars);
		
	return data;
} // build_data()

RooDataSet *build_data_split(RooWorkspace *wspace, double nsObsB, double nsObsE, double ndObsB, double ndObsE, double nbObsB, double nbObsE)
{
	RooDataSet *data = new RooDataSet("data","",*wspace->set("obs"));
	RooArgSet vars;
	
	vars.addClone(*wspace->set("obs"));
	((RooRealVar&)vars["NbObs_B"]).setVal(nbObsB);
	((RooRealVar&)vars["NbObs_E"]).setVal(nbObsE);
	((RooRealVar&)vars["NsObs_B"]).setVal(nsObsB);
	((RooRealVar&)vars["NsObs_E"]).setVal(nsObsE);
	((RooRealVar&)vars["NdObs_B"]).setVal(ndObsB);
	((RooRealVar&)vars["NdObs_E"]).setVal(ndObsE);
	
	data->add(vars);
	
	return data;
} // build_data_split()

/* Start values are set the following way:
 *	nu_b = NbObs
 *	( Pss	Psd ) (mu_s NuS) + nu_b	(TauS) = (NsObs)
 *	( Pds	Pdd ) (mu_d NuD)		(TauD) = (NdObs)
 */
void estimate_start_values(RooWorkspace *wspace, RooDataSet *data, set<int> *channels)
{
	double p[2][2];
	double tau[2];
	double munu[2];
	double N[2];
	double nu_b;
	double det;
	RooStats::ModelConfig *splusbConfig = dynamic_cast<RooStats::ModelConfig*> (wspace->obj("splusbConfig"));
	RooStats::ModelConfig *bConfig = dynamic_cast<RooStats::ModelConfig*> (wspace->obj("bConfig"));
	
	for (set<int>::const_iterator ch = channels->begin(); ch != channels->end(); ++ch) {
		
		// set the constants
		p[0][0] = wspace->var(Form("Pss_%d",*ch))->getVal();
		p[0][1] = wspace->var(Form("Psd_%d",*ch))->getVal();
		p[1][0] = wspace->var(Form("Pds_%d",*ch))->getVal();
		p[1][1] = wspace->var(Form("Pdd_%d",*ch))->getVal();
		det = p[0][0]*p[1][1] - p[1][0]*p[0][1];
		
		tau[0] = wspace->var(Form("TauS_%d",*ch))->getVal();
		tau[1] = wspace->var(Form("TauD_%d",*ch))->getVal();
		
		N[0] = ((RooRealVar&)(*data->get(0))[Form("NsObs_%d",*ch)]).getVal();
		N[1] = ((RooRealVar&)(*data->get(0))[Form("NdObs_%d",*ch)]).getVal();
		
		// we need to set nui & poi
		nu_b = ((RooRealVar&)(*data->get(0))[Form("NbObs_%d",*ch)]).getVal();
		munu[0] = (p[1][1]*(N[0] - nu_b*tau[0]) - p[0][1]*(N[1] - nu_b*tau[1])) / det;
		munu[1] = (-p[1][0]*(N[0] - nu_b*tau[0]) + p[0][0]*(N[1]-  nu_b*tau[1])) / det;
		
		wspace->var(Form("nu_b_%d",*ch))->setVal(nu_b);
		wspace->var(Form("nu_b_%d",*ch))->setConstant(kFALSE);
		
		wspace->var("mu_s")->setVal(munu[0] / wspace->var(Form("NuS_%d",*ch))->getVal());
		wspace->var("mu_s")->setConstant(kFALSE);
		wspace->var("mu_d")->setVal(munu[1] / wspace->var(Form("NuD_%d",*ch))->getVal());
		wspace->var("mu_d")->setConstant(kFALSE);
	}
	
	// do a likelihood fit to the data to get the real values...
	wspace->pdf("total_pdf")->fitTo(*data);
	splusbConfig->SetSnapshot(*wspace->set("poi"));
	
	wspace->var("mu_s")->setVal(0.0);
	wspace->var("mu_d")->setVal(0.0);
	wspace->var("mu_s")->setConstant(kTRUE);
	wspace->var("mu_d")->setConstant(kTRUE);
	wspace->pdf("total_pdf")->fitTo(*data);
	wspace->var("mu_s")->setConstant(kFALSE);
	wspace->var("mu_d")->setConstant(kFALSE);
	bConfig->SetSnapshot(*wspace->set("poi"));
} // estimate_start_values()

// FIXME: Feldman-Cousins broken
RooStats::ConfInterval *est_ul_fc(RooWorkspace *wspace, RooDataSet *data, double cLevel, double *ulLimit, double *cpuUsed)
{
	using namespace RooStats;
	FeldmanCousins fc(*data,*(dynamic_cast<ModelConfig*>(wspace->obj("splusbConfig"))));
	PointSetInterval *psInterval = NULL;
	TStopwatch swatch;
	
	swatch.Start(kTRUE);
	
	// configure Feldman Cousins
	fc.SetTestSize(1.-cLevel);
	fc.UseAdaptiveSampling(true); // adaptive sampling (disable later on)
	fc.AdditionalNToysFactor(5.0);
	fc.SetNBins(25);
	fc.FluctuateNumDataEntries(false);
	
	// estimate_start_values(wspace,data);
	psInterval = fc.GetInterval();
	
	if (ulLimit)
		*ulLimit = psInterval->UpperLimit(*wspace->var("mu_s"));
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	return psInterval;
} // est_ul()

RooStats::ConfInterval *est_ul_bc(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, double cLevel, double *ulLimit, double *cpuUsed)
{
	using namespace RooStats;
	BayesianCalculator bc(*data,*(dynamic_cast<ModelConfig*>(wspace->obj("splusbConfig"))));
	SimpleInterval *simpleInt = NULL;
	TStopwatch swatch;
	
	swatch.Start(kTRUE);
	
	// configure BayesianCalculator
	bc.SetConfidenceLevel(cLevel);
	bc.SetLeftSideTailFraction(0.0); // compute upper limit
	
	estimate_start_values(wspace, data, channels);
	simpleInt = bc.GetInterval();
	
	if (ulLimit)
		*ulLimit = simpleInt->UpperLimit();
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	return simpleInt;
} // est_ul_bc()

// FIXME: Markov Chain MC Broken
RooStats::ConfInterval *est_ul_mc(RooWorkspace *wspace, RooDataSet *data, double cLevel, double *ulLimit, bool splitModel, double *cpuUsed)
{
	using namespace RooStats;
	MCMCCalculator mc(*data,*(dynamic_cast<ModelConfig*> (wspace->obj("splusbConfig"))));
	MCMCInterval *mcInt = NULL;
	TStopwatch swatch;
	
	swatch.Start(kTRUE);
	
	// configure MCMCCalculator
	mc.SetNumBins(1000);
	mc.SetNumIters(100000);
	mc.SetNumBurnInSteps(10);
	mc.SetConfidenceLevel(cLevel);
	mc.SetLeftSideTailFraction(0.0); // compute upper limit
	
//	estimate_start_values(wspace, data);
	mcInt = mc.GetInterval();

	if (ulLimit)
		*ulLimit = mcInt->UpperLimit(*wspace->var("mu_s"));
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	return mcInt;
} // est_ul_mc()

// FIXME: CLs not yet implemented
RooStats::ConfInterval *est_ul_cls(RooWorkspace *wspace, RooDataSet *data, double cLevel, double *ulLimit, bool splitModel, double *cpuUsed)
{
	using namespace RooStats;
	ModelConfig *bModel = dynamic_cast<ModelConfig*> (wspace->obj("bConfig"));
	ModelConfig *sbModel = dynamic_cast<ModelConfig*> (wspace->obj("splusbConfig"));
	HybridCalculator hc(*data,*sbModel,*bModel);
	ToyMCSampler *sampler = (ToyMCSampler*)hc.GetTestStatSampler();
	
	
	sampler->SetNEventsPerToy(1); // number counting
	hc.SetToys(1000, 1000);
	
	
	return NULL;
} // est_ul_cls()

void compute_vars(map<bmm_param,double> *bmm, bool bstomumu)
{
	double tot_bplus,nu;
	set<int> channels;
	set<int>::const_iterator it;
	double c = bstomumu ? c_s_theory : c_d_theory;
	
	add_channels(bmm, &channels);
	
	for (it = channels.begin(); it != channels.end(); ++it) {
		
		if (bmm->count(make_pair(kTot_bplus, *it)) > 0)
			tot_bplus = (*bmm)[make_pair(kTot_bplus, *it)];
		else
			tot_bplus = ((*bmm)[make_pair(kObs_bplus, *it)]) / compute_efftot_bplus(bmm,*it);
		
		nu = c * compute_efftot_bmm(bmm,*it) * tot_bplus;
		
		(*bmm)[make_pair(kExp_bmm, *it)] = nu;
	}
} // compute_vars()

double compute_efftot_bplus(map<bmm_param,double> *bmm, int channel)
{
	return (*bmm)[make_pair(kAcc_bplus,channel)] * (*bmm)[make_pair(kEff_mu_bplus,channel)] * (*bmm)[make_pair(kEff_trig_bplus,channel)] * (*bmm)[make_pair(kEff_cand_bplus,channel)] * (*bmm)[make_pair(kEff_ana_bplus,channel)];
} // compute_efftot_bplus()

double compute_efftot_bmm(map<bmm_param,double> *bmm, int channel)
{
	double eff;
	
	if(bmm->count(make_pair(kEff_total_bmm, channel)) > 0)
		eff = (*bmm)[make_pair(kEff_total_bmm, channel)];
	else
		eff = (*bmm)[make_pair(kAcc_bmm,channel)] * (*bmm)[make_pair(kEff_mu_bmm,channel)] * (*bmm)[make_pair(kEff_trig_bmm,channel)] * (*bmm)[make_pair(kEff_cand_bmm,channel)] * (*bmm)[make_pair(kEff_ana_bmm,channel)];
	
	return eff;
} // compute_efftot_bmm()
