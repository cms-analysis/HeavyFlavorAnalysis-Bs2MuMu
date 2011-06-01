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
void add_channels(map<bmm_param,measurement_t> *bmm, set<int> *channels)
{
	map<bmm_param,measurement_t>::const_iterator it;
	for (it = bmm->begin(); it != bmm->end(); ++it)
		channels->insert(it->first.second);
} // add_channels()

RooWorkspace *build_model_nchannel(map<bmm_param,measurement_t> *bsmm, map<bmm_param,measurement_t> *bdmm, bool no_errors, int verbosity) {
	RooStats::ModelConfig *splusbModel = NULL;
	RooStats::ModelConfig *bModel = NULL;
	RooWorkspace *wspace = new RooWorkspace("wspace");
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
		
		/////////////////////////
		// Construction of NuS //
		/////////////////////////
		if ( ((*bsmm)[make_pair(kExp_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bsmm)[make_pair(kExp_bmm, *chan)];
			wspace->factory(Form("NuS0_%d[%f]", *chan, m.getVal())); // fixed mean variable
			wspace->factory(Form("NuSErr_%d[%f]", *chan, m.getErr())); // fixed error variable
			wspace->factory(Form("NuS_%d[%f,%f,%f]", *chan, m.getVal(), 0.0, m.getVal() + 100*m.getErr()));
			wspace->factory(Form("Gaussian::NuS_Gauss_%d(NuS_%d,NuS0_%d,NuSErr_%d)",*chan,*chan,*chan,*chan)); // error gaussian.
		} else {
			// no error associated to this variable. just create the default one
			wspace->factory(Form("NuS_%d[%f]", *chan, ((*bsmm)[make_pair(kExp_bmm, *chan)]).getVal()));
		}
		
		/////////////////////////
		// Construction of NuD //
		/////////////////////////
		if ( ((*bdmm)[make_pair(kExp_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bdmm)[make_pair(kExp_bmm, *chan)];
			wspace->factory(Form("NuD0_%d[%f]", *chan, m.getVal())); // fixed mean variable
			wspace->factory(Form("NuDErr_%d[%f]", *chan, m.getErr())); // fixed error variable
			wspace->factory(Form("NuD_%d[%f,%f,%f]", *chan, m.getVal(), 0.0, m.getVal() + 100*m.getErr()));
			wspace->factory(Form("Gaussian::NuD_Gauss_%d(NuD_%d,NuD0_%d,NuDErr_%d)",*chan,*chan,*chan,*chan)); // error gaussian
		} else {
			// no error assigned, just make a constant
			wspace->factory(Form("NuD_%d[%f]", *chan, ((*bdmm)[make_pair(kExp_bmm, *chan)]).getVal()));
		}
		
		/////////////////////////
		// Construction of Pss //
		/////////////////////////
		if ( ((*bsmm)[make_pair(kProb_swind_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bsmm)[make_pair(kProb_swind_bmm, *chan)];
			wspace->factory(Form("Pss0_%d[%f]", *chan, m.getVal()));
			wspace->factory(Form("PssErr_%d[%f]", *chan, m.getErr()));
			wspace->factory(Form("Pss_%d[%f,%f,%f]", *chan, m.getVal(), 0.0, 1.0));
			wspace->factory(Form("Gaussian::Pss_Gauss_%d(Pss_%d,Pss0_%d,PssErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			// no error assigned, just make a constant.
			wspace->factory(Form("Pss_%d[%f]", *chan, ((*bsmm)[make_pair(kProb_swind_bmm, *chan)]).getVal()));
		}
		
		/////////////////////////
		// Construction of Psd //
		/////////////////////////
		if ( ((*bdmm)[make_pair(kProb_swind_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bdmm)[make_pair(kProb_swind_bmm, *chan)];
			wspace->factory(Form("Psd0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PsdErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("Psd_%d[%f,%f,%f]", *chan, m.getVal(), 0.0, 1.0));
			wspace->factory(Form("Gaussian::Psd_Gauss_%d(Psd_%d,Psd0_%d,PsdErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			// no error assigned, just make a constant.
			wspace->factory(Form("Psd_%d[%f]", *chan, ((*bdmm)[make_pair(kProb_swind_bmm, *chan)]).getVal()));
		}
		
		/////////////////////////
		// Construction of Pds //
		/////////////////////////
		if ( ((*bsmm)[make_pair(kProb_dwind_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bsmm)[make_pair(kProb_dwind_bmm, *chan)];
			wspace->factory(Form("Pds0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PdsErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("Pds_%d[%f,%f,%f]",*chan,m.getVal(),0.0,1.0));
			wspace->factory(Form("Gaussian::Pds_Gauss_%d(Pds_%d,Pds0_%d,PdsErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			// no error assigned, just make a constant.
			wspace->factory(Form("Pds_%d[%f]", *chan, ((*bsmm)[make_pair(kProb_dwind_bmm, *chan)]).getVal()));
		}
		
		/////////////////////////
		// Construction of Pdd //
		/////////////////////////
		if ( ((*bdmm)[make_pair(kProb_dwind_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bdmm)[make_pair(kProb_dwind_bmm, *chan)];
			wspace->factory(Form("Pdd0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PddErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("Pdd_%d[%f,%f,%f]",*chan,m.getVal(),0.0,1.0));
			wspace->factory(Form("Gaussian::Pdd_Gauss_%d(Pdd_%d,Pdd0_%d,PddErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("Pdd_%d[%f]", *chan, ((*bdmm)[make_pair(kProb_dwind_bmm, *chan)]).getVal()));
		}
		
		wspace->factory(Form("Poisson::bkg_window_%d(NbObs_%d,nu_b_%d)",*chan,*chan,*chan));
		wspace->factory(Form("Poisson::bs_window_%d(NsObs_%d,FormulaVar::bs_mean_%d(\"TauS_%d*nu_b_%d + Pss_%d*NuS_%d*mu_s + Psd_%d*NuD_%d*mu_d\",{TauS_%d,nu_b_%d,Pss_%d,NuS_%d,mu_s,Psd_%d,NuD_%d,mu_d}))",*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan));
		wspace->factory(Form("Poisson::bd_window_%d(NdObs_%d,FormulaVar::bd_mean_%d(\"TauD_%d*nu_b_%d + Pds_%d*NuS_%d*mu_s + Pdd_%d*NuD_%d*mu_d\",{TauD_%d,nu_b_%d,Pds_%d,NuS_%d,mu_s,Pdd_%d,NuD_%d,mu_d}))",*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan));
	}
	
	// build the product
	totalPdf = new RooProdPdf("total_pdf","Total Model PDF",RooArgList(wspace->allPdfs()));
	wspace->import(*totalPdf, ((verbosity > 0) ? RooCmdArg::none() : RooFit::Silence(kTRUE)) );
	
	// uniform prior in case of bayesian code
	wspace->factory("Uniform::prior_mus(mu_s)");
	
	// define the sets
	wspace->defineSet("obs", observables);
	wspace->defineSet("poi", "mu_s"); // Parameter of interest (at the moment, just mu_s)
	nuisanceParams.addClone(wspace->allVars());
	nuisanceParams.remove(*wspace->set("obs"), kTRUE, kTRUE);
	nuisanceParams.remove(*wspace->set("poi"), kTRUE, kTRUE);
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
	
	if (verbosity > 0) {
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

/* Start values are set the following way:
 *	nu_b = NbObs
 *	( Pss	Psd ) (mu_s NuS) + nu_b	(TauS) = (NsObs)
 *	( Pds	Pdd ) (mu_d NuD)		(TauD) = (NdObs)
 */
void estimate_start_values(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, int verbosity)
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
	wspace->pdf("total_pdf")->fitTo(*data, ((verbosity > 0) ? RooCmdArg::none() : RooFit::PrintLevel(-1)));
	splusbConfig->SetSnapshot(*wspace->set("poi"));
	
	wspace->var("mu_s")->setVal(0.0);
	wspace->var("mu_d")->setVal(0.0);
	wspace->var("mu_s")->setConstant(kTRUE);
	wspace->var("mu_d")->setConstant(kTRUE);
	wspace->pdf("total_pdf")->fitTo(*data, ((verbosity > 0) ? RooCmdArg::none() : RooFit::PrintLevel(-1)));
	wspace->var("mu_s")->setConstant(kFALSE);
	wspace->var("mu_d")->setConstant(kFALSE);
	bConfig->SetSnapshot(*wspace->set("poi"));
} // estimate_start_values()

RooStats::ConfInterval *est_ul_fc(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, double cLevel, int verbosity, uint32_t nbins, pair<double,double> *rg, double err, double *ulLimit, double *cpuUsed)
{
	using namespace RooStats;
	FeldmanCousins fc(*data,*(dynamic_cast<ModelConfig*>(wspace->obj("splusbConfig"))));
	PointSetInterval *psInterval = NULL;
	ConfidenceBelt *belt;
	TStopwatch swatch;
	RooDataSet *custom_points;
	RooDataSet *pointsInInterval;
	RooAbsData *paramPts;
	Double_t v;
	Int_t j;
	
	swatch.Start(kTRUE);
	
	// configure Feldman Cousins
	fc.CreateConfBelt(true);
	fc.SetTestSize(1.-cLevel);
//	fc.UseAdaptiveSampling(true); // adaptive sampling (disable later on)
	fc.AdditionalNToysFactor(6.0);
	fc.SetNBins(nbins);
	fc.FluctuateNumDataEntries(false);
	if (rg) {
		RooArgSet poi;
		poi.addClone(*wspace->set("poi"));
		custom_points = new RooDataSet("","",poi);
		
		for (j = 0; j < (Int_t)nbins; j++) {
			v = rg->first + (rg->second - rg->first) * j / (nbins - 1);
			((RooRealVar*)poi.first())->setVal(v); // Currenlty only one parameter of interest supported
			custom_points->add(poi);
		}
		fc.SetPOIPointsToTest(*custom_points); // note, this is owned by Feldman-Cousins.
	}
	
	estimate_start_values(wspace, data, channels, verbosity);
	psInterval = fc.GetInterval();
	
	// if 'err' specified, correct for numerical inaccurancies...
	if (err > 0) {
		// construct our own PointSetInterval
		pointsInInterval = new RooDataSet("pointsInInterval","points in interval", *(fc.GetPointsToScan()->get(0)));
		delete psInterval; // delete the old one (we are replacing it!)
		
		belt = fc.GetConfidenceBelt();
		paramPts = fc.GetPointsToScan();
		for (j = 0; j < paramPts->numEntries(); j++) {
			RooArgSet params;
			params.add(*paramPts->get(j));
			v = fc.GetTestStatSampler()->EvaluateTestStatistic(*data, params);
			if (v - err < belt->GetAcceptanceRegionMax(params))
				pointsInInterval->add(params);
		}
		
		psInterval = new PointSetInterval("ClassicalConfidenceInterval",*pointsInInterval);
		psInterval->SetConfidenceLevel(cLevel);
	}
	
	if (ulLimit)
		*ulLimit = psInterval->UpperLimit(*wspace->var("mu_s"));
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	return psInterval;
} // est_ul()

RooStats::ConfInterval *est_ul_bc(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, double cLevel, int verbosity, double *ulLimit, double *cpuUsed)
{
	using namespace RooStats;
	BayesianCalculator bc(*data,*(dynamic_cast<ModelConfig*>(wspace->obj("splusbConfig"))));
	SimpleInterval *simpleInt = NULL;
	TStopwatch swatch;
	
	swatch.Start(kTRUE);
	
	// configure BayesianCalculator
	bc.SetConfidenceLevel(cLevel);
	bc.SetLeftSideTailFraction(0.0); // compute upper limit
	
	//bc.SetIntegrationType("TOYMC");
	// FIXME: Multidimensional Integration funktioniert nicht mehr!
	//	bc.SetIntegrationType("PLAIN");
	
	estimate_start_values(wspace, data, channels, verbosity);
	simpleInt = bc.GetInterval();
	
	if (ulLimit)
		*ulLimit = simpleInt->UpperLimit();
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	return simpleInt;
} // est_ul_bc()

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

void compute_vars(map<bmm_param,measurement_t> *bmm, bool bstomumu)
{
	set<int> channels;
	set<int>::const_iterator it;
	measurement_t tot_bplus,nu;
	measurement_t c = bstomumu ? c_s_theory() : c_d_theory();
	
	add_channels(bmm, &channels);
	
	for (it = channels.begin(); it != channels.end(); ++it) {
		
		if (bmm->count(make_pair(kTot_bplus, *it)) > 0)
			tot_bplus = (*bmm)[make_pair(kTot_bplus, *it)];
		else
			tot_bplus = (*bmm)[make_pair(kObs_bplus, *it)] / compute_efftot_bplus(bmm,*it);
		
		nu = c * compute_efftot_bmm(bmm,*it) * tot_bplus;
		
		(*bmm)[make_pair(kExp_bmm, *it)] = nu;
	}
} // compute_vars()

measurement_t compute_efftot_bplus(map<bmm_param,measurement_t> *bmm, int channel)
{
	measurement_t eff = (*bmm)[make_pair(kAcc_bplus,channel)] * (*bmm)[make_pair(kEff_mu_bplus,channel)] * (*bmm)[make_pair(kEff_trig_bplus,channel)] * (*bmm)[make_pair(kEff_cand_bplus,channel)] * (*bmm)[make_pair(kEff_ana_bplus,channel)];
	
	return eff;
} // compute_efftot_bplus()

measurement_t compute_efftot_bmm(map<bmm_param,measurement_t> *bmm, int channel)
{
	measurement_t eff;
	
	if(bmm->count(make_pair(kEff_total_bmm, channel)) > 0)
		eff = (*bmm)[make_pair(kEff_total_bmm, channel)];
	else
		eff = (*bmm)[make_pair(kAcc_bmm,channel)] * (*bmm)[make_pair(kEff_mu_bmm,channel)] * (*bmm)[make_pair(kEff_trig_bmm,channel)] * (*bmm)[make_pair(kEff_cand_bmm,channel)] * (*bmm)[make_pair(kEff_ana_bmm,channel)];
	
	return eff;
} // compute_efftot_bmm()
