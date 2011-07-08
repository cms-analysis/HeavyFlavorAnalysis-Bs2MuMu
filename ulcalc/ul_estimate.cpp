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
#include <RooStats/RooStatsUtils.h>
#include <RooStats/FrequentistCalculator.h>
#include <RooStats/HypoTestInverter.h>
#include <RooStats/ProfileLikelihoodTestStat.h>
#include <RooStats/RatioOfProfiledLikelihoodsTestStat.h>


/* Add all the channels present in bmm to channels */
void add_channels(map<bmm_param,measurement_t> *bmm, set<int> *channels)
{
	map<bmm_param,measurement_t>::const_iterator it;
	for (it = bmm->begin(); it != bmm->end(); ++it)
		channels->insert(it->first.second);
} // add_channels()

RooWorkspace *build_model_nchannel(map<bmm_param,measurement_t> *bsmm, map<bmm_param,measurement_t> *bdmm, bool no_errors, int verbosity, bool compute_bd_ul)
{
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
	wspace->factory("mu_s[1,0,20]");	// initialize to standard model
	wspace->factory("mu_d[1,0,200]");	// initialize to standard model
	
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
		
		///////////////////////////////////
		// build the background ratio Bs //
		///////////////////////////////////
		if ( ((*bsmm)[make_pair(kTau_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bsmm)[make_pair(kTau_bmm, *chan)];
			wspace->factory(Form("TauS0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("TauSErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("TauS_%d[%f,%f,%f]",*chan,m.getVal(),0.0,1.0));
			wspace->factory(Form("Gaussian::TauS_Gauss_%d(TauS_%d,TauS0_%d,TauSErr_%d)",*chan,*chan,*chan,*chan));
		}
		else {
			wspace->factory(Form("TauS_%d[%f]",*chan, ((*bsmm)[make_pair(kTau_bmm, *chan)]).getVal()));
		}
		
		////////////////////////////////////
		// build the background ration Bd //
		////////////////////////////////////
		if ( ((*bdmm)[make_pair(kTau_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bdmm)[make_pair(kTau_bmm, *chan)];
			wspace->factory(Form("TauD0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("TauDErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("TauD_%d[%f,%f,%f]",*chan,m.getVal(),0.0,1.0));
			wspace->factory(Form("Gaussian::TauD_Gauss_%d(TauD_%d,TauD0_%d,TauDErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("TauD_%d[%f]",*chan,((*bdmm)[make_pair(kTau_bmm, *chan)]).getVal()));
		}
		
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
		
		//////////////////////////////////////
		// Construction of rare backgrounds //
		//////////////////////////////////////
		if ( ((*bsmm)[make_pair(kPeakBkgOff_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bsmm)[make_pair(kPeakBkgOff_bmm, *chan)];
			wspace->factory(Form("PeakBkgSB0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PeakBkgSBErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("PeakBkgSB_%d[%f,%f,%f]",*chan,m.getVal(),0.0,m.getVal() + 10*m.getErr()));
			wspace->factory(Form("Gaussian::PeakBkgSB_Gauss_%d(PeakBkgSB_%d,PeakBkgSB0_%d,PeakBkgSBErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("PeakBkgSB_%d[%f]", *chan, ((*bsmm)[make_pair(kPeakBkgOff_bmm, *chan)]).getVal()));
		}
		
		if ( ((*bsmm)[make_pair(kPeakBkgOn_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bsmm)[make_pair(kPeakBkgOn_bmm, *chan)];
			wspace->factory(Form("PeakBkgBs0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PeakBkgBsErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("PeakBkgBs_%d[%f,%f,%f]",*chan,m.getVal(),0.0,m.getVal() + 10*m.getErr()));
			wspace->factory(Form("Gaussian::PeakBkgBs_Gauss_%d(PeakBkgBs_%d,PeakBkgBs0_%d,PeakBkgBsErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("PeakBkgBs_%d[%f]",*chan,((*bsmm)[make_pair(kPeakBkgOn_bmm, *chan)]).getVal()));
		}
		
		if ( ((*bdmm)[make_pair(kPeakBkgOn_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			measurement_t m = (*bdmm)[make_pair(kPeakBkgOn_bmm, *chan)];
			wspace->factory(Form("PeakBkgBd0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PeakBkgBdErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("PeakBkgBd_%d[%f,%f,%f]",*chan, m.getVal(), 0.0 ,m.getVal() + 10*m.getErr()));
			wspace->factory(Form("Gaussian::PeakBkgBd_Gauss_%d(PeakBkgBd_%d,PeakBkgBd0_%d,PeakBkgBdErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("PeakBkgBd_%d[%f]",*chan,((*bdmm)[make_pair(kPeakBkgOn_bmm, *chan)]).getVal()));
		}
		
		wspace->factory(Form("Poisson::bkg_window_%d(NbObs_%d,FormulaVar::bkg_mean_%d(\"nu_b_%d + PeakBkgSB_%d\",{nu_b_%d,PeakBkgSB_%d}))",*chan,*chan,*chan,*chan,*chan,*chan,*chan));
		
		wspace->factory(Form("Poisson::bs_window_%d(NsObs_%d,FormulaVar::bs_mean_%d(\"TauS_%d*nu_b_%d + PeakBkgBs_%d + Pss_%d*NuS_%d*mu_s + Psd_%d*NuD_%d*mu_d\",{TauS_%d,nu_b_%d,PeakBkgBs_%d,Pss_%d,NuS_%d,mu_s,Psd_%d,NuD_%d,mu_d}))",*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan));
		
		wspace->factory(Form("Poisson::bd_window_%d(NdObs_%d,FormulaVar::bd_mean_%d(\"TauD_%d*nu_b_%d + PeakBkgBd_%d + Pds_%d*NuS_%d*mu_s + Pdd_%d*NuD_%d*mu_d\",{TauD_%d,nu_b_%d,PeakBkgBd_%d,Pds_%d,NuS_%d,mu_s,Pdd_%d,NuD_%d,mu_d}))",*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan));
	}
	
	// build the product
	totalPdf = new RooProdPdf("total_pdf","Total Model PDF",RooArgList(wspace->allPdfs()));
	wspace->import(*totalPdf, ((verbosity > 0) ? RooCmdArg::none() : RooFit::Silence(kTRUE)) );
	
	// uniform prior in case of bayesian code
	wspace->factory(Form("Uniform::prior_mu(mu_%c)",(compute_bd_ul ? 'd' : 's')));
	
	// define the sets
	wspace->defineSet("obs", observables);
	if (compute_bd_ul)	wspace->defineSet("poi", "mu_d");
	else				wspace->defineSet("poi", "mu_s");
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
	splusbModel->SetPriorPdf(*wspace->pdf("prior_mu"));
	wspace->import(*splusbModel);
	
	bModel = new RooStats::ModelConfig("bConfig");
	bModel->SetWorkspace(*wspace); // set the workspace
	bModel->SetPdf(*wspace->pdf("total_pdf"));
	bModel->SetParametersOfInterest(*wspace->set("poi"));
	bModel->SetObservables(*wspace->set("obs"));
	bModel->SetNuisanceParameters(*wspace->set("nui"));
	bModel->SetPriorPdf(*wspace->pdf("prior_mu"));
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

// FIXME: introduce gamma prior for the background, too!
// FIXME: Merge with the above one.
RooWorkspace *build_model_light(map<bmm_param,measurement_t> *bsmm, int verbosity)
{
	RooStats::ModelConfig *splusbModel = NULL;
	RooWorkspace *wspace = new RooWorkspace("wspace");
	RooProdPdf *totalPdf = NULL;
	RooArgSet observables;
	RooArgSet nuisanceParams;
	
	// make sure we cover the entire physical range
	wspace->factory("mu_s[1,0,1000]");	// initialize to standard model
	wspace->factory("mu_d[1,0,1000]");	// initialize to standard model
	
	// Observables
	wspace->factory("NbObs[0,1000]"); // Observed Background
	wspace->factory("NsObs[0,1000]"); // Observed Evts in Bs window
	observables.add(*wspace->var("NbObs"));
	observables.add(*wspace->var("NsObs"));

	// nuisance parameter
	wspace->factory("nu_b[0,0,1000]"); // background strength
	
	// build the constants
	wspace->factory(Form("TauS[%f]", ((*bsmm)[make_pair(kTau_bmm, 0)]).getVal()));
	wspace->factory(Form("NuS[%f]", ((*bsmm)[make_pair(kExp_bmm, 0)]).getVal()));
	wspace->factory(Form("Pss[%f]", ((*bsmm)[make_pair(kProb_swind_bmm, 0)]).getVal()));
	
	// pdfs
	wspace->factory("Poisson::bkg_window(NbObs,nu_b)");
	wspace->factory("Poisson::bs_window(NsObs,FormulaVar::bs_mean(\"TauS*nu_b + Pss*NuS*mu_s\",{TauS,nu_b,Pss,NuS,mu_s}))");
	totalPdf = new RooProdPdf("total_pdf","Total Model PDF",RooArgList(wspace->allPdfs()));
	wspace->import(*totalPdf, ((verbosity > 0) ? RooCmdArg::none() : RooFit::Silence(kTRUE)) );
	wspace->factory("Uniform::prior_mu(mu_s)");
	
	wspace->defineSet("obs", observables);
	wspace->defineSet("poi", "mu_s"); // Parameter of interest (at the moment, just mu_s)
	nuisanceParams.addClone(wspace->allVars());
	nuisanceParams.remove(*wspace->set("obs"), kTRUE, kTRUE);
	nuisanceParams.remove(*wspace->set("poi"), kTRUE, kTRUE);
	RooStats::RemoveConstantParameters(&nuisanceParams);
	wspace->defineSet("nui", nuisanceParams);
	
	splusbModel = new RooStats::ModelConfig("splusbConfig");
	splusbModel->SetWorkspace(*wspace); // set the workspace
	splusbModel->SetPdf(*wspace->pdf("total_pdf"));
	splusbModel->SetParametersOfInterest(*wspace->set("poi"));
	splusbModel->SetObservables(*wspace->set("obs"));
	splusbModel->SetNuisanceParameters(*wspace->set("nui"));
	splusbModel->SetPriorPdf(*wspace->pdf("prior_mu"));
	wspace->import(*splusbModel);
	
	if (verbosity > 0) {
		cout << "-------------------------------------" << endl;
		cout << "Workspace Configuration:" << endl;
		wspace->Print();
		cout << "-------------------------------------" << endl;
		cout << "S+B ModelConfig configuration:" << endl;
		splusbModel->Print();
		cout << "-------------------------------------" << endl;
		cout << "Variables:" << endl;
		wspace->allVars().Print("v");
		cout << "-------------------------------------" << endl;
	}
	
	delete splusbModel;
	delete totalPdf;
	
	return wspace;
} // build_model_light()

/* Start values are set the following way:
 *	nu_b = NbObs
 *	( Pss	Psd ) (mu_s NuS) + nu_b	(TauS) = (NsObs)
 *	( Pds	Pdd ) (mu_d NuD)		(TauD) = (NdObs)
 */
void measure_params(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, int verbosity)
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
} // measure_params()

RooStats::ConfInterval *est_ul_fc(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, double cLevel, int verbosity, uint32_t nbins, pair<double,double> *rg, double *ulLimit, double *loLimit, double *cpuUsed)
{
	using namespace RooStats;
	ModelConfig *splusbConfig = dynamic_cast<ModelConfig*>(wspace->obj("splusbConfig"));
	FeldmanCousins fc(*data,*splusbConfig);
	PointSetInterval *psInterval = NULL;
	TStopwatch swatch;
	RooDataSet *custom_points;
	Double_t v;
	Int_t j;
	
	swatch.Start(kTRUE);
	
	// configure Feldman Cousins
	fc.CreateConfBelt(true);
	fc.SetTestSize(1.-cLevel);
//	fc.UseAdaptiveSampling(true); // adaptive sampling (disable later on)
	fc.AdditionalNToysFactor(10.0);
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
	
	measure_params(wspace, data, channels, verbosity);
	splusbConfig->LoadSnapshot();
	psInterval = fc.GetInterval();
	
	if (ulLimit)
		*ulLimit = psInterval->UpperLimit((RooRealVar&)*wspace->set("poi")->first());
	if (loLimit)
		*loLimit = psInterval->LowerLimit((RooRealVar&)*wspace->set("poi")->first());
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	return psInterval;
} // est_ul()

RooStats::ConfInterval *est_ul_bc(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, double cLevel, int verbosity, double *ulLimit, double *cpuUsed)
{
	using namespace RooStats;
	BayesianCalculator bc(*data,*(dynamic_cast<ModelConfig*>(wspace->obj("splusbConfig"))));
	RooStats::ModelConfig *splusbConfig;
	SimpleInterval *simpleInt = NULL;
	TStopwatch swatch;
	
	swatch.Start(kTRUE);
	
	// configure BayesianCalculator
	bc.SetConfidenceLevel(cLevel);
	bc.SetLeftSideTailFraction(0.0); // compute upper limit
	
	measure_params(wspace, data, channels, verbosity);
	
	// Load the Snapshot of the s+b configuration
	splusbConfig = dynamic_cast<RooStats::ModelConfig*> (wspace->obj("splusbConfig"));
	splusbConfig->LoadSnapshot();
	simpleInt = bc.GetInterval();
	
	if (ulLimit)
		*ulLimit = simpleInt->UpperLimit();
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	return simpleInt;
} // est_ul_bc()

// FIXME: Incorporate this one to the 'standard' bayesian estimate
RooStats::ConfInterval *est_ul_bc_light(RooWorkspace *wspace, RooDataSet *data, double cLevel, int verbosity, double *upperLimit, double *cpuUsed)
{
	using namespace RooStats;
	ModelConfig *splusbConfig = dynamic_cast<ModelConfig*>(wspace->obj("splusbConfig"));
	BayesianCalculator bc(*data,*splusbConfig);
	SimpleInterval *simpleInt = NULL;
	TStopwatch swatch;
	double nu_b,mu_s;
	
	swatch.Start(kTRUE);
	
	// configure
	bc.SetConfidenceLevel(cLevel);
	bc.SetLeftSideTailFraction(0.0); // compute upper limit
	
	// set variables...
	nu_b = ((RooRealVar&)(*data->get(0))["NbObs"]).getVal();
	mu_s = (((RooRealVar&)(*data->get(0))["NsObs"]).getVal() - wspace->var("TauS")->getVal()* nu_b) / (wspace->var("NuS")->getVal() * wspace->var("Pss")->getVal());
	
	wspace->var("nu_b")->setVal(nu_b);
	wspace->var("mu_s")->setVal(mu_s);
	
	wspace->pdf("total_pdf")->fitTo(*data, ((verbosity > 0) ? RooCmdArg::none() : RooFit::PrintLevel(-1)));
	splusbConfig->SetSnapshot(*wspace->set("poi"));
	
	// get interval...
	simpleInt = bc.GetInterval();
	if (upperLimit) *upperLimit = simpleInt->UpperLimit();
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	return simpleInt;
} // est_ul_bc_light()

RooStats::ConfInterval *est_ul_cls(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, double cLevel, int verbosity, double err, double *ulLimit, double *cpuUsed)
{
	using namespace RooStats;
	ModelConfig *bModel = dynamic_cast<ModelConfig*> (wspace->obj("bConfig"));
	ModelConfig *sbModel = dynamic_cast<ModelConfig*> (wspace->obj("splusbConfig"));
	ProfileLikelihoodTestStat testStat(*wspace->pdf("total_pdf"));
	ToyMCSampler *mcSampler = new ToyMCSampler(testStat,1000);
	FrequentistCalculator frequCalc(*data,*bModel,*sbModel,mcSampler); // Note null = sb, alt = b
	HypoTestInverter hypoInv(frequCalc, (RooRealVar*)sbModel->GetParametersOfInterest()->first(), 1.0 - cLevel);
	double limitErr = 0.1; // FIXME: Understand this parameter
	TStopwatch swatch;
	bool ok;
	
	swatch.Start(kTRUE);
	
	mcSampler->SetNEventsPerToy(1);
	hypoInv.SetAutoScan();
	hypoInv.UseCLs();
	hypoInv.SetTestSize(1.0 - cLevel);
	hypoInv.SetVerbose(verbosity);
		
	measure_params(wspace, data, channels, verbosity);
	sbModel->LoadSnapshot();
	ok = hypoInv.RunLimit(*ulLimit,limitErr);
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	delete mcSampler;
	return NULL;
} // est_ul_cls()

void est_ul_clb(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity, double err, double *pvalue)
{
	using namespace RooStats;
	ModelConfig *bModel = dynamic_cast<ModelConfig*> (wspace->obj("bConfig"));
	ModelConfig *sbModel = dynamic_cast<ModelConfig*> (wspace->obj("splusbConfig"));
	ProfileLikelihoodTestStat testStat(*wspace->pdf("total_pdf"));
	ToyMCSampler *mcSampler = new ToyMCSampler(testStat,1000);
	FrequentistCalculator frequCalc(*data,*sbModel,*bModel,mcSampler); // null = bModel interpreted as signal, alt = s+b interpreted as bkg
	HypoTestResult *result;
	
	// congigure
	mcSampler->SetNEventsPerToy(1);
	measure_params(wspace, data, channels, verbosity);
	
	// set the signal strengths to constant (=0)
	wspace->var("mu_s")->setConstant(kTRUE);
	wspace->var("mu_d")->setConstant(kTRUE);
	
	result = frequCalc.GetHypoTest();
	result->SetBackgroundAsAlt(kTRUE);
	*pvalue = result->CLsplusb();
	
	delete mcSampler;
	delete result;
} // est_ul_clb()

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
