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
#include <RooStats/HybridCalculator.h>
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

RooWorkspace *build_model_nchannel(map<bmm_param,measurement_t> *bsmm, map<bmm_param,measurement_t> *bdmm, bool no_errors, int verbosity, bool compute_bd_ul, bool create_gammas, bool fixed_bkg)
{
	RooStats::ModelConfig *splusbModel = NULL;
	RooStats::ModelConfig *bModel = NULL;
	RooWorkspace *wspace = new RooWorkspace("wspace");
	RooProdPdf *totalPdf = NULL;
	RooProdPdf *priorPdf = NULL;
	RooArgList priorList;
	RooArgList poissonList;
	RooArgSet observables;
	RooArgSet nuisanceParams;
	set<int> channels;
	set<int>::const_iterator chan;
	measurement_t m;
	
	add_channels(bsmm,&channels);
	add_channels(bdmm,&channels);
	
	// make sure we cover the entire physical range
	wspace->factory("mu_s[1,0,20]");	// initialize to standard model
	wspace->factory("mu_d[1,0,200]");	// initialize to standard model
	
	// global correlated variables of all channels
	m = f_ratio();
	if (m.getErr() > 0 && !no_errors) {
		wspace->factory(Form("fratio0[%f]",m.getVal()));
		wspace->factory(Form("fratioErr[%f]",m.getErr()));
		wspace->factory(Form("fratio[%f,0,10]",m.getVal()));
		wspace->factory("Gaussian::fratio_Gauss(fratio,fratio0,fratioErr)");
	} else {
		wspace->factory(Form("fratio[%f,0,1]",m.getVal()));
		wspace->var("fratio")->setConstant(kTRUE);
	}
	
	// Create channel specific variables
	for (chan = channels.begin(); chan != channels.end(); ++chan) {
		
		// Observables
		wspace->factory(Form("NsObs_%d[0,1000]",*chan)); // Observed Evts in Bs window
		wspace->factory(Form("NdObs_%d[0,1000]",*chan)); // Observed Evts in Bd window
		wspace->factory(Form("NbObs_%d[0,1000]",*chan)); // Observed Background
		observables.add(*wspace->var(Form("NsObs_%d",*chan)));
		observables.add(*wspace->var(Form("NdObs_%d",*chan)));
		if (fixed_bkg)
			wspace->var(Form("NbObs_%d",*chan))->setConstant(kTRUE);
		else
			observables.add(*wspace->var(Form("NbObs_%d",*chan)));
		
		// nuisance parameter
		wspace->factory(Form("nu_b_%d[1,0,1000]",*chan)); // background strength
		
		///////////////////////////////////
		// build the background ratio Bs //
		///////////////////////////////////
		if ( ((*bsmm)[make_pair(kTau_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bsmm)[make_pair(kTau_bmm, *chan)];
			wspace->factory(Form("TauS0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("TauSErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("TauS_%d[%f,%f,%f]",*chan,m.getVal(),0.0,1.0));
			wspace->factory(Form("Gaussian::TauS_Gauss_%d(TauS_%d,TauS0_%d,TauSErr_%d)",*chan,*chan,*chan,*chan));
		}
		else {
			wspace->factory(Form("TauS_%d[%f,0,10]",*chan, ((*bsmm)[make_pair(kTau_bmm, *chan)]).getVal()));
			wspace->var(Form("TauS_%d",*chan))->setConstant(kTRUE);
		}
		
		////////////////////////////////////
		// build the background ration Bd //
		////////////////////////////////////
		if ( ((*bdmm)[make_pair(kTau_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bdmm)[make_pair(kTau_bmm, *chan)];
			wspace->factory(Form("TauD0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("TauDErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("TauD_%d[%f,%f,%f]",*chan,m.getVal(),0.0,1.0));
			wspace->factory(Form("Gaussian::TauD_Gauss_%d(TauD_%d,TauD0_%d,TauDErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("TauD_%d[%f,0,10]",*chan,((*bdmm)[make_pair(kTau_bmm, *chan)]).getVal()));
			wspace->var(Form("TauD_%d",*chan))->setConstant(kTRUE);
		}
		
		/////////////////////////
		// Construction of NuS //
		/////////////////////////
		if ( ((*bsmm)[make_pair(kExpUncor_bmm, *chan)]).getErr() > 0 && !no_errors) {
			m = (*bsmm)[make_pair(kExpUncor_bmm, *chan)];
			wspace->factory(Form("NuSUncor0_%d[%f]", *chan, m.getVal())); // fixed mean variable
			wspace->factory(Form("NuSUncorErr_%d[%f]", *chan, m.getErr())); // fixed error variable
			wspace->factory(Form("NuSUncor_%d[%f,%f,%f]", *chan, m.getVal(), 0.0, m.getVal() + 100*m.getErr()));
			wspace->factory(Form("Gaussian::NuSUncor_Gauss_%d(NuSUncor_%d,NuSUncor0_%d,NuSUncorErr_%d)",*chan,*chan,*chan,*chan)); // error gaussian
		} else { // no error associated to this variable
			wspace->factory(Form("NuSUncor_%d[%f,0,100]", *chan, ((*bsmm)[make_pair(kExpUncor_bmm, *chan)]).getVal()));
			wspace->var(Form("NuSUncor_%d", *chan))->setConstant(kTRUE);
		}
		
		/////////////////////////
		// Construction of NuD //
		/////////////////////////
		if ( ((*bdmm)[make_pair(kExp_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bdmm)[make_pair(kExp_bmm, *chan)];
			wspace->factory(Form("NuD0_%d[%f]", *chan, m.getVal())); // fixed mean variable
			wspace->factory(Form("NuDErr_%d[%f]", *chan, m.getErr())); // fixed error variable
			wspace->factory(Form("NuD_%d[%f,%f,%f]", *chan, m.getVal(), 0.0, m.getVal() + 100*m.getErr()));
			wspace->factory(Form("Gaussian::NuD_Gauss_%d(NuD_%d,NuD0_%d,NuDErr_%d)",*chan,*chan,*chan,*chan)); // error gaussian
		} else {
			// no error assigned, just make a constant
			wspace->factory(Form("NuD_%d[%f,0,100]", *chan, ((*bdmm)[make_pair(kExp_bmm, *chan)]).getVal()));
			wspace->var(Form("NuD_%d",*chan))->setConstant(kTRUE);
		}
		
		/////////////////////////
		// Construction of Pss //
		/////////////////////////
		if ( ((*bsmm)[make_pair(kProb_swind_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bsmm)[make_pair(kProb_swind_bmm, *chan)];
			wspace->factory(Form("Pss0_%d[%f]", *chan, m.getVal()));
			wspace->factory(Form("PssErr_%d[%f]", *chan, m.getErr()));
			wspace->factory(Form("Pss_%d[%f,%f,%f]", *chan, m.getVal(), 0.0, 1.0));
			wspace->factory(Form("Gaussian::Pss_Gauss_%d(Pss_%d,Pss0_%d,PssErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			// no error assigned, just make a constant.
			wspace->factory(Form("Pss_%d[%f,0,1]", *chan, ((*bsmm)[make_pair(kProb_swind_bmm, *chan)]).getVal()));
			wspace->var(Form("Pss_%d",*chan))->setConstant(kTRUE);
		}
		
		/////////////////////////
		// Construction of Psd //
		/////////////////////////
		if ( ((*bdmm)[make_pair(kProb_swind_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bdmm)[make_pair(kProb_swind_bmm, *chan)];
			wspace->factory(Form("Psd0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PsdErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("Psd_%d[%f,%f,%f]", *chan, m.getVal(), 0.0, 1.0));
			wspace->factory(Form("Gaussian::Psd_Gauss_%d(Psd_%d,Psd0_%d,PsdErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			// no error assigned, just make a constant.
			wspace->factory(Form("Psd_%d[%f,0,1]", *chan, ((*bdmm)[make_pair(kProb_swind_bmm, *chan)]).getVal()));
			wspace->var(Form("Psd_%d",*chan))->setConstant(kTRUE);
		}
		
		/////////////////////////
		// Construction of Pds //
		/////////////////////////
		if ( ((*bsmm)[make_pair(kProb_dwind_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bsmm)[make_pair(kProb_dwind_bmm, *chan)];
			wspace->factory(Form("Pds0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PdsErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("Pds_%d[%f,%f,%f]",*chan,m.getVal(),0.0,1.0));
			wspace->factory(Form("Gaussian::Pds_Gauss_%d(Pds_%d,Pds0_%d,PdsErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			// no error assigned, just make a constant.
			wspace->factory(Form("Pds_%d[%f,0,1]", *chan, ((*bsmm)[make_pair(kProb_dwind_bmm, *chan)]).getVal()));
			wspace->var(Form("Pds_%d",*chan))->setConstant(kTRUE);
		}
		
		/////////////////////////
		// Construction of Pdd //
		/////////////////////////
		if ( ((*bdmm)[make_pair(kProb_dwind_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bdmm)[make_pair(kProb_dwind_bmm, *chan)];
			wspace->factory(Form("Pdd0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PddErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("Pdd_%d[%f,%f,%f]",*chan,m.getVal(),0.0,1.0));
			wspace->factory(Form("Gaussian::Pdd_Gauss_%d(Pdd_%d,Pdd0_%d,PddErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("Pdd_%d[%f,0,1]", *chan, ((*bdmm)[make_pair(kProb_dwind_bmm, *chan)]).getVal()));
			wspace->var(Form("Pdd_%d",*chan))->setConstant(kTRUE);
		}
		
		//////////////////////////////////////
		// Construction of rare backgrounds //
		//////////////////////////////////////
		if ( ((*bsmm)[make_pair(kPeakBkgOff_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bsmm)[make_pair(kPeakBkgOff_bmm, *chan)];
			wspace->factory(Form("PeakBkgSB0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PeakBkgSBErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("PeakBkgSB_%d[%f,%f,%f]",*chan,m.getVal(),0.0,m.getVal() + 10*m.getErr()));
			wspace->factory(Form("Gaussian::PeakBkgSB_Gauss_%d(PeakBkgSB_%d,PeakBkgSB0_%d,PeakBkgSBErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("PeakBkgSB_%d[%f,0,50]", *chan, ((*bsmm)[make_pair(kPeakBkgOff_bmm, *chan)]).getVal()));
			wspace->var(Form("PeakBkgSB_%d",*chan))->setConstant(kTRUE);
		}
		
		if ( ((*bsmm)[make_pair(kPeakBkgOn_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bsmm)[make_pair(kPeakBkgOn_bmm, *chan)];
			wspace->factory(Form("PeakBkgBs0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PeakBkgBsErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("PeakBkgBs_%d[%f,%f,%f]",*chan,m.getVal(),0.0,m.getVal() + 10*m.getErr()));
			wspace->factory(Form("Gaussian::PeakBkgBs_Gauss_%d(PeakBkgBs_%d,PeakBkgBs0_%d,PeakBkgBsErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("PeakBkgBs_%d[%f,0,50]",*chan,((*bsmm)[make_pair(kPeakBkgOn_bmm, *chan)]).getVal()));
			wspace->var(Form("PeakBkgBs_%d",*chan))->setConstant(kTRUE);
		}
		
		if ( ((*bdmm)[make_pair(kPeakBkgOn_bmm, *chan)]).getErr() > 0 && !no_errors ) {
			m = (*bdmm)[make_pair(kPeakBkgOn_bmm, *chan)];
			wspace->factory(Form("PeakBkgBd0_%d[%f]",*chan,m.getVal()));
			wspace->factory(Form("PeakBkgBdErr_%d[%f]",*chan,m.getErr()));
			wspace->factory(Form("PeakBkgBd_%d[%f,%f,%f]",*chan, m.getVal(), 0.0 ,m.getVal() + 10*m.getErr()));
			wspace->factory(Form("Gaussian::PeakBkgBd_Gauss_%d(PeakBkgBd_%d,PeakBkgBd0_%d,PeakBkgBdErr_%d)",*chan,*chan,*chan,*chan));
		} else {
			wspace->factory(Form("PeakBkgBd_%d[%f,0,50]",*chan,((*bdmm)[make_pair(kPeakBkgOn_bmm, *chan)]).getVal()));
			wspace->var(Form("PeakBkgBd_%d",*chan))->setConstant(kTRUE);
		}
		
		//////////////////////////////
		// Construction of formulas //
		//////////////////////////////
		wspace->factory(Form("FormulaVar::NuS_%d(\"fratio*NuSUncor_%d\",{fratio,NuSUncor_%d})",*chan,*chan,*chan));
		wspace->factory(Form("FormulaVar::bkg_mean_%d(\"nu_b_%d + PeakBkgSB_%d\",{nu_b_%d,PeakBkgSB_%d})",*chan,*chan,*chan,*chan,*chan));
		wspace->factory(Form("FormulaVar::bs_mean_%d(\"TauS_%d*nu_b_%d + PeakBkgBs_%d + Pss_%d*NuS_%d*mu_s + Psd_%d*NuD_%d*mu_d\",{TauS_%d,nu_b_%d,PeakBkgBs_%d,Pss_%d,NuS_%d,mu_s,Psd_%d,NuD_%d,mu_d})",*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan));
		wspace->factory(Form("FormulaVar::bd_mean_%d(\"TauD_%d*nu_b_%d + PeakBkgBd_%d + Pds_%d*NuS_%d*mu_s + Pdd_%d*NuD_%d*mu_d\",{TauD_%d,nu_b_%d,PeakBkgBd_%d,Pds_%d,NuS_%d,mu_s,Pdd_%d,NuD_%d,mu_d})",*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan,*chan));
		
		//////////////////////
		// Main Poissonians //
		//////////////////////		
		wspace->factory(Form("Poisson::bs_window_%d(NsObs_%d,bs_mean_%d)",*chan,*chan,*chan));
		wspace->factory(Form("Poisson::bd_window_%d(NdObs_%d,bd_mean_%d)",*chan,*chan,*chan));
		wspace->factory(Form("Poisson::bkg_window_%d(NbObs_%d,bkg_mean_%d)",*chan,*chan,*chan));
		
		// add the poissonians of this channel
		poissonList.add(*wspace->pdf(Form("bs_window_%d",*chan)));
		poissonList.add(*wspace->pdf(Form("bd_window_%d",*chan)));
		if (!fixed_bkg) poissonList.add(*wspace->pdf(Form("bkg_window_%d",*chan)));
	}
	
	// define the sets
	wspace->defineSet("obs", observables);
	nuisanceParams.addClone(wspace->allVars());
	// FIXME: Think if other signal strength should be handled as nuisance parameter or not!
	if (compute_bd_ul)
		wspace->defineSet("poi", "mu_d");
	else
		wspace->defineSet("poi", "mu_s");
	nuisanceParams.remove(*wspace->set("obs"), kTRUE, kTRUE);
	nuisanceParams.remove(*wspace->set("poi"), kTRUE, kTRUE);
	RooStats::RemoveConstantParameters(&nuisanceParams);
	wspace->defineSet("nui", nuisanceParams);
	
	wspace->factory("beta[1]");
	wspace->factory("mu[0]");
	
	// add gamma distributions if requested!
	if (create_gammas) {
		for (chan = channels.begin(); chan != channels.end(); ++chan) {
			// background measurement constraint
			wspace->factory(Form("Gamma::bkg_prior_%d(bkg_mean_%d,gamma_%d[1],beta,mu)",*chan,*chan,*chan));
		}
	}
	
	// build the prior
	priorList.add(wspace->allPdfs());
	priorList.remove(poissonList, kTRUE, kTRUE);
	if (priorList.getSize() > 0) {
		priorPdf = new RooProdPdf("prior_pdf","Total Prior PDF",priorList);
		poissonList.add(*priorPdf);
	}
	
	totalPdf = new RooProdPdf("total_pdf","Total Model PDF",poissonList);
	wspace->import(*totalPdf, ((verbosity > 0) ? RooCmdArg::none() : RooFit::Silence(kTRUE)) );	
	
	// uniform prior in case of bayesian code
	wspace->factory(Form("Uniform::prior_mu(mu_%c)",(compute_bd_ul ? 'd' : 's')));

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
	delete priorPdf;
	delete totalPdf;
	
	return wspace;
} // build_model_nchannel()

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
	RooStats::ModelConfig *splusbConfig = dynamic_cast<RooStats::ModelConfig*> (wspace->obj("splusbConfig"));
	RooStats::ModelConfig *bConfig = dynamic_cast<RooStats::ModelConfig*> (wspace->obj("bConfig"));
	wspace->allVars() = *data->get(0);
	
	((RooRealVar*)wspace->set("poi")->first())->setVal(0.0);
	((RooRealVar*)wspace->set("poi")->first())->setConstant(kTRUE);
	wspace->pdf("total_pdf")->fitTo(*data, ((verbosity > 0) ? RooCmdArg::none() : RooFit::PrintLevel(-1)));
	((RooRealVar*)wspace->set("poi")->first())->setConstant(kFALSE);
	bConfig->SetSnapshot(*wspace->set("poi"));
	
	// do a likelihood fit to the data to get the real values...
	wspace->pdf("total_pdf")->fitTo(*data, ((verbosity > 0) ? RooCmdArg::none() : RooFit::PrintLevel(-1)));
	splusbConfig->SetSnapshot(*wspace->set("poi"));
} // measure_params()

RooStats::ConfInterval *est_ul_fc(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, double cLevel, int verbosity, double err, double *ulLimit, double *loLimit, pair<double,double> *rg, uint32_t *inBins, double *cpuUsed, uint32_t nbrProof, int nToys)
{
	using namespace RooStats;
	ModelConfig *splusbConfig = dynamic_cast<ModelConfig*>(wspace->obj("splusbConfig"));
	FeldmanCousins fc(*data,*splusbConfig);
	ToyMCSampler *toySampler;
	PointSetInterval *psInterval = NULL;
	TStopwatch swatch;
	RooDataSet *custom_points;
	Double_t v;
	Int_t j;
	ProofConfig *pc = NULL;
	string proofString = Form("workers=%u",nbrProof);
	
	swatch.Start(kTRUE);
	
	// configure Feldman Cousins
	fc.CreateConfBelt(true);
	fc.SetTestSize(1.-cLevel);
	fc.SetNBins(*inBins);
	fc.FluctuateNumDataEntries(false);
	toySampler = dynamic_cast<ToyMCSampler*> (fc.GetTestStatSampler());
	toySampler->SetNToys(nToys);
	if (nbrProof > 1) {
		uint32_t nPackages = (((nToys + 999)/1000+(nbrProof-1))/nbrProof)*nbrProof;
		pc = new ProofConfig(*wspace, nPackages, proofString.c_str(), kFALSE);
		toySampler->SetProofConfig(pc);
	}
	if (rg) {
		RooArgSet poi;
		poi.addClone(*wspace->set("poi"));
		custom_points = new RooDataSet("","",poi);
		
		for (j = 0; j < (Int_t)*inBins; j++) {
			v = rg->first + (rg->second - rg->first) * j / (*inBins - 1);
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
	
	delete pc;
	
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

// try with hybrid approach and ratioofprofiled
RooStats::ConfInterval *est_ul_hybrid(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, double cLevel, int verbosity, double err, double *ulLimit, pair<double,double> *rg, uint32_t* inBins, double *cpuUsed, uint32_t nbrProof, int nToys, bool bdmm)
 {
	using namespace RooStats;
	ModelConfig *bModel = dynamic_cast<ModelConfig*> (wspace->obj("bConfig"));
	ModelConfig *sbModel = dynamic_cast<ModelConfig*> (wspace->obj("splusbConfig"));
	((RooRealVar*)wspace->set("poi")->first())->setVal(0); // for background
	RatioOfProfiledLikelihoodsTestStat testStat(*sbModel->GetPdf(),*bModel->GetPdf(),wspace->set("poi"));
	testStat.SetSubtractMLE(false);
	ToyMCSampler *mcSampler = new ToyMCSampler(testStat,nToys);
	HybridCalculator hybCalc(*data,*bModel,*sbModel,mcSampler);
	HypoTestInverter *hypoInv = NULL;
	HypoTestInverterResult *result = NULL;
	double limitErr = 0.1;
	double obs;
	TStopwatch swatch;
	bool ok;
	RooArgSet mu;
	uint32_t nBins = inBins ? *inBins : 10; // default 10 bins
	ProofConfig *pc =  NULL;
	string proofString = Form("workers=%u",nbrProof);
	double beta = 0;
	
	if (nbrProof > 1) {
		uint32_t nPackages = (((nToys + 999)/1000+(nbrProof-1))/nbrProof)*nbrProof;
		pc = new ProofConfig(*wspace, nPackages, proofString.c_str(), kFALSE);
		mcSampler->SetProofConfig(pc);
	}
	
	swatch.Start(kTRUE);

	wspace->var("gamma_0")->setVal(((RooRealVar&)((*data->get(0))["NbObs_0"])).getVal()+1);
	wspace->var("gamma_1")->setVal(((RooRealVar&)((*data->get(0))["NbObs_1"])).getVal()+1);
	measure_params(wspace, data, channels, verbosity);
	sbModel->LoadSnapshot();
	
	if (bdmm) {
		for (set<int>::const_iterator it = channels->begin(); it != channels->end(); ++it)
			beta += wspace->var(Form("Pss_%d",*it))->getVal() * wspace->var(Form("NuS_%d",*it))->getVal();
		
		obs = beta * wspace->var("mu_s")->getVal();
		
		beta = 1./beta; // beta is actually the inverse thereof in RooFit
		wspace->factory(Form("Gamma::mu_prior_gamma(mu_s,gamma_b[%f],beta_b[%f],mu)",1.0 + obs,beta));
	} else {
		for (set<int>::const_iterator it = channels->begin(); it != channels->end(); ++it)
			beta += wspace->var(Form("Pdd_%d",*it))->getVal() * wspace->var(Form("NuD_%d",*it))->getVal();
		
		obs = beta * wspace->var("mu_d")->getVal();
		beta = 1./beta; // beta is actually the inverse thereof in RooFit
		wspace->factory(Form("Gamma::mu_prior_gamma(mu_d,gamma_b[%f],beta_b[%f],mu)",1.0 + obs, beta));
	}
	wspace->factory("PROD::nui_sampling(prior_pdf,mu_prior_gamma)");
	hybCalc.ForcePriorNuisanceAlt(*wspace->pdf("nui_sampling"));
	hybCalc.ForcePriorNuisanceNull(*wspace->pdf("nui_sampling"));
	mcSampler->SetNEventsPerToy(1);
	hypoInv = new HypoTestInverter(hybCalc, (RooRealVar*)sbModel->GetParametersOfInterest()->first(), 1.0 - cLevel);
	hypoInv->SetAutoScan();
	hypoInv->UseCLs(true);
	hypoInv->SetTestSize(1.0 - cLevel);
	hypoInv->SetVerbose(verbosity);
	
	if (rg)	ok = hypoInv->RunFixedScan(nBins, rg->first, rg->second);
	else	ok = hypoInv->RunLimit(*ulLimit,limitErr);

	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	result = hypoInv->GetInterval();
	*ulLimit = result->UpperLimit();

	delete hypoInv;
	delete pc;
	return result;
} // est_ul_hybrid()

RooStats::ConfInterval *est_ul_cls(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, double cLevel, int verbosity, double err, double *ulLimit, pair<double,double> *rg, uint32_t *npts, double *cpuUsed, uint32_t nbrProof, int nToys)
{
	using namespace RooStats;
	ModelConfig *bModel = dynamic_cast<ModelConfig*> (wspace->obj("bConfig"));
	ModelConfig *sbModel = dynamic_cast<ModelConfig*> (wspace->obj("splusbConfig"));
	ProfileLikelihoodTestStat testStat(*wspace->pdf("total_pdf"));
	testStat.SetOneSided(true);
	ToyMCSampler *mcSampler = new ToyMCSampler(testStat,nToys);
	FrequentistCalculator frequCalc(*data,*bModel,*sbModel,mcSampler); // Note null = sb, alt = b
	HypoTestInverter *hypoInv = NULL;
	HypoTestInverterResult *result = NULL;
	double limitErr = 0.1;
	TStopwatch swatch;
	bool ok;
	int nbr = npts ? *npts : 10; // defaults to 10 points if not otherwise specified.
	ProofConfig *pc = NULL;
	string proofString = Form("workers=%u",nbrProof);
	
	if (nbrProof > 1) {
		uint32_t nPackages = (((nToys + 999)/1000+(nbrProof-1))/nbrProof)*nbrProof;
		pc = new ProofConfig(*wspace, nPackages, proofString.c_str(), kFALSE);
		mcSampler->SetProofConfig(pc);
	}
	
	swatch.Start(kTRUE);
	
	measure_params(wspace, data, channels, verbosity);
	sbModel->LoadSnapshot();
	
	mcSampler->SetNEventsPerToy(1);
	hypoInv = new HypoTestInverter(frequCalc, (RooRealVar*)sbModel->GetParametersOfInterest()->first(), 1.0 - cLevel);
	hypoInv->SetAutoScan();
	hypoInv->UseCLs();
	hypoInv->SetTestSize(1.0 - cLevel);
	hypoInv->SetVerbose(verbosity);
	
	if (rg)	ok = hypoInv->RunFixedScan(nbr, rg->first, rg->second);
	else	ok = hypoInv->RunLimit(*ulLimit,limitErr);
	
	swatch.Stop();
	if (cpuUsed) *cpuUsed = swatch.CpuTime();
	
	result = hypoInv->GetInterval();
	*ulLimit = result->UpperLimit();
	
	delete pc;
	delete hypoInv;
	return result;
} // est_ul_cls()

RooStats::HypoTestResult *est_ul_clb_hybrid(RooWorkspace *wspace, RooDataSet *data, std::set<int> *channels, int verbosity, double err, double *pvalue, uint32_t nbrProof, int nToys, bool bdmm)
{
	using namespace RooStats;
	ModelConfig *bModel = dynamic_cast<ModelConfig*> (wspace->obj("bConfig"));
	ModelConfig *sbModel = dynamic_cast<ModelConfig*> (wspace->obj("splusbConfig"));
	ProfileLikelihoodTestStat testStat(*wspace->pdf("total_pdf"));
	ToyMCSampler *mcSampler = new ToyMCSampler(testStat,nToys);
	HybridCalculator hybCalc(*data,*sbModel,*bModel,mcSampler); // null = bModel interpreted as signal, alt = s+b interpreted as bkg
	HypoTestResult *result;
	RooArgSet mu;
	ProofConfig *pc = NULL;
	string proofString = Form("workers=%u",nbrProof);
	double beta = 0;
	double obs;
	
	if (nbrProof > 1) {
		uint32_t nPackages = (((nToys + 999)/1000+(nbrProof-1))/nbrProof)*nbrProof;
		pc = new ProofConfig(*wspace,nPackages,proofString.c_str(),kFALSE);
		mcSampler->SetProofConfig(pc);
	}
	
	// set background prior to measurement
	wspace->var("gamma_0")->setVal(((RooRealVar&)((*data->get(0))["NbObs_0"])).getVal()+1);
	wspace->var("gamma_1")->setVal(((RooRealVar&)((*data->get(0))["NbObs_1"])).getVal()+1);
	measure_params(wspace, data, channels, verbosity);
	bModel->LoadSnapshot(); // FIXME: adjust measure_params s.t. nuisance mu is also in background model
	
	// nuisance signal
	if (bdmm) {
		for (set<int>::const_iterator it = channels->begin(); it != channels->end(); ++it)
			beta += wspace->var(Form("Pss_%d",*it))->getVal() * wspace->var(Form("NuS_%d",*it))->getVal();
		
		obs = beta * wspace->var("mu_s")->getVal();
		
		beta = 1./beta; // beta is actually the inverse thereof in RooFit
		wspace->factory(Form("Gamma::mu_prior_gamma(mu_s,gamma_b[%f],beta_b[%f],mu)",1.0 + obs,beta));
	} else {
		for (set<int>::const_iterator it = channels->begin(); it != channels->end(); ++it)
			beta += wspace->var(Form("Pdd_%d",*it))->getVal() * wspace->var(Form("NuD_%d",*it))->getVal();
		
		obs = beta * wspace->var("mu_d")->getVal();
		beta = 1./beta; // beta is actually the inverse thereof in RooFit
		wspace->factory(Form("Gamma::mu_prior_gamma(mu_d,gamma_b[%f],beta_b[%f],mu)",1.0 + obs, beta));
	}
	wspace->factory("PROD::nui_sampling(prior_pdf,mu_prior_gamma)");
	hybCalc.ForcePriorNuisanceAlt(*wspace->pdf("nui_sampling"));
	hybCalc.ForcePriorNuisanceNull(*wspace->pdf("nui_sampling"));
	mcSampler->SetNEventsPerToy(1);
	result = hybCalc.GetHypoTest();
	result->SetBackgroundAsAlt(kTRUE);
	
	*pvalue = result->CLsplusb();
	
	delete pc;
	delete mcSampler;
	
	return result;
} // est_ul_clb()

RooStats::HypoTestResult *est_ul_clb(RooWorkspace *wspace, RooDataSet *data, set<int> *channels, int verbosity, double err, double *pvalue, uint32_t nbrProof, int nToys)
{
	using namespace RooStats;
	ModelConfig *bModel = dynamic_cast<ModelConfig*> (wspace->obj("bConfig"));
	ModelConfig *sbModel = dynamic_cast<ModelConfig*> (wspace->obj("splusbConfig"));
	ProfileLikelihoodTestStat testStat(*wspace->pdf("total_pdf"));
	ToyMCSampler *mcSampler = new ToyMCSampler(testStat,nToys);
	FrequentistCalculator frequCalc(*data,*sbModel,*bModel,mcSampler); // null = bModel interpreted as signal, alt = s+b interpreted as bkg
	HypoTestResult *result = NULL;
	ProofConfig *pc = NULL;
	string proofString = Form("workers=%u",nbrProof);
	RooArgList mu;
	
	if (nbrProof > 1) {
		uint32_t nPackages = (((nToys + 999)/1000+(nbrProof-1))/nbrProof)*nbrProof;
		pc = new ProofConfig(*wspace,nPackages,proofString.c_str(),kFALSE);
		mcSampler->SetProofConfig(pc);
	}
	
	// set background prior to measurement
	measure_params(wspace, data, channels, verbosity);
	bModel->LoadSnapshot();
	mcSampler->SetNEventsPerToy(1);
	
	result = frequCalc.GetHypoTest();
	result->SetBackgroundAsAlt(kTRUE);
	*pvalue = result->CLsplusb();
	
	delete mcSampler;
	delete pc;
	
	return result;
} // est_ul_clb()

void compute_vars(map<bmm_param,measurement_t> *bmm, bool bstomumu)
{
	set<int> channels;
	set<int>::const_iterator it;
	measurement_t tot_bplus,nu;
	measurement_t c = bstomumu ? c_s_theory() : c_d_theory();
	measurement_t f_rat = f_ratio();
	
	add_channels(bmm, &channels);
	
	for (it = channels.begin(); it != channels.end(); ++it) {
		
		if (bmm->count(make_pair(kTot_bplus, *it)) > 0)
			tot_bplus = (*bmm)[make_pair(kTot_bplus, *it)];
		else
			tot_bplus = (*bmm)[make_pair(kObs_bplus, *it)] / compute_efftot_bplus(bmm,*it);
		
		nu = c * compute_efftot_bmm(bmm,*it) * tot_bplus;
		(*bmm)[make_pair(kExpUncor_bmm, *it)] = nu; // uncorrelated error from measurements
		
		if (bstomumu) nu = nu * f_rat;
		(*bmm)[make_pair(kExp_bmm, *it)] = nu; // for convenience compute expectation as well.
	}
} // compute_vars()
