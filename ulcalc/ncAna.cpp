/*
 *  ncAna.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 28.06.12.
 *
 */

// my headers
#include "ncAna.h"
#include "ncTree.h"
#include "ncFormula.h"

#include "external_constants.h"
#include "../rootutils/NCRootUtils.h"

#include <iostream>

// ROOT headers
#include <TMath.h>
#include <TTree.h>

using std::cout;
using std::endl;


ncAna::ncAna() : dataFile(NULL), mcFile(NULL), fNbrMassBins(50), fMassRange(4.9,5.9), fNbrSidebandBins(50), fLowSidebandRange(5.05,5.15), fHighSidebandRange(5.4,5.5), fPlotDir("plots")
{
	typedef ncCut::ncCutRange range_t;
	
	// initialize the files
	initFiles();
	
	// FIXME: enable the other ones again
	fSidebandCuts.push_back(ncCut("pt_mu1",		range_t(0.0,20.0),	range_t(4.5,1e30)));
	fSidebandCuts.push_back(ncCut("pt_mu2",		range_t(0.0,10.0),	range_t(4.0,1e30)));
	fSidebandCuts.push_back(ncCut("pt",			range_t(0.0,40.0),	range_t(6.5,1e30)));
	fSidebandCuts.push_back(ncCut("ip",			range_t(0.0,0.02),	range_t(-1e30,0.008)));
	fSidebandCuts.push_back(ncCut("ip/ipe",		range_t(0.0,5.0),	range_t(-1e30,2.0),"sig_ip"));
	fSidebandCuts.push_back(ncCut("alpha",		range_t(0.0,0.2),	range_t(-1e30,0.05)));
	fSidebandCuts.push_back(ncCut("chi2/Ndof",	range_t(0.0,6.0),	range_t(-1e30,2.2),"norm_chi2"));
	fSidebandCuts.push_back(ncCut("d3/d3e",		range_t(0.0,100.0),	range_t(13.0,1e30),"sig_d3"));
	fSidebandCuts.push_back(ncCut("d3",			range_t(0.0,1.0),	range_t(-1e30,1.5)));
	fSidebandCuts.push_back(ncCut("iso_mor12",	range_t(0.0,1.0),	range_t(0.8,1e30)));
	fSidebandCuts.push_back(ncCut("doca0",		range_t(0.0,0.2),	range_t(0.015,1e30)));
	fSidebandCuts.push_back(ncCut("ntrk",		range_t(0.0,20.0),	range_t(-1e30,2.0)));
} // ncAna()

ncAna::~ncAna()
{
	if (dataFile) delete dataFile;
	if (mcFile) delete mcFile;
} // ~ncAna

void ncAna::initFiles()
{
	dataFile = new TFile("/Users/cn/CMSData/Reduced/data-2011.root");
	mcFile = new TFile("/Users/cn/CMSData/Reduced/production-mix-general.root");
} // initFiles()

TH1D *ncAna::getDefaultMassHistogram(const char *name)
{
	return new TH1D(name,"",fNbrMassBins,fMassRange.first, fMassRange.second);
} // getDefaultMassHistogram()

void ncAna::showSidebandSubtraction(bool massConstraint)
{
	TTree *dataTree = (TTree*)dataFile->Get("T");
	TTree *mcTree = (TTree*)mcFile->Get("T");
	TH1D *hsg,*hlo,*hhi,*hmc;
	TH1D *mHisto;
	TCut cut;
	size_t j,k;
	double mu,res,err;
	fit_t fit;
	TCanvas *c = new TCanvas;
	std::string mHistoName;
	
	TCut signalCut;
	TCut lowSBCut;
	TCut highSBCut;
	
	measurement_t sig;
	measurement_t bkg;
	
	init_gauss_p2(&fit, 5.05, 5.45);
	
	// all cuts stored in the sideband stuff
	for (j = 0; j < fSidebandCuts.size(); j++) {
		
		// plotte variable j mit allen cuts applied ausser cut j
		cout << "ncAna: Processing variable '" << fSidebandCuts[j].getName() << "'" << endl;
		
		// create the various histograms
		mHistoName = std::string(Form("mass_%s",fSidebandCuts[j].getName()));
		mHisto = getDefaultMassHistogram(mHistoName.c_str());
		hsg = new TH1D("hsg","Signal Region",fNbrSidebandBins,fSidebandCuts[j].getDomain().first,fSidebandCuts[j].getDomain().second);
		hlo = new TH1D("hlo","Low Sideband",fNbrSidebandBins,fSidebandCuts[j].getDomain().first,fSidebandCuts[j].getDomain().second);
		hhi = new TH1D("hhi","High Sideband",fNbrSidebandBins,fSidebandCuts[j].getDomain().first,fSidebandCuts[j].getDomain().second);
		hmc = new TH1D("hmc","MC Histogram",fNbrSidebandBins,fSidebandCuts[j].getDomain().first,fSidebandCuts[j].getDomain().second);
		
		hmc->SetFillColor(kRed);
		hmc->SetLineColor(kRed);
		hmc->SetFillStyle(3001);
		
		hsg->GetXaxis()->SetTitle(fSidebandCuts[j].getFormula());
		hlo->GetXaxis()->SetTitle(fSidebandCuts[j].getFormula());
		hhi->GetXaxis()->SetTitle(fSidebandCuts[j].getFormula());
		hmc->GetXaxis()->SetTitle(fSidebandCuts[j].getFormula());
		
		// enable sumw2
		hsg->Sumw2();
		hlo->Sumw2();
		hhi->Sumw2();
		
		// cut konstruieren
		cut = cutNormCand() && cutMuon() && cutTrigger(false);
		for (k = 0; k < fSidebandCuts.size(); k++) {
			if (j == k) continue;
			cut = cut && TCut(Form("%f < %s && %s < %f",fSidebandCuts[k].getCut().first,fSidebandCuts[k].getFormula(),fSidebandCuts[k].getFormula(),fSidebandCuts[k].getCut().second));
		}
		
		// mass control plot zeichnen & speichern
		cout << "Applying cut '" << cut.GetTitle() << "'" << endl;
		dataTree->Draw(Form("%s >> %s", (massConstraint ? "mass_c" : "mass"), mHistoName.c_str()),cut);
		adjust_parameter_gauss_p2(mHisto, fit.fit_fct);
		mHisto->Fit(fit.fit_fct,"R"); // fit using chi2 (very stable to get values)
		mHisto->Fit(fit.fit_fct,"LR"); // fit using binned likelihood
		// zeichnen
		mHisto->SetTitle(Form("Mass control plot for variable %s",fSidebandCuts[j].getName()));
		mHisto->Draw("E1");
		c->SaveAs(Form("%s/%s.eps",fPlotDir.c_str(),mHistoName.c_str()));
		
		mu = fit.fit_fct->GetParameter(fit.mean_param);
		res = fit.fit_fct->GetParameter(fit.sigma_param);
		
		// get the number of events
		sig.setVal(signal_events_gauss_p2(fit.fit_fct, mHisto->GetBinWidth(1), &err)); // 2 sigma
		sig.setErr(err);
		
		bkg.setVal(background_events_gauss_p2(fit.fit_fct, mHisto->GetBinWidth(1))); // 2 sigma
		bkg.setErr(TMath::Sqrt(bkg.getVal())); // FIXME: no correct error here
		
		// histogram einteilen
		signalCut = TCut(Form("%f < mass && mass < %f",mu - 2*res, mu + 2*res));
		lowSBCut = TCut(Form("%f < mass && mass < %f",fLowSidebandRange.first, fLowSidebandRange.second));
		highSBCut = TCut(Form("%f < mass && mass < %f",fHighSidebandRange.first, fHighSidebandRange.second));
		
		// variable zeichnen und speichern
		dataTree->Draw(Form("%s >> hsg",fSidebandCuts[j].getFormula()),cut && signalCut);
		hsg->Draw("E1");
		c->SaveAs(Form("%s/%s_s+b.eps",fPlotDir.c_str(),fSidebandCuts[j].getName()));
		// low sideband
		dataTree->Draw(Form("%s >> hlo",fSidebandCuts[j].getFormula()),cut && lowSBCut);
		hlo->Draw("E1");
		c->SaveAs(Form("%s/%s_lo.eps",fPlotDir.c_str(),fSidebandCuts[j].getName()));
		// high sideband
		dataTree->Draw(Form("%s >> hhi",fSidebandCuts[j].getFormula()), cut && highSBCut);
		hhi->Draw("E1");
		c->SaveAs(Form("%s/%s_hi.eps",fPlotDir.c_str(),fSidebandCuts[j].getName()));
		
		// add sidebands together
		hlo->Add(hhi);
		hsg->Add(hlo, -bkg.getVal()/hlo->Integral());
		// plot the mc truth signal
		mcTree->Draw(Form("%s >> hmc",fSidebandCuts[j].getFormula()),cut && cutNormTruth());
		c->SaveAs(Form("%s/%s_mc.eps",fPlotDir.c_str(),fSidebandCuts[j].getName()));
		// rescale the mc to the same integral as the signal
		hmc->Scale(hsg->Integral() / hmc->Integral());
		hmc->SetTitle("Sideband subtracted signal");
		hmc->Draw("B");
		hsg->Draw("E1same");
		c->SaveAs(Form("%s/%s_sig.eps",fPlotDir.c_str(),fSidebandCuts[j].getName()));
	}
	
	delete c;
} // showSidebandSubtraction()

// returns the systematic uncertainty of the cut
// FIXME: Intertwine with sideband subtraction
double ncAna::getSystematicsEffAna(TCut cut, bool massConstraint)
{
	using std::string;
	string str = string(cut.GetTitle());
	string::iterator it;
	double result = 0.0;
	ncTree<ncFormula> *tree;
	std::set<string> cuts;
	
	// parse the string & extract all the cuts
	tree = read_formula(string(cut.GetTitle()));
	cout << "Cut read: '" << tree->toString(false) << "'" << endl;
	cuts = get_cuts(tree);
	delete tree;
	
	// dump the cuts
	for (std::set<string>::const_iterator it = cuts.begin(); it != cuts.end(); ++it)
		cout << "Found the cut '" << *it << "'" << endl;
	
	// test all the cuts
	result = getSystematicsEffAna(cuts,massConstraint);
	
	return result;
} // getSystematicsEffAna()

double ncAna::getSystematicsEffAna(std::set<std::string> cuts, bool massConstraint)
{
	using std::set; using std::string;
	set<string>::const_iterator itSig,itBkg;
	set<string> sigDeps,bkgDeps;
	ncTree<ncFormula> *tree;
	TCut bkgCut,sgCut; // signal == variable we are processing, background == all other variables
	TTree *dataTree = (TTree*)dataFile->Get("T");
	TTree *mcTree = (TTree*)mcFile->Get("T");
	double total_cands,passed_cands;
	double eff_mc,eff_data,rel;
	double systemErr2 = 0;
	int counter = 0;
	TH1D *mHisto = getDefaultMassHistogram("mass");
	fit_t fit;
	TCanvas *c = new TCanvas;
	
	init_gauss_p2(&fit, 5.05, 5.45);
	
	// process all cuts
	for (itSig = cuts.begin(); itSig != cuts.end(); ++itSig) {
		
		tree = read_formula(*itSig); // read the current formula we are processing
		sigDeps = get_dependencies(tree);
		delete tree;
		
		// set up the signal cut
		sgCut = TCut(itSig->c_str());
		bkgCut = TCut();
		for (itBkg = cuts.begin(); itBkg != cuts.end(); ++itBkg) {
			if (itBkg == itSig) continue;
			tree = read_formula(*itBkg);
			bkgDeps = get_dependencies(tree);
			delete tree;
			
			if(find_first_of(bkgDeps.begin(), bkgDeps.end(), sigDeps.begin(), sigDeps.end()) == bkgDeps.end())
				bkgCut = bkgCut && TCut(itBkg->c_str());
		}
		
		std::cout << "--------------------" << std::endl;
		std::cout << "Testing signal Cut '" << sgCut.GetTitle() << "'" << std::endl;
		std::cout << "With background Cut '" << bkgCut.GetTitle() << "'" << std::endl;
		std::cout << "--------------------" << std::endl;
		
		// get efficiency from data
		bkgCut = bkgCut && ncAna::cutTrigger(false); // normalization trigger
		
		dataTree->Draw(Form("%s >> mass", (massConstraint ? "mass_c" : "mass")),bkgCut);
		adjust_parameter_gauss_p2(mHisto, fit.fit_fct);
		mHisto->Fit(fit.fit_fct,"R");
		mHisto->Fit(fit.fit_fct,"LR");
		mHisto->SetTitle(Form("Bkg Fit for cut '%s'",itSig->c_str()));
		mHisto->Draw("E1");
		c->SaveAs(Form("%s/SideBand_%d.eps",fPlotDir.c_str(),counter++));
		total_cands = signal_events_gauss_p2(fit.fit_fct, mHisto->GetBinWidth(1), NULL);
		
		dataTree->Draw(Form("%s >> mass", (massConstraint ? "mass_c" : "mass")), bkgCut && sgCut);
		adjust_parameter_gauss_p2(mHisto, fit.fit_fct);
		mHisto->Fit(fit.fit_fct,"R");
		mHisto->Fit(fit.fit_fct,"LR");
		mHisto->SetTitle(Form("Sig Fit for cut '%s'",itSig->c_str()));
		mHisto->Draw("E1");
		c->SaveAs(Form("%s/SideBand_%d.eps",fPlotDir.c_str(),counter++));
		passed_cands = signal_events_gauss_p2(fit.fit_fct, mHisto->GetBinWidth(1), NULL);
		
		eff_data = passed_cands / total_cands;
		
		// get efficiency from MC
		bkgCut = bkgCut && ncAna::cutNormTruth(); // truth in MC
		total_cands = (double)mcTree->Draw("", bkgCut);
		passed_cands = (double)mcTree->Draw("", bkgCut && sgCut);
		eff_mc = passed_cands / total_cands;
		
		rel = TMath::Abs(eff_data - eff_mc) / eff_mc;
		systemErr2 += rel*rel;
		
		std::cout << "--------------------" << std::endl;
		cout << "Data efficiency for Cut '" << *itSig << "': " << eff_data << endl;
		cout << "MC efficiency for Cut '" << *itSig << "': " << eff_mc << endl;
		cout << "Difference: " << TMath::Abs(eff_data - eff_mc) << endl;
		cout << "Rel Difference: " << TMath::Abs(eff_data - eff_mc) / eff_mc << endl;
		std::cout << "--------------------" << std::endl;
		cout << "Total Rel Uncertainty: " << TMath::Sqrt(systemErr2) << endl;
		std::cout << "--------------------" << std::endl;
	}
	
	delete mHisto;
	delete c;
	
	return TMath::Sqrt(systemErr2);
} // getSystematicsEffAna()
