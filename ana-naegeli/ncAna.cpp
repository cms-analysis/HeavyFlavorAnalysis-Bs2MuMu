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
#include "ncVarReader.h";

#include "../rootutils/NCRootUtils.h"

#include <iostream>

// import the utility functions from the ulcalc directory
#include "../ulcalc/external_constants.cpp"

// ROOT headers
#include <TMath.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <TEventList.h>
#include <TTreeFormula.h>

// RooFit
#include <RooDataSet.h>
#include <RooExtendPdf.h>
#include <RooRealVar.h>
#include <RooPlot.h>

// TMVA
#include <TMVA/Factory.h>

using std::cout;
using std::endl;

ncAna::ncAna() : fNormFileName(""), fSignalDataFileName("/Users/cn/CMSData/Reduced/cms-Run2012-sig.root"), fDataFileName("/Users/cn/CMSData/Reduced/data-newvars.root"), fMCFileName("/Users/cn/CMSData/Reduced/mc-newvars.root"), fPeakFileName("/Users/cn/CMSData/Reduced/production-2e33-general.root"), fAccFileName("/Users/cn/CMSData/Reduced/production-2e33-acceptance.root"), fNbrMassBins(50), fMassRange(4.9,5.9), fNbrSidebandBins(50), fLowNoSidebandRange(5.05,5.15), fHighNoSidebandRange(5.4,5.5), fLowCoSidebandRange(5.1,5.2), fHighCoSidebandRange(5.5,5.8), fBlindedRegion(5.2,5.45), fBsWindowRegion(5.3,5.45),fBdWindowRegion(5.2,5.3), fNoSignalRegion(5.2,5.35), fFitRangeNorm(4.95,5.6), fPlotDir("plots"), fSystematicsTable(NULL), fMisIDKaonPion(0.001,0.0002,0.0002), fMisIDProton(0.0005,0.0001,0.0001), fVarFileName("cuts/newvars.def")
{
	using std::string;
	using std::make_pair;
	typedef ncCut::ncCutRange range_t;
	
	// channel initialization
	fChannels.push_back(make_pair(range_t(0,1.4),string("pt_mu1 > 4.5 && pt_mu2 > 4.0 && pt > 6.5 && ip < 0.008 && ip/ipe < 2 && alpha < 0.05 && chi2/Ndof < 2.2 && d3/d3e > 13 && d3 < 1.5 && iso_mor12 > 0.8 && doca0 > 0.015 && ntrk < 2 && trk_weight > 0.6")));
	fChannels.push_back(make_pair(range_t(1.4,10),string("pt_mu1 > 4.5 && pt_mu2 > 4.2 && pt > 8.5 && ip < 0.008 && ip/ipe < 2 && alpha < 0.03 && chi2/Ndof < 1.8 && d3/d3e > 15 && d3 < 1.5 && iso_mor12 > 0.8 && doca0 > 0.015 && ntrk < 2 && trk_weight > 0.6")));
	
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
	
	loadSystematics();
} // ncAna()

ncAna::~ncAna()
{
//	if (fSystematicsTable) delete fSystematicsTable;
} // ~ncAna

TH1D *ncAna::getDefaultMassHistogram(const char *name)
{
	return new TH1D(name,"",fNbrMassBins,fMassRange.first, fMassRange.second);
} // getDefaultMassHistogram()

void ncAna::showSidebandSubtraction(bool massConstraint, bool cs)
{
	TFile dataFile(fSignalDataFileName.c_str());
	TFile mcFile(fMCFileName.c_str());
	TTree *dataTree = (TTree*)dataFile.Get("T");
	TTree *mcTree = (TTree*)mcFile.Get("T");
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
		sig.setErrHi(err);
		sig.setErrLo(err);
		
		bkg.setVal(background_events_gauss_p2(fit.fit_fct, mHisto->GetBinWidth(1))); // 2 sigma
		bkg.setErrHi(TMath::Sqrt(bkg.getVal()));
		bkg.setErrLo(TMath::Sqrt(bkg.getVal()));
		
		// histogram einteilen
		signalCut = TCut(Form("%f < mass && mass < %f",mu - 2*res, mu + 2*res));
		lowSBCut = (cs ? TCut(Form("%f < mass && mass < %f",fLowCoSidebandRange.first, fLowCoSidebandRange.second))
					: TCut(Form("%f < mass && mass < %f",fLowNoSidebandRange.first, fLowNoSidebandRange.second)));
		highSBCut = (cs ? TCut(Form("%f < mass && mass < %f",fHighCoSidebandRange.first, fHighCoSidebandRange.second))
					 : TCut(Form("%f < mass && mass < %f",fHighNoSidebandRange.first, fHighNoSidebandRange.second)));
		
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
	TFile dataFile(fDataFileName.c_str());
	TFile mcFile(fMCFileName.c_str());
	TTree *dataTree = (TTree*)dataFile.Get("T");
	TTree *mcTree = (TTree*)mcFile.Get("T");
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

// FIXME: Check Systematics here
void ncAna::loadSystematics(bool def)
{
	using std::make_pair;
	if (fSystematicsTable) delete fSystematicsTable;
	
	fSystematicsTable = new std::map<systematics_t,double>;
	
	fSystematicsTable->insert(make_pair(g_sys_acc_efftrack, 0.04));
	fSystematicsTable->insert(make_pair(g_sys_acc_ppro_barrel, 0.035));
	fSystematicsTable->insert(make_pair(g_sys_acc_ppro_endcap, 0.05));
	if(def) fSystematicsTable->insert(make_pair(g_sys_effana, 0.04));
	else	fSystematicsTable->insert(make_pair(g_sys_effana, this->getSystematicsEffAna(TCut("pt_mu1 > 4.5 && pt_mu2 > 4.0 && pt > 6.5 && ip < 0.008 && ip/ipe < 2 && alpha < 0.05 && chi2/Ndof < 2.2 && d3/d3e > 13 && d3 < 1.5 && iso_mor12 > 0.8 && doca0 > 0.015 && ntrk < 2 && trk_weight > 0.6"), true))); // FIXME: correct cut here (!!!)
	fSystematicsTable->insert(make_pair(g_sys_massscale, 0.03));
	fSystematicsTable->insert(make_pair(g_sys_effmu_barrel, 0.04));
	fSystematicsTable->insert(make_pair(g_sys_effmu_endcap, 0.08));
	fSystematicsTable->insert(make_pair(g_sys_efftrig_barrel, 0.03));
	fSystematicsTable->insert(make_pair(g_sys_efftrig_endcap, 0.06));
	fSystematicsTable->insert(make_pair(g_sys_normfit, 0.05));
	fSystematicsTable->insert(make_pair(g_sys_shapecombbkg, 0.04));
} // loadSystematics()

void ncAna::varAna()
{
	using std::cout; using std::endl;
	std::set<ncCut> *varSet;
	std::set<ncCut>::const_iterator it;
	std::string name;
	ncVarReader reader;
	size_t j;
	
	try {
		reader.loadFile(fVarFileName.c_str());
	} catch (std::string err) {
		cout << "Error loading the Variables..." << endl;
		cout << "Caught string '" << err << "'" << endl;
	}
	
	for (j = 0; j < reader.getNbr(); j++) {
		varSet = reader.getVars(j);
		
		cout << "====================" << endl;
		cout << "set[j] = " << j << endl;
		for (it = varSet->begin(); it != varSet->end(); ++it) it->dump();
		cout << "====================" << endl;
		
		// process the files...
		name = std::string(Form("variable-%d",(int)j));
		processVarSet(varSet,name.c_str());
	}
} // varAna()

void ncAna::processVarSet(std::set<ncCut> *varSet, const char *name)
{
	std::set<ncCut>::iterator it;
	TFile *mcFile = TFile::Open(fMCFileName.c_str());
	TFile *dataFile = TFile::Open(fDataFileName.c_str());
	TFile tmvaFile(Form("TMVA-%s.root",name),"recreate");
	TTree *mcTree = (TTree*)mcFile->Get("T");
	TTree *dataTree = (TTree*)dataFile->Get("T");
	TString factOptions("Transformations=I");
	TString prepOptions("");
	TMVA::Factory *factory;
	TFile *tmpFile;
	char *tmp = tempnam(".","training-");
	string *trainFilename;
	TCut cutPresel,cut;
	std::vector<std::pair<double,std::string> > ranking;
	size_t j;
	
	// build the cut
	for (it = varSet->begin(); it != varSet->end(); ++it) {
		ncCut c = *it;
		cutPresel = cutPresel && TCut(Form("%e < %s && %s < %e", c.getCut().first, c.getFormula(), c.getFormula(), c.getCut().second));
	}
	
	trainFilename = new string(Form("%s.root",tmp));
	tmpFile = new TFile(trainFilename->c_str(),"recreate");
	tmpFile->cd();
	
	// extracting signal
	cut = cutPresel && TCut("candidate == 1313") && cutSigTruth(true) && (cutTrigger(true, true) || cutTrigger(true, false));
	cout << "Selecting signal..." << flush;
	mcTree = mcTree->CopyTree(cut.GetTitle());
	cout << "	done" << endl;
	
	// extracting background
	cut = cutPresel && TCut("candidate == 1313") && (TCut(Form("%e < mass && mass < %e", fLowNoSidebandRange.first, fLowNoSidebandRange.second)) || TCut(Form("%e < mass && mass < %e", fHighCoSidebandRange.first, fHighCoSidebandRange.second))) && (cutTrigger(true, true) || cutTrigger(true, false));
	cout << "Selecting background..." << flush;
	dataTree = dataTree->CopyTree(cut.GetTitle());
	cout << "	done" << endl;
	
	mcTree->Write("sigT");
	dataTree->Write("bkgT");
	
	factory = new TMVA::Factory("varAna",&tmvaFile,factOptions);
	factory->AddSignalTree(mcTree,1.0);
	factory->AddBackgroundTree(dataTree,1.0);
	
	for (it = varSet->begin(); it != varSet->end(); ++it) {
		ncCut c = *it;
		factory->AddVariable(Form("%s := %s", c.getName(), c.getFormula()), 'F');
		ranking.push_back(make_pair(computeRanking(&c, mcTree, dataTree), std::string(c.getName())));
	}
	
	factory->PrepareTrainingAndTestTree("",prepOptions);
	factory->BookMethod(TMVA::Types::kCuts,"cuts","FitMethod=MC:SampleSize=100");
	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();
	
	// sort the variable ranking
	std::sort(ranking.begin(),ranking.end());
	cout << "Ranking of Variables" << endl;
	for (j = 0; j < ranking.size(); j++)
		cout << '\t' << ranking[j].second << ":		" << ranking[j].first << endl;
	
	unlink(trainFilename->c_str());
	free(tmp);
	delete tmpFile;
	delete trainFilename;
	delete factory;
	delete mcFile;
	delete dataFile;
} // processVarSet()

double ncAna::computeRanking(ncCut *nc, TTree *sigTree, TTree *bkgTree)
{
	TTreeFormula sigForm("sigForm",nc->getFormula(),sigTree);
	TTreeFormula bkgForm("bkgForm",nc->getFormula(),bkgTree);
	Long64_t j;
	double value;
	double sigMean = 0,sigSigma = 0;
	double bkgMean = 0,bkgSigma = 0;
	double nsig,nbkg;
	
	// process signal tree
	for (j = 0; j < sigTree->GetEntries(); j++) {
		
		sigTree->GetEntry(j);
		value = sigForm.EvalInstance();
		sigMean += value;
		sigSigma += value*value;
	}
	
	// process bkg tree
	for (j = 0; j < bkgTree->GetEntries(); j++) {
		
		bkgTree->GetEntry(j);
		value = bkgForm.EvalInstance();
		bkgMean += value;
		bkgSigma += value*value;
	}
	
	// compute mean and std dev
	nsig = (double)sigTree->GetEntries();
	nbkg = (double)bkgTree->GetEntries();
	sigMean /= nsig;
	bkgMean /= nbkg;
	
	sigSigma = sigSigma/(nsig-1.) - nsig/(nsig-1.)*(sigMean*sigMean); // sigma^2
	bkgSigma = bkgSigma/(nbkg-1.) - nbkg/(nbkg-1.)*(bkgMean*bkgMean); // sigma^2
	
	// add quadratically the sigmas
	sigSigma = TMath::Sqrt(sigSigma + bkgSigma);
	
	value = TMath::Abs(sigMean-bkgMean)/sigSigma;
	
	return value;
} // computeRanking()

#pragma mark -


TCut ncAna::cutPreselection(int nKaons)
{
	TCut result = TCut("d3 < 2.0 && d3e > 0 && ipe > 0") && cutHisto() && cutAcceptanceData(nKaons);
	
	switch (nKaons) {
		case 2:
			result = result && TCut("0.995 < mass_dikaon && mass_dikaon < 1.045 && deltaR_kaons < 0.25");
		case 1:
			result = result && TCut("pt > 7.0 && 3.0 < fMassJPsi && fMassJPsi < 3.2");
			break;
		default:
			result = TCut("pt > 7.5");
			break;
	}
	
	return result;
}

TCut ncAna::cutAcceptanceMC(int nKaons)
{
	TCut accCut("TMath::Abs(eta_mu1_gen) < 2.5 && TMath::Abs(eta_mu2_gen) < 2.5 && pt_mu1_gen > 1. && pt_mu2_gen > 1.");
	for (int n = 1; n <= nKaons; n++)
		accCut = accCut && TCut(Form("pt_kp%d_gen > 0.4 && TMath::Abs(eta_kp%d_gen) < 2.5",n,n));
	
	return accCut;
} // cutAcceptanceMC()

TCut ncAna::cutAcceptanceMCHard(int nKaons)
{
	TCut accCut("TMath::Abs(eta_mu1_gen) < 2.5 && TMath::Abs(eta_mu2_gen) < 2.5 && pt_mu1_gen > 3.5 && pt_mu2_gen > 3.5 && (track_qual_mu1 & 4) && (track_qual_mu2 & 4)");
	for (int n = 1; n <= nKaons; n++)
		accCut = accCut && TCut(Form("pt_kp%d_gen > 0.4 && TMath::Abs(eta_kp%d_gen) < 2.5",n,n));
	
	return accCut;
} // cutAcceptanceMCHard()

TCut ncAna::cutAcceptanceData(int nKaons)
{
	TCut accCut("(q_mu1*q_mu2 < 0) && TMath::Abs(eta_mu1) < 2.4 && TMath::Abs(eta_mu2) < 2.4 && pt_mu1 > 1. && pt_mu2 > 1. && (track_qual_mu1 & 4) && (track_qual_mu2 & 4) && (q_kp1 * q_kp2 != 1)");
	for (int n = 1; n <= nKaons; n++)
		accCut = accCut && TCut(Form("pt_kp%d > 0.5 && TMath::Abs(eta_kp%d) < 2.4 && (track_qual_kp%d & 4)",n,n,n));
	
	return accCut;
} // cutAcceptanceData()

TCut ncAna::cutSigCandGen(bool bsmm, bool reco)
{
	int cand = bsmm ? 80 : 90;
	if (!reco) cand = -cand;
	return TCut(Form("candidate == %d",cand));
} // cutSigCandGen()

TCut ncAna::cutChannel(unsigned ix)
{
	TCut up,low;
	
	if (ix >= fChannels.size()) throw std::string("Out of Channel exception!!!");
	
	low = TCut(Form("TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f", fChannels[ix].first.first, fChannels[ix].first.first));
	up = TCut(Form("TMath::Abs(eta_mu1) < %f && TMath::Abs(eta_mu2) < %f", fChannels[ix].first.second, fChannels[ix].first.second));
	
	return (up && !low);
} // cutChannel()

TCut ncAna::cutAna(unsigned ix, bool jpsi)
{
	TCut result;
	
	if (ix >= fChannels.size()) throw std::string("Out of Channel exception!!!");
	
	result = TCut(fChannels[ix].second.c_str());
	if (jpsi)
		result = result && TCut("3.0 < mass_dimuon && mass_dimuon < 3.2");
	
	return result;
} // cutAna()

#pragma mark -

void ncAna::computeBplus(unsigned channelIx, bool reload)
{
	using std::flush;
	using std::string;
	using std::make_pair;
	TEventList *elist;
	TEfficiency eff("eff","",1,0,1); // we only need one bin
	TFile *file;
	TTree *tree;
	string listName;
	ulc_t bplus;
	TCut cut;
	std::map<systematics_t,double>::iterator syst_it;
	measurement_t mes;
	
	double nbr_gens,nbr_acc,nbr_acc2;			// acceptance file
	double nbr_ana,nbr_reco,nbr_mu,nbr_trig;	// mc file
	
	cout << "===> ncAna::compute_bplus(" << channelIx << ")" << endl;
	
	// configure efficiency calculation
	eff.SetConfidenceLevel(0.68);
	eff.SetStatisticOption(TEfficiency::kFCP);
	
	/*******************
	 * Acceptance file *
	 *******************/
	file = new TFile(fAccFileName.c_str(),"update");
	tree = (TTree*)file->Get("T");
	
	// generated
	cout << "\tCounting generator canidates..." << flush;
	elist = NULL;
	listName = string("normGens");
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		// create the event list
		tree->Draw(Form(">>%s",listName.c_str()),cutNormCandGen(false));
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_gens = (double)elist->GetN();
	cout << '\t' << nbr_gens << endl;
	
	// accepted
	cout << "\tCounting accepted candidates..." << flush;
	elist = NULL;
	listName = string(Form("normGensAcc_%d",channelIx));
	if(!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		// create the event list
		cut = cutNormCandGen(true) && cutAcceptanceMC(1) && cutAcceptanceData(1) && cutChannel(channelIx);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_acc = (double)elist->GetN();
	cout << '\t' << nbr_acc << flush;
	
	elist = NULL;
	listName = string(Form("normGensAccHard_%d",channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		// create the event list
		cut = cutNormCandGen(true) && cutAcceptanceMCHard(1) && cutAcceptanceData(1) && cutChannel(channelIx);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_acc2 = (double)elist->GetN();
	cout << " (" << nbr_acc2 << " / " << flush;
	
	// close the file
	delete file;
	
	/***********
	 * MC File *
	 ***********/
	
	file = new TFile(fMCFileName.c_str(),"update");
	tree = (TTree*)file->Get("T");
	
	elist = NULL;
	listName = string(Form("normGensAccHard_%d",channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		// create the event list
		cut = cutNormCandGen(true) && cutAcceptanceMCHard(1) && cutAcceptanceData(1) && cutChannel(channelIx);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_reco = (double)elist->GetN();
	cout << nbr_reco << ")" << endl;
	
	// analysis
	cout << "\tCounting analysis candidates..." << flush;
	elist = NULL;
	listName = string(Form("normAna_%d",channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		// create the event list
		cut = cutNormCand() && cutAcceptanceMCHard(1) && cutAcceptanceData(1) && cutChannel(channelIx) && cutAna(channelIx,true) && cutNormTruth();
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_ana = (double)elist->GetN();
	cout << '\t' << nbr_ana << endl;
	
	// muon
	cout << "\tCounting muon candidates..." << flush;
	elist = NULL;
	listName = string(Form("normMuon_%d",channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		// create the event list
		cut = cutNormCand() && cutAcceptanceMCHard(1) && cutAcceptanceData(1) && cutChannel(channelIx) && cutAna(channelIx, true) && cutMuon() && cutNormTruth();
		tree->Draw(Form(">>%s", listName.c_str()), cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_mu = (double)elist->GetN();
	cout << '\t' << nbr_mu << endl;
	
	// trigger
	cout << "\tCounting trigger candidates..." << flush;
	elist = NULL;
	listName = string(Form("normTrig_%d",channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutNormCand() && cutAcceptanceMCHard(1) && cutAcceptanceData(1) && cutChannel(channelIx) && cutAna(channelIx, true) && cutMuon() && cutTrigger(false) && cutNormTruth();
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_trig = (double)elist->GetN();
	cout << '\t' << nbr_trig << endl;
	
	delete file;
	
	/***********
	 * Compute *
	 ***********/
	
	cout << "	Computing efficiencies..." << flush;
	// acceptance
	eff.SetTotalEvents(0, nbr_gens);
	eff.SetPassedEvents(0, nbr_acc);
	bplus[make_pair(kAcc_bplus, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	// systematics tracking efficiency
	if (!fSystematicsTable) loadSystematics(true);
	if( (syst_it = fSystematicsTable->find(g_sys_acc_efftrack)) != fSystematicsTable->end() ) {
		mes = bplus[make_pair(kAcc_bplus, channelIx)];
		mes = measurement_t(0, mes.getVal() * syst_it->second, mes.getVal() * syst_it->second);
		bplus[make_pair(kAcc_bplus, channelIx)] = mes + bplus[make_pair(kAcc_bplus, channelIx)];
	}
	
	// Candidate efficiency
	bplus[make_pair(kEff_cand_bplus, channelIx)] = measurement_t(1., 0., 0.);
	
	// Analysis efficiency
	eff.SetPassedEvents(0, nbr_ana);
	eff.SetTotalEvents(0, nbr_reco*nbr_acc/nbr_acc2);
	bplus[make_pair(kEff_ana_bplus, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	// muon efficiency
	eff.SetTotalEvents(0, nbr_ana);
	eff.SetPassedEvents(0, nbr_mu);
	bplus[make_pair(kEff_mu_bplus, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	// trigger efficiency
	eff.SetTotalEvents(0, nbr_mu);
	eff.SetPassedEvents(0, nbr_trig);
	bplus[make_pair(kEff_trig_bplus, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	cout << "\tdone" << endl;
	
	/***************
	 * Observation *
	 ***************/
	
	cout << "\tMeasuring number of B+ -> J/psi K+ decays..." << endl;
	bplus[make_pair(kObs_bplus,channelIx)] = fitBplus(channelIx,reload);
	if ( (syst_it = fSystematicsTable->find(g_sys_normfit)) != fSystematicsTable->end() ) {
		mes = bplus[make_pair(kObs_bplus,channelIx)];
		mes = measurement_t(0, mes.getVal() * syst_it->second, mes.getVal() * syst_it->second);
		bplus[make_pair(kObs_bplus,channelIx)] = mes + bplus[make_pair(kObs_bplus,channelIx)];
	}
	cout << "\tdone" << endl;
	
	// save locally
	fUlcBplus = bplus;
} // compute_bplus()

void ncAna::computeBmm(unsigned channelIx, bool bsmm, bool reload)
{
	using std::flush;
	using std::string;
	using std::make_pair;
	map<systematics_t,double>::const_iterator syst_it;
	measurement_t m;
	TEfficiency eff("eff","",1,0,1); // we only need one bin
	TEventList *elist;
	string listName;
	ulc_t bmm;
	TFile *file;
	TTree *tree;
	TCut cut;
	const char *prefix = bsmm ? "bsmm" : "bdmm";
	
	double nbr_gens,nbr_reco,nbr_ana,nbr_cand;
	double nbr_mu,nbr_trig, nbr_bs, nbr_bd;
	
	cout << "===> ncAna::computeB" << (bsmm ? 's' : 'd') << "mm(" << channelIx << ")" << endl;
	
	// configure the efficiency calculation
	eff.SetConfidenceLevel(0.68);
	eff.SetStatisticOption(TEfficiency::kFCP);
	
	/*********
	 * COUNT *
	 *********/
	
	file = new TFile(fMCFileName.c_str(),"update");
	tree = (TTree*)file->Get("T");
	
	// FIXME: Merge the subsequent block into a function
	// generated
	cout << "\tCounting generator candidates..." << flush;
	elist = NULL;
	listName = string(Form("%sGens",prefix));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		tree->Draw(Form(">>%s",listName.c_str()),cutSigCandGen(bsmm, false));
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_gens = (double)elist->GetN();
	cout << '\t' << nbr_gens << endl;
	
	// accepted
	cout << "\tCounting reco candidates..." << flush;
	elist = NULL;
	listName = string(Form("%sReco_%u",prefix,channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutSigCandGen(bsmm, true) && cutAcceptanceMC(0) && cutAcceptanceData(0) && cutChannel(channelIx);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_reco = (double)elist->GetN();
	cout << '\t' << nbr_reco << endl;
	
	// candidate efficiency
	cout << "\tCounting cand candidates..." << flush;
	elist = NULL;
	listName = string(Form("%sCand_%u",prefix,channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	else {
		cut = cutSigCand() && cutAcceptanceMC(0) && cutAcceptanceData(0) && cutChannel(channelIx) && cutSigTruth(bsmm);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_cand = (double)elist->GetN();
	cout << '\t' << nbr_cand << endl;
	
	// analysis
	cout << "\tCounting ana candiates..." << flush;
	elist = NULL;
	listName = string(Form("%sAna_%u",prefix,channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutSigCand() && cutAcceptanceMC(0) && cutAcceptanceData(0) && cutChannel(channelIx) && cutAna(channelIx, false) && cutSigTruth(bsmm);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_ana = (double)elist->GetN();
	cout << '\t' << nbr_ana << endl;
	
	// muon
	cout << "\tCounting muon candidates..." << flush;
	elist = NULL;
	listName = string(Form("%sMuon_%u",prefix,channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutSigCand() && cutAcceptanceMC(0) && cutAcceptanceData(0) && cutChannel(channelIx) && cutAna(channelIx, false) && cutMuon() && cutSigTruth(bsmm);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_mu = (double)elist->GetN();
	cout << '\t' << nbr_mu << endl;
	
	// Trigger
	cout << "\tCounting trigger candidates..." << flush;
	elist = NULL;
	listName = string(Form("%sTrig_%u",prefix,channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutSigCand() && cutAcceptanceMC(0) && cutAcceptanceData(0) && cutChannel(channelIx) && cutAna(channelIx, false) && cutMuon() && cutTrigger(true) && cutSigTruth(bsmm);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_trig = (double)elist->GetN();
	cout << '\t' << nbr_trig << endl;
	
	// cross feed
	cout << "\tComputing cross feed..." << flush;
	elist = NULL;
	listName = string(Form("%sBsWind_%u",prefix,channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutSigCand() && cutAcceptanceMC(0) && cutAcceptanceData(0) && cutChannel(channelIx) && cutAna(channelIx, false) && cutMuon() && cutTrigger(true) && cutSigTruth(bsmm) && cutSigWindow(true);
		tree->Draw(Form(">>%s",listName.c_str()), cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_bs = (double)elist->GetN();
	
	elist = NULL;
	listName = string(Form("%sBdWind_%u",prefix,channelIx));
	if (!reload) elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutSigCand() && cutAcceptanceMC(0) && cutAcceptanceData(0) && cutChannel(channelIx) && cutAna(channelIx, false) && cutMuon() && cutTrigger(true) && cutSigTruth(bsmm) && cutSigWindow(false);
		tree->Draw(Form(">>%s",listName.c_str()), cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_bd = (double)elist->GetN();
	cout << "\tdone" << endl;
	
	delete file;
	
	/***********
	 * COMPUTE *
	 ***********/
	
	cout << "\tComputing efficiencies..." << flush;
	// acceptance
	eff.SetTotalEvents(0, nbr_gens);
	eff.SetPassedEvents(0, nbr_reco);
	bmm[make_pair(kAcc_bmm, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	// candidate efficiency
	eff.SetTotalEvents(0,nbr_reco);
	eff.SetTotalEvents(0,nbr_cand);
	bmm[make_pair(kEff_cand_bmm, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	// analysis efficiency
	eff.SetTotalEvents(0, nbr_cand);
	eff.SetPassedEvents(0, nbr_ana);
	bmm[make_pair(kEff_ana_bmm, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	if( (syst_it = fSystematicsTable->find(g_sys_effana)) != fSystematicsTable->end() ) {
		m = bmm[make_pair(kEff_ana_bmm, channelIx)];
		m = measurement_t(0, m.getVal() * syst_it->second,m.getVal() * syst_it->second);
		bmm[make_pair(kEff_ana_bmm, channelIx)] = m + bmm[make_pair(kEff_ana_bmm, channelIx)];
	}
	
	// muon efficiency
	eff.SetTotalEvents(0, nbr_ana);
	eff.SetPassedEvents(0, nbr_mu);
	bmm[make_pair(kEff_mu_bmm, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	// trigger efficiency
	eff.SetTotalEvents(0, nbr_mu);
	eff.SetPassedEvents(0, nbr_trig);
	bmm[make_pair(kEff_trig_bmm, channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	// cross feed
	eff.SetTotalEvents(0,nbr_trig);
	eff.SetPassedEvents(0,nbr_bs);
	bmm[make_pair(kProb_swind_bmm,channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	
	eff.SetTotalEvents(0,nbr_trig);
	eff.SetPassedEvents(0,nbr_bd);
	bmm[make_pair(kProb_dwind_bmm,channelIx)] = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	cout << "\tdone" << endl;
	
	/***************
	 * OBSERVATION *
	 ***************/
	
	file = new TFile(fDataFileName.c_str(),"update");
	tree = (TTree*)file->Get("T");
	
	cout << "\tComputing background observation..." << flush;
	elist = NULL;
	listName = string(Form("bmmObsBkg_%u",channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutSigCand() && cutAcceptanceData(0) && cutChannel(channelIx) && cutAna(channelIx, false) && cutMuon() && cutTrigger(true) && cutHisto() && !cutSigWindow(true) && !cutSigWindow(false);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	bmm[make_pair(kObsBkg_bmm, channelIx)] = measurement_t((double)elist->GetN(),0.0,0.0);
	cout << "\tdone" << endl;
	
	cout << "\tComputing signal observation..." << flush;
	elist = NULL;
	listName = string(Form("%sObs_%u",prefix,channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		cut = cutSigCand() && cutAcceptanceData(0) && cutChannel(channelIx) && cutAna(channelIx, false) && cutMuon() && cutTrigger(true) && cutSigWindow(bsmm);
		tree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	bmm[make_pair(kObsB_bmm, channelIx)] = measurement_t((double)elist->GetN(),0.0,0.0);
	cout << "\tdone" << endl;
	
	delete file;
	
	if (bsmm)	fUlcBsmm = bmm;
	else		fUlcBdmm = bmm;
} // computeBmm()

void ncAna::computePeaking(unsigned channelIx, bool reload)
{
	TEfficiency eff("eff","",2,0,1);
	TFile *file;
	TTree *tree;
	double nbr_gens,nbr_acc;
	measurement_t filterHad,filterSemi;
	measurement_t effSingleMu,m;
	
	cout << "===> ncAna::computePeaking(" << channelIx << ")" << endl;
	
	// configure efficiency calculation
	eff.SetConfidenceLevel(0.68);
	eff.SetStatisticOption(TEfficiency::kFCP);
	
	// compute hadronic filter efficiency based on bsmm
	file = new TFile(fMCFileName.c_str());
	tree = (TTree*)gDirectory->Get("T");
	
	nbr_gens = (double)tree->Draw("",cutSigCandGen(true, false));
	nbr_acc = (double)tree->Draw("",cutSigCandGen(true, true) && cutAcceptanceMCHard(0));
	eff.SetTotalEvents(0,nbr_gens);
	eff.SetPassedEvents(0,nbr_acc);
	filterHad = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
	delete file;
	
	// compute semileptonic filter efficiency based on Bd -> pi mu nu
	file = new TFile(fAccFileName.c_str());
	tree = (TTree*)gDirectory->Get("T");
	
	nbr_gens = (double)tree->Draw("", "candidate == -95");
	nbr_acc = (double)tree->Draw("", TCut("candidate == 1000095") && cutAcceptanceMCHard(0));
	eff.SetTotalEvents(1,nbr_gens);
	eff.SetPassedEvents(1,nbr_acc);
	filterSemi = measurement_t(eff.GetEfficiency(1), eff.GetEfficiencyErrorUp(1), eff.GetEfficiencyErrorLow(1));
	delete file;
	
	// compute single muon efficiency
	if (fUlcBsmm.count(make_pair(kEff_mu_bmm,channelIx)) == 0) computeBmm(channelIx, true, reload);
	effSingleMu = fUlcBsmm[make_pair(kEff_mu_bmm,channelIx)];
	effSingleMu = measurement_t(TMath::Sqrt(effSingleMu.getVal()), effSingleMu.getErrHi() / TMath::Sqrt(2*effSingleMu.getVal()), effSingleMu.getErrLo() / TMath::Sqrt(2*effSingleMu.getVal()));
	
	// init value
	fUlcBsmm[make_pair(kPeakBkgOff_bmm,channelIx)] = measurement_t(0,0,0);
	fUlcBsmm[make_pair(kPeakBkgOn_bmm,channelIx)] = measurement_t(0,0,0);
	fUlcBdmm[make_pair(kPeakBkgOff_bmm,channelIx)] = measurement_t(0,0,0);
	fUlcBdmm[make_pair(kPeakBkgOn_bmm,channelIx)] = measurement_t(0,0,0);
	
	// compute the peaking background efficiency of different channels
	appendPeakingChannel(82, fMisIDKaonPion, fMisIDKaonPion, filterHad, bf_BsToKK()*f_ratio(), channelIx, reload);
	appendPeakingChannel(83, fMisIDKaonPion, fMisIDKaonPion, filterHad, bf_BsToKPi()*f_ratio(), channelIx, reload);
	appendPeakingChannel(84, fMisIDKaonPion, fMisIDKaonPion, filterHad, bf_BsToPiPi()*f_ratio(), channelIx, reload);
	appendPeakingChannel(86, fMisIDKaonPion, effSingleMu, filterSemi, bf_BsToKMuNu()*f_ratio(), channelIx, reload);
	appendPeakingChannel(91, fMisIDKaonPion, fMisIDKaonPion, filterHad, bf_BdToPiPi(), channelIx, reload);
	appendPeakingChannel(92, fMisIDKaonPion, fMisIDKaonPion, filterHad, bf_BdToKPi(), channelIx, reload);
	appendPeakingChannel(93, fMisIDKaonPion, fMisIDKaonPion, filterHad, bf_BdToKK(), channelIx, reload);
	appendPeakingChannel(95, fMisIDKaonPion, effSingleMu, filterSemi, bf_BdToPiMuNu(), channelIx, reload);
	appendPeakingChannel(60, fMisIDProton, fMisIDKaonPion, filterHad, bf_LambdaBToPPi()*f_ratio_lb(), channelIx, reload);
	appendPeakingChannel(61, fMisIDProton, fMisIDKaonPion, filterHad, bf_LambdaBToPK()*f_ratio_lb(), channelIx, reload);
	appendPeakingChannel(62, fMisIDProton, effSingleMu, filterSemi, bf_LambdaBToPMuNu()*f_ratio_lb(), channelIx, reload);
	
	m = fUlcBsmm[make_pair(kPeakBkgOff_bmm,channelIx)];
	m = fUlcBsmm[make_pair(kPeakBkgOn_bmm,channelIx)];
	m = fUlcBdmm[make_pair(kPeakBkgOff_bmm,channelIx)];
	m = fUlcBdmm[make_pair(kPeakBkgOn_bmm,channelIx)];
} // computePeaking()

void ncAna::appendPeakingChannel(int trueCand, measurement_t muMis1, measurement_t muMis2, measurement_t effFilter, measurement_t bf, unsigned channelIx, bool reload)
{
	using std::endl; using std::cout;
	TEfficiency eff("eff","",3,0,1); // bkg, bd, bs
	TFile peakFile(fPeakFileName.c_str(),"update");
	double nbr_reco,nbr_off,nbr_bs,nbr_bd;
	TTree *tree = (TTree*)peakFile.Get("T");
	TEventList *elist;
	string listName;
	TCut candCut(Form("candidate == %d",1000000+trueCand));
	TCut cut = cutHisto() && cutChannel(channelIx) && cutAna(channelIx, false);
	measurement_t effTrigger = fUlcBsmm[make_pair(kEff_trig_bmm,channelIx)];
	measurement_t m,efficiencies = effFilter * muMis1 * muMis2 * effTrigger * bf;
	measurement_t tot_bu;
	
	if (fUlcBplus.count(make_pair(kTot_bplus,channelIx)) == 0) {
		computeBplus(channelIx, reload);
		fUlcBplus[make_pair(kTot_bplus, channelIx)] = fUlcBplus[make_pair(kObs_bplus, channelIx)] / compute_efftot_bplus(&fUlcBplus, channelIx);
	}
	tot_bu = fUlcBplus[make_pair(kTot_bplus,channelIx)] / (bf_Bu2JpsiKp() * bf_PsiToMuMu());
	
	// maybe the fit bplus has changed current file
	peakFile.cd();
	
	eff.SetConfidenceLevel(0.68);
	eff.SetStatisticOption(TEfficiency::kFCP);
	
	cout << "	appendPeakingChannel(" << trueCand << ")" << endl;
	
	// reco
	elist = NULL;
	listName = string(Form("peak%dReco_%u",trueCand,channelIx));
	if (!reload)
		elist = (TEventList*)gDirectory->Get(listName.c_str());
	if (!elist) {
		tree->Draw(Form(">>%s",listName.c_str()),candCut && cutAcceptanceMCHard(0));
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(), TObject::kOverwrite);
	}
	nbr_reco = (double)elist->GetN();
	
	// fix this as base
	tree->SetEventList(elist);
	
	// nbr of passed candidates
	nbr_off = (double)tree->Draw("", cut && !cutSigWindow(true) && !cutSigWindow(false));
	nbr_bd = (double)tree->Draw("", cut && cutSigWindow(false));
	nbr_bs = (double)tree->Draw("", cut && cutSigWindow(true));
	
	if (nbr_reco > 0) {
		
		// bkg
		eff.SetTotalEvents(0, nbr_reco);
		eff.SetPassedEvents(0, nbr_off);
		m = measurement_t(eff.GetEfficiency(0),eff.GetEfficiencyErrorUp(0),eff.GetEfficiencyErrorLow(0));
		m = m * efficiencies;
		m = m * tot_bu;
		fUlcBsmm[make_pair(kPeakBkgOff_bmm,channelIx)] = m + fUlcBsmm[make_pair(kPeakBkgOff_bmm,channelIx)];
		fUlcBdmm[make_pair(kPeakBkgOff_bmm,channelIx)] = m + fUlcBdmm[make_pair(kPeakBkgOff_bmm,channelIx)];
		
		// bd
		eff.SetTotalEvents(1, nbr_reco);
		eff.SetPassedEvents(1, nbr_bd);
		m = measurement_t(eff.GetEfficiency(1),eff.GetEfficiencyErrorUp(1),eff.GetEfficiencyErrorLow(1));
		fUlcBdmm[make_pair(kPeakBkgOn_bmm,channelIx)] = m*efficiencies*tot_bu + fUlcBdmm[make_pair(kPeakBkgOn_bmm,channelIx)];
		
		
		eff.SetTotalEvents(2, nbr_reco);
		eff.SetTotalEvents(2, nbr_bs);
		m = measurement_t(eff.GetEfficiency(2),eff.GetEfficiencyErrorUp(2),eff.GetEfficiencyErrorLow(2));
		fUlcBsmm[make_pair(kPeakBkgOn_bmm,channelIx)] = m*efficiencies*tot_bu + fUlcBsmm[make_pair(kPeakBkgOn_bmm,channelIx)];
	} else
		std::cerr << "===> no reco candidate for rare background " << trueCand << "." << endl;
} // appendPeakingChannel()

measurement_t ncAna::fitBplus(unsigned channelIx, bool reload, RooWorkspace **wout)
{
	using std::string;
	RooWorkspace *w = new RooWorkspace("w");
	RooDataSet *data;
	string listName(Form("normChannel_%d",channelIx));
	TEventList *elist = NULL;
	TFile *dataFile = NULL;
	TTree *sourceTree,*copyTree;
	char *treeFile = NULL;
	TFile *smallFile;
	RooExtendPdf *ex;
	measurement_t m1,m2;
	TCut cut;
	TCanvas *c = NULL;
	
	std::cout << "===> fitBplus(" << channelIx << ")" << std::endl;
	
	// get the event list
	std::cout << "		Getting candidate list..." << std::flush;
	dataFile = new TFile(fDataFileName.c_str(),"update");
	sourceTree = (TTree*)dataFile->Get("T");
	if (!reload)
		elist = (TEventList*)dataFile->Get(listName.c_str());
	if (!elist) {
		cut = cutNormCand() && cutAcceptanceData(1) && cutChannel(channelIx) && cutAna(channelIx, true) && cutMuon() && cutTrigger(false);
		sourceTree->Draw(Form(">>%s",listName.c_str()),cut);
		elist = (TEventList*)gDirectory->Get(listName.c_str());
		elist->Write(listName.c_str(),TObject::kOverwrite);
	}
	std::cout << "	done" << std::endl;
	
	// get the working tree
	std::cout << "		Extracting Tree..." << std::flush;
	sourceTree->SetEventList(elist);
	treeFile = tempnam(".", "fitting-");
	smallFile = new TFile(treeFile, "recreate");
	copyTree = sourceTree->CopyTree("");
	std::cout << "	done" << std::endl;
	
	// import data to roofit
	w->factory(Form("mass_c[5.28,%f,%f]",fMassRange.first,fMassRange.second));
	w->defineSet("obs","mass_c");
	
	data = new RooDataSet("data","Norm mass distribution",copyTree,*w->set("obs"));
	w->import(*data);
	delete data;
	
	// setup model
	w->factory("Gaussian::sig(mass_c,mu[5.28,4.9,5.9],sigma[0.01,0,0.02])");
	w->factory("Gaussian::sig2(mass_c,mu,sigma2[0.03,0.02,0.4])");
	w->factory("Gaussian::bump(mass_c,mu_bump[5.1,5.0,5.2],sigma_bump[0.01,0.0,0.1])");
	w->factory("Gaussian::bump2(mass_c,mu_bump2[4.95,4.8,5.05],sigma_bump[0.01,0.0,0.1])");
	w->factory("Exponential::bkg_comb(mass_c,c[0,-1e30,1e30])");
	w->factory("GenericPdf::bkg_peak(\"TMath::Erf(sc*(mass_c - th))/2.0 + 0.5\",{sc[-10,-1000,1000],mass_c,th[5.15,4.9,5.2]})");
	w->factory(Form("SUM::model(nsig[%f,0,1e7]*sig,nsig2[0,0,1e7]*sig2,ncomb[%f,0,1e7]*bkg_comb,npeak[%f,0,1e7]*bkg_peak,nbump[0,0,1e7]*bump,nbump2[0,0,1e7]*bump2)",(double)data->numEntries(),(double)data->numEntries(),(double)data->numEntries()));
	
	// fit
	
	// initial value setup
	ex = new RooExtendPdf((const char *)"",(const char *)"",*w->pdf("bkg_peak"),*w->var("npeak"));
	ex->fitTo(*w->data("data"),RooFit::Range(5.0,5.15));
	delete ex;
	
	ex = new RooExtendPdf("","",*w->pdf("bkg_comb"),*w->var("ncomb"));
	ex->fitTo(*w->data("data"),RooFit::Range(5.4,5.6));
	delete ex;
	
	ex = new RooExtendPdf("","",*w->pdf("sig"),*w->var("nsig"));
	ex->fitTo(*w->data("data"),RooFit::Range(5.25,5.35));
	delete ex;
	
	// final fit
	w->pdf("model")->fitTo(*w->data("data"), RooFit::Range(fFitRangeNorm.first,fFitRangeNorm.second));
	
	// extract result
	m1 = measurement_t(w->var("nsig")->getVal(),w->var("nsig")->getError(),w->var("nsig")->getError());
	m2 = measurement_t(w->var("nsig2")->getVal(),w->var("nsig2")->getError(),w->var("nsig2")->getError());
	
	w->writeToFile(Form("%s/bplus_fit_ch%u.root",fPlotDir.c_str(),channelIx));
	
	{
		RooPlot *p = w->var("mass_c")->frame();
		c = new TCanvas;
		w->data("data")->plotOn(p,RooFit::Binning(50));
		w->pdf("model")->plotOn(p);
		p->Draw();
		c->SaveAs(Form("%s/bplus_fit_ch%u.pdf",fPlotDir.c_str(),channelIx));
	}
	
	// clean up
	if (wout)	*wout = w;
	else		delete w;
	delete smallFile;
	if(c) delete c;
	unlink(treeFile);
	free(treeFile);
	
	return (m1 + m2);
} // fitBplus()

void ncAna::writeULC(const char *ulcname, bool reload)
{
	using std::make_pair;
	unsigned channelIx;
	measurement_t m;
	FILE *outputFile = fopen(ulcname, "w");
	
	for (channelIx = 0; channelIx < fChannels.size(); channelIx++) {
		
		/*************************
		 * NORMALIZATION CHANNEL *
		 *************************/
		
		// Normalization stuff
		computeBplus(channelIx,reload);
		
		// add the numbers to the file
		fprintf(outputFile, "####################################################\n");
		fprintf(outputFile, "# B+ -> J/psi K+ (channel = %u, %.2f < eta < %.2f) #\n", channelIx, fChannels[channelIx].first.first, fChannels[channelIx].first.second);
		fprintf(outputFile, "####################################################\n");
		
		// acceptance
		m = fUlcBplus[make_pair(kAcc_bplus, channelIx)];
		fprintf(outputFile, "ACC_BPLUS\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		// muon efficiency
		m = fUlcBplus[make_pair(kEff_mu_bplus, channelIx)];
		fprintf(outputFile, "EFF_MU_BPLUS\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		// trigger efficiency
		m = fUlcBplus[make_pair(kEff_trig_bplus, channelIx)];
		fprintf(outputFile, "EFF_TRIG_BPLUS\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		// cand efficiency
		m = fUlcBplus[make_pair(kEff_cand_bplus, channelIx)];
		fprintf(outputFile, "EFF_CAND_BPLUS\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		// analysis efficiency
		m = fUlcBplus[make_pair(kEff_ana_bplus, channelIx)];
		fprintf(outputFile, "EFF_ANA_BPLUS\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		// observed bpluses
		m = fUlcBplus[make_pair(kObs_bplus, channelIx)];
		fprintf(outputFile, "OBS_BPLUS\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		
		// for convenience, compute eff_tot and tot_bplus in comment line...
		m = compute_efftot_bplus(&fUlcBplus, channelIx);
		fprintf(outputFile, "# EFF_TOT_BPLUS\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(), m.getErrLo());
		m = fUlcBplus[make_pair(kObs_bplus, channelIx)] / m;
		fprintf(outputFile, "# TOT_BPLUS\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		
		// save the total bplus -> j/psi k+
		fUlcBplus[make_pair(kTot_bplus, channelIx)] = m;
		
		/******************
		 * SIGNAL CHANNEL *
		 ******************/
		
		computeBmm(channelIx, true, reload);
		computeBmm(channelIx, false, reload);
		
		// add the numbers to the file
		fprintf(outputFile, "#######################################\n");
		fprintf(outputFile, "# B -> mumu (channel = %u, %.2f < eta < %.2f) #\n", channelIx, fChannels[channelIx].first.first, fChannels[channelIx].first.second);
		fprintf(outputFile, "#######################################\n");
		fprintf(outputFile, "LOW_BD\t%u\t%f\n", channelIx, fBdWindowRegion.first);
		fprintf(outputFile, "HIGH_BD\t%u\t%f\n", channelIx, fBdWindowRegion.second);
		fprintf(outputFile, "LOW_BS\t%u\t%f\n", channelIx, fBsWindowRegion.first);
		fprintf(outputFile, "HIGH_BS\t%u\t%f\n", channelIx, fBsWindowRegion.second);
		
		// crossfeed
		m = fUlcBsmm[make_pair(kProb_swind_bmm, channelIx)];
		fprintf(outputFile, "PSS\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		m = fUlcBdmm[make_pair(kProb_swind_bmm, channelIx)];
		fprintf(outputFile, "PSD\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		m = fUlcBsmm[make_pair(kProb_dwind_bmm, channelIx)];
		fprintf(outputFile, "PDS\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		m = fUlcBdmm[make_pair(kProb_dwind_bmm, channelIx)];
		fprintf(outputFile, "PDD\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		
		/********
		 * BSMM *
		 ********/
		// acceptance
		m = fUlcBsmm[make_pair(kAcc_bmm, channelIx)];
		fprintf(outputFile, "ACC_BSMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// muon efficiency
		m = fUlcBsmm[make_pair(kEff_mu_bmm, channelIx)];
		fprintf(outputFile, "EFF_MU_BSMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// trigger efficiency
		m = fUlcBsmm[make_pair(kEff_trig_bmm, channelIx)];
		fprintf(outputFile, "EFF_TRIG_BSMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// cand efficiency
		m = fUlcBsmm[make_pair(kEff_cand_bmm, channelIx)];
		fprintf(outputFile, "EFF_CAND_BSMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// ana effeciency
		m = fUlcBsmm[make_pair(kEff_ana_bmm, channelIx)];
		fprintf(outputFile, "EFF_ANA_BSMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// total efficiency for convenience
		m = compute_efftot_bmm(&fUlcBsmm, channelIx);
		fprintf(outputFile, "# EFF_TOT_BSMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		
		/********
		 * BDMM *
		 ********/
		// acceptance
		m = fUlcBdmm[make_pair(kAcc_bmm, channelIx)];
		fprintf(outputFile, "ACC_BDMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// muon efficiency
		m = fUlcBdmm[make_pair(kEff_mu_bmm, channelIx)];
		fprintf(outputFile, "EFF_MU_BDMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// trigger efficiency
		m = fUlcBdmm[make_pair(kEff_trig_bmm, channelIx)];
		fprintf(outputFile, "EFF_TRIG_BDMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// cand efficiency
		m = fUlcBdmm[make_pair(kEff_cand_bmm, channelIx)];
		fprintf(outputFile, "EFF_CAND_BDMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// ana efficiency
		m = fUlcBdmm[make_pair(kEff_ana_bmm, channelIx)];
		fprintf(outputFile, "EFF_ANA_BDMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		// total efficiency for convenience
		m = compute_efftot_bmm(&fUlcBdmm, channelIx);
		fprintf(outputFile, "# EFF_TOT_BDMM\t%u\t%f\t%f\t%f\n", channelIx, m.getVal(), m.getErrHi(),m.getErrLo());
		
		/***************
		 * OBSERVATION *
		 ***************/
		fprintf(outputFile, "OBS_BKG\t%u\t%f\n",  channelIx, fUlcBsmm[make_pair(kObsBkg_bmm, channelIx)].getVal());
		fprintf(outputFile, "OBS_BSMM\t%u\t%f\n", channelIx, fUlcBsmm[make_pair(kObsB_bmm, channelIx)].getVal());
		fprintf(outputFile, "OBS_BDMM\t%u\t%f\n", channelIx, fUlcBdmm[make_pair(kObsB_bmm, channelIx)].getVal());		
		
		/***********
		 * PEAKING *
		 ***********/
		computePeaking(channelIx, reload);
		m = fUlcBsmm[make_pair(kPeakBkgOff_bmm, channelIx)];
		fprintf(outputFile, "PEAK_BKG_OFF\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		m = fUlcBsmm[make_pair(kPeakBkgOn_bmm, channelIx)];
		fprintf(outputFile, "PEAK_BKG_BS\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		m = fUlcBdmm[make_pair(kPeakBkgOn_bmm, channelIx)];
		fprintf(outputFile, "PEAK_BKG_BD\t%u\t%f\t%f\t%f\n",channelIx,m.getVal(),m.getErrHi(),m.getErrLo());
		
		// Print tau
		fprintf(outputFile, "TAU_BS\t%u\t%f\t%f\t%f\n", channelIx, fUlcBsmm[make_pair(kTau_bmm, channelIx)].getVal(),fUlcBsmm[make_pair(kTau_bmm, channelIx)].getErrHi(),fUlcBsmm[make_pair(kTau_bmm, channelIx)].getErrLo());
		fprintf(outputFile, "TAU_BD\t%u\t%f\t%f\t%f\n", channelIx, fUlcBdmm[make_pair(kTau_bmm, channelIx)].getVal(),fUlcBdmm[make_pair(kTau_bmm, channelIx)].getErrHi(),fUlcBdmm[make_pair(kTau_bmm, channelIx)].getErrLo());
		
		// Print extra stuff as comment
		// needed for mva training
		{
			measurement_t tot_bpjpsik;
			measurement_t exp_bmm;
			measurement_t gen_bmm;
			TFile file(fMCFileName.c_str());
			TEventList *elist;
			
			fprintf(outputFile, "# Scaling factor for Background given by TAU_BS and TAU_BD\n");
			tot_bpjpsik = fUlcBplus[make_pair(kTot_bplus, channelIx)];
			
			// scaling bsmm
			elist = (TEventList*)file.Get("bsmmGens");
			exp_bmm = f_ratio() * tot_bpjpsik * bf_ratio_bsmm() / (bf_Bu2JpsiKp() * bf_PsiToMuMu());
			gen_bmm = measurement_t((double)elist->GetN(),0,0);
			fprintf(outputFile, "# SCALE_BS = %f\n", exp_bmm.getVal()/gen_bmm.getVal());
			
			// scaling bdmm
			elist = (TEventList*)file.Get("bdmmGens");
			exp_bmm = tot_bpjpsik * bf_ratio_bdmm() / (bf_Bu2JpsiKp() * bf_PsiToMuMu());
			gen_bmm = measurement_t((double)elist->GetN(),0,0);
			fprintf(outputFile, "# SCALE_BD = %f\n", exp_bmm.getVal()/gen_bmm.getVal());			
			
			// dump output, all B+
			exp_bmm = tot_bpjpsik / (bf_Bu2JpsiKp() * bf_PsiToMuMu());
			fprintf(outputFile, "# TOTAL PRODUCED B+ (not B+ -> J/psi K+) = %.2f + %.2f - %.2f\n", exp_bmm.getVal(), exp_bmm.getErrHi(),exp_bmm.getErrLo());
		}
	}
	
	fclose(outputFile);
} // writeULC()

void ncAna::showBsmmPlots()
{
	TFile *mcFile = TFile::Open(fSignalMCFileName.c_str());
	TFile *dataFile = TFile::Open(fDataFileName.c_str());
	TTree *mcTree = (TTree*)mcFile->Get("T");
	TTree *dataTree = (TTree*)dataFile->Get("T");
	set<ncCut> variables;
	set<ncCut>::const_iterator it;
	ncVarReader reader;
	TCut presel = cutSigCand() && cutAcceptanceData(0) && (cutTrigger(true, true) || cutTrigger(true, false)) && cutMuon() && cutHisto();
	TCut mcCut = cutSigTruth(true) && cutAcceptanceMC(0);
	TCanvas *c;
	
	try {
		reader.loadFile(fVarFileName.c_str());
	}
	catch (std::string err) {
		cout << Form("ncAna::showBsmmPlots(): error reading variable file. '%s'",err.c_str()) << endl;
		goto bail;;
	}
	
	if (reader.getNbr() == 0)
		goto bail;
	
	variables = *reader.getVars(0);
	
	// show the histograms.
	for (it = variables.begin(); it != variables.end(); ++it) {
		TH1D *histoMC = new TH1D(Form("%s_mc",it->getName()),"",50,it->getCut().first,it->getCut().second);
		TH1D *histoData = new TH1D(Form("%s_data",it->getName()),"",50,it->getCut().first,it->getCut().second);
		
		setHistoStyle(histoMC, kHistoStyle_Norm); // same as normalization
		
		dataTree->Draw(Form("%s >> %s_data",it->getFormula(),it->getName()), presel);
		mcTree->Draw(Form("%s >> %s_mc",it->getFormula(),it->getName()), presel && mcCut);
		
		c = new TCanvas;
		histoMC->Draw();
		histoData->Draw("sameE1");
	}
	
bail:
	delete mcFile;
	delete dataFile;
} // showVarPlots()

void ncAna::showNormPlots()
{
	// load the production mc
	TFile *mcFile = TFile::Open(fNormFileName.c_str());
	TFile *dataFile = TFile::Open(fDataFileName.c_str());
	TTree *mcTree = (TTree*)mcFile->Get("T");
	TTree *dataTree = (TTree*)dataFile->Get("T");
	TH1D *nPVData = new TH1D("nPVData","",50,0,40);
	TH1D *nPVMC = new TH1D("nPVMC","",50,0,40);
	TCanvas *c = new TCanvas("pvCan","Nbr of PV");
	TCut presel = cutNormCand() && cutAcceptanceData(1) && cutTrigger(false) && cutMuon() && cutHisto();
	TCut mcCut = cutNormTruth() && cutAcceptanceMC(1);
	int j;
	
	setHistoStyle(nPVMC, kHistoStyle_Norm);
	
	for (j = 0; j <= 1; j++) { // iterate through channels
		// FIXME: choose correct cut here
		dataTree->Draw("nbr_pv >> nPVData", presel && cutChannel(j));
		mcTree->Draw("nbr_pv >> nPVMC", presel && mcCut && cutChannel(j));
		c->cd();
		nPVMC->Draw();
		nPVData->Draw("sameE1");
		
		c->SaveAs(Form("plots/systematics-npv_NoData_NoMc_chan%d.pdf",j));
	}
	
	delete mcFile;
	delete dataFile;
} // showNormPlots()

void ncAna::setHistoStyle(TH1D *h, nc_histostyle style)
{
	h->SetLineColor(kBlue);
	h->SetFillColor(kBlue);
	h->SetFillStyle(3005);
} // setHistoStyle()
