#include <stddef.h>
#include <iostream>
#include <string>
#include <fstream>
#include "TFile.h"
#include "TROOT.h"
#include "mclimit_csm.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TStopwatch.h"
// Modified on Oct 18, 2012
// Notes: Upgrading for the new analysis. (more to come)
//

using namespace std;

int main(int argc, char **argv)
{
	TStopwatch t;
	t.Start();

	delete gRandom;
	gRandom = (TRandom*) new TRandom3;
	
	if (argc<4||argc>4) {
		cout << endl;
		cout << "Usage: bin/multibin_cls option1 option2" << endl;
		cout << "Set option1 to: bsbarrel, bdbarrel, bsendcap, bdendcap, bscomb, or bdcomb" << endl; 
		cout << "Set option2 to any int (it will be used to set the random seed" << endl;
		cout << "NOTE: In TRandom3, SetSeed(0) is always random so the results are NOT reproducible." << endl;
		cout << "Option3 is a float used to scale the signal histogram." << endl;
		cout << "***NOTE: Should be set to 1.0 for normal calculations or range from 0.5-2.0 for CLs scans.***" << endl;
		cout << endl;
		return 0;		
	}
	
	char bschn0[] = "bsbarrel"; char bdchn0[] = "bdbarrel"; 
	char bschn1[] = "bsendcap"; char bdchn1[] = "bdendcap"; 
	char combs[] = "bscomb"; char combd[] = "bdcomb";
	bool bsmm = false; bool bdmm = false;
	bool barrel = false; bool endcap = false; 
	bool bscombined = false; bool bdcombined = false;

	if (!strcmp(bschn0,argv[1])) {
		bsmm = true;
		barrel = true;
	}
	else if (!strcmp(bdchn0,argv[1])) {
		bdmm = true;
		barrel = true;
	}
	else if (!strcmp(bschn1,argv[1])) {
		bsmm = true;
		endcap = true;
	}
	else if (!strcmp(bdchn1,argv[1])) {
		bdmm = true;
		endcap = true;
	}
	else if (!strcmp(combs,argv[1])) {
		bsmm = true;
		barrel = true;
		endcap = true;
		bscombined = true;
	}
	else if (!strcmp(combd,argv[1])) {
		bdmm = true;
		barrel = true;
		endcap = true;
		bdcombined = true;
	}
	
	int i_numb = atoi(argv[2]);
	double d_scale = atof(argv[3]);

	cout << "Analyzing signal " << argv[1] << ", for channel " << argv[2] 
		 << ", with random seed " << i_numb << endl;
	gRandom->SetSeed(i_numb);

	TString s_outputfile;
	if (bscombined) s_outputfile = Form("CMS_W12_Bs_Comb_ShpReslts%f_%i.txt",d_scale,i_numb);
	else if (bdcombined) s_outputfile = Form("CMS_W12_Bd_Comb_ShpReslts%f_%i.txt",d_scale,i_numb);
	else if(barrel && bsmm) s_outputfile = Form("CMS_W12_Bs_brrl_ShpReslts%f_%i.txt",d_scale,i_numb);
	else if(endcap && bsmm) s_outputfile = Form("CMS_W12_Bs_endc_ShpReslts%f_%i.txt",d_scale,i_numb);
	else if(barrel && bdmm) s_outputfile = Form("CMS_W12_Bd_brrl_ShpReslts%f_%i.txt",d_scale,i_numb);
	else if(endcap && bdmm) s_outputfile = Form("CMS_W12_Bd_endc_ShpReslts%f_%i.txt",d_scale,i_numb);
	
	ofstream f_outputfile;
	f_outputfile.open(s_outputfile);

//#include "my_prepare.h"
#include "prepshape.h"
  
	cout << "Back from my_prepare.h" << endl;
	
	//========================================
	// Have a visualization of how fitting looks like
	// It has nothing to do with limit calculation
	TCanvas * mycanvas; 
	if (bscombined || bdcombined) {
		mycanvas = (TCanvas *) new TCanvas("Canvas1","Canvas1",0,0,1000,800);
		mycanvas->Divide(2,2);
	}
	else {
		mycanvas = (TCanvas *) new TCanvas("Canvas1","Canvas1",0,0,500,800);
		mycanvas->Divide(1,2);
	}
	mycanvas->cd(1);
	
	cout << ">>>>>>> Begin calculating!!! <<<<<<<<" << endl;
	csm* mycsm = new csm();
	if (bscombined || bdcombined) {
		mycsm->set_htofit(h_dataB,channameB);
		mycsm->set_htofit(h_dataE,channameE);
	}
	else if (barrel) mycsm->set_htofit(h_dataB,channameB);
	else if (endcap) mycsm->set_htofit(h_dataE,channameE);
	mycsm->set_modeltofit(testhyp);
	double chisq = mycsm->chisquared();
	
	csm_model* bestnullfit = mycsm->getbestmodel();
	if (bscombined || bdcombined) {
		bestnullfit->plotwithdata(channameB,h_dataB);
		mycanvas->cd(2);
		bestnullfit->plotwithdata(channameE,h_dataE);
	}
	else if (barrel) bestnullfit->plotwithdata(channameB,h_dataB);
	else if (endcap) bestnullfit->plotwithdata(channameE,h_dataE);

	cout << "chisq from fitter " << chisq << endl; 
	delete mycsm;
	
	
	//======================================
	//Construct limit calculator - mclimit_csm
	//  Sensitivity, Significance calculation
	
	mclimit_csm* mymclimit = (mclimit_csm*) new mclimit_csm();		
	//print out pseudo-experiments details 
	mymclimit->setpxprintflag(0);	
	cout << ">>>>>> setting up hypotheses <<<<<<<" << endl;
	mymclimit->set_null_hypothesis(nullhyp);
	mymclimit->set_test_hypothesis(testhyp);
	mymclimit->set_null_hypothesis_pe(nullhyp_pe);
	mymclimit->set_test_hypothesis_pe(testhyp_pe);
	if (bscombined || bdcombined) {
		mymclimit->set_datahist(h_dataB,channameB);
		mymclimit->set_datahist(h_dataE,channameE);
	}
	else if (barrel) mymclimit->set_datahist(h_dataB,channameB);
	else if (endcap) mymclimit->set_datahist(h_dataE,channameE);

	cout << ">>>>>>> PRINTING Pseudo-exp (debug purposes) <<<<<<<<" << endl;
	testhyp_pe->print();
	mymclimit->set_npe(10000);
	mymclimit->run_pseudoexperiments();
	
	
	cout << "<<<<<<<< Getting Results >>>>>>>>" << endl;
	Double_t tsobs = mymclimit->ts(); Double_t d_tsbmed = mymclimit->tsbmed(); Double_t d_tssmed = mymclimit->tssmed();
	cout << "Test statistics of observed data --> " << "ts: " << tsobs << endl;
	cout << "Test statistic in null hypothesis (H0), 2sig low edge --> " << "tsbm2: "   << mymclimit->tsbm2() << endl;
	cout << "Test statistic in null hypothesis (H0), 1sig low edge --> " << "tsbm1: "   << mymclimit->tsbm1() << endl;
	cout << "Test statistic in null hypothesis (H0), medium --> "        << "tsbmed: "  << d_tsbmed << endl;
	cout << "Test statistic in null hypothesis (H0), 1sig upper edge --> " << "tsbp1: " << mymclimit->tsbp1() << endl;
	cout << "Test statistic in null hypothesis (H0), 2sig upper edge --> " << "tsbp2: " << mymclimit->tsbp2() << endl;
	cout << "Test statistic in test hypothesis (H1), 2sig low edge --> "   << "tssm2: " << mymclimit->tssm2() << endl;
	cout << "Test statistic in test hypothesis (H1), 1sig low edge --> "   << "tssm1: " << mymclimit->tssm1() << endl;
	cout << "Test statistic in test hypothesis (H1), medium --> "         << "tssmed: " << d_tssmed << endl;
	cout << "Test statistic in test hypothesis (H1), 1sig upper edge --> " << "tssp1: " << mymclimit->tssp1() << endl;
	cout << "Test statistic in test hypothesis (H1), 2sig upper edge --> " << "tssp2: " << mymclimit->tssp2() << endl;
	
	f_outputfile << "CLb = " << mymclimit->clb() << endl;
	f_outputfile << "CLsb = " << mymclimit->clsb() << endl;
	f_outputfile << "CLs (CLsb/CLb) = " << mymclimit->cls() << endl << endl;
	double pval = mymclimit->omclb();

	cout << "<<<<<<<< Getting Expected CLs Results for Null Hyp >>>>>>>>" << endl;
	Double_t d_omclbexpbmed = mymclimit->omclbexpbmed();
	f_outputfile << "Expected CLs in null hypothesis: -2 sigma = " << mymclimit->clsexpbm2() << endl;
	f_outputfile << "Expected CLs in null hypothesis: -1 sigma = " << mymclimit->clsexpbm1() << endl;
	f_outputfile << "Expected CLs in null hypothesis: median = " << mymclimit->clsexpbmed() << endl;
	f_outputfile << "Expected CLs in null hypothesis: +1 sigma = " << mymclimit->clsexpbp1() << endl;
	f_outputfile << "Expected CLs in null hypothesis: +2 sigma = " << mymclimit->clsexpbp2() << endl << endl;
	cout << "<<<<<<<< Could background only explain the observed data? >>>>>>>>" << endl;	
	f_outputfile << "Expected p-value in null hypothesis: -2 sigma = " << mymclimit->omclbexpbm2() << endl;
	f_outputfile << "Expected p-value in null hypothesis: -1 sigma = " << mymclimit->omclbexpbm1() << endl;
	f_outputfile << "Expected p-value in null hypothesis: median = " << d_omclbexpbmed << endl;
	f_outputfile << "Expected p-value in null hypothesis: +1 sigma = " << mymclimit->omclbexpbp1() << endl;
	f_outputfile << "Expected p-value in null hypothesis: +1 sigma = " << mymclimit->omclbexpbp2() << endl << endl;
	
	cout << "<<<<<<<< Getting Expected CLs Results for Test Hyp >>>>>>>>" << endl;
	Double_t d_omclbexpsmed = mymclimit->omclbexpsmed();
	f_outputfile << "Expected CLs in test hypothesis: -2 sigma = " << mymclimit->clsexpsm2() << endl;
	f_outputfile << "Expected CLs in test hypothesis: -1 sigma = " << mymclimit->clsexpsm1() << endl;
	f_outputfile << "Expected CLs in test hypothesis: median = " << mymclimit->clsexpsmed() << endl;
	f_outputfile << "Expected CLs in test hypothesis: +1 sigma = " << mymclimit->clsexpsp1() << endl;
	f_outputfile << "Expected CLs in test hypothesis: +2 sigma = " << mymclimit->clsexpsp2() << endl << endl;
	cout << "<<<<<<<< Could signal+background explain the observed data? >>>>>>>>" << endl;
	f_outputfile << "Expected p-value in test hypothesis: -2 sigma = " << mymclimit->omclbexpsm2() << endl;
	f_outputfile << "Expected p-value in test hypothesis: -1 sigma = " << mymclimit->omclbexpsm1() << endl;
	f_outputfile << "Expected p-value in test hypothesis: median = " << d_omclbexpsmed << endl;
	f_outputfile << "Expected p-value in test hypothesis: +1 sigma = " << mymclimit->omclbexpsp1() << endl;
	f_outputfile << "Expected p-value in test hypothesis: +2 sigma = " << mymclimit->omclbexpsp2() << endl << endl;
	
	cout << "<<<<<<<< Getting Probabilities assuming test Hyp is TRUE >>>>>>>>" << endl;
	// Sensitivity of test hypothesis. Probability of a x-sigma evidence assuming test hyp. is true
	f_outputfile << "Probability of a 2 sigma evidence assuming test hyp. is true: Prob = " << mymclimit->p2sigmat() << endl; 
	f_outputfile << "Probability of a 3 sigma evidence assuming test hyp. is true: Prob = " << mymclimit->p3sigmat() << endl; 
	f_outputfile << "Probability of a 5 sigma discovery assuming test hyp. is true: Prob = " << mymclimit->p5sigmat() << endl 
																													  << endl; 
	// Sensitivity of null hypothesis. Probability of a x-sigma evidence assuming null hyp. is true
	f_outputfile << "Probability of a 2 sigma evidence assuming null hyp. is true: Prob = " << mymclimit->p2sigman() << endl; 
	f_outputfile << "Probability of a 3 sigma evidence assuming null hyp. is true: Prob = " << mymclimit->p3sigman() << endl; 
	f_outputfile << "Probability of a 5 sigma discovery assuming null hyp. is true: Prob = " << mymclimit->p5sigman() << endl 
																													 << endl; 
	
	cout << "<<<<<<<< Calculating scale factors at 95% CL" << endl;
	double d_sf95 = mymclimit->s95(); double d_s95m2 = mymclimit->s95m2();
	double d_s95m1 = mymclimit->s95m1(); double d_s95med = mymclimit->s95med();
	double d_s95p1 = mymclimit->s95p1(); double d_s95p2 = mymclimit->s95p2();
	double d_lumi3s = mymclimit->lumi3s(); double d_lumi5s = mymclimit->lumi5s(); //double d_lumi95 = mymclimit->lumi95();
	cout << "<<<<<<<< Finshed calculating scale factors at 95% CL >>>>>>>>" << endl;
	f_outputfile  << "Observed scale factor of 95% CL excluded signal: sc_f = " << d_sf95 << endl;
	f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(-2sig) = " << d_s95m2 << endl;
	f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(-1sig) = " << d_s95m1 << endl;
	f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(median) = " << d_s95med << endl;
	f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(+1sig) = " << d_s95p1 << endl;
	f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(+2sig) = " << d_s95p2 << endl << endl;
	//	f_outputfile << "Lumi scale factor of 95% CL excluded signal: lumiSf = " << d_lumi95 << endl;
	f_outputfile << "Lumi needed for 3 sigma discovery: lumi3sig = " << d_lumi3s*cmslumiw12 << " fb-1" << endl;
	f_outputfile << "Lumi needed for 5 sigma discovery: lumi5sig = " << d_lumi5s*cmslumiw12 << " fb-1" << endl << endl;
	cout << "<<<<<<<< Finished calculating luminosity scale factors for discovery 95% CL >>>>>>>>" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "<<<<<<<< Calculating scale factors at 90% CL >>>>>>>>" << endl;
	double d_sf90 = mymclimit->s90(); double d_s90m2 = mymclimit->s90m2();
	double d_s90m1 = mymclimit->s90m1(); double d_s90med = mymclimit->s90med();
	double d_s90p1 = mymclimit->s90p1(); double d_s90p2 = mymclimit->s90p2();	
	cout << "<<<<<<<< Finshed calculating scale factors at 90% CL >>>>>>>>" << endl;
	f_outputfile  << "Observed scale factor of 90% CL excluded signal: sc_f = " << d_sf90 << endl;
	f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(-2sig) = " << d_s90m2 << endl;
	f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(-1sig) = " << d_s90m1 << endl;
	f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(median) = " << d_s90med << endl;
	f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(+1sig) = " << d_s90p1 << endl;
	f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(+2sig) = " << d_s90p2 << endl << endl;
	
	if (bscombined || bdcombined) {
		if (bsmm) {
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (combined): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;

			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (combined): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		if (bdmm) {
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (combined): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (combined): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
	}
	else if (barrel) {
		if (bsmm) {
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (barrel): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (barrel): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		if (bdmm) {
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (barrel): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (barrel): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
	}
	else if (endcap) {
		if (bsmm) {
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (endcap): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (endcap): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		if (bdmm) {
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (endcap): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (endcap): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
	}
		
	//Null hypothesis 
	f_outputfile << "Null Hypothesis: Median of test statistics = " << d_tsbmed 
				 << ", p-value = " << d_omclbexpbmed << ", significance = " 
		         << TMath::ErfcInverse(d_omclbexpbmed*2)*sqrt(2) << endl;
	
	
	//Expected significance
	f_outputfile << "Test Hypothesis: Median of test statistics = " << d_tssmed 
				 << ", p-value = " << d_omclbexpsmed << ", significance = " 
				 << TMath::ErfcInverse(d_omclbexpsmed*2)*sqrt(2) << endl;
	
	// Significance of input data 
	double significance = TMath::ErfcInverse(pval*2)*sqrt(2); //convert pval to one-side gaussian significance	
	f_outputfile << "Input Data: test statistics = " << tsobs << ", p-value = " << pval 
				 << ", significance = " << significance << endl;
	
	cout << "<<<<<<<< Getting Bayes sf Results >>>>>>>>" << endl;
//	double conf_lv = 0.95;
//	double sf, sfE; 
//	mymclimit->bayes_heinrich(conf_lv, &sf, &sfE);
//	cout << "Bayes scale factor at 95% CL = " << sf << " +/- " << sfE << endl;
//	if (bscombined || bdcombined) {
//		if (bsmm) cout << "Bayes 95% CL upper limit (combined): Br(Bs->MuMu) = " << sf*Brbsmm << endl;
//		if (bdmm) cout << "Bayes 95% CL upper limit (combined): Br(Bd->MuMu) = " << sf*Brbdmm << endl;
//	}
//	else if (barrel) {
//		if (bsmm) cout << "Bayes 95% CL upper limit (barrel): Br(Bs->MuMu) = " << sf*Brbsmm << endl;
//		if (bdmm) cout << "Bayes 95% CL upper limit (barrel): Br(Bd->MuMu) = " << sf*Brbdmm << endl;
//	}
//	else if (endcap) {
//		if (bsmm) cout << "Bayes 95% CL upper limit (endcap): Br(Bs->MuMu) = " << sf*Brbsmm << endl;
//		if (bdmm) cout << "Bayes 95% CL upper limit (endcap): Br(Bd->MuMu) = " << sf*Brbdmm << endl;
//	}
//	
//	double sf90, sf90E;
//	mymclimit->bayes_heinrich(0.90, &sf90, &sf90E);
//	cout << "Bayes scale factor at 90% CL = " << sf90 << " +/- " << sf90E << endl;
//	if (bscombined || bdcombined) {
//		if (bsmm) cout << "Bayes 90% CL upper limit (combined): Br(Bs->MuMu) = " << sf90*Brbsmm << endl;
//		if (bdmm) cout << "Bayes 90% CL upper limit (combined): Br(Bd->MuMu) = " << sf90*Brbdmm << endl;
//	}
//	else if (barrel) {
//		if (bsmm) cout << "Bayes 90% CL upper limit (barrel): Br(Bs->MuMu) = " << sf90*Brbsmm << endl;
//		if (bdmm) cout << "Bayes 90% CL upper limit (barrel): Br(Bd->MuMu) = " << sf90*Brbdmm << endl;
//	}
//	else if (endcap) {
//		if (bsmm) cout << "Bayes 90% CL upper limit (endcap): Br(Bs->MuMu) = " << sf90*Brbsmm << endl;
//		if (bdmm) cout << "Bayes 90% CL upper limit (endcap): Br(Bd->MuMu) = " << sf90*Brbdmm << endl;
//	}
//	
//	cout << "<<<<<<<< Getting Bayes with expect sf Results >>>>>>>>" << endl;
//	double sm2, sm1, smed, sp1, sp2; 
//	int npx(10000); 
//	mymclimit->bayes_heinrich_withexpect(conf_lv, &sf, &sfE, npx, &sm2, &sm1, &smed, &sp1, &sp2);
//	cout << "Bayes scale factor (sf) = " << sf << " +/- " << sfE  << endl;
//	cout << "-2sig (exp limit) = " << sm2 << ", -1sig (exp limit) = " << sm1 << endl;
//	cout << "Medium (exp limit) = " << smed << endl; 
//	cout << "+1sig (exp limit) = " << sp1 << ", +2sig (exp limit) = " << sp2 << endl;
	
	//==============================================
	//Draw the test statistics distribution 
	//test hypothesis - we can calculate probability of observering more than N sigma 
	//null hypothesis - we can calculate p-value of a measured test statistics  
	TH1F* ts_test = new TH1F("ts_test","",100,-50,50);
	TH1F* ts_null = new TH1F("ts_null","",100,-50,50);
	////////////////////////////////////////////////////////////////////////
	TH1F* mcb_hist = new TH1F("mcb_hist","",nbins,0,2);
	TH1F* mcs_hist = new TH1F("mcs_hist","",nbins,0,2);
	TH1F* data_hist = new TH1F("data_hist","",nbins,0,2);
	////////////////////////////////////////////////////////////////////////
	mymclimit->tshists(ts_test,ts_null); 
	if (bscombined || bdcombined) {
		mycanvas->cd(3);
		mymclimit->plotlnsb(mcb_hist,mcs_hist,data_hist);
		mycanvas->cd(4);
	}
	else mycanvas->cd(2);
	mycanvas->SetLogy(1);
	ts_null->SetLineColor(2);//red
	ts_test->SetLineColor(4);//blue
	ts_null->Draw();
	ts_test->Draw("same");
	
	TArrow* arrow = new TArrow(tsobs,0.3*ts_null->GetMaximum(),tsobs,0);
	arrow->Draw();
	
	TString s_pdffilename;
	if (bscombined ) s_pdffilename = Form("CMS_W12_Bs_Comb_ShpReslts%f_%i.pdf",d_scale,i_numb);
	else if (bdcombined ) s_pdffilename = Form("CMS_W12_Bd_Comb_ShpReslts%f_%i.pdf",d_scale,i_numb);
	else if (bsmm && barrel) s_pdffilename = Form("CMS_W12_Bs_brrl_ShpReslts%f_%i.pdf",d_scale,i_numb);
	else if (bsmm && endcap) s_pdffilename = Form("CMS_W12_Bs_endc_ShpReslts%f_%i.pdf",d_scale,i_numb);
	else if (bdmm && barrel) s_pdffilename = Form("CMS_W12_Bd_brrl_ShpReslts%f_%i.pdf",d_scale,i_numb);
	else if (bdmm && endcap) s_pdffilename = Form("CMS_W12_Bd_endc_ShpReslts%f_%i.pdf",d_scale,i_numb);
	mycanvas->Print(s_pdffilename);
	
	delete mymclimit;
	f_outfile.Write();
	f_outfile.Close();
	f_outputfile.close();

	t.Stop();
	t.Print();
}
