#include <stddef.h>
#include <iostream>
#include <string>
#include <fstream>
#include "TFile.h"
#include "TROOT.h"
#include "mclimit_csm.hh"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TStopwatch.h"

using namespace std;

int main(int argc, char **argv)
{
	TStopwatch t;
	t.Start();

	delete gRandom;
	gRandom = (TRandom*) new TRandom3;
	gRandom->SetSeed(0);
	
	if (argc<4||argc>4) {
		cout << endl;
		cout << "Usage: bin/multibin_cls option1 option2 option3" << endl;
		cout << "Set option1 to bs or bd." << endl;
		cout << "Set option2 to chan0, chan1 or comb for barrel, encap or combined, respectively." << endl;
		cout << "Set option3 to any int (it will be used to set the random seed" << endl;
		cout << "NOTE: In TRandom3, SetSeed(0) is always random so the results are NOT reproducible." << endl;
		cout << endl;
		return 0;		
	}
	
	char bsmumu[] = "bs"; char bdmumu[] = "bd"; 
	char chann0[] = "chan0"; char chann1[] = "chan1"; char combd[] = "comb";
	bool bsmm = false; bool bdmm = false;
	bool barrel = false; bool endcap = false; bool combined = false;
	int i_numb = atoi(argv[3]);

	if (!strcmp(bsmumu,argv[1])) {
		bsmm = true;
	}
	else if (!strcmp(bdmumu,argv[1])) {
		bdmm = true;
	}
	if (!strcmp(chann0,argv[2])) {
		barrel = true;
	}
	else if (!strcmp(chann1,argv[2])) {
		endcap = true;
	}
	else if (!strcmp(combd,argv[2])) {
		barrel = true;
		endcap = true;
		combined = true;
	}

	cout << "Analyzing signal " << argv[1] << ", for channel " << argv[2] 
		 << ", with random seed " << i_numb << endl;

	TString s_outputfile;
	if (combined) s_outputfile = Form("CMS_Comb_Results%i.txt",i_numb);
	else if(barrel) s_outputfile = Form("CMS_brrl_Results%i.txt",i_numb);
	else if(endcap) s_outputfile = Form("CMS_endc_Results%i.txt",i_numb);
	
	ofstream f_outputfile;
	f_outputfile.open(s_outputfile);

#include "my_prepare.h"
  
	cout << "Back from my_prepare.h" << endl;
	
	//========================================
	// Have a visualization of how fitting looks like
	// It has nothing to do with limit calculation
	TCanvas * mycanvas; 
	if (combined) {
		mycanvas = (TCanvas *) new TCanvas("Canvas1","Canvas1",0,0,1000,800);
		mycanvas->Divide(2,2);
	}
	else {
		mycanvas = (TCanvas *) new TCanvas("Canvas1","Canvas1",0,0,600,800);
		mycanvas->Divide(1,2);
	}
	mycanvas->cd(1);
	
	cout << ">>>>>>> Begin calculating!!! <<<<<<<<" << endl;
	csm* mycsm = new csm();
	if (combined) {
		mycsm->set_htofit(h_dataB,channameB);
		mycsm->set_htofit(h_dataE,channameE);
		barrel = false;
		endcap = false;
	}
	if (barrel) mycsm->set_htofit(h_dataB,channameB);
	if (endcap) mycsm->set_htofit(h_dataE,channameE);
	mycsm->set_modeltofit(testhyp);
	double chisq = mycsm->chisquared();
	
	csm_model* bestnullfit = mycsm->getbestmodel();
	if (combined) {
		bestnullfit->plotwithdata(channameB,h_dataB);
		mycanvas->cd(2);
		bestnullfit->plotwithdata(channameE,h_dataE);
	}
	if (barrel) bestnullfit->plotwithdata(channameB,h_dataB);
	if (endcap) bestnullfit->plotwithdata(channameE,h_dataE);

	cout << "chisq from fitter " << chisq << endl; 
	delete mycsm;
	
	
	//======================================
	//Construct limit calculator - mclimit_csm
	//  Sensitivity, Significance calculation
	
	mclimit_csm* mymclimit = (mclimit_csm*) new mclimit_csm();		
	//print out pseudo-experiments details 
	//mymclimit->setpxprintflag(1);	
	cout << ">>>>>> setting up hypotheses <<<<<<<" << endl;
	mymclimit->set_null_hypothesis(nullhyp);
	mymclimit->set_test_hypothesis(testhyp);
	mymclimit->set_null_hypothesis_pe(nullhyp_pe);
	mymclimit->set_test_hypothesis_pe(testhyp_pe);
	if (combined) {
		mymclimit->set_datahist(h_dataB,channameB);
		mymclimit->set_datahist(h_dataE,channameE);
	}
	if (barrel) mymclimit->set_datahist(h_dataB,channameB);
	if (endcap) mymclimit->set_datahist(h_dataE,channameE);

	cout << ">>>>>>> PRINTING Pseudo-exp (debug purposes) <<<<<<<<" << endl;
	testhyp_pe->print();
	mymclimit->set_npe(10000);
	mymclimit->run_pseudoexperiments();
	
	
	cout << "<<<<<<<< Getting Results >>>>>>>>" << endl;
	Double_t tsobs = mymclimit->ts();// test statistics of observed data
	cout << "Test statistics of observed data --> " << "ts: " << tsobs << endl;
	cout << "Test statistic in null hypothesis (H0), 2sig low edge --> " << "tsbm2: "   << mymclimit->tsbm2() << endl;
	cout << "Test statistic in null hypothesis (H0), 1sig low edge --> " << "tsbm1: "   << mymclimit->tsbm1() << endl;
	cout << "Test statistic in null hypothesis (H0), medium --> "        << "tsbmed: "  << mymclimit->tsbmed() << endl;
	cout << "Test statistic in null hypothesis (H0), 1sig upper edge --> " << "tsbp1: " << mymclimit->tsbp1() << endl;
	cout << "Test statistic in null hypothesis (H0), 2sig upper edge --> " << "tsbp2: " << mymclimit->tsbp2() << endl;
	cout << "Test statistic in null hypothesis (H1), 2sig low edge --> "   << "tssm2: " << mymclimit->tssm2() << endl;
	cout << "Test statistic in null hypothesis (H1), 1sig low edge --> "   << "tssm1: " << mymclimit->tssm1() << endl;
	cout << "Test statistic in null hypothesis (H0), medium --> "         << "tssmed: " << mymclimit->tssmed() << endl;
	cout << "Test statistic in null hypothesis (H1), 1sig upper edge --> " << "tssp1: " << mymclimit->tssp1() << endl;
	cout << "Test statistic in null hypothesis (H1), 2sig upper edge --> " << "tssp2: " << mymclimit->tssp2() << endl;
	
	f_outputfile << "CLb = " << mymclimit->clb() << endl;
	f_outputfile << "CLsb = " << mymclimit->clsb() << endl;
	f_outputfile << "CLs (CLsb/CLb) = " << mymclimit->cls() << endl;

	cout << "<<<<<<<< Getting Expected CLs Results for Null Hyp >>>>>>>>" << endl;
	f_outputfile << "CLs -2sigma (bkg): " << mymclimit->clsexpbm2() << endl;
	f_outputfile << "CLs -1sigma (bkg): " << mymclimit->clsexpbm1() << endl;
	f_outputfile << "CLs median  (bkg): " << mymclimit->clsexpbmed() << endl;
	f_outputfile << "CLs +1sigma (bkg): " << mymclimit->clsexpbp1() << endl;
	f_outputfile << "CLs +2sigma (bkg): " << mymclimit->clsexpbp2() << endl;
	cout << "<<<<<<<< Could background only explain the observed data? >>>>>>>>" << endl;
	f_outputfile << "1-CLb -2sigma (sig): " << mymclimit->omclbexpsm2() << endl;
	f_outputfile << "1-CLb -1sigma (sig): " << mymclimit->omclbexpsm1() << endl;
	f_outputfile << "1-CLb median  (sig): " << mymclimit->omclbexpsmed() << endl;
	f_outputfile << "1-CLb +1sigma (sig): " << mymclimit->omclbexpsp1() << endl;
	f_outputfile << "1-CLb +2sigma (sig): " << mymclimit->omclbexpsp2() << endl;
	
	cout << "<<<<<<<< Getting Probabilities assuming test Hyp is TRUE >>>>>>>>" << endl;
	// Sensitivity of test hypothesis. Probability of a x-sigma evidence assuming test hyp. is true
	f_outputfile << "2-sigma prob (H1 is true): " << mymclimit->p2sigmat() << endl; 
	f_outputfile << "3-sigma prob (H1 is true): " << mymclimit->p3sigmat() << endl; 
	f_outputfile << "5-sigma prob (H1 is true): " << mymclimit->p5sigmat() << endl; 

	cout << "<<<<<<<< Getting results from Rate calculations >>>>>>>>" << endl;
	double d_sf95 = mymclimit->s95();
	f_outputfile << "Scale factor of 95% CL excluded signal: s95 = " << d_sf95 << endl;
	if (combined) {
		if (bsmm) f_outputfile << "95% CL upper limit (combined): Br(Bs->MuMu) = " << d_sf95*(3.2e-9) << endl;
		if (bdmm) f_outputfile << "95% CL upper limit (combined): Br(Bd->MuMu) = " << d_sf95*(1.0e-10) << endl;
	}
	if (barrel) {
		if (bsmm) f_outputfile << "95% CL upper limit (barrel): Br(Bs->MuMu) = " << d_sf95*(3.2e-9) << endl;
		if (bdmm) f_outputfile << "95% CL upper limit (barrel): Br(Bd->MuMu) = " << d_sf95*(1.0e-10) << endl;
	}
	if (endcap) {
		if (bsmm) f_outputfile << "95% CL upper limit (endcap): Br(Bs->MuMu) = " << d_sf95*(3.2e-9) << endl;
		if (bdmm) f_outputfile << "95% CL upper limit (endcap): Br(Bd->MuMu) = " << d_sf95*(1.0e-10) << endl;
	}
	
	f_outputfile << "s95 -2sigma (sig): " << mymclimit->s95m2() << endl;
	f_outputfile << "s95 -1sigma (sig): " << mymclimit->s95m1() << endl;
	f_outputfile << "s95 median  (sig): " << mymclimit->s95med() << endl;
	f_outputfile << "s95 +1sigma (sig): " << mymclimit->s95p1() << endl;
	f_outputfile << "s95 +2sigma (sig): " << mymclimit->s95p2() << endl;
		
	//Null hypothesis 
	f_outputfile << "Null Hypothesis:"
		 << " Median of test statistics = " << mymclimit->tsbmed() 
		 << ", p-value = " << mymclimit->omclbexpbmed()
		 << ", significance = "<< TMath::ErfcInverse(mymclimit->omclbexpbmed()*2)*sqrt(2)
		 << endl;
	
	
	//Expected significance
	f_outputfile << "Test Hypothesis:"
		 << " Median of test statistics = " << mymclimit->tssmed() 
		 << ", p-value = " << mymclimit->omclbexpsmed() 
		 << ", significance = " << TMath::ErfcInverse(mymclimit->omclbexpsmed()*2)*sqrt(2)
		 << endl;
	
	// Significance of input data 
	double ts_data = mymclimit->ts();
	double pval = mymclimit->omclb();
	double significance = TMath::ErfcInverse(pval*2)*sqrt(2); //convert pval to one-side gaussian significance
	
	f_outputfile << "Input Data: test statistics = " << ts_data << ", p-val = " << pval 
		 << ", significance = " << significance << endl;
	
//	cout << "<<<<<<<< Getting Bayes sf Results >>>>>>>>" << endl;
//	double conf_lv = 0.95;
//	double sf, sfE; 
//	mymclimit->bayes_heinrich(conf_lv, &sf, &sfE);
//	cout << "Bayes scale factor at 95% CL = " << sf << " +/- " << sfE << endl;
//	if (combined) {
//		if (bsmm) cout << "Bayes 95% CL upper limit (combined): Br(Bs->MuMu) = " << sf*(3.2e-9) << endl;
//		if (bdmm) cout << "Bayes 95% CL upper limit (combined): Br(Bd->MuMu) = " << sf*(1.0e-10) << endl;
//	}
//	if (barrel) {
//		if (bsmm) cout << "Bayes 95% CL upper limit (barrel): Br(Bs->MuMu) = " << sf*(3.2e-9) << endl;
//		if (bdmm) cout << "Bayes 95% CL upper limit (barrel): Br(Bd->MuMu) = " << sf*(1.0e-10) << endl;
//	}
//	if (endcap) {
//		if (bsmm) cout << "Bayes 95% CL upper limit (endcap): Br(Bs->MuMu) = " << sf*(3.2e-9) << endl;
//		if (bdmm) cout << "Bayes 95% CL upper limit (endcap): Br(Bd->MuMu) = " << sf*(1.0e-10) << endl;
//	}
//	
//	double sf90, sf90E;
//	mymclimit->bayes_heinrich(0.90, &sf90, &sf90E);
//	cout << "Bayes scale factor at 90% CL = " << sf90 << " +/- " << sf90E << endl;
//	if (combined) {
//		if (bsmm) cout << "Bayes 90% CL upper limit (combined): Br(Bs->MuMu) = " << sf90*(3.2e-9) << endl;
//		if (bdmm) cout << "Bayes 90% CL upper limit (combined): Br(Bd->MuMu) = " << sf90*(1.0e-10) << endl;
//	}
//	if (barrel) {
//		if (bsmm) cout << "Bayes 90% CL upper limit (barrel): Br(Bs->MuMu) = " << sf90*(3.2e-9) << endl;
//		if (bdmm) cout << "Bayes 90% CL upper limit (barrel): Br(Bd->MuMu) = " << sf90*(1.0e-10) << endl;
//	}
//	if (endcap) {
//		if (bsmm) cout << "Bayes 90% CL upper limit (endcap): Br(Bs->MuMu) = " << sf90*(3.2e-9) << endl;
//		if (bdmm) cout << "Bayes 90% CL upper limit (endcap): Br(Bd->MuMu) = " << sf90*(1.0e-10) << endl;
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
	mymclimit->tshists(ts_test,ts_null); 
	if (combined) mycanvas->cd(3);
	else mycanvas->cd(2);
	mycanvas->SetLogy(1);
	ts_null->SetLineColor(2);//red
	ts_test->SetLineColor(4);//blue
	ts_null->Draw();
	ts_test->Draw("same");
	
	TArrow* arrow = new TArrow(ts_data,0.3*ts_null->GetMaximum(),ts_data,0);
	arrow->Draw();
	
	if (combined) {
		if (bsmm) mycanvas->Print("multibin_clsBsmmComb.pdf");
		if (bdmm) mycanvas->Print("multibin_clsBdmmComb.pdf");
	}
	if (barrel) {
		if (bsmm) mycanvas->Print("multibin_clsBsmmBrrl.pdf");
		if (bdmm) mycanvas->Print("multibin_clsBdmmBrrl.pdf");
	}
	if (endcap) {
		if (bsmm) mycanvas->Print("multibin_clsBsmmEndc.pdf");
		if (bdmm) mycanvas->Print("multibin_clsBdmmEndc.pdf");
	}
	
	delete mymclimit;
	f_outfile.Write();
	f_outfile.Close();
	f_outputfile.close();

	t.Stop();
	t.Print();
}
