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
#include "TString.h"

#define N_MASS_BINS 6
#define N_BDT_BINS 4
#define N_MASS_BINS12 9
#define N_BDT_BINS12 8

using namespace std;

int main(int argc, char **argv)
{
	TStopwatch t;
	t.Start();
	if (argc<5||argc>5) {
		cout << endl;
		cout << "Usage: bin/upperlimitCalc option1 option2 option3 option4" << endl;
		cout << "Option1 can be:" << endl;
		cout << "---> cms11bs or cmb11bd for summer 11 Bs or Bd" << endl; 
		cout << "---> cms12bs or cms12bd for Spring 12 Bs or Bd" << endl; 
		cout << "---> lhcb the combined 10 and 11 lhcb Bs results" << endl;
		cout << "---> lhcb10, lhcb11 or lhcb12 for either bs data set" << endl;
		cout << "---> lhcb11bd, lhcb12bd for lhcb Bd 11 or 12 results" << endl;
		cout << "---> comb11 or comb12 for CMS-LHCb combined bs result" << endl;
		cout << "---> comb11bd or comb12bd for CMS-LHCb combined Bd result" << endl;
		cout << "---> atlas12bs or atlas12bd for winter Bs or Bd" << endl; 
		cout << "Option2 is a number (int) that is used to set the random seed." << endl;
		cout << "NOTE: In TRandom3, SetSeed(0) is always random so the results are NOT reproducible." << endl;
		cout << "Option3 is a float used to scale the signal histogram." << endl;
		cout << "***NOTE: Should be set to 1.0 for normal calculations or range from 0.5-2.0 for Bs CLs scans.***" << endl;
		cout << "***NOTE: Should be set to 1.0 for normal calculations or range from 0.5-80.0 for Bd CLs scans.***" << endl;
		cout << "Option4 can be: " << endl;
		cout << "95 or 90 and is used to set the confidence level for the calculation" << endl;
		cout << endl;
		return 0;		
	}
	
	delete gRandom;
	gRandom = (TRandom*) new TRandom3;
	
	char cms11bschan[] = "cms11bs"; char cms11bdchan[] = "cms11bd";
	char cms12bschan[] = "cms12bs"; char cms12bdchan[] = "cms12bd"; 
	char lhcbchan[] = "lhcb"; char lhcb10chan[] = "lhcb10"; 
	char lhcb11chan[] = "lhcb11"; char lhcb11bdchan[] = "lhcb11bd";
	char lhcb12chan[] = "lhcb12"; char lhcb12bdchan[] = "lhcb12bd";
	char combned11[] = "comb11"; char combned11bd[] = "comb11bd";
	char combned12[] = "comb12"; char combned12bd[] = "comb12bd";
	char atlas12bschan[] = "atlas12bs";
	bool combined11 = false; bool combined11bd = false;
	bool combined12 = false; bool combined12bd = false;
	bool cms11bs = false; bool cms11bd = false;
	bool cms12bs = false; bool cms12bd = false;
	bool lhcbs = false; bool lhcbs10 = false; 
	bool lhcbs11 = false; bool lhcbd11 = false;
	bool lhcbs12 = false; bool lhcbd12 = false;
	bool atlasbs12 = false;
	bool barrel = false; bool endcap = false; 
	char conf_level95[] = "95"; char conf_level90[] = "90"; 
	bool conflevel95 = false; bool conflevel90 = false; 

	if (!strcmp(cms11bschan,argv[1])) {
		cms11bs = true;
		barrel = true;
		endcap = true;
	}
	else if (!strcmp(cms11bdchan,argv[1])) {
		cms11bd = true;
		barrel = true;
		endcap = true;
	}
	else if (!strcmp(cms12bschan,argv[1])) {
		cms12bs = true;
//		lhcbs12 = true;
		barrel = true;
		endcap = true;
	}	
	else if (!strcmp(cms12bdchan,argv[1])) {
		cms12bd = true;
		barrel = true;
		endcap = true;
	}	
	else if (!strcmp(lhcbchan,argv[1])) {
		lhcbs = true;
		lhcbs10 = true;
		lhcbs11 = true;
	}
	else if	(!strcmp(combned11,argv[1])) {
		combined11 = true;
		cms11bs = true;
		barrel = true;
		endcap = true;
		
		lhcbs = true;
		lhcbs10 = true;
		lhcbs11 = true;

	}
///////////////////
	else if	(!strcmp(combned11bd,argv[1])) {
		combined11bd = true;
		cms11bd = true;
		barrel = true;
		endcap = true;
		
//		lhcbd = true;
//		lhcbd10 = true;
		lhcbd11 = true;
		
	}
//////////////////
	else if	(!strcmp(combned12,argv[1])) {
		combined12 = true;
		cms12bs = true;
		barrel = true;
		endcap = true;
		
		lhcbs = true;
		//lhcbs10 = true;
		//lhcbs11 = true;
		lhcbs12 = true;
		atlasbs12 = true;

	}
	else if	(!strcmp(combned12bd,argv[1])) {
		combined12bd = true;
		cms12bd = true;
		barrel = true;
		endcap = true;
		
		lhcbd12 = true;
//		atlasbd12 = true;
		
	}
	else if (!strcmp(lhcb10chan,argv[1])) {
		lhcbs10 = true;
		
	}
	else if (!strcmp(lhcb11chan,argv[1])) {
		lhcbs11 = true;
		
	}	
	else if (!strcmp(lhcb11bdchan,argv[1])) {
		lhcbd11 = true;
		
	}	
	else if (!strcmp(lhcb12chan,argv[1])) {
		lhcbs12 = true;
		
	}	
	else if (!strcmp(lhcb12bdchan,argv[1])) {
		lhcbd12 = true;
		
	}		
	else if (!strcmp(atlas12bschan,argv[1])) {
		atlasbs12 = true;
		
	}	
	else {
		cout << "Did not recognise your option, Please try again." << endl;
		return 0;		
	}

	
	int i_numb = atoi(argv[2]);

	double d_scale = atof(argv[3]);
		
	if (!strcmp(conf_level95,argv[4])) {
		conflevel95 = true;
	}
	else if (!strcmp(conf_level90,argv[4])) {
		conflevel90 = true;
	}
	gRandom->SetSeed(i_numb);

	TString s_outputfile;
	if (combined11) s_outputfile = Form("CMS_LHCb_S11_Comb_Bs_Results%f_%i.txt",d_scale,i_numb);
	else if (combined12) s_outputfile = Form("CMS_LHCb_W12_Comb_Bs_Results%f_%i.txt",d_scale,i_numb);
	else if (combined12bd) s_outputfile = Form("CMS_LHCb_W12_Comb_Bd_Results%f_%i.txt",d_scale,i_numb);
	else if(cms11bs) s_outputfile = Form("CMS_S11_Bs_Results%f_%i.txt",d_scale,i_numb);
	else if(cms11bd) s_outputfile = Form("CMS_S11_Bd_Results%f_%i.txt",d_scale,i_numb);
	else if(cms12bs) s_outputfile = Form("CMS_W12_Bs_Results%f_%i.txt",d_scale,i_numb);
	else if(cms12bd) s_outputfile = Form("CMS_W12_Bd_Results%f_%i.txt",d_scale,i_numb);
	else if(lhcbs) s_outputfile = Form("LHCb_10_11_Comb_Bs_Results%f_%i.txt",d_scale,i_numb);
	else if(lhcbs10) s_outputfile = Form("LHCb_10_Bs_Results%f_%i.txt",d_scale,i_numb);
	else if(lhcbs11) s_outputfile = Form("LHCb_11_Bs_Results%f_%i.txt",d_scale,i_numb);
	else if(lhcbd11) s_outputfile = Form("LHCb_11_Bd_Results%f_%i.txt",d_scale,i_numb);
	else if(lhcbs12) s_outputfile = Form("LHCb_12_Bs_Results%f_%i.txt",d_scale,i_numb);
	else if(lhcbd12) s_outputfile = Form("LHCb_12_Bd_Results%f_%i.txt",d_scale,i_numb);
	else if(atlasbs12) s_outputfile = Form("ATLAS_12_Bs_Results%f_%i.txt",d_scale,i_numb);
	
	ofstream f_outputfile;
	f_outputfile.open(s_outputfile);
	
	f_outputfile << "Setting the random seed to = " << i_numb << endl;
	
#include "prepUprLimCalc.h"
  
	cout << "Back from the prepare header file" << endl;
	
	//========================================
	// Have a visualization of how fitting looks like
	// It has nothing to do with limit calculation
	TCanvas *mycanvas = (TCanvas *) new TCanvas("Canvas1","Canvas1",0,0,600,600);
	if ((cms11bs || cms11bd || cms12bs || cms12bd) && !(combined11 || combined11bd || combined12 || combined12bd)) {
		mycanvas->Divide(2,2);
		mycanvas->cd(1);
	}
		
	cout << ">>>>>>> Begin calculating!!! <<<<<<<<" << endl;
	csm* mycsm = new csm();
	if (cms11bs || cms11bd || cms12bs || cms12bd) {
		mycsm->set_htofit(h_dataB,channameB);
		mycsm->set_htofit(h_dataE,channameE);
	}
	if (lhcbs10) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {				
				char * c_lhcbchan1 = new char[s_lhcbchan1[j][l].size() + 1];
				copy(s_lhcbchan1[j][l].begin(), s_lhcbchan1[j][l].end(), c_lhcbchan1);
				c_lhcbchan1[s_lhcbchan1[j][l].size()] = '\0'; 			
				mycsm->set_htofit(h_lhcbdat1[j][l],c_lhcbchan1);
				
				delete[] c_lhcbchan1;
			}
		}
	}
	if (lhcbs11) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {
				char * c_lhcbchan = new char[s_lhcbchan[j][l].size() + 1];
				copy(s_lhcbchan[j][l].begin(), s_lhcbchan[j][l].end(), c_lhcbchan);
				c_lhcbchan[s_lhcbchan[j][l].size()] = '\0'; 
				mycsm->set_htofit(h_lhcbdat[j][l],c_lhcbchan);				
				
				delete[] c_lhcbchan;
			}
		}
	}
	if (lhcbd11) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {
				char * c_lhcbchan4 = new char[s_lhcbchan4[j][l].size() + 1];
				copy(s_lhcbchan4[j][l].begin(), s_lhcbchan4[j][l].end(), c_lhcbchan4);
				c_lhcbchan4[s_lhcbchan4[j][l].size()] = '\0'; 
				mycsm->set_htofit(h_lhcbdat4[j][l],c_lhcbchan4);				
				
				delete[] c_lhcbchan4;
			}
		}
	}
	if (lhcbs12) {
		for (int j = 0; j<N_MASS_BINS12; j++) {
			for ( int l = 0; l<N_BDT_BINS12; l++) {				
				char * c_lhcbchan2 = new char[s_lhcbchan2[j][l].size() + 1];
				copy(s_lhcbchan2[j][l].begin(), s_lhcbchan2[j][l].end(), c_lhcbchan2);
				c_lhcbchan2[s_lhcbchan2[j][l].size()] = '\0'; 			
				mycsm->set_htofit(h_lhcbdat2[j][l],c_lhcbchan2);
				
				delete[] c_lhcbchan2;
			}
		}
	}
	if (lhcbd12) {
		for (int j = 0; j<N_MASS_BINS12; j++) {
			for ( int l = 0; l<N_BDT_BINS12; l++) {				
				char * c_lhcbchan3 = new char[s_lhcbchan3[j][l].size() + 1];
				copy(s_lhcbchan3[j][l].begin(), s_lhcbchan3[j][l].end(), c_lhcbchan3);
				c_lhcbchan3[s_lhcbchan3[j][l].size()] = '\0'; 			
				mycsm->set_htofit(h_lhcbdat3[j][l],c_lhcbchan3);
				
				delete[] c_lhcbchan3;
			}
		}
	}
	if (atlasbs12) {
		mycsm->set_htofit(h_atlasdata1,atchnname1);
		mycsm->set_htofit(h_atlasdata2,atchnname2);
		mycsm->set_htofit(h_atlasdata3,atchnname3);
	}
	
	mycsm->set_modeltofit(testhyp);
	double chisq = mycsm->chisquared();
	
	csm_model* bestnullfit = mycsm->getbestmodel();
	if (cms11bs || cms11bd || cms12bs || cms12bd) { //////////////
		bestnullfit->plotwithdata(channameB,h_dataB);
		if ((cms11bs || cms11bd || cms12bs || cms12bd) && !(combined11 || combined11bd || combined12 || combined12bd)) mycanvas->cd(2);//
		bestnullfit->plotwithdata(channameE,h_dataE);
	}
	if (lhcbs10) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {
				char * c_lhcbchan1 = new char[s_lhcbchan1[j][l].size() + 1];
				copy(s_lhcbchan1[j][l].begin(), s_lhcbchan1[j][l].end(), c_lhcbchan1);
				c_lhcbchan1[s_lhcbchan1[j][l].size()] = '\0'; 			
				bestnullfit->plotwithdata(c_lhcbchan1,h_lhcbdat1[j][l]);
				
				delete[] c_lhcbchan1;
			}
		}
	}
	if (lhcbs11) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {
				char * c_lhcbchan = new char[s_lhcbchan[j][l].size() + 1];
				copy(s_lhcbchan[j][l].begin(), s_lhcbchan[j][l].end(), c_lhcbchan);
				c_lhcbchan[s_lhcbchan[j][l].size()] = '\0'; 
				bestnullfit->plotwithdata(c_lhcbchan,h_lhcbdat[j][l]);
								
				delete[] c_lhcbchan;
			}
		}
	}
	if (lhcbd11) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {
				char * c_lhcbchan4 = new char[s_lhcbchan4[j][l].size() + 1];
				copy(s_lhcbchan4[j][l].begin(), s_lhcbchan4[j][l].end(), c_lhcbchan4);
				c_lhcbchan4[s_lhcbchan4[j][l].size()] = '\0'; 
				bestnullfit->plotwithdata(c_lhcbchan4,h_lhcbdat4[j][l]);
				
				delete[] c_lhcbchan4;
			}
		}
	}
	if (lhcbs12) {
		for (int j = 0; j<N_MASS_BINS12; j++) {
			for ( int l = 0; l<N_BDT_BINS12; l++) {
				char * c_lhcbchan2 = new char[s_lhcbchan2[j][l].size() + 1];
				copy(s_lhcbchan2[j][l].begin(), s_lhcbchan2[j][l].end(), c_lhcbchan2);
				c_lhcbchan2[s_lhcbchan2[j][l].size()] = '\0'; 			
				bestnullfit->plotwithdata(c_lhcbchan2,h_lhcbdat2[j][l]);
				
				delete[] c_lhcbchan2;
			}
		}
	}
	if (lhcbd12) {
		for (int j = 0; j<N_MASS_BINS12; j++) {
			for ( int l = 0; l<N_BDT_BINS12; l++) {
				char * c_lhcbchan3 = new char[s_lhcbchan3[j][l].size() + 1];
				copy(s_lhcbchan3[j][l].begin(), s_lhcbchan3[j][l].end(), c_lhcbchan3);
				c_lhcbchan3[s_lhcbchan3[j][l].size()] = '\0'; 			
				bestnullfit->plotwithdata(c_lhcbchan3,h_lhcbdat3[j][l]);
				
				delete[] c_lhcbchan3;
			}
		}
	}
	if (atlasbs12) { 
		bestnullfit->plotwithdata(atchnname1,h_atlasdata1);
		bestnullfit->plotwithdata(atchnname2,h_atlasdata2);
		bestnullfit->plotwithdata(atchnname3,h_atlasdata3);
	}
	
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
	if (cms11bs || cms11bd || cms12bs || cms12bd) { ////////////////
		mymclimit->set_datahist(h_dataB,channameB);
		mymclimit->set_datahist(h_dataE,channameE);
	}
	if (lhcbs10) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {
				char * c_lhcbchan1 = new char[s_lhcbchan1[j][l].size() + 1];
				copy(s_lhcbchan1[j][l].begin(), s_lhcbchan1[j][l].end(), c_lhcbchan1);
				c_lhcbchan1[s_lhcbchan1[j][l].size()] = '\0'; 			
				mymclimit->set_datahist(h_lhcbdat1[j][l],c_lhcbchan1);
				
				delete[] c_lhcbchan1;
			}
		}
	}
	if (lhcbs11) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {
				char * c_lhcbchan = new char[s_lhcbchan[j][l].size() + 1];
				copy(s_lhcbchan[j][l].begin(), s_lhcbchan[j][l].end(), c_lhcbchan);
				c_lhcbchan[s_lhcbchan[j][l].size()] = '\0'; 
				mymclimit->set_datahist(h_lhcbdat[j][l],c_lhcbchan);
				
				delete[] c_lhcbchan;
			}
		}
	}
	if (lhcbd11) {
		for (int j = 0; j<N_MASS_BINS; j++) {
			for ( int l = 0; l<N_BDT_BINS; l++) {
				char * c_lhcbchan4 = new char[s_lhcbchan4[j][l].size() + 1];
				copy(s_lhcbchan4[j][l].begin(), s_lhcbchan4[j][l].end(), c_lhcbchan4);
				c_lhcbchan4[s_lhcbchan4[j][l].size()] = '\0'; 
				mymclimit->set_datahist(h_lhcbdat4[j][l],c_lhcbchan4);
				
				delete[] c_lhcbchan4;
			}
		}
	}
	if (lhcbs12) {
		for (int j = 0; j<N_MASS_BINS12; j++) {
			for ( int l = 0; l<N_BDT_BINS12; l++) {
				char * c_lhcbchan2 = new char[s_lhcbchan2[j][l].size() + 1];
				copy(s_lhcbchan2[j][l].begin(), s_lhcbchan2[j][l].end(), c_lhcbchan2);
				c_lhcbchan2[s_lhcbchan2[j][l].size()] = '\0'; 			
				mymclimit->set_datahist(h_lhcbdat2[j][l],c_lhcbchan2);
				
				delete[] c_lhcbchan2;
			}
		}
	}
	if (lhcbd12) {
		for (int j = 0; j<N_MASS_BINS12; j++) {
			for ( int l = 0; l<N_BDT_BINS12; l++) {
				char * c_lhcbchan3 = new char[s_lhcbchan3[j][l].size() + 1];
				copy(s_lhcbchan3[j][l].begin(), s_lhcbchan3[j][l].end(), c_lhcbchan3);
				c_lhcbchan3[s_lhcbchan3[j][l].size()] = '\0'; 			
				mymclimit->set_datahist(h_lhcbdat3[j][l],c_lhcbchan3);
				
				delete[] c_lhcbchan3;
			}
		}
	}
	if (atlasbs12) { ////////////////
		mymclimit->set_datahist(h_atlasdata1,atchnname1);
		mymclimit->set_datahist(h_atlasdata2,atchnname2);
		mymclimit->set_datahist(h_atlasdata3,atchnname3);
	}
	
	cout << ">>>>>>> PRINTING Pseudo-exp (debug purposes) <<<<<<<<" << endl;
	testhyp_pe->print();
	mymclimit->set_npe(50000);
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
	
	f_outputfile << "The signal scale factor for 95% calculation is = " << d_scale << endl;
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
	f_outputfile << "Expected 1-CLb in null hypothesis: -2 sigma = " << mymclimit->omclbexpbm2() << endl;
	f_outputfile << "Expected 1-CLb in null hypothesis: -1 sigma = " << mymclimit->omclbexpbm1() << endl;
	f_outputfile << "Expected 1-CLb in null hypothesis: median = " << d_omclbexpbmed << endl;
	f_outputfile << "Expected 1-CLb in null hypothesis: +1 sigma = " << mymclimit->omclbexpbp1() << endl;
	f_outputfile << "Expected 1-CLb in null hypothesis: +1 sigma = " << mymclimit->omclbexpbp2() << endl << endl;
	cout << "<<<<<<<< Could signal+background explain the observed data? >>>>>>>>" << endl;	
	f_outputfile << "Expected CLsb in null hypothesis: -2 sigma = " << mymclimit->clsbexpbm2() << endl;
	f_outputfile << "Expected CLsb in null hypothesis: -1 sigma = " << mymclimit->clsbexpbm1() << endl;
	f_outputfile << "Expected CLsb in null hypothesis: median = " << mymclimit->clsbexpbmed() << endl;
	f_outputfile << "Expected CLsb in null hypothesis: +1 sigma = " << mymclimit->clsbexpbp1() << endl;
	f_outputfile << "Expected CLsb in null hypothesis: +2 sigma = " << mymclimit->clsbexpbp2() << endl << endl;

	cout << "<<<<<<<< Getting Expected CLs Results for Test Hyp >>>>>>>>" << endl;
	Double_t d_omclbexpsmed = mymclimit->omclbexpsmed();
	f_outputfile << "Expected CLs in test hypothesis: -2 sigma = " << mymclimit->clsexpsm2() << endl;
	f_outputfile << "Expected CLs in test hypothesis: -1 sigma = " << mymclimit->clsexpsm1() << endl;
	f_outputfile << "Expected CLs in test hypothesis: median = " << mymclimit->clsexpsmed() << endl;
	f_outputfile << "Expected CLs in test hypothesis: +1 sigma = " << mymclimit->clsexpsp1() << endl;
	f_outputfile << "Expected CLs in test hypothesis: +2 sigma = " << mymclimit->clsexpsp2() << endl << endl;
	cout << "<<<<<<<< Could signal+background explain the observed data? >>>>>>>>" << endl;
	f_outputfile << "Expected 1-CLb in test hypothesis: -2 sigma = " << mymclimit->omclbexpsm2() << endl;
	f_outputfile << "Expected 1-CLb in test hypothesis: -1 sigma = " << mymclimit->omclbexpsm1() << endl;
	f_outputfile << "Expected 1-CLb in test hypothesis: median = " << d_omclbexpsmed << endl;
	f_outputfile << "Expected 1-CLb in test hypothesis: +1 sigma = " << mymclimit->omclbexpsp1() << endl;
	f_outputfile << "Expected 1-CLb in test hypothesis: +2 sigma = " << mymclimit->omclbexpsp2() << endl << endl;
	
	cout << "<<<<<<<< Getting Probabilities assuming test Hyp is TRUE >>>>>>>>" << endl;
	// Sensitivity of test hypothesis. Probability of a x-sigma evidence assuming test hyp. is true
	f_outputfile << "Probability of a 2 sigma evidence assuming test hyp. is true: Prob = " << mymclimit->p2sigmat() << endl; 
	f_outputfile << "Probability of a 3 sigma evidence assuming test hyp. is true: Prob = " << mymclimit->p3sigmat() << endl; 
	f_outputfile << "Probability of a 5 sigma discovery assuming test hyp. is true: Prob = " << mymclimit->p5sigmat() << endl 
	<< endl; 
		
	double d_totlumi = 0.0; 
	if (combined11) d_totlumi = cmslumis11 + lhcblumis10 + lhcblumis11;
	else if (combined12)  d_totlumi = cmslumiw12 + lhcblumiw12;
	else if (combined12bd)  d_totlumi = cmslumiw12 + lhcblumiw12;
	else if (cms11bs || cms11bd) d_totlumi = cmslumis11;
	else if (cms12bs || cms12bd) d_totlumi = cmslumiw12;
	else if (lhcbs) d_totlumi = lhcblumis10 + lhcblumis11;
	else if (lhcbs10) d_totlumi = lhcblumis10;
	else if (lhcbs11) d_totlumi = lhcblumis11;
	else if (lhcbs12) d_totlumi = lhcblumiw12 + lhcblumis10;
	else if (lhcbd12) d_totlumi = lhcblumiw12;
	else if (atlasbs12) d_totlumi = atlaslumiw12;
	
	if (conflevel95) {
		cout << "<<<<<<<< Calculating scale factors at 95% CL" << endl;
		double d_sf95 = mymclimit->s95(); double d_s95m2 = mymclimit->s95m2();
		double d_s95m1 = mymclimit->s95m1(); double d_s95med = mymclimit->s95med();
		double d_s95p1 = mymclimit->s95p1(); double d_s95p2 = mymclimit->s95p2();	
		double d_lumi3s = mymclimit->lumi3s(); double d_lumi5s = mymclimit->lumi5s(); 
		cout << "<<<<<<<< Finshed calculating scale factors at 95% CL >>>>>>>>" << endl;
		f_outputfile  << "Observed scale factor of 95% CL excluded signal: sc_f = " << d_sf95 << endl;
		f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(-2sig) = " << d_s95m2 << endl;
		f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(-1sig) = " << d_s95m1 << endl;
		f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(median) = " << d_s95med << endl;
		f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(+1sig) = " << d_s95p1 << endl;
		f_outputfile  << "Expected scale factor at 95% CL in the null hyp.: sc_f(+2sig) = " << d_s95p2 << endl << endl;
		//	f_outputfile << "Lumi scale factor of 95% CL excluded signal: lumiSf = " << d_lumi95 << endl;
		f_outputfile << "Lumi needed for 3 sigma discovery: lumi3sig = " << d_lumi3s*d_totlumi << " fb-1" << endl;
		f_outputfile << "Lumi needed for 5 sigma discovery: lumi5sig = " << d_lumi5s*d_totlumi << " fb-1" << endl << endl;
		cout << "<<<<<<<< Finished calculating luminosity scale factors for discovery 95% CL >>>>>>>>" << endl;
		
		if (combined11) {
			f_outputfile << "CMS-LHCb summer 11 combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (CERN combined summer 11): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		else if (combined11bd) {
			f_outputfile << "CMS-LHCb summer 11 combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (CERN combined summer 11): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
		}
		else if (combined12) {
			f_outputfile << "CMS-LHCb winter 11-12 combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (CERN combined winter 11-12): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		else if (combined12bd) {
			f_outputfile << "CMS-LHCb winter 11-12 combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (CERN combined winter 11-12): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
		}
		else if (cms11bs) {
			f_outputfile << "CMS summer 11 Bs->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (CMS summer 11): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		else if (cms11bd) {
			f_outputfile << "CMS summer 11 Bd->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (CMS summer 11): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
		}
		else if (cms12bs) {
			f_outputfile << "CMS winter 11-12 Bs->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (CMS winter 11-12): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		else if (cms12bd) {
			f_outputfile << "CMS winter 11-12 Bd->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (CMS winter 11-12): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
		}
		else if (lhcbs) {
			f_outputfile << "LHCb 10-11 Bs->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (LHCb): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		else if (lhcbs10) {
			f_outputfile << "LHCb 10 Bs->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (LHCb 10): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		else if (lhcbs11) {
			f_outputfile << "LHCb 11 Bs->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (LHCb 11): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		else if (lhcbd11) {
			f_outputfile << "LHCb summer 11 Bd->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (LHCb summer 11): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
		}
		else if (lhcbs12) {
			f_outputfile << "LHCb 11-12 Bs->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (LHCb 11-12): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		else if (lhcbd12) {
			f_outputfile << "LHCb winter 11-12 Bd->MuMu combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bd->MuMu) = " << d_s95m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bd->MuMu) = " << d_s95m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bd->MuMu) = " << d_s95med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bd->MuMu) = " << d_s95p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bd->MuMu) = " << d_s95p2*Brbdmm << endl;
			f_outputfile << "95% CL upper limit (LHCb winter 11-12): Br(Bd->MuMu) = " << d_sf95*Brbdmm << endl << endl;
			
		}
		else if (atlasbs12) {
			f_outputfile << "ATLAS 12 Bs->MuMu result:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -2sig: Br(Bs->MuMu) = " << d_s95m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL -1sig: Br(Bs->MuMu) = " << d_s95m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL median: Br(Bs->MuMu) = " << d_s95med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +1sig: Br(Bs->MuMu) = " << d_s95p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 95% CL +2sig: Br(Bs->MuMu) = " << d_s95p2*Brbsmm << endl;
			f_outputfile << "95% CL upper limit (ATLAS 12): Br(Bs->MuMu) = " << d_sf95*Brbsmm << endl << endl;
			
		}
		
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	if (conflevel90) {
		cout << "<<<<<<<< Calculating scale factors at 90% CL >>>>>>>>" << endl;
		double d_sf90 = mymclimit->s90();double d_s90m2 = mymclimit->s90m2();
		double d_s90m1 = mymclimit->s90m1(); double d_s90med = mymclimit->s90med();
		double d_s90p1 = mymclimit->s90p1(); double d_s90p2 = mymclimit->s90p2();	
		cout << "<<<<<<<< Finshed calculating scale factors at 90% CL >>>>>>>>" << endl;
		f_outputfile  << "Observed scale factor of 90% CL excluded signal: sc_f = " << d_sf90 << endl;
		f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(-2sig) = " << d_s90m2 << endl;
		f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(-1sig) = " << d_s90m1 << endl;
		f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(median) = " << d_s90med << endl;
		f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(+1sig) = " << d_s90p1 << endl;
		f_outputfile  << "Expected scale factor at 90% CL in the null hyp.: sc_f(+2sig) = " << d_s90p2 << endl << endl;
		
		if (combined11) {
			f_outputfile << "CMS-LHCb summer 11 combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (CERN combined summer 11): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl;
		}
		else if (combined11bd) {
			f_outputfile << "CMS-LHCb summer 11 combined results:" << endl << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (CERN combined summer 11): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
		else if (combined12) {
			f_outputfile << "CMS-LHCb winter 11-12 combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (CERN combined winter 11-12): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		else if (combined12bd) {
			f_outputfile << "CMS-LHCb winter 11-12 combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (CERN combined winter 11-12): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
		else if (cms11bs) {
			f_outputfile << "CMS summer 11 Bs->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (CMS summer 11): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		else if (cms11bd) {
			f_outputfile << "CMS summer 11 Bd->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (CMS summer 11): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
		else if (cms12bs) {
			f_outputfile << "CMS winter 11-12 Bs->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (CMS winter 11-12): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		else if (cms12bd) {
			f_outputfile << "CMS winter 11-12 Bd->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (CMS winter 11-12): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
		else if (lhcbs) {
			f_outputfile << "LHCb 10-11 Bs->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (LHCb): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		else if (lhcbs10) {
			f_outputfile << "LHCb 10 Bs->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (LHCb 10): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		else if (lhcbs11) {
			f_outputfile << "LHCb 11 Bs->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (LHCb 11): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		else if (lhcbd11) {
			f_outputfile << "LHCb summer 11 Bd->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (LHCb summer 11): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
		else if (lhcbs12) {
			f_outputfile << "LHCb 11-12 Bs->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (LHCb 11-12): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
		}
		else if (lhcbd12) {
			f_outputfile << "LHCb winter 11-12 Bd->MuMu combined results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bd->MuMu) = " << d_s90m2*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bd->MuMu) = " << d_s90m1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bd->MuMu) = " << d_s90med*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bd->MuMu) = " << d_s90p1*Brbdmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bd->MuMu) = " << d_s90p2*Brbdmm << endl;
			f_outputfile << "90% CL upper limit (LHCb winter 11-12): Br(Bd->MuMu) = " << d_sf90*Brbdmm << endl << endl;
		}
		else if (atlasbs12) {
			f_outputfile << "ATLAS 12 Bs->MuMu results:" << endl << endl;
			
			f_outputfile << "Exp null hyp. upper limit at 90% CL -2sig: Br(Bs->MuMu) = " << d_s90m2*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL -1sig: Br(Bs->MuMu) = " << d_s90m1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL median: Br(Bs->MuMu) = " << d_s90med*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +1sig: Br(Bs->MuMu) = " << d_s90p1*Brbsmm << endl;
			f_outputfile << "Exp null hyp. upper limit at 90% CL +2sig: Br(Bs->MuMu) = " << d_s90p2*Brbsmm << endl;
			f_outputfile << "90% CL upper limit (ATLAS 12): Br(Bs->MuMu) = " << d_sf90*Brbsmm << endl << endl;
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
	
	
	//==============================================
	//Draw the test statistics distribution 
	//test hypothesis - we can calculate probability of observering more than N sigma 
	//null hypothesis - we can calculate p-value of a measured test statistics  
	TH1F* ts_test = new TH1F("ts_test","",200,-100,100);
	TH1F* ts_null = new TH1F("ts_null","",200,-100,100);
	////////////////////////////////////////////////////////////////////////
	TH1F* mcb_hist = new TH1F("mcb_hist","",20,0,2);
	TH1F* mcs_hist = new TH1F("mcs_hist","",20,0,2);
	TH1F* data_hist = new TH1F("data_hist","",20,0,2);
	////////////////////////////////////////////////////////////////////////
	mymclimit->tshists(ts_test,ts_null); 
	if ((cms11bs || cms11bd || cms12bs || cms12bd) && !(combined11 || combined12)) {
		mycanvas->cd(3);//////
		mymclimit->plotlnsb(mcb_hist,mcs_hist,data_hist);
		mycanvas->cd(4);
	}
	ts_null->SetLineColor(2);//red
	ts_test->SetLineColor(4);//blue
	ts_null->Draw();
	ts_test->Draw("same");
	
	TArrow* arrow = new TArrow(tsobs,0.3*ts_null->GetMaximum(),tsobs,0);
	arrow->Draw();
	TString s_pdffilename;
	if (combined11) s_pdffilename = Form("CMS_LHCb_S11_Comb_Bs_Results%f_%i.pdf",d_scale,i_numb);
	else if (combined11bd) s_pdffilename = Form("CMS_LHCb_S11_Comb_Bd_Results%f_%i.pdf",d_scale,i_numb);
	else if (combined12) s_pdffilename = Form("CMS_LHCb_W12_Comb_Bs_Results%f_%i.pdf",d_scale,i_numb);
	else if (combined12bd) s_pdffilename = Form("CMS_LHCb_W12_Comb_Bd_Results%f_%i.pdf",d_scale,i_numb);
	else if (cms11bs) s_pdffilename = Form("CMS_S11_Bs_Results%f_%i.pdf",d_scale,i_numb);
	else if (cms11bd) s_pdffilename = Form("CMS_S11_Bd_Results%f_%i.pdf",d_scale,i_numb);
	else if (cms12bs) s_pdffilename = Form("CMS_W12_Bs_Results%f_%i.pdf",d_scale,i_numb);
	else if (cms12bd) s_pdffilename = Form("CMS_W12_Bd_Results%f_%i.pdf",d_scale,i_numb);
	else if (lhcbs) s_pdffilename = Form("LHCb_10_11_Comb_Bs_Results%f_%i.pdf",d_scale,i_numb);
	else if (lhcbs10) s_pdffilename = Form("LHCb_10_Bs_Results%f_%i.pdf",d_scale,i_numb);
	else if (lhcbs11) s_pdffilename = Form("LHCb_11_Bs_Results%f_%i.pdf",d_scale,i_numb);
	else if (lhcbd11) s_pdffilename = Form("LHCb_11_Bd_Results%f_%i.pdf",d_scale,i_numb);
	else if (lhcbs12) s_pdffilename = Form("LHCb_12_Bs_Results%f_%i.pdf",d_scale,i_numb);
	else if (lhcbd12) s_pdffilename = Form("LHCb_12_Bd_Results%f_%i.pdf",d_scale,i_numb);
	else if (atlasbs12) s_pdffilename = Form("ATLAS_12_Bs_Results%f_%i.pdf",d_scale,i_numb);
	mycanvas->Print(s_pdffilename);
	
	
	delete mymclimit;
	f_outfile.Write();
	f_outfile.Close();
	f_outputfile.close();

	t.Stop();
	t.Print();
}
