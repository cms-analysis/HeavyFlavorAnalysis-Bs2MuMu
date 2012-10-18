// By: Jose Lazo-Flores
// Modified on: Oct 18, 2012
// common source for preparing t-channel likelihood function
// limit inputs -- using SM with SM single top
// as the null hypothesis and SM with SM single top and extra anomalous
// single top as the test hypothesis

gROOT->SetStyle("Plain");
gROOT->ForceStyle();

// declare some constant variables
static const double Brbsmm = 3.2e-9;
static const double Brbdmm = 1.0e-10;
static const double cmslumiw12 = 5.0;

//declare nusiance parameters
const int NSYS = 10;
char *ename[NSYS];
char accname[]="ACC";// acceptance
char toteffname[]="TOTEFF";
char pdfhadron[]="HPDF";
char pssname[]="PSS";
//char psdname[]="PSD";
//char pdsname[]="PDS";
char pddname[]="PDD";
char bkgerrname[]="BKGERR";
char experrname[]="EXPERR";
char peakerrname[]="PCKERR";
char unconstrned[]="UNCONSTRAINED";

char channameB[] = "bmmB";
char channameE[] = "bmmE";

//declare the shape functions
double nps_low[NSYS];
double nps_high[NSYS];
double lowsigma[NSYS];
double highsigma[NSYS];

double sfact; // scale factor
int nps_count; // number of nuisance paramaters
int pssnflg;// Poisson flag -- 1 if Poisson, 0 if not.
int sclflg;// Scale flag -- set to 1 if you want this template to be scaled in the s95 calculation, 0 if you don't.
Double_t bsmsswin = 5.45 - 5.3; Double_t bdmsswin = 5.3 - 5.2; Double_t msswin = 5.9 - 4.9 - 0.25;

TH1 *lowshape[NSYS];
TH1 *highshape[NSYS];


//Double_t xbins[3] = {0, 1, 2};

// Opent the data file
TFile *f_data = TFile::Open("anaBmm.default-11.root");
TString s_outfilename;
if (bscombined) s_outfilename = Form("CMS_W12_Bs_Comb_ShpReslts%f_%i.root",d_scale,i_numb);
else if (bdcombined) s_outfilename = Form("CMS_W12_Bd_Comb_ShpReslts%f_%i.root",d_scale,i_numb);
else if(barrel && bsmm) s_outfilename = Form("CMS_W12_Bs_brrl_ShpReslts%f_%i.root",d_scale,i_numb);
else if(endcap && bsmm) s_outfilename = Form("CMS_W12_Bs_endc_ShpReslts%f_%i.root",d_scale,i_numb);
else if(barrel && bdmm) s_outfilename = Form("CMS_W12_Bd_brrl_ShpReslts%f_%i.root",d_scale,i_numb);
else if(endcap && bdmm) s_outfilename = Form("CMS_W12_Bd_endc_ShpReslts%f_%i.root",d_scale,i_numb);
cout << s_outfilename << endl;
TFile f_outfile(s_outfilename,"RECREATE");

//Templates and input Data
TH1D *h_sigbsB = (TH1D*)f_data->Get("Bs_cnc_chan0");
TH1D *h_sigbdB = (TH1D*)f_data->Get("Bd_cnc_chan0");
TH1D *h_bkgB = (TH1D*)f_data->Get("bRare_cnc");
TH1D *h_obsB = (TH1D*)f_data->Get("hMassWithAllCuts_cnc_5_chan0");
TH1D *h_sigbsE = (TH1D*)f_data->Get("Bs_cnc_chan1");
TH1D *h_sigbdE = (TH1D*)f_data->Get("Bd_cnc_chan1");
TH1D *h_bkgE = (TH1D*)f_data->Get("eRare_cnc");
TH1D *h_obsE = (TH1D*)f_data->Get("hMassWithAllCuts_cnc_5_chan1");

int lowmassbin = h_sigbsB->FindBin(4.9); int himassbin = h_sigbsB->FindBin(5.9);
int lobldbin = h_sigbsB->FindBin(5.2); int hibldbin = h_sigbsB->FindBin(5.45); 
int lobsmssbin = h_sigbsB->FindBin(5.30);
Int_t nbins = himassbin - lowmassbin;
cout << "number of bins = " << nbins << endl;

TH1D *h_rarebkgB = new TH1D("h_rarebkgB","rare background barrel",nbins,4.9,5.9); 
TH1D *h_combbkgB = new TH1D("h_combbkgB","combinatorial background barrel",nbins,4.9,5.9); 
TH1D *h_signalbsB = new TH1D("h_signalbsB","signal bs barrel",nbins,4.9,5.9); 
TH1D *h_signalbdB = new TH1D("h_signalbdB","signal bd barrel",nbins,4.9,5.9); 
TH1D *h_dataB = new TH1D("h_dataB","Observed events in barrel",nbins,4.9,5.9);
TH1D *h_rarebkgE	= new TH1D("h_rarebkgE","rare background endcaps",nbins,4.9,5.9); 
TH1D *h_combbkgE = new TH1D("h_combbkgE","combinatorial background endcaps",nbins,4.9,5.9); 
TH1D *h_signalbsE = new TH1D("h_signalbsE","signal bs endcaps",nbins,4.9,5.9); 
TH1D *h_signalbdE = new TH1D("h_signalbdE","signal bd endcaps",nbins,4.9,5.9); 
TH1D *h_dataE = new TH1D("h_dataE","Observed events in endcaps",nbins,4.9,5.9);

cout << "<<<<<< Done reading data >>>>>>>" << endl;

TF1 *myft = new TF1("myft","[0]+[1]*x",0,1);
myft->SetParameters(1,0);
//h_test->FillRandom("myft",10000);

//return 0;
//========================================
//Construct test/null hypothesis for fitting
csm_model* nullhyp = new csm_model();
csm_model* testhyp = new csm_model();
csm_model* nullhyp_pe = new csm_model();
csm_model* testhyp_pe = new csm_model();



if (barrel) {
	Double_t obsrvtot = h_obsB->Integral(lowmassbin,himassbin);
	Double_t rarebkgtot = h_bkgB->Integral(lowmassbin,himassbin);
	cout << "Observed = " << obsrvtot << ", rareBkg = " << rarebkgtot << endl << endl;

	Double_t obsrvleft = h_obsB->Integral(lowmassbin,lobldbin); Double_t obsrvright = h_obsB->Integral(hibldbin,himassbin);
	Double_t obsrvshldrs = obsrvleft + obsrvright;
	cout << "Observed (blind) = " << obsrvshldrs << endl << endl;
	Double_t obsrvbswin = h_obsB->Integral(lobsmssbin, hibldbin); Double_t obsrvbdwin = h_obsB->Integral(lobldbin,lobsmssbin);
	cout << "Observed in Bs window = " << obsrvbswin << ", in Bd window = " << obsrvbdwin << endl << endl;

	Double_t rarebkgbdwin = h_bkgB->Integral(lobldbin,lobsmssbin); Double_t rarebkgbdwinerr = 0.070;
	Double_t rarebkgbswin = h_bkgB->Integral(lobsmssbin,hibldbin); Double_t rarebkgbswinerr = 0.060;
	Double_t rarebkgleft = h_bkgB->Integral(lowmassbin,lobldbin); Double_t rarebkgright = h_bkgB->Integral(hibldbin,himassbin);
	Double_t rarebkgshldrs = rarebkgleft + rarebkgright; Double_t rarebkgshldrserr = 0.632963;
	cout << "Rare bkgd (shoulders) = " << rarebkgshldrs << endl;
	cout << "Rare bkgd (peaking in Bd signal window) = " << rarebkgbdwin << endl;
	cout << "Rare bkgd (peaking in Bs signal window) = " << rarebkgbswin << endl << endl;
	
	Double_t combkg = (obsrvshldrs - rarebkgshldrs);
	Double_t combkgbd = (obsrvshldrs - rarebkgleft)*bdmsswin/msswin; Double_t combkgbderr = 0.34;
	Double_t combkgbs = (obsrvshldrs - rarebkgleft)*bsmsswin/msswin; Double_t combkgbserr = 0.50;
	Double_t combkgtot = combkg + combkgbd + combkgbs;
	cout << "CombBkgd (shoulders) = " << combkg << endl;
	cout << "CombBkgd (in Bd blind window) = " << combkgbd << " +- " << combkgbderr << endl;
	cout << "CombBkgd (in Bs blind window) = " << combkgbs << " +- " << combkgbserr  << endl << endl;
	cout << "CombBkgd (total) = " << combkgtot << endl;

	Double_t exptbdevnts = h_sigbdB->Integral(lobldbin,lobsmssbin); Double_t exptbdevntserr = 0.020;
	Double_t exptbdevntsbswin = h_sigbdB->Integral(lobsmssbin,hibldbin);
	Double_t exptbsevnts = h_sigbsB->Integral(lobsmssbin,hibldbin); Double_t exptbsevntserr = 0.41;
	Double_t exptbsevntsbdwin = h_sigbsB->Integral(lobldbin,lobsmssbin);
	cout << "Expected Bd signal events = " << exptbdevnts << " +- " << exptbdevntserr << endl;
	cout << "Expected Bd signal events in Bs window = " << exptbdevntsbswin << endl;	
	cout << "Expected Bs signal events = " << exptbsevnts << " +- " << exptbsevntserr  << endl;
	cout << "Expected Bs signal events in Bd window = " << exptbsevntsbdwin << endl << endl;

	Double_t totbdwin = exptbdevnts + combkgbd + rarebkgbdwin + exptbsevntsbdwin;// including exptbsevntsbdwin
	Double_t totbdwinerr = 0.35;
	Double_t totbswin = exptbsevnts + combkgbs + rarebkgbswin + exptbdevntsbswin;// including exptbsevntsbswin
	Double_t totbswinerr = 0.65;
	Double_t totblndwin = totbswin + totbdwin;
	cout << "Total (expected in Bd blind window) = " << totbdwin << " ~ " << floor(totbdwin + 0.5) << endl;
	cout << "Total (expected in Bs blind window) = " << totbswin << " ~ " << floor(totbswin + 0.5) << endl;
	cout << "Total (expected in blind window) = " << totblndwin << " ~ " << floor(totblndwin + 0.5) << endl << endl;

	Double_t toteff = 0.0029; Double_t totefferr = 0.0002;
	
	// Setup the Histos
	for (int i = 0; i < nbins; i++) {
		//Rare bkgd Histogram
		h_rarebkgB->SetBinContent(i,h_bkgB->GetBinContent(lowmassbin+i));
		
		//Signal Histograms
		h_signalbsB->SetBinContent(i,h_sigbsB->GetBinContent(lowmassbin+i));
		h_signalbdB->SetBinContent(i,h_sigbdB->GetBinContent(lowmassbin+i));
		
		// Data histogram
		h_dataB->SetBinContent(i,h_obsB->GetBinContent(lowmassbin+i));
	}
	if (bsmm) h_signalbsB->Scale(d_scale);
	if (bdmm) h_signalbdB->Scale(d_scale);

	// Combinatorial background histo
	h_combbkgB->FillRandom("myft",10000);
	Double_t scale = combkgtot/h_combbkgB->Integral();
	h_combbkgB->Scale(scale);
	
	cout << "<<<<<< Done Setting up Histos for Barrel Channel>>>>>>>" << endl;


	// Add combinatorial background templates
	for(int i=0;i<NSYS;i++)
	{
		nps_low[i]  = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i]= 0;
		lowshape[i] = 0;
		highshape[i]= 0;
	}
	
	sfact = 1;
	if (bsmm) {
		ename[0] = bkgerrname;
		nps_low[0] = 0.;
		nps_high[0] = 0.;
		
		nps_count=0;
	}	
	if (bdmm) {
		ename[0] = bkgerrname;
		nps_low[0] = 0;
		nps_high[0] = 0;
		
		nps_count=0;
	}

	pssnflg = 0;
	sclflg = 0;

	// Construct test/null hypothesis for pseudo-experiments.
	nullhyp_pe->add_template(h_combbkgB,sfact,nps_count,ename,nps_low,nps_high,
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	testhyp_pe->add_template(h_combbkgB,sfact,nps_count,ename,nps_low,nps_high,
						 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	nps_count=0;
	nullhyp->add_template(h_combbkgB,sfact,nps_count,ename,nps_low,nps_high,
					  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	testhyp->add_template(h_combbkgB,sfact,nps_count,ename,nps_low,nps_high,
					  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	cout << "Finished setting up nullhyp, testhyp, nullhyp_pe, testhyp_pe for combinatorial background" << endl;


	// Add rare background templates
	for(int i=0;i<NSYS;i++)
	{
		nps_low[i]  = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i]= 0;
		lowshape[i] = 0;
		highshape[i]= 0;
	}
	
	sfact = 1;
	if (bsmm) {
		ename[0] = bkgerrname;
		nps_low[0] = -rarebkgshldrserr/rarebkgshldrs;
		nps_high[0] = rarebkgshldrserr/rarebkgshldrs;
		
//		ename[1] = psdname;
//		nps_low[1] = -0.0224451/0.291961;
//		nps_high[1] = 0.0224451/0.291961;
		
		nps_count=1;	
	}	
	if (bdmm) {
		ename[0] = bkgerrname;
		nps_low[0] = -rarebkgshldrserr/rarebkgshldrs;
		nps_high[0] = rarebkgshldrserr/rarebkgshldrs;
		
//		ename[1] = pdsname;
//		nps_low[1] = -0.00741162/0.0712015;
//		nps_high[1] = 0.00741162/0.0712015;		
		
		nps_count=1;
	}
	
	pssnflg = 0;
	sclflg = 0;
	
	// Construct test/null hypothesis for pseudo-experiments.
	nullhyp_pe->add_template(h_rarebkgB,sfact,nps_count,ename,nps_low,nps_high,
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	testhyp_pe->add_template(h_rarebkgB,sfact,nps_count,ename,nps_low,nps_high,
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	nps_count=0;
	nullhyp->add_template(h_rarebkgB,sfact,nps_count,ename,nps_low,nps_high,
						  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	testhyp->add_template(h_rarebkgB,sfact,nps_count,ename,nps_low,nps_high,
						  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	cout << "Finished setting up nullhyp, testhyp, nullhyp_pe, testhyp_pe for rare background" << endl;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Add Bd/bs background templates
	for(int i=0;i<NSYS;i++)
	{
		nps_low[i]  = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i]= 0;
		lowshape[i] = 0;
		highshape[i]= 0;
	}
	
	sfact = 1;
	if (bsmm) {
		ename[0] = unconstrned;
		nps_low[0] = -0.0224451/0.291961;
		nps_high[0] = 0.0224451/0.291961;
		
//		ename[1] = psdname;
//		nps_low[1] = -0.0224451/0.291961;
//		nps_high[1] = 0.0224451/0.291961;
		
		nps_count=1;	
		pssnflg = 0;
		sclflg = 0;
		
		// Construct test/null hypothesis for pseudo-experiments.
		nullhyp_pe->add_template(h_signalbdB,sfact,nps_count,ename,nps_low,nps_high,
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		testhyp_pe->add_template(h_signalbdB,sfact,nps_count,ename,nps_low,nps_high,
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		nps_count=0;
		nullhyp->add_template(h_signalbdB,sfact,nps_count,ename,nps_low,nps_high,
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		testhyp->add_template(h_signalbdB,sfact,nps_count,ename,nps_low,nps_high,
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		cout << "Finished setting up nullhyp, testhyp, nullhyp_pe, testhyp_pe for rare background" << endl;
	}	
	if (bdmm) {
		ename[0] = unconstrned;
		nps_low[0] = -0.00741162/0.0712015;
		nps_high[0] = 0.00741162/0.0712015;
		
//		ename[1] = pdsname;
//		nps_low[1] = -0.00741162/0.0712015;
//		nps_high[1] = 0.00741162/0.0712015;		
		
		nps_count=1;
		pssnflg = 0;
		sclflg = 0;
		
		// Construct test/null hypothesis for pseudo-experiments.
		nullhyp_pe->add_template(h_signalbsB,sfact,nps_count,ename,nps_low,nps_high,
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		testhyp_pe->add_template(h_signalbsB,sfact,nps_count,ename,nps_low,nps_high,
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		nps_count=0;
		nullhyp->add_template(h_signalbsB,sfact,nps_count,ename,nps_low,nps_high,
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		testhyp->add_template(h_signalbsB,sfact,nps_count,ename,nps_low,nps_high,
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		cout << "Finished setting up nullhyp, testhyp, nullhyp_pe, testhyp_pe for rare background" << endl;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////

	// Add signal templates
	for(int i=0;i<NSYS;i++) {
		nps_low[i] = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i] = 0;
		lowshape[i] = 0;
		highshape[i] = 0;
	}

	if (bsmm) {
		
		ename[0] = pdfhadron;
		nps_low[0] = -0.021/0.267;
		nps_high[0] = 0.021/0.267;

		ename[1] = toteffname;
		nps_low[1] = -totefferr/toteff;
		nps_high[1] = totefferr/toteff;
		
		ename[2] = experrname;
		nps_low[2] = -totbswinerr/totbswin;
		nps_high[2] = totbswinerr/totbswin;
		
		ename[3] = accname;
		nps_low[3] = -0.009/0.248;
		nps_high[3] = 0.009/0.248;
		
		ename[4] = pssname;
		nps_low[4] = -0.0439992/0.862683;
		nps_high[4] = 0.0439992/0.862683;

		ename[5] = peakerrname;
		nps_low[5] = -rarebkgbswinerr/rarebkgbswin;
		nps_high[5] = rarebkgbswinerr/rarebkgbswin;
				
//		ename[6] = combbgerrname;
//		nps_low[6] = -combkgbserr/combkgbs;
//		nps_high[6] = combkgbserr/combkgbs;
		
		nps_count = 6;

		sfact = 1.;
		
		pssnflg = 0;// 
		sclflg = 1;// this is set to 1 if signal, 
		testhyp_pe->add_template(h_signalbsB,sfact,nps_count,ename,nps_low,nps_high, 
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		nps_count = 0;
		testhyp->add_template(h_signalbsB,sfact,nps_count,ename,nps_low,nps_high, 
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		
	}
	if (bdmm) {

		ename[0] = toteffname;
		nps_low[0] = -totefferr/toteff;
		nps_high[0] = totefferr/toteff;
		
		ename[1] = experrname;
		nps_low[1] = -totbdwinerr/totbdwin;
		nps_high[1] = totbdwinerr/totbdwin;
		
		ename[2] = accname;
		nps_low[2] = -0.009/0.247;
		nps_high[2] = 0.009/0.247;

		ename[3] = pddname;
		nps_low[3] = -0.0366707/0.638928;
		nps_high[3] = 0.0366707/0.638928;

		ename[4] = peakerrname;
		nps_low[4] = -rarebkgbdwinerr/rarebkgbdwin;
		nps_high[4] = rarebkgbdwinerr/rarebkgbdwin;

//		ename[5] = combbgerrname;
//		nps_low[5] = -combkgbderr/combkgbd;
//		nps_high[5] = combkgbderr/combkgbd;		

		nps_count = 5;

		sfact = 1.;
		
		pssnflg = 0;// 
		sclflg = 1;// this is set to 1 if signal, 
		testhyp_pe->add_template(h_signalbdB,sfact,nps_count,ename,nps_low,nps_high, 
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
		nps_count = 0;
		testhyp->add_template(h_signalbdB,sfact,nps_count,ename,nps_low,nps_high, 
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	}
	
	cout << "Finished setting up testhyp, testhyp_pe for signal in the barrel" << endl;	
}

if (endcap) {
	Double_t obsrvtot = h_obsE->Integral(lowmassbin,himassbin);
	Double_t rarebkgtot = h_bkgE->Integral(lowmassbin,himassbin);
	cout << "Observed = " << obsrvtot << ", rareBkg = " << rarebkgtot << endl << endl;
	
	Double_t obsrvleft = h_obsE->Integral(lowmassbin,lobldbin); Double_t obsrvright = h_obsE->Integral(hibldbin,himassbin);
	Double_t obsrvshldrs = obsrvleft + obsrvright;
	cout << "Observed (blind) = " << obsrvshldrs << endl << endl;
	Double_t obsrvbswin = h_obsE->Integral(lobsmssbin, hibldbin); Double_t obsrvbdwin = h_obsE->Integral(lobldbin,lobsmssbin);
	cout << "Observed in Bs window = " << obsrvbswin << ", in Bd window = " << obsrvbdwin << endl << endl;
	
	Double_t rarebkgbdwin = h_bkgE->Integral(lobldbin,lobsmssbin); Double_t rarebkgbdwinerr = 0.030;
	Double_t rarebkgbswin = h_bkgE->Integral(lobsmssbin,hibldbin); Double_t rarebkgbswinerr = 0.020;
	Double_t rarebkgleft = h_bkgE->Integral(lowmassbin,lobldbin); Double_t rarebkgright = h_bkgE->Integral(hibldbin,himassbin);
	Double_t rarebkgshldrs = rarebkgleft + rarebkgright; Double_t rarebkgshldrserr = 0.238872;
	cout << "Rare bkgd (shoulders) = " << rarebkgshldrs << endl;
	cout << "Rare bkgd (peaking in Bd signal window) = " << rarebkgbdwin << endl;
	cout << "Rare bkgd (peaking in Bs signal window) = " << rarebkgbswin << endl << endl;

	Double_t combkg = (obsrvshldrs - rarebkgshldrs);
	Double_t combkgbd = (obsrvshldrs - rarebkgleft)*bdmsswin/msswin; Double_t combkgbderr = 0.35;
	Double_t combkgbs = (obsrvshldrs - rarebkgleft)*bsmsswin/msswin; Double_t combkgbserr = 0.53;
	Double_t combkgtot = combkg + combkgbd + combkgbs;
	cout << "CombBkgd (shoulders) = " << combkg << endl;
	cout << "CombBkgd (in Bd blind window) = " << combkgbd << " +- " << combkgbderr << endl;
	cout << "CombBkgd (in Bs blind window) = " << combkgbs << "+- " << combkgbserr << endl << endl;
	cout << "CombBkgd (total) = " << combkgtot << endl;

	Double_t exptbdevnts = h_sigbdE->Integral(lobldbin,lobsmssbin); Double_t exptbdevntserr = 0.010;
	Double_t exptbdevntsbswin = h_sigbdE->Integral(lobsmssbin,hibldbin);
	Double_t exptbsevnts = h_sigbsE->Integral(lobsmssbin,hibldbin); Double_t exptbsevntserr = 0.18;
	Double_t exptbsevntsbdwin = h_sigbsE->Integral(lobldbin,lobsmssbin);
	cout << "Expected Bd signal events = " << exptbdevnts << " +- " << exptbdevntserr << endl;
	cout << "Expected Bd signal events in Bs window = " << exptbdevntsbswin << endl;	
	cout << "Expected Bs signal events = " << exptbsevnts << " +- " << exptbsevntserr<< endl;
	cout << "Expected Bs signal events in Bd window = " << exptbsevntsbdwin << endl << endl;

	Double_t totbdwin = exptbdevnts + combkgbd + rarebkgbdwin + exptbsevntsbdwin;// including exptbsevntsbdwin
	Double_t totbdwinerr = 0.35;
	Double_t totbswin = exptbsevnts + combkgbs + rarebkgbswin + exptbdevntsbswin;// including exptbsevntsbswin
	Double_t totbswinerr = 0.56;
	Double_t totblndwin = totbswin + totbdwin;
	cout << "Total (expected in Bd blind window) = " << totbdwin << " ~ " << floor(totbdwin + 0.5) << endl;
	cout << "Total (expected in Bs blind window) = " << totbswin << " ~ " << floor(totbswin + 0.5) << endl;
	cout << "Total (expected in blind window) = " << totblndwin << " ~ " << floor(totblndwin + 0.5) << endl << endl;
	
	Double_t toteff = 0.0016; Double_t totefferr = 0.0002;

	// Setup the Histos
	// Setup the Histos
	for (int i = 0; i < nbins; i++) {
		//Rare bkgd Histogram
		h_rarebkgE->SetBinContent(i,h_bkgE->GetBinContent(lowmassbin+i));
		
		//Signal Histograms
		h_signalbsE->SetBinContent(i,h_sigbsE->GetBinContent(lowmassbin+i));
		h_signalbdE->SetBinContent(i,h_sigbdE->GetBinContent(lowmassbin+i));
		
		// Data histogram
		h_dataE->SetBinContent(i,h_obsE->GetBinContent(lowmassbin+i));
	}
	if (bsmm) h_signalbsE->Scale(d_scale);
	if (bdmm) h_signalbdE->Scale(d_scale);
	
	// Combinatorial background histo
	h_combbkgE->FillRandom("myft",10000);
	Double_t scale = combkgtot/h_combbkgE->Integral();
	h_combbkgE->Scale(scale);
	
	cout << "<<<<<< Done Setting up Histos for Endcap Channel>>>>>>>" << endl;
		
	
	// Add combinatorial background templates
	for(int i=0;i<NSYS;i++)
	{
		nps_low[i]  = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i]= 0;
		lowshape[i] = 0;
		highshape[i]= 0;
	}
	
	sfact = 1;
	if (bsmm) {
		ename[0] = bkgerrname;
		nps_low[0] = 0.;
		nps_high[0] = 0.;
		
		nps_count=0;
	}	
	if (bdmm) {
		ename[0] = bkgerrname;
		nps_low[0] = 0.;
		nps_high[0] = 0.;
		
		nps_count=0;
	}
	
	pssnflg = 0;
	sclflg = 0;
	
	// Construct test/null hypothesis for pseudo-experiments.
	nullhyp_pe->add_template(h_combbkgE,sfact,nps_count,ename,nps_low,nps_high,
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	testhyp_pe->add_template(h_combbkgE,sfact,nps_count,ename,nps_low,nps_high,
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	nps_count=0;
	nullhyp->add_template(h_combbkgE,sfact,nps_count,ename,nps_low,nps_high,
						  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	testhyp->add_template(h_combbkgE,sfact,nps_count,ename,nps_low,nps_high,
						  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	cout << "Finished setting up nullhyp, testhyp, nullhyp_pe, testhyp_pe for combinatorial background" << endl;
	
	// Add rare background templates
	// Initialize everything
	for(int i=0;i<NSYS;i++)
	{
		nps_low[i]  = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i]= 0;
		lowshape[i] = 0;
		highshape[i]= 0;
	}
	
	sfact = 1;
	if (bsmm) {
		ename[0] = bkgerrname;
		nps_low[0] = -rarebkgshldrserr/rarebkgshldrs;
		nps_high[0] = rarebkgshldrserr/rarebkgshldrs;
		
//		ename[1] = psdname;
//		nps_low[1] = -0.0285105/0.318653;
//		nps_high[1] = 0.0285105/0.318653;
		
		nps_count=1;	
	}	
	if (bdmm) {
		ename[0] = bkgerrname;
		nps_low[0] = -rarebkgshldrserr/rarebkgshldrs;
		nps_high[0] = rarebkgshldrserr/rarebkgshldrs;
		
//		ename[1] = pdsname;
//		nps_low[1] = -0.0148983/0.163842;
//		nps_high[1] = 0.0148983/0.163842;
		
		nps_count=1;
	}
	
	pssnflg = 0;
	sclflg = 0;
	
	// Construct test/null hypothesis for pseudo-experiments.
	nullhyp_pe->add_template(h_rarebkgE,sfact,nps_count,ename,nps_low,nps_high,
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	testhyp_pe->add_template(h_rarebkgE,sfact,nps_count,ename,nps_low,nps_high,
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	nps_count=0;
	nullhyp->add_template(h_rarebkgE,sfact,nps_count,ename,nps_low,nps_high,
						  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	testhyp->add_template(h_rarebkgE,sfact,nps_count,ename,nps_low,nps_high,
						  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	cout << "Finished setting up nullhyp, testhyp, nullhyp_pe, testhyp_pe for rare background" << endl;

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Add Bs/Bd background templates
	// Initialize everything
	for(int i=0;i<NSYS;i++)
	{
		nps_low[i]  = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i]= 0;
		lowshape[i] = 0;
		highshape[i]= 0;
	}
	
	sfact = 1;
	if (bsmm) {
		ename[0] = unconstrned;
		nps_low[0] = -0.0285105/0.318653;
		nps_high[0] = 0.0285105/0.318653;
		
		nps_count=1;	
		pssnflg = 0;
		sclflg = 0;
		
		// Construct test/null hypothesis for pseudo-experiments.
		nullhyp_pe->add_template(h_signalbdE,sfact,nps_count,ename,nps_low,nps_high,
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		testhyp_pe->add_template(h_signalbdE,sfact,nps_count,ename,nps_low,nps_high,
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		nps_count=0;
		nullhyp->add_template(h_signalbdE,sfact,nps_count,ename,nps_low,nps_high,
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		testhyp->add_template(h_signalbdE,sfact,nps_count,ename,nps_low,nps_high,
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		cout << "Finished setting up nullhyp, testhyp, nullhyp_pe, testhyp_pe for rare background" << endl;
	}	
	if (bdmm) {
		ename[0] = unconstrned;
		nps_low[0] = -0.0148983/0.163842;
		nps_high[0] = 0.0148983/0.163842;
		
		nps_count=2;
		pssnflg = 0;
		sclflg = 0;
		
		// Construct test/null hypothesis for pseudo-experiments.
		nullhyp_pe->add_template(h_signalbsE,sfact,nps_count,ename,nps_low,nps_high,
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		testhyp_pe->add_template(h_signalbsE,sfact,nps_count,ename,nps_low,nps_high,
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		nps_count=0;
		nullhyp->add_template(h_signalbsE,sfact,nps_count,ename,nps_low,nps_high,
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		testhyp->add_template(h_signalbsE,sfact,nps_count,ename,nps_low,nps_high,
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		cout << "Finished setting up nullhyp, testhyp, nullhyp_pe, testhyp_pe for rare background" << endl;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////

	// Add signal templates
	for(int i=0;i<NSYS;i++) {
		nps_low[i] = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i] = 0;
		lowshape[i] = 0;
		highshape[i] = 0;
	}
	
	if (bsmm) {
		ename[0] = pdfhadron;
		nps_low[0] = -0.021/0.267;
		nps_high[0] = 0.021/0.267;
		
		ename[1] = toteffname;
		nps_low[1] = -totefferr/toteff;
		nps_high[1] = totefferr/toteff;
		
		ename[2] = experrname;
		nps_low[2] = -totbswinerr/totbswin;
		nps_high[2] = totbswinerr/totbswin;

		ename[3] = accname;
		nps_low[3] = -0.011/0.229;
		nps_high[3] = 0.011/0.229;

		ename[4] = pssname;
		nps_low[4] = -0.0387952/0.714124;
		nps_high[4] = 0.0387952/0.714124;
		
		ename[5] = peakerrname;
		nps_low[5] = -rarebkgbswinerr/rarebkgbswin;
		nps_high[5] = rarebkgbswinerr/rarebkgbswin;
		
//		ename[6] = combbgerrname;
//		nps_low[6] = -combkgbserr/combkgbs;
//		nps_high[6] = combkgbserr/combkgbs;
		
		nps_count = 6;
		
		sfact = 1.;
		
		pssnflg = 0;// 
		sclflg = 1;// this is set to 1 if signal, 
		testhyp_pe->add_template(h_signalbsE,sfact,nps_count,ename,nps_low,nps_high, 
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		nps_count = 0;
		testhyp->add_template(h_signalbsE,sfact,nps_count,ename,nps_low,nps_high, 
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	}
	if (bdmm) {
		
		ename[0] = toteffname;
		nps_low[0] = -0.0002/0.0016;
		nps_high[0] = 0.0002/0.0016;
		
		ename[1] = experrname;
		nps_low[1] = -totbdwinerr/totbdwin;
		nps_high[1] = totbdwinerr/totbdwin;

		ename[2] = accname;
		nps_low[2] = -0.011/0.226;
		nps_high[2] = 0.011/0.226;

		ename[3] = pddname;
		nps_low[3] = -0.036679/0.531088;
		nps_high[3] = 0.036679/0.531088;
		
		ename[4] = peakerrname;
		nps_low[4] = -rarebkgbdwinerr/rarebkgbdwin;
		nps_high[4] = rarebkgbdwinerr/rarebkgbdwin;
		
//		ename[5] = combbgerrname;
//		nps_low[5] = -combkgbderr/combkgbd;
//		nps_high[5] = combkgbderr/combkgbd;		
		
		nps_count = 5;

		sfact = 1.;
		
		pssnflg = 0;// 
		sclflg = 1;// this is set to 1 if signal, 
		testhyp_pe->add_template(h_signalbdE,sfact,nps_count,ename,nps_low,nps_high, 
								 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
		nps_count = 0;
		testhyp->add_template(h_signalbdE,sfact,nps_count,ename,nps_low,nps_high, 
							  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	}	
	
	cout << "Finished setting up testhyp, testhyp_pe for signal in the endcaps" << endl;			
}

cout << ">>>>>> End of templates <<<<<<<" << endl;
//return 0;
