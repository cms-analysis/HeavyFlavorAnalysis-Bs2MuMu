// common source for preparing t-channel likelihood function
// limit inputs -- using SM with SM single top
// as the null hypothesis and SM with SM single top and extra anomalous
// single top as the test hypothesis

struct expectation_t {
	double value;
	double hi_err;
	double lo_err;
public:
	expectation_t() : value(0.0),hi_err(0.0),lo_err(0.0) {}
	expectation_t(double v, double hi, double lo) : value(v),hi_err(hi),lo_err(lo) {}
};
//expectation_t e;
//cout << e.value << " + " << e.hi_err << " - " << e.lo_err << endl;

//declare nusiance parameters
const int NSYS = 10;
char *ename[NSYS];
//char accname[]="ACC";// acceptance
char pdfhadron[]="HPDF";
char toteffname[]="TOTEFF";
char bkgerrname[]="BKGERR";
char combbkgerrname[]="CBKGERR";
char experrname[]="EXPERR";
//LHCb nuisance parameters
char lhcbbkgerrna[]="LHCBBKGERR";
char lhcbmisiderna[]="MISIDERR";
char lhcbsigerrna[]="SIGERR";
// Posible additional CMS parameters
//char rarename[]="RARE";
//char pssname[]="PSS";
//char psdname[]="PSD";
//char pdsname[]="PDS";
//char pddname[]="PDD";
//char peakerrname[]="PCKERR";

char channameB[] = "bsmmB";
char channameE[] = "bsmmE";

//declare the shape functions
double nps_low[NSYS];
double nps_high[NSYS];
double lowsigma[NSYS];
double highsigma[NSYS];

double sfact; // scale factor
int nps_count; // number of nuisance paramaters
int pssnflg;// Poisson flag -- 1 if Poisson, 0 if not.
int sclflg;// Scale flag -- set to 1 if you want this template to be scaled in the s95 calculation, 0 if you don't.

TH1 *lowshape[NSYS];
TH1 *highshape[NSYS];

Int_t nbins = 1;
Double_t xbins[2] = {0, 1};

// Opent the data file
//TFile *f_data = TFile::Open("anaBmm.default-11.root");
TString s_outfilename;
if (combined) s_outfilename = Form("CMS_LHCb_Comb_Results%i.root",i_numb);
else if(cmsbs) s_outfilename = Form("CMS_Results%i.root",i_numb);
else if(lhcbs) s_outfilename = Form("LHCb_Comb_Results%i.root",i_numb);
else if(lhcbs10) s_outfilename = Form("LHCb_10_Results%i.root",i_numb);
else if(lhcbs11) s_outfilename = Form("LHCb_11_Results%i.root",i_numb);
cout << s_outfilename << endl;
TFile f_outfile(s_outfilename,"RECREATE");

TH1D *h_rarebkgB = new TH1D("h_rarebkgB","Rare bkgd brrl",nbins,xbins); 
TH1D *h_combbkgB = new TH1D("h_combbkgB","Comb bkgd brrl",nbins,xbins); 
TH1D *h_signalbsB = new TH1D("h_signalbsB","Sigl bs brrl",nbins,xbins); 
TH1D *h_dataB = new TH1D("h_dataB","Obsd evts brrl",nbins,xbins);
TH1D *h_rarebkgE	= new TH1D("h_rarebkgE","Rare bkgd endc",nbins,xbins); 
TH1D *h_combbkgE = new TH1D("h_combbkgE","Comb bkgd endc",nbins,xbins); 
TH1D *h_signalbsE = new TH1D("h_signalbsE","Sigl bs endc",nbins,xbins); 
TH1D *h_dataE = new TH1D("h_dataE","Obsd evts in endc",nbins,xbins);
// Setting up LHCb histos for GL bin1

cout << "<<<<<< Done reading data >>>>>>>" << endl;

//========================================
//Construct test/null hypothesis for fitting
csm_model* nullhyp = new csm_model();
csm_model* testhyp = new csm_model();
csm_model* nullhyp_pe = new csm_model();
csm_model* testhyp_pe = new csm_model();


if (barrel) {
	Double_t toteff = 0.0036;
	Double_t totefferr = 0.0004;
	
	Double_t rarebkgbswin = 0.07; 
	Double_t rarebkgbswinerr = 0.02;
	cout << "Rare bkgd (peaking in Bs signal window) = " << rarebkgbswin << endl << endl;
	
	Double_t combkgbs = 0.60; 
	Double_t combkgbserr = 0.35;
	cout << "CombBkgd (in Bs blind window) = " << combkgbs << endl << endl;
	
	Double_t exptbsevnts = 0.80; 
	Double_t exptbsevntserr = 0.16;
	cout << "Expected Bs signal events = " << exptbsevnts << endl;
	
	Double_t obsrvbswin = 2; 
	cout << "Observed in Bs window = " << obsrvbswin << endl << endl;
		
	// Setup the Histos
	//Rare bkgd Histogram
	h_rarebkgB->SetBinContent(1,rarebkgbswin); 
	
	//Combinatorial bkgd Histogram
	h_combbkgB->SetBinContent(1,combkgbs); 
	
	//Signal Histograms		
	h_signalbsB->SetBinContent(1,exptbsevnts);
	
	// Data histogram
	h_dataB->SetBinContent(1,obsrvbswin);
	
	cout << "<<<<<< Done Setting up Histos for Barrel Channel>>>>>>>" << endl;
	
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
	
	// Add rare background templates
	sfact = 1;
	
	ename[0] = bkgerrname;
	nps_low[0] = -rarebkgbswinerr/rarebkgbswin;
	nps_high[0] = rarebkgbswinerr/rarebkgbswin;
	
//	ename[1] = psdname;
//	nps_low[1] = -0.0144844/0.289687;
//	nps_high[1] = 0.0144844/0.289687;
	
	nps_count=1;	
		
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
	
	ename[0] = combbkgerrname;
	nps_low[0] = -combkgbserr/combkgbs;
	nps_high[0] = combkgbserr/combkgbs;
//	nps_low[0] = 0;
//	nps_high[0] = 0;
	
	nps_count=1;	
	
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
	
	
	// Add signal templates
	for(int i=0;i<NSYS;i++) {
		nps_low[i] = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i] = 0;
		lowshape[i] = 0;
		highshape[i] = 0;
	}
	
	
	ename[0] = pdfhadron;
	nps_low[0] = -0.021/0.267;
	nps_high[0] = 0.021/0.267;
	//		nps_low[0] = -0.1192/3.5487;
	//		nps_high[0] = 0.1192/3.5487;
	
	ename[1] = toteffname;
	nps_low[1] = -totefferr/toteff;
	nps_high[1] = totefferr/toteff;
	
	ename[2] = experrname;
	nps_low[2] = -exptbsevntserr/exptbsevnts;
	nps_high[2] = exptbsevntserr/exptbsevnts;
	
//	ename[3] = pssname;
//	nps_low[3] = -0.043758/0.875161;
//	nps_high[3] = 0.043758/0.875161;
	
//	ename[5] = peakerrname;
//	nps_low[5] = -0.0569132/0.2438;
//	nps_high[5] = 0.0569132/0.2438;
	
//	ename[3] = accname;
//	nps_low[3] = -0.009/0.248;
//	nps_high[3] = 0.009/0.248;
	
	nps_count = 3;
	
	sfact = 1.;
	
	pssnflg = 0;// 
	sclflg = 1;// this is set to 1 if signal, 
	testhyp_pe->add_template(h_signalbsB,sfact,nps_count,ename,nps_low,nps_high, 
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	nps_count = 0;
	testhyp->add_template(h_signalbsB,sfact,nps_count,ename,nps_low,nps_high, 
						  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameB);
	
	
	cout << "Finished setting up testhyp, testhyp_pe for signal in the barrel" << endl;	
}

if (endcap) {
	
	Double_t toteff = 0.0021;
	Double_t totefferr = 0.0002;
	
	Double_t rarebkgbswin = 0.04; 
	Double_t rarebkgbswinerr = 0.01;
	cout << "Rare bkgd (peaking in Bs signal window) = " << rarebkgbswin << endl << endl;
	
	Double_t combkgbs = 0.8; 
	Double_t combkgbserr = 0.4;
	cout << "CombBkgd (in Bs blind window) = " << combkgbs << endl << endl;
	
	Double_t exptbsevnts = 0.36; 
	Double_t exptbsevntserr = 0.07;
	cout << "Expected Bs signal events = " << exptbsevnts << endl;
	
	Double_t obsrvbswin = 1; 
	cout << "Observed in Bs window = " << obsrvbswin << endl << endl;
	
	
	// Setup the Histos
	//Rare bkgd Histogram
	h_rarebkgE->SetBinContent(1,rarebkgbswin); 
	
	//Combinatorial bkgd Histogram
	h_combbkgE->SetBinContent(1,combkgbs); 
	
	//Signal Histograms
	h_signalbsE->SetBinContent(1,exptbsevnts);
	
	// Data histogram
	h_dataE->SetBinContent(1,obsrvbswin);
	
	cout << "<<<<<< Done Setting up Histos for Endcap Channel>>>>>>>" << endl;
	
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
	
	// Add rare background templates
	sfact = 1;
	
	ename[0] = bkgerrname;
	nps_low[0] = -rarebkgbswinerr/rarebkgbswin;
	nps_high[0] = rarebkgbswinerr/rarebkgbswin;
	
	
//	ename[1] = psdname;
//	nps_low[1] = -0.0180617/0.361233;
//	nps_high[1] = 0.0180617/0.361233;
	
	nps_count=1;	
	
	
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
	
	ename[0] = combbkgerrname;
	nps_low[0] = -combkgbserr/combkgbs;
	nps_high[0] = combkgbserr/combkgbs;
	
	nps_count=1;
	
	
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
	
	// Add signal templates
	for(int i=0;i<NSYS;i++) {
		nps_low[i] = 0;
		nps_high[i] = 0;
		lowsigma[i] = 0;
		highsigma[i] = 0;
		lowshape[i] = 0;
		highshape[i] = 0;
	}
	
	
	ename[0] = pdfhadron;
	nps_low[0] = -0.021/0.267;
	nps_high[0] = 0.021/0.267;
	//		nps_low[0] = -0.1192/3.5487;
	//		nps_high[0] = 0.1192/3.5487;
	
	ename[1] = toteffname;
	nps_low[1] = -totefferr/toteff;
	nps_high[1] = totefferr/toteff;
	
	ename[2] = experrname;
	nps_low[2] = -exptbsevntserr/exptbsevnts;
	nps_high[2] = exptbsevntserr/exptbsevnts;
	
//	ename[3] = pssname;
//	nps_low[3] = -0.0348352/0.696703;
//	nps_high[3] = 0.0348352/0.696703;
	
//	ename[5] = peakerrname;
//	nps_low[5] = -0.0228963/0.1253;
//	nps_high[5] = 0.0228963/0.1253;
//	
//	ename[3] = accname;
//	nps_low[3] = -0.011/0.229;
//	nps_high[3] = 0.011/0.229;
	
	nps_count = 3;
	
	sfact = 1.;
	
	pssnflg = 0;// 
	sclflg = 1;// this is set to 1 if signal, 
	testhyp_pe->add_template(h_signalbsE,sfact,nps_count,ename,nps_low,nps_high, 
							 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	nps_count = 0;
	testhyp->add_template(h_signalbsE,sfact,nps_count,ename,nps_low,nps_high, 
						  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,channameE);
	
	
	cout << "Finished setting up testhyp, testhyp_pe for signal in the endcaps" << endl;			
}

//Setting up the LHCb 2010 data
const static expectation_t bkg1[N_MASS_BINS][N_BDT_BINS] = {
	{expectation_t(56.9,1.1,1.1), expectation_t(1.31,0.19,0.17), expectation_t(0.282,0.076,0.065), expectation_t(0.016,0.021,0.010)},
	{expectation_t(56.1,1.1,1.1), expectation_t(1.28,0.18,0.17), expectation_t(0.269,0.072,0.062), expectation_t(0.0151,0.0195,0.0094)},
	{expectation_t(55.3,1.1,1.1), expectation_t(1.24,0.17,0.16), expectation_t(0.257,0.069,0.059), expectation_t(0.0139,0.0179,0.0086)},
	{expectation_t(54.4,1.1,1.1), expectation_t(1.21,0.17,0.16), expectation_t(0.246,0.066,0.057), expectation_t(0.0128,0.0165,0.0080)},
	{expectation_t(53.6,1.1,1.0), expectation_t(1.18,0.17,0.15), expectation_t(0.235,0.063,0.054), expectation_t(0.0118,0.0152,0.0073)},
	{expectation_t(52.8,1.0,1.0), expectation_t(1.14,0.16,0.15), expectation_t(0.224,0.060,0.052), expectation_t(0.0108,0.0140,0.0068)}
};

const static expectation_t sig1[N_MASS_BINS][N_BDT_BINS] = {
	{expectation_t(0.0076,0.0034,0.0030),expectation_t(0.0050,0.0027,0.0020),expectation_t(0.0037,0.0015,0.0011),expectation_t(0.0047,0.0015,0.0010)},
	{expectation_t(0.0220,0.0084,0.0081),expectation_t(0.0146,0.0067,0.0054),expectation_t(0.0107,0.0036,0.0027),expectation_t(0.0138,0.0035,0.0025)},
	{expectation_t(0.0380,0.0150,0.0150),expectation_t(0.0250,0.0120,0.0100),expectation_t(0.0183,0.0063,0.0047),expectation_t(0.0235,0.0060,0.0044)},
	{expectation_t(0.0380,0.0150,0.0150),expectation_t(0.0250,0.0120,0.0100),expectation_t(0.0183,0.0063,0.0047),expectation_t(0.0235,0.0060,0.0044)},
	{expectation_t(0.0220,0.0084,0.0081),expectation_t(0.0146,0.0067,0.0054),expectation_t(0.0107,0.0036,0.0027),expectation_t(0.0138,0.0035,0.0025)},
	{expectation_t(0.0076,0.0031,0.0027),expectation_t(0.0050,0.0025,0.0019),expectation_t(0.0037,0.0013,0.0010),expectation_t(0.0047,0.0013,0.0010)}
};

const static int obs1[N_MASS_BINS][N_BDT_BINS] = {
	{39,2,1,0},
	{55,2,0,0},
	{73,0,0,0},
	{60,0,0,0},
	{53,2,0,0},
	{55,1,0,0}
};
// Histos to hold the 2010 data 
TH1D *h_lhcbbkg1[N_MASS_BINS][N_BDT_BINS];
TH1D *h_lhcbsig1[N_MASS_BINS][N_BDT_BINS];
TH1D *h_lhcbdat1[N_MASS_BINS][N_BDT_BINS];
string s_lhcbchan1[N_MASS_BINS][N_BDT_BINS];

//Setting up the LHCb summer11 data
const static expectation_t bkg[N_MASS_BINS][N_BDT_BINS] = {
	{expectation_t(514,12,11), expectation_t(4.32,.39,.39), expectation_t(0.504,0.158,0.095), expectation_t(0.118,0.078,0.039)},
	{expectation_t(506,12,11), expectation_t(4.25,0.38,0.38), expectation_t(0.502,0.157,0.094), expectation_t(0.115,0.076,0.038)},
	{expectation_t(499,11,11), expectation_t(4.19,0.38,0.38), expectation_t(0.499,0.156,0.094), expectation_t(0.112,0.074,0.037)},
	{expectation_t(491,11,11), expectation_t(4.13,0.37,0.37), expectation_t(0.496,0.155,0.093), expectation_t(0.109,0.072,0.036)},
	{expectation_t(483,11,10), expectation_t(4.07,0.37,0.37), expectation_t(0.494,0.154,0.093), expectation_t(0.106,0.070,0.035)},
	{expectation_t(476,11,10), expectation_t(4.01,0.36,0.36), expectation_t(0.491,0.154,0.092), expectation_t(0.103,0.069,0.034)}
};

const static expectation_t misid[N_MASS_BINS][N_BDT_BINS] = {
	{expectation_t(0.52,0.066,0.045), expectation_t(0.052,0.065,0.045), expectation_t(0.050,0.065,-0.045), expectation_t(0.052,0.066,0.046)},
	{expectation_t(0.029,0.028,0.020), expectation_t(0.028,0.027,0.018), expectation_t(0.028,0.027,0.019), expectation_t(0.028,0.028,0.020)},
	{expectation_t(0.019,0.023,0.021), expectation_t(0.019,0.022,0.020), expectation_t(0.020,0.023,0.022), expectation_t(0.019,0.022,0.020)},
	{expectation_t(0.014,0.017,0.016), expectation_t(0.014,0.017,0.016), expectation_t(0.014,0.017,0.016), expectation_t(0.015,0.019,0.017)},
	{expectation_t(0.011,0.014,0.012), expectation_t(0.012,0.015,0.013), expectation_t(0.011,0.014,0.013), expectation_t(0.011,0.014,0.013)},
	{expectation_t(0.0085,0.0109,0.0098), expectation_t(0.0084,0.0108,0.0097), expectation_t(0.0082,0.0108,0.0096), expectation_t(0.009,0.011,0.010)}
};

const static expectation_t sig[N_MASS_BINS][N_BDT_BINS] = {
	{expectation_t(0.058,0.016,0.014),expectation_t(0.0280,0.0096,0.0075),expectation_t(0.0306,0.0074,0.0057),expectation_t(0.0332,0.0079,0.0061)},
	{expectation_t(0.199,0.046,0.044),expectation_t(0.097,0.030,0.024),expectation_t(0.106,0.021,0.016),expectation_t(0.114,0.023,0.018)},
	{expectation_t(0.371,0.084,0.081),expectation_t(0.181,0.056,0.044),expectation_t(0.197,0.039,0.029),expectation_t(0.214,0.043,0.032)},
	{expectation_t(0.371,0.085,0.080),expectation_t(0.181,0.056,0.045),expectation_t(0.197,0.039,0.029),expectation_t(0.214,0.043,0.032)},
	{expectation_t(0.199,0.047,0.044),expectation_t(0.097,0.030,0.024),expectation_t(0.106,0.021,0.016),expectation_t(0.114,0.023,0.017)},
	{expectation_t(0.057,0.017,0.014),expectation_t(0.0276,0.0095,0.0074),expectation_t(0.0302,0.0077,0.0058),expectation_t(0.0327,0.0083,0.0064)}
};

const static int obs[N_MASS_BINS][N_BDT_BINS] = {
	{486,5,1,0},
	{483,3,0,1},
	{511,6,1,1},
	{472,3,0,0},
	{484,4,1,0},
	{436,5,0,0}
};
// Histos to hold the summer 2011 data 
TH1D *h_lhcbbkg[N_MASS_BINS][N_BDT_BINS];
TH1D *h_lhcbmisid[N_MASS_BINS][N_BDT_BINS];
TH1D *h_lhcbsig[N_MASS_BINS][N_BDT_BINS];
TH1D *h_lhcbdat[N_MASS_BINS][N_BDT_BINS];
string s_lhcbchan[N_MASS_BINS][N_BDT_BINS];	


int counter1 = 0;
if (lhcbs10) {
	for (int j = 0; j<N_MASS_BINS; j++) {
		for (int l = 0; l<N_BDT_BINS; l++) {
			s_lhcbchan1[j][l] = Form("lhcbChanA%i%i",j,l);
			// Make a char* from the string
			char * c_lhcbchan1 = new char[s_lhcbchan1[j][l].size() + 1];
			copy(s_lhcbchan1[j][l].begin(), s_lhcbchan1[j][l].end(), c_lhcbchan1);
			c_lhcbchan1[s_lhcbchan1[j][l].size()] = '\0'; 
			
			h_lhcbbkg1[j][l] = new TH1D(Form("h_lhcbbkg1%i%i",j,l),Form("LHCb1 Bkgd M_bn%i G_bn%i",j,l),nbins,xbins); 
			h_lhcbbkg1[j][l]->SetBinContent(1,bkg1[j][l].value);
			
			h_lhcbsig1[j][l] = new TH1D(Form("h_lhcbsig1%i%i",j,l),Form("LHCb1 Sgnl M_bn%i G_bn%i",j,l),nbins,xbins);
			h_lhcbsig1[j][l]->SetBinContent(1,sig1[j][l].value);
			
			h_lhcbdat1[j][l] = new TH1D(Form("h_lhcbdat1%i%i",j,l),Form("LHCb1 Obsd M_bn%i G_bn%i",j,l),nbins,xbins);
			h_lhcbdat1[j][l]->SetBinContent(1,obs1[j][l]);
			
			// Add background templates
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
			// The error may be correlated with the 2011 data set
//			ename[0] = lhcbbkgerrna1; 
			ename[0] = lhcbbkgerrna;
			nps_low[0] = -bkg1[j][l].lo_err/bkg1[j][l].value;
			nps_high[0] = bkg1[j][l].hi_err/bkg1[j][l].value;
			
			nps_count=1;		
			pssnflg = 0;
			sclflg = 0;
			
			// Construct test/null hypothesis for pseudo-experiments.
			nullhyp_pe->add_template(h_lhcbbkg1[j][l],sfact,nps_count,ename,nps_low,nps_high,
									 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan1);
			testhyp_pe->add_template(h_lhcbbkg1[j][l],sfact,nps_count,ename,nps_low,nps_high,
									 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan1);
			nps_count=0;
			nullhyp->add_template(h_lhcbbkg1[j][l],sfact,nps_count,ename,nps_low,nps_high,
								  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan1);
			testhyp->add_template(h_lhcbbkg1[j][l],sfact,nps_count,ename,nps_low,nps_high,
								  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan1);
			cout << "Finished setting up the background nullhyp, testhyp, nullhyp_pe, and testhyp_pe for channel " << counter1 << endl;
			
			// Add signal templates
			for(int i=0;i<NSYS;i++) {
				nps_low[i] = 0;
				nps_high[i] = 0;
				lowsigma[i] = 0;
				highsigma[i] = 0;
				lowshape[i] = 0;
				highshape[i] = 0;
			}
			// Add GLbin1 misid background templates
			sfact = 1;	
			
			ename[0] = pdfhadron;
			nps_low[0] = -0.021/0.267;
			nps_high[0] = 0.021/0.267;
			// The error may be correlated with 2011 data set
			ename[1] = lhcbsigerrna;
			nps_low[1] = -sig1[j][l].lo_err/sig1[j][l].value;
			nps_high[1] = sig1[j][l].hi_err/sig1[j][l].value;
			
			nps_count = 2;
			
			pssnflg = 0;// 
			sclflg = 1;// this is set to 1 if signal, 
			testhyp_pe->add_template(h_lhcbsig1[j][l],sfact,nps_count,ename,nps_low,nps_high, 
									 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan1);
			nps_count = 0;
			testhyp->add_template(h_lhcbsig1[j][l],sfact,nps_count,ename,nps_low,nps_high, 
								  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan1);
			
			cout << "Finished setting up the signal testhyp, testhyp_pe for channel " << counter1 << endl;	
			
			delete[] c_lhcbchan1;
			counter1++;
			
		}
	}
	cout << "Fnished setting up LHCb 2010 results, total number of channels = " << counter1 << endl;
}

int counter = 0;
if (lhcbs11) {
	for (int j = 0; j<N_MASS_BINS; j++) {
		for (int l = 0; l<N_BDT_BINS; l++) {
			//		cout << "FUBAR1" << endl;
			s_lhcbchan[j][l] = Form("lhcbChan%i%i",j,l);
			// Make a char* from the string
			char * c_lhcbchan = new char[s_lhcbchan[j][l].size() + 1];
			copy(s_lhcbchan[j][l].begin(), s_lhcbchan[j][l].end(), c_lhcbchan);
			c_lhcbchan[s_lhcbchan[j][l].size()] = '\0'; 
			
			// Make the histograms
			h_lhcbbkg[j][l] = new TH1D(Form("h_lhcbbkg%i%i",j,l),Form("LHCb Bkgd M_bn%i G_bn%i",j,l),nbins,xbins); 
			h_lhcbbkg[j][l]->SetBinContent(1,bkg[j][l].value);
			
			h_lhcbmisid[j][l] = new TH1D(Form("h_lhcbmisid%i%i",j,l),Form("LHCb msID M_bn%i G_bn%i",j,l),nbins,xbins);
			h_lhcbmisid[j][l]->SetBinContent(1,misid[j][l].value);
			
			h_lhcbsig[j][l] = new TH1D(Form("h_lhcbsig%i%i",j,l),Form("LHCb Sgnl M_bn%i G_bn%i",j,l),nbins,xbins);
			h_lhcbsig[j][l]->SetBinContent(1,sig[j][l].value);
			
			h_lhcbdat[j][l] = new TH1D(Form("h_lhcbdat%i%i",j,l),Form("LHCb Obsd M_bn%i G_bn%i",j,l),nbins,xbins);
			h_lhcbdat[j][l]->SetBinContent(1,obs[j][l]);
			
			// Add background templates
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
			
			ename[0] = lhcbbkgerrna;
			nps_low[0] = -bkg[j][l].lo_err/bkg[j][l].value;
			nps_high[0] = bkg[j][l].hi_err/bkg[j][l].value;
			
			nps_count=1;		
			pssnflg = 0;
			sclflg = 0;
			
			// Construct test/null hypothesis for pseudo-experiments.
			nullhyp_pe->add_template(h_lhcbbkg[j][l],sfact,nps_count,ename,nps_low,nps_high,
									 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			testhyp_pe->add_template(h_lhcbbkg[j][l],sfact,nps_count,ename,nps_low,nps_high,
									 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			nps_count=0;
			nullhyp->add_template(h_lhcbbkg[j][l],sfact,nps_count,ename,nps_low,nps_high,
								  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			testhyp->add_template(h_lhcbbkg[j][l],sfact,nps_count,ename,nps_low,nps_high,
								  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			cout << "Finished setting up the background nullhyp, testhyp, nullhyp_pe, and testhyp_pe for channel " << counter << endl;
			
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
			
			// Add GLbin1 misid background templates
			sfact = 1;
			
			ename[0] = lhcbmisiderna;
			nps_low[0] = -misid[j][l].lo_err/misid[j][l].value;
			nps_high[0] = misid[j][l].hi_err/misid[j][l].value;
			
			nps_count=1;		
			pssnflg = 0;
			sclflg = 0;
			
			// Construct test/null hypothesis for pseudo-experiments.
			nullhyp_pe->add_template(h_lhcbmisid[j][l],sfact,nps_count,ename,nps_low,nps_high,
									 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			testhyp_pe->add_template(h_lhcbmisid[j][l],sfact,nps_count,ename,nps_low,nps_high,
									 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			nps_count=0;
			nullhyp->add_template(h_lhcbmisid[j][l],sfact,nps_count,ename,nps_low,nps_high,
								  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			testhyp->add_template(h_lhcbmisid[j][l],sfact,nps_count,ename,nps_low,nps_high,
								  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			cout << "Finished setting up the misID bkgnd nullhyp, testhyp, nullhyp_pe, testhyp_pe for channel " << counter << endl;
			
			// Add signal templates
			for(int i=0;i<NSYS;i++) {
				nps_low[i] = 0;
				nps_high[i] = 0;
				lowsigma[i] = 0;
				highsigma[i] = 0;
				lowshape[i] = 0;
				highshape[i] = 0;
			}
			// Add GLbin1 misid background templates
			sfact = 1;	
			
			ename[0] = pdfhadron;
			nps_low[0] = -0.021/0.267;
			nps_high[0] = 0.021/0.267;
			
			ename[1] = lhcbsigerrna;
			nps_low[1] = -sig[j][l].lo_err/sig[j][l].value;
			nps_high[1] = sig[j][l].hi_err/sig[j][l].value;
			
			nps_count = 2;
			
			pssnflg = 0;// 
			sclflg = 1;// this is set to 1 if signal, 
			testhyp_pe->add_template(h_lhcbsig[j][l],sfact,nps_count,ename,nps_low,nps_high, 
									 lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			nps_count = 0;
			testhyp->add_template(h_lhcbsig[j][l],sfact,nps_count,ename,nps_low,nps_high, 
								  lowshape,lowsigma,highshape,highsigma,pssnflg,sclflg,c_lhcbchan);
			
			cout << "Finished setting up the signal testhyp, testhyp_pe for channel " << counter << endl;	
			
			delete[] c_lhcbchan;
			counter++;
		}
	}
	cout << "Fnished setting up LHCb summer 2011 results, total number of channels = " << counter << endl;
}

//return 0;

cout << ">>>>>> End of templates <<<<<<<" << endl;
