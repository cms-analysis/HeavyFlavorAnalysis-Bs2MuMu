#include "anaBmm.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "AnalysisDistribution.hh"
#include "initFunc.hh"
#include "mclimit_csm.hh"
#include "bayesianlimit.hh"

#include "TF1.h"
#include "THStack.h"
#include "TKey.h"
#include "TTree.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVirtualPad.h"  // access to gPad
#include "TMinuit.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"

#include <iomanip>
#include <string>
#include <list>

// Usage: root[0] anaBmm(

using namespace std; 
using std::string; 

ClassImp(anaBmm)

// ----------------------------------------------------------------------
anaBmm::anaBmm(const char *files, const char *cuts, const char *dir, int mode) { 

  fDoPrint = true; 

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

  mclimit_csm *p = new mclimit_csm();
  p->print_version();

  double beta  = 0.9;
  double e0    = 1.; 
  double esig  = 0.25; 
  double alpha = 1.; 
  int nobs = 0;
  double b0 = 0.2;
  double bsig = 0.2;

  double ulmean =  blimit(beta, nobs, e0, esig, b0, bsig, alpha);
  cout << "Bayesian blimit: " << ulmean << endl;

  init(files, cuts, dir, mode);
}

// ----------------------------------------------------------------------
anaBmm::~anaBmm() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void anaBmm::init(const char *files, const char *cuts, const char *dir, int mode) {

  double confLevel = 0.9;
  
  fBF = 0.0593*1.014e-3;
  // PDG 2010:
  fu  = 0.402;
  fs  = 0.105;
  // LHCb:
  //   fu = 0.404; 
  //   fs = 0.109;

  fpFunc  = new initFunc(); 
  fpRolke = new TRolke(confLevel);

  fMassLo = 4.5; 
  fMassHi = 6.5;

  fSigLo = 5.27;
  fSigHi = 5.47;

  fNormLo = 5.1;
  fNormHi = 5.4;

  fCsLo = 5.27;
  fCsHi = 5.47;

  fBgLo = 4.8;
  fBgHi = 6.0;

  fNData = fNMc = 0; 
  fSgData = fSgMc = fNoData = fNoMc = -1; 
  fCsData = fCsMc = -1; 

  // -- initialize cuts
  cout << "Reading cuts from " << Form("anaBmm.%s.cuts", cuts) << endl;
  readCuts(Form("anaBmm.%s.cuts", cuts)); 

  printCuts(cout); 

  int HBINS(15); 
  double HLO(0.), HHI(45.); 

  numbers *a(0); 
  fNchan = fCuts.size(); 
  fChan = -1; 
  int NBINS = (fMassHi - fMassLo)/0.025+1;
  TH1D *h; 
  for (unsigned int i = 0; i < fNchan; ++i) {
    // -- signal Bs2MuMu
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("signal Bs2MuMu %i", i); 
    //    a->effGenFilter  = 0.63;
    //    a->effGenFilter  = 1.0;
    //    a-> effGenFilterE = 0.03;
    fNumbersBs.push_back(a); 

    // -- signal Bd2MuMu
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("signal Bd2MuMu %i", i); 
    //    a->effGenFilter  = 0.63;
    //    a->effGenFilter  = 1.0;
    //    a-> effGenFilterE = 0.03;
    fNumbersBd.push_back(a); 

    // --  normalization
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name = Form("normalization %i", i); 
    //    a->effGenFilter  = 0.24;
    //    a->effGenFilter  = 1.0;
    //    a->effGenFilterE = 0.013;
    fNumbersNorm.push_back(a); 

    // --  control sample
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("control sample %i", i); 
    //    a->effGenFilter  = 0.22;
    //    a->effGenFilter  = 1.0;
    //    a-> effGenFilterE = 0.015;
    fNumbersCS.push_back(a); 

    h = new TH1D(Form("hMassWithMassCuts%d", i), Form("hMassWithMassCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMassCuts.push_back(h); 
    h = new TH1D(Form("fhMassWithMassCutsManyBins%d", i), Form("fhMassWithMassCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithMassCutsManyBins.push_back(h); 
  
    h = new TH1D(Form("hMassWithAllCuts%d", i), Form("hMassWithAllCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCuts.push_back(h); 
    h = new TH1D(Form("hMassWithAllCutsManyBins%d", i), Form("hMassWithAllCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithAllCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassWithTriggerCuts%d", i), Form("hMassWithTriggerCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithTriggerCuts.push_back(h); 
    h = new TH1D(Form("hMassWithTriggerCutsManyBins%d", i), Form("hMassWithTriggerCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithTriggerCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassWithMuonCuts%d", i), Form("hMassWithMuonCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMuonCuts.push_back(h); 
    h = new TH1D(Form("hMassWithMuonCutsManyBins%d", i), Form("hMassWithMuonCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithMuonCutsManyBins.push_back(h); 

    h = new TH1D(Form("fhMassWithAnaCuts%d", i), Form("hMassChan%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAnaCuts.push_back(h); 
    h = new TH1D(Form("fhMassWithAnaCutsManyBins%d", i), Form("hMassChan%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithAnaCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassNoCuts%d", i), Form("hMassNoCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassNoCuts.push_back(h); 

    h = new TH1D(Form("hMassAbsNoCuts%d", i), Form("hMassAbsNoCuts%d", i), 100, 0, 10);
    fhMassAbsNoCuts.push_back(h); 

    h = new TH1D(Form("hMuId%d", i), Form("hMuId%d", i), 100, 0., 1.);
    fhMuId.push_back(h); 

    h = new TH1D(Form("hMuTr%d", i), Form("hMuTr%d", i), 100, 0., 1.);
    fhMuTr.push_back(h); 

    h = new TH1D(Form("h0PidTrigger%d", i), Form("hPidTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0PidTrigger.push_back(h); 
    h = new TH1D(Form("h1PidTrigger%d", i), Form("hPidTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1PidTrigger.push_back(h); 
    h = new TH1D(Form("h0PidMuID%d", i), Form("hPidMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0PidMuID.push_back(h); 
    h = new TH1D(Form("h1PidMuID%d", i), Form("hPidMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1PidMuID.push_back(h); 

    h = new TH1D(Form("h0MCTrigger%d", i), Form("hMCTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0MCTrigger.push_back(h); 
    h = new TH1D(Form("h1MCTrigger%d", i), Form("hMCTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1MCTrigger.push_back(h); 
    h = new TH1D(Form("h0MCMuID%d", i), Form("hMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0MCMuID.push_back(h); 
    h = new TH1D(Form("h1MCMuID%d", i), Form("hMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1MCMuID.push_back(h); 
  }

  cout << "----------------------------" << endl;
  cout << "fNchan = " << fNchan << endl;
  cout << "----------------------------" << endl;

  fFont = 42; 
  fMode = mode;  

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (c0 == 0) {
    cout << "TCanvas c0 not found. Creating my own version." << endl;
    c0 = new TCanvas("c0","--c0--",356,0,656,700);
  }
  tl  = new TLatex();
  tl->SetNDC(kTRUE); 
  tl->SetTextSize(0.07);
  tl->SetTextFont(fFont);

  pl  = new TLine();
  pa  = new TArrow();
  box = new TBox();

  f0 = new TF1("f0", f_p1, 0., 7., 2); 

  f1 = new TF1("f1", f_expo, 0., 7.0, 2);

  f2 = new TF1("f2", f_p1aG, 0., 7., 5);
  f2->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f3 = new TF1("f3", f_eaG, 0., 7., 5);
  f3->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f4 = new TF1("f4", f_ea2G, 0., 7., 8);
  f4->SetParNames("Area", "Peak", "Sigma", "Fraction2", "Peak2", "Sigma2", "Offset", "Slope"); 

  fDirectory = dir; 
  fSuffix    = cuts; 
  if (fMode > 0)  fSuffix += Form("-%d", fMode); 
  cout << "--> Dumping output into " << fDirectory << endl;
  fNumbersFileName = fDirectory + "/anaBmm." + fSuffix + ".txt";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fOUT.open(fNumbersFileName.c_str(), ios::app);

  fNumbersFileName = fDirectory + "/anaBmm." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);
  

  fUlcalcFileName  = fDirectory + "/anaBmm." + fSuffix + ".ulc";
  system(Form("/bin/rm -f %s", fUlcalcFileName.c_str()));

  loadFiles(files);
  string hfname  = fDirectory + "/anaBmm." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

  printCuts(fOUT); 

  dumpSamples();
  fF["SgMc"]->cd();
  dumpCutNames();

}


// ----------------------------------------------------------------------
void anaBmm::loadFiles(const char *files) {

  char buffer[1000];
  ifstream is(files);
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    
    string sbuffer = string(buffer); 

    string::size_type m1 = sbuffer.find("dset="); 
    string stype = sbuffer.substr(5, m1-6); 

    string::size_type m2 = sbuffer.find("lumi="); 
    string sdset = sbuffer.substr(m1+5, m2-m1-6); 
    m1 = m2; 

    m2 = sbuffer.find("file="); 
    string slumi = sbuffer.substr(m1+5, m2-m1-6); 
    string sfile = sbuffer.substr(m2+5); 
    string sname; 
    
    // -- DATA
    if (string::npos != sdset.find("data")) {
      fpData[fNData] = loadFile(sfile, stype); 
      fDataLumi[fNData] = atof(slumi.c_str()); 
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	fSgData = fNData;
	sname = "SgData"; 
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-}")); 
      }
      if (string::npos != stype.find("2010") && string::npos != stype.find("sg")) {
	sname = "SgData10";
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-} (2010)")); 
      }
      if (string::npos != stype.find("2011") && string::npos != stype.find("sg")) {
	sname = "SgData11";
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-} (2011)")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("no")) {
	fNoData = fNData;
	sname = "NoData";
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-}K^{+} (2011)")); 
      }
      if (string::npos != stype.find("2010") && string::npos != stype.find("no")) {
	sname = "NoData10"; 
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-}K^{+} (2010)")); 
      }
      if (string::npos != stype.find("2011") && string::npos != stype.find("no")) {
	sname = "NoData11"; 
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-}K^{+} (2011)")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	fCsData = fNData;
	sname = "CsData"; 
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-}K^{+}K^{-} (2011)")); 
      }
      if (string::npos != stype.find("2010") && string::npos != stype.find("cs")) {
	sname = "CsData10"; 
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-}K^{+}K^{-} (2010)")); 
      }
      if (string::npos != stype.find("2011") && string::npos != stype.find("cs")) {
	sname = "CsData11"; 
	fF.insert(make_pair(sname, fpData[fNData])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#mu^{+}#mu^{-}K^{+}K^{-} (2011)")); 
      }
      cout << "open data " << fNData << " " << sfile << " as " << sname << " (" << stype << ") with lumi = " << fDataLumi[fNData] << endl;
      ++fNData;
    } else {
      string sfilter = sdset; 
      replaceAll(sfilter, "mc,", ""); 
      double effFilter = atof(sfilter.c_str());
      // -- MC
      fpMc[fNMc] = loadFile(sfile, stype); 
      fMcLumi[fNMc] = atof(slumi.c_str()); 
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	fSgMc = fNMc;
	sname = "SgMc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 	
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("sg")) {
	sname = "SgMcAcc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (Spring11)")); 
      }
      if (string::npos != stype.find("Summer11") && string::npos != stype.find("sg")) {
	sname = "SgMcSummer11";
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (Summer11)")); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("sg")) {
	sname = "SgMc1e33"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC 1e33)")); 
      }
      if (string::npos != stype.find("5e32") && string::npos != stype.find("sg")) {
	sname = "SgMc5e32"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC 5e32)")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("bd")) {
	fBdMc = fNMc;
	sname = "BdMc"; 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("bd")) {
	sname = "BdMcAcc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (Spring11)")); 
      }
      if (string::npos != stype.find("Summer11") && string::npos != stype.find("bd")) {
	sname = "BdMcSummer11"; 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (Summer11)")); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("bd")) {
	sname = "BdMc1e33"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (MC 1e33)")); 
      }
      if (string::npos != stype.find("5e32") && string::npos != stype.find("bd")) {
	sname = "BdMc5e32"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (MC 5e32)")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("no")) {
	fNoMc = fNMc;
	sname = "NoMc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #mu^{+}#mu^{-}K^{+} (MC)")); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("no")) {
	sname = "NoMc1e33"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #mu^{+}#mu^{-}K^{+} (MC 1e33)")); 
      }
      if (string::npos != stype.find("5e32") && string::npos != stype.find("no")) {
	sname = "NoMc5e32"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #mu^{+}#mu^{-}K^{+} (MC 5e32)")); 
      }
      if (string::npos != stype.find("Summer11") && string::npos != stype.find("no")) {
	sname = "NoMcSummer11"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #mu^{+}#mu^{-}K^{+} (Summer11)")); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("no")) {
	sname = "NoMcAcc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #mu^{+}#mu^{-}K^{+} (Spring11)")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	fCsMc = fNMc;
	sname = "CsMc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-}K^{+}K^{-} (MC)")); 
      }	
      if (string::npos != stype.find("1e33") && string::npos != stype.find("cs")) {
	sname = "CsMc1e33"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-}K^{+}K^{-} (MC 1e33)")); 
      }	
      if (string::npos != stype.find("5e32") && string::npos != stype.find("cs")) {
	sname = "CsMc5e32"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-}K^{+}K^{-} (MC 5e32)")); 
      }	
      if (string::npos != stype.find("Summer11") && string::npos != stype.find("cs")) {
	sname = "CsMcSummer11"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-}K^{+}K^{-} (Summer11)")); 
      }	
      if (string::npos != stype.find("acc") && string::npos != stype.find("cs")) {
	sname = "CsMcAcc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0}_{s} #rightarrow #mu^{+}#mu^{-}K^{+}K^{-} (Spring11)")); 
      }



      if (string::npos != stype.find("bg,82")) {
	sname = "bg82"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{+}K^{-}")); 
      }	
      if (string::npos != stype.find("bg,83")) {
	sname = "bg83"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}K^{-}")); 
      }	
      if (string::npos != stype.find("bg,84")) {
	sname = "bg84"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}#pi^{-}")); 
      }	
      if (string::npos != stype.find("bg,93")) {
	sname = "bg93"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}K^{-}")); 
      }	
      if (string::npos != stype.find("bg,92")) {
	sname = "bg92"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}#pi^{-}")); 
      }	
      if (string::npos != stype.find("bg,91")) {
	sname = "bg91"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{+}#pi^{-}")); 
      }	
      if (string::npos != stype.find("bg,60")) {
	sname = "bg60"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p K^{-}")); 
      }	
      if (string::npos != stype.find("bg,61")) {
	sname = "bg61"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p #pi^{-}")); 
      }	
      cout << "open MC " << fNMc << " " << sfile  << " as " << sname << " (" << stype << ") with lumi = " << fMcLumi[fNMc] << endl;
      ++fNMc;
    }
  }

  cout << "fSgMc = " << fSgMc << endl;
  cout << "fBdMc = " << fBdMc << endl;
  cout << "fNoMc = " << fNoMc << endl;
  cout << "fCsMc = " << fCsMc << endl;

  cout << "fSgData = " << fSgData << endl;
  cout << "fNoData = " << fNoData << endl;
  cout << "fCsData = " << fCsData << endl;

}


// ----------------------------------------------------------------------
TFile* anaBmm::loadFile(string file, string type) {
  TFile *f = TFile::Open(file.c_str());
  return f; 
}


// ----------------------------------------------------------------------
void anaBmm::makeAll(int channel) {

  // -- some special plots
  if ((0x1<<0) & channel) {
    allInvertedIso();

    triggerSignal();    
    triggerNorm();    

    plotWithCut("fl3d", "", 2.0, "l_{3d} [cm]", 0., 5);
    plotWithCut("fl3dE", "", -1.0, "#sigma(l_{3d}) [cm]", 0., 0.5);
    
    // -- barrel
    puEff("B_iso4", 0.75, "#epsilon(I>0.75)", "NoData", "Ao");
    puEff("B_iso4", 0.75, "#epsilon(I>0.75)", "CsData", "Ao");
    puEff("B_iso5", 0.75, "#epsilon(I>0.75)", "NoData", "Ao");
    puEff("B_iso5", 0.75, "#epsilon(I>0.75)", "CsData", "Ao");
    puEff("B_iso", 0.75, "#epsilon(I>0.75)", "CsData", "Ao");

    puEff("B_fls3d", 10, "#epsilon(l_{3d}/#sigma(l_{3d})) > 10)", "NoData", "Ao");
    puEff("B_fls3d", 10, "#epsilon(l_{3d}/#sigma(l_{3d})) > 10)", "CsData", "Ao");

    puEff("B_flsxy", 10, "#epsilon(l_{xy}/#sigma(l_{xy})) > 10)", "NoData", "Ao");
    puEff("B_flsxy", 10, "#epsilon(l_{xy}/#sigma(l_{xy})) > 10)", "CsData", "Ao");

    // -- endcap
    puEff("E_iso4", 0.75, "#epsilon(I>0.75)", "NoData", "Ao");
    puEff("E_iso4", 0.75, "#epsilon(I>0.75)", "CsData", "Ao");
    puEff("E_iso5", 0.75, "#epsilon(I>0.75)", "NoData", "Ao");
    puEff("E_iso5", 0.75, "#epsilon(I>0.75)", "CsData", "Ao");
    puEff("E_iso", 0.75, "#epsilon(I>0.75)", "NoData", "Ao");
    puEff("E_iso", 0.75, "#epsilon(I>0.75)", "CsData", "Ao");

    puEff("E_fls3d", 10, "#epsilon(l_{3d}/#sigma(l_{3d})) > 10)", "NoData", "Ao");
    puEff("E_fls3d", 10, "#epsilon(l_{3d}/#sigma(l_{3d})) > 10)", "CsData", "Ao");

    puEff("E_flsxy", 10, "#epsilon(l_{xy}/#sigma(l_{xy})) > 10)", "NoData", "Ao");
    puEff("E_flsxy", 10, "#epsilon(l_{xy}/#sigma(l_{xy})) > 10)", "CsData", "Ao");


    vector<string> todo;
    todo.push_back("iso5");
    todo.push_back("fls3d");
    todo.push_back("chi2dof");
    todo.push_back("alpha");
    
    for (unsigned int i = 0; i < todo.size(); ++i) {
      sbsDistributionOverlaySameFile("NoData", Form("B_%s", todo[i].c_str()), Form("E_%s", todo[i].c_str()), "Ao", "Barrel", "Endcap");
      sbsDistributionOverlaySameFile("NoMc", Form("B_%s", todo[i].c_str()), Form("E_%s", todo[i].c_str()), "Ao", "Barrel", "Endcap");
      
      sbsDistributionOverlaySameFile("CsData", Form("B_%s", todo[i].c_str()), Form("E_%s", todo[i].c_str()), "Ao", "Barrel", "Endcap");
      sbsDistributionOverlaySameFile("CsMc", Form("B_%s", todo[i].c_str()), Form("E_%s", todo[i].c_str()), "Ao", "Barrel", "Endcap");
    }
  }

  
  sbsDistributionOverlaySameFile("NoData", "APV0_iso", "APV1_iso", "Ao", "low pileup", "high pileup");
  sbsDistributionOverlaySameFile("NoMcSummer11", "APV0_iso", "APV1_iso", "Ao", "low pileup", "high pileup");

  sbsDistributionOverlaySameFile("CsData", "APV0_iso", "APV1_iso", "Ao", "low pileup", "high pileup");
  sbsDistributionOverlaySameFile("CsMcSummer11", "APV0_iso", "APV1_iso", "Ao", "low pileup", "high pileup");



  // -- 2
  if ((0x1<<1) & channel) {
    vector<string> regions;
    regions.push_back("A");
    regions.push_back("B");
    regions.push_back("E");
    for (unsigned int i = 0; i < regions.size(); ++i) {
      sbsDistributionOverlay("SgData", "SgMc", "Ao", regions[i].c_str());
      sbsDistributionOverlay("NoData", "NoMcSummer11", "Ao", regions[i].c_str());
      sbsDistributionOverlay("CsData", "CsMcSummer11", "Ao", regions[i].c_str());
      sbsDistributionOverlay("NoData", "NoMc", "Ao", regions[i].c_str());
      sbsDistributionOverlay("CsData", "CsMc", "Ao", regions[i].c_str());
      //       sbsDistributionOverlay("NoData", "CsData", "Ao", regions[i].c_str());
    }

    //    sbsDistributionOverlay("NoData", "NoMc1e33", "Ao", regions[i].c_str());
  }

  // -- 4
  if ((0x1<<2) & channel) {
    allEffTables();
  }

  // -- 8
  if ((0x1<<3) & channel) {
    rareBg();
  }

  // -- 16 
  if ((0x1<<4) & channel) {
    computeNormUL();
    computeCsBF();
  }

}


// ----------------------------------------------------------------------
void anaBmm::dumpSamples() {

  std::map<string, double> ngen;
  ngen.insert(make_pair("bg60", 20.e6)); 
  ngen.insert(make_pair("bg61", 20.e6)); 

  ngen.insert(make_pair("bg84", 20.e6)); 
  ngen.insert(make_pair("bg83", 20.e6)); 
  ngen.insert(make_pair("bg82", 20.e6)); 

  ngen.insert(make_pair("bg91", 20.e6)); 
  ngen.insert(make_pair("bg92", 20.e6)); 
  ngen.insert(make_pair("bg93", 20.e6)); 

  ngen.insert(make_pair("SgMcAcc", 100.e6)); 

  ngen.insert(make_pair("BdMcAcc", 40.e6)); 

  ngen.insert(make_pair("NoMcAcc", 160.e6)); 
  ngen.insert(make_pair("CsMcAcc", 40.e6)); 
  
  //  ofstream OUT(fNumbersFileName.c_str(), ios::app);

  fTEX << "% ----------------------------------------------------------------------" << endl;
  string name; 
  double lumi, n, f; 
  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    cout << "===> " << imap->first;
    cout << ": " << fF[imap->first]->GetName();
    cout << " -> " << imap->second << endl;
    name = imap->second; 
    lumi = fLumi[imap->first];
    n = ngen[imap->first];
    f = ((TH1D*)fF[imap->first]->Get("monEvents"))->GetBinContent(1);
    replaceAll(name, "#", "\\"); 
    //    cout <<  Form("\\vdef{%s:sampleName:%s}   {\\ensuremath{{%s } } }", fSuffix.c_str(), imap->first.c_str(), name.c_str()) << endl;
    fTEX <<  Form("\\vdef{%s:sampleName:%s}   {\\ensuremath{{%s } } }", fSuffix.c_str(), imap->first.c_str(), name.c_str()) << endl;
    //    cout <<  Form("\\vdef{%s:lumi:%s}   {\\ensuremath{{%4.1f } } }", fSuffix.c_str(), imap->first.c_str(), lumi) << endl;
    fTEX <<  Form("\\vdef{%s:lumi:%s}   {\\ensuremath{{%4.1f } } }", fSuffix.c_str(), imap->first.c_str(), lumi) << endl;
    if (n>0) {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f } } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
      fTEX <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f} } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
    } else {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str() ) << endl;
      fTEX <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
    if (f>0) {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f } } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
      fTEX <<  Form("\\vdef{%s:nfilt:%s}   {\\ensuremath{{%4.0f} } }", fSuffix.c_str(), imap->first.c_str(), f) << endl;
    } else {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str() ) << endl;
      fTEX <<  Form("\\vdef{%s:nfilt:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
  }
  fTEX.flush();

}



// ----------------------------------------------------------------------
void anaBmm::dumpCutNames() {

  //  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  
  TH1D *h = (TH1D*)gFile->Get("hcuts"); 
  string cut, empty(""), cutstring, cutline, cutvalue; 
  double value; 
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i));
    value = h->GetBinContent(i);
    if (cut == empty) continue;
    string::size_type n1 = cut.find_first_of(":"); 
    string::size_type n2 = cut.find_last_of(":")-1; 
    
    cutstring = cut.substr(0, n1-1); 
    cutline   = cut.substr(n1+3, n2-n1-4);
    replaceAll(cutline, "#", "\\"); 
    cutvalue  = cut.substr(n2+2, cut.size()-n2);

    //    cout << cut << " cutstring==>" << cutstring << "<== cutvalue==>" << cutvalue << " cutline==>" << cutline << "<==" << endl;
    if (string::npos != cutstring.find("JSON")) {
      replaceAll(cutvalue, "_", "\\_"); 
      replaceAll(cutvalue, "cuts/", ""); 
      fTEX <<  Form("\\vdef{%s:%s:cutLine}   {{%s } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
      fTEX <<  Form("\\vdef{%s:%s:cutValue}  {{%s } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
    } else if (string::npos != cutstring.find("CANDCOSALPHA")) {
      fTEX <<  Form("\\vdef{%s:CANDALPHA:cutLine}   {\\ensuremath{{\\alpha } } }", fSuffix.c_str()) << endl;
      fTEX <<  Form("\\vdef{%s:CANDALPHA:cutValue}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), TMath::ACos(atof(cutvalue.c_str()))) 
	   << endl;
    } else {
      fTEX <<  Form("\\vdef{%s:%s:cutLine}   {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
      fTEX <<  Form("\\vdef{%s:%s:cutValue}  {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
    }

  }
  fTEX.flush();
}


// ----------------------------------------------------------------------
void anaBmm::allEffTables() {
  effTable("SgMc", "A");
  effTable("SgData", "A");
  effTable("NoData", "A");
  effTable("NoMc", "A");
  effTable("CsData", "A");
  effTable("CsMc", "A");

  effTable("BdMcSummer11", "APV0");
  effTable("BdMcSummer11", "APV1");

  effTable("NoData", "APV0");
  effTable("NoData", "APV1");

  effTable("NoMcSummer11", "APV0");
  effTable("NoMcSummer11", "APV1");

  effTable("CsData", "APV0");
  effTable("CsData", "APV1");

  effTable("CsMcSummer11", "APV0");
  effTable("CsMcSummer11", "APV1");

  fTEX.flush();
  // -- Compute difference
  double deltaIso     = computeDelta(Form("{%s:A_iso5eff:NoData}", fSuffix.c_str()), Form("{%s:A_iso5eff:NoMc}", fSuffix.c_str()));
  double deltaAlpha   = computeDelta(Form("{%s:A_alphaeff:NoData}", fSuffix.c_str()), Form("{%s:A_alphaeff:NoMc}", fSuffix.c_str()));
  double deltaFls     = computeDelta(Form("{%s:A_fls3deff:NoData}", fSuffix.c_str()), Form("{%s:A_fls3deff:NoMc}", fSuffix.c_str()));
  double deltaChi2dof = computeDelta(Form("{%s:A_chi2dofeff:NoData}", fSuffix.c_str()), Form("{%s:A_chi2dofeff:NoMc}", fSuffix.c_str()));
  double deltaTot     = TMath::Sqrt(deltaIso*deltaIso + deltaAlpha*deltaAlpha + deltaFls*deltaFls + deltaChi2dof*deltaChi2dof); 
  double deltaTotNoIso= TMath::Sqrt(deltaAlpha*deltaAlpha + deltaFls*deltaFls + deltaChi2dof*deltaChi2dof); 

  fTEX << Form("%s", (formatTex(deltaIso, fSuffix+":A_deltaIso", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaAlpha, fSuffix+":A_deltaAlpha", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaFls, fSuffix+":A_deltaFls", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaChi2dof, fSuffix+":A_deltaChi2dof", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(100*deltaTot, fSuffix+":A_deltaTot", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(100*deltaTotNoIso, fSuffix+":A_deltaTotNoIso", "No")).c_str()) << endl;


  deltaIso     = computeDelta(Form("{%s:A_iso5eff:CsData}", fSuffix.c_str()), Form("{%s:A_iso5eff:CsMc}", fSuffix.c_str()));
  deltaAlpha   = computeDelta(Form("{%s:A_alphaeff:CsData}", fSuffix.c_str()), Form("{%s:A_alphaeff:CsMc}", fSuffix.c_str()));
  deltaFls     = computeDelta(Form("{%s:A_fls3deff:CsData}", fSuffix.c_str()), Form("{%s:A_fls3deff:CsMc}", fSuffix.c_str()));
  deltaChi2dof = computeDelta(Form("{%s:A_chi2dofeff:CsData}", fSuffix.c_str()), Form("{%s:A_chi2dofeff:CsMc}", fSuffix.c_str()));
  deltaTot     = TMath::Sqrt(deltaIso*deltaIso + deltaAlpha*deltaAlpha + deltaFls*deltaFls + deltaChi2dof*deltaChi2dof); 
  deltaTotNoIso= TMath::Sqrt(deltaAlpha*deltaAlpha + deltaFls*deltaFls + deltaChi2dof*deltaChi2dof); 
  
  fTEX << Form("%s", (formatTex(deltaIso, fSuffix+":A_deltaIso", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaAlpha, fSuffix+":A_deltaAlpha", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaFls, fSuffix+":A_deltaFls", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaChi2dof, fSuffix+":A_deltaChi2dof", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(100*deltaTot, fSuffix+":A_deltaTot", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(100*deltaTotNoIso, fSuffix+":A_deltaTotNoIso", "Cs")).c_str()) << endl;
  
  fTEX.flush();
}


// ----------------------------------------------------------------------
double anaBmm::computeDelta(const char *s1, const char *s2, int relative) {

  char  buffer[200];
  string line;
  ifstream is(fNumbersFileName.c_str());

  float f1(-1.), f2(-1.);
  while (is.getline(buffer, 200, '\n')) {
    line = buffer; 
    if (string::npos != line.find(s1)) {
      string::size_type m1 = line.find("ensuremath{{"); 
      string number = line.substr(m1+12, 6); 
      f1 = atof(number.c_str());
      cout << line << " --> " << number << " => " << f1 << endl;
    }
    if (string::npos != line.find(s2)) {
      string::size_type m1 = line.find("ensuremath{{"); 
      string number = line.substr(m1+12, 6); 
      f2 = atof(number.c_str());
      cout << line << " --> " << number << " => " << f2 << endl;
    }
  }

  double delta = TMath::Abs(f1 - f2);
  is.close();

  if (relative) {
    delta = 2.*delta/(f1+f2);
  }
		
  return delta;
  
}


// ----------------------------------------------------------------------
void anaBmm::testEff(const char *additionalCuts, const char *basicCuts) {
  
  //  TH1D *h1 = new TH1D("h1", "h1", 40, 4.8, 6.0); 

  TTree *td = (TTree*)(fpMc[fSgMc]->Get("events"));
  double n2 = td->Draw("m", Form("%s", basicCuts));
  double n1 = td->Draw("m", Form("%s&&%s", basicCuts, additionalCuts));
  
  cout << additionalCuts << endl;
  cout << n1 << "/" << n2 << " = " << n1/n2 << endl;

}


// ----------------------------------------------------------------------
void anaBmm::effTable(string smode, const char *region) {
  gStyle->SetOptTitle(1); 
  gStyle->SetOptFit(1112); 
  tl->SetTextSize(0.05); 
  fF[smode]->cd();
  cout << "effTable for file: " << gFile->GetName() << " and region: " << region << endl;
  TH1D *h = (TH1D*)gFile->Get("analysisDistributions"); 
  if (0 == h) {
    cout << "no histogram analysisDistributions found, returning" << endl;
    return;
  }

  double massLo(5.1), massHi(5.6), massPeak(-1.), massSigma(-1.);
  if (string::npos != smode.find("No")) {
    massLo = 5.0; 
    massHi = 5.5;
    massPeak = 5.28;
    massSigma = 0.035;
  }

  if (string::npos != smode.find("Cs")) {
    massLo = 5.15; 
    massHi = 5.6;
    massPeak = 5.37;
    massSigma = 0.035;
  }
  
  //  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << " % --- " << smode << endl;
  int mode(0); 
  if (string::npos != smode.find("SgMc")) {
    cout << "==> This is signal MC, counting events in histogram" << endl;
    mode = 0; 
  }
  if (string::npos != smode.find("BdMc")) {
    cout << "==> This is Bd2MuMu signal MC, counting events in histogram" << endl;
    mode = 0; 
  }
  if (string::npos != smode.find("SgData")) {
    cout << "==> This is signal Data, counting events in sidebands" << endl;
    mode = 0; 
  }
  
  if (string::npos != smode.find("No") || string::npos != smode.find("Cs")) {
    if (0 == fMode) {
      cout << "==> This is normalization/control sample, fitting mode = 11 = expo+Gaus" << endl;
      mode = 11; 
    } else {
      cout << "==> This is normalization/control sample, fitting mode = " << fMode << endl;
      mode = fMode; 
    }
  }

  string cut, pdfname;
  double n(0.), nE(0.); 
  double norm(0.),  normE(0.), eff(0.), effE(0.);

  // -- normalization
  AnalysisDistribution *an = new AnalysisDistribution(Form("%s_docatrk", region));
  an->fMassLo    = massLo; 
  an->fMassHi    = massHi; 
  an->fMassPeak  = massPeak; 
  an->fMassSigma = massSigma; 
  an->hMassCu->SetMinimum(0.);
  norm = an->fitMass(an->hMassCu, normE, mode); 
  fTEX << Form("%s", (formatTex(norm, fSuffix+":"+cut+"Norm", smode)).c_str()) << endl;
  fTEX << Form("%s", (formatTex(normE, fSuffix+":"+cut+"NormE", smode)).c_str()) << endl;
  
  tl->DrawLatex(0.22, 0.75, Form("%4.2f+/-%4.2f", norm, normE)); 
  pdfname = Form("%s/%s_%s_hMassNorm.pdf", fDirectory.c_str(), fSuffix.c_str(), smode.c_str());
  cout << "AD for " << cut << " results in " << pdfname << endl;
  if (fDoPrint) c0->SaveAs(pdfname.c_str(), "Portrait");

  delete an;
  
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i)); 

    if (string::npos == cut.find(Form("%s_", region))) {
      //      cout << "region " << region << " not part of name, skipping" << endl;
      continue;
    }

    if (cut == string("")) {
      cout << "empty string found, break ..." << endl;
      //OUT.close();
      break;
    }
    pdfname = Form("%s/%s_%s.pdf", fDirectory.c_str(), smode.c_str(), cut.c_str());
    AnalysisDistribution *a = new AnalysisDistribution(cut.c_str());
    a->fMassLo    = massLo; 
    a->fMassHi    = massHi; 
    a->fMassPeak  = massPeak; 
    a->fMassSigma = massSigma; 
    a->hMassAo->SetMinimum(0.);

    n    = a->fitMass(a->hMassAo, nE, mode); 
    if (nE > n) nE = 1.2*TMath::Sqrt(n); 
    tl->DrawLatex(0.22, 0.75, Form("%4.2f+/-%4.2f", n, nE)); 
    pdfname = Form("%s/%s_%s_%s_hMassAo.pdf", fDirectory.c_str(), fSuffix.c_str(), smode.c_str(), cut.c_str());
    cout << "AD for " << cut << " results in " << pdfname << endl;
    // c0->SaveAs(pdfname.c_str(), "Portrait");

    delete a; 

    if ((string::npos != cut.find("tracks")) || (string::npos != cut.find("muons"))) {
      n *=0.5; 
      nE *=0.5; 
    }

    if (nE > n) nE = 0.02*n;
    if (nE < normE) nE = normE;
    eff    = norm/n;
    if (eff > 1.0) eff = 1.0; 
    effE   = dEff(norm, normE, n, nE); 

    cout << cut << " norm = " << norm << "  n = " << n << "+/-" << nE 
	 << " eff = " << eff << "+/-" << effE
	 << endl;

    fTEX << Form("%s", (formatTex(n, fSuffix+":"+cut+"N", smode)).c_str()) << endl;
    fTEX << Form("%s", (formatTex(nE, fSuffix+":"+cut+"NE", smode)).c_str()) << endl;
    fTEX << Form("%s", (formatTex(eff, fSuffix+":"+cut+"eff", smode)).c_str()) << endl;
    fTEX << Form("%s", (formatTex(effE, fSuffix+":"+cut+"effE", smode)).c_str()) << endl;
    //     OUT << Form("%s", (formatTex(relEff, fSuffix+":"+cut+"eRel", smode)).c_str()) << endl;
    //     OUT << Form("%s", (formatTex(relEffE, fSuffix+":"+cut+"eRelE", smode)).c_str()) << endl;
    //     OUT << Form("%s", (formatTex(cumEff, fSuffix+":"+cut+"eCum", smode)).c_str()) << endl;
    //     OUT << Form("%s", (formatTex(cumEffE, fSuffix+":"+cut+"eCumE", smode)).c_str()) << endl;

    //     nprev = n; 
    //     nprevE= nE; 
  }
  fTEX.flush();
}


// ----------------------------------------------------------------------
void anaBmm::isoMean(const char *var, const char *cuts, int maxPvN) {

  cout << "normalize to " << fDataLumi[fSgData] << endl;

  gROOT->cd();

  double eps(0.1); 
  TH1D *h1 = (TH1D*)gFile->Get("isoMeanSgMc");  
  if (0 == h1) {h1 = new TH1D("isoMeanSgMc", "", maxPvN+1, 0., maxPvN+1); h1->Sumw2();} else  h1->Reset();

  TH1D *h2 = (TH1D*)gFile->Get("isoMeanSgDa");  
  if (0 == h2) {h2 = new TH1D("isoMeanSgDa", "", maxPvN+1, 0.+eps, maxPvN+1+eps); h2->Sumw2();} else  h2->Reset();

  TH1D *h3 = (TH1D*)gFile->Get("isoMeanNoMc");  
  if (0 == h3) {h3 = new TH1D("isoMeanNoMc", "", maxPvN+1, 0.+eps*2, maxPvN+1+eps*2); h3->Sumw2();} else  h3->Reset();

  TH1D *h4 = (TH1D*)gFile->Get("isoMeanNoDa");  
  if (0 == h4) {h4 = new TH1D("isoMeanNoDa", "", maxPvN+1, 0.+eps*3, maxPvN+1+eps*3); h4->Sumw2();} else  h4->Reset();
  

  TTree *t; 

  // -- signal MC
  fpMc[fSgMc]->cd(); 
  TH1D *h = (TH1D*)gFile->Get("isoMeanDist");  
  if (0 == h) {
    h = new TH1D("isoMeanDist", "", 101, 0., 1.01);
  } else{
    h->Reset();
  }

  t = (TTree*)gFile->Get("events"); 
  for (int i = 1; i <= maxPvN; ++i) {
    h->Reset();
    cout << gFile->GetName() << " npv = " << i << " cuts: " << Form("%s&&pvn==%d", cuts, i) << endl;
    t->Draw(Form("%s>>isoMeanDist", var), Form("%s&&pvn==%d", cuts, i));
    h->Draw();
    cout << "mean: "  << h->GetMean() << endl;
    c0->Modified();
    c0->Update();
    if (h->GetSumOfWeights() > 20) {
      h1->SetBinContent(i, h->GetMean()); 
      h1->SetBinError(i, h->GetRMS()); 
    } else {
      //      break;
    }
  } 


  // -- signal Data
  fpData[fSgData]->cd(); 
  h = (TH1D*)gFile->Get("isoMeanDist");  
  if (0 == h) {
    h = new TH1D("isoMeanDist", "", 101, 0., 1.01);
  } else{
    h->Reset();
  }

  t = (TTree*)gFile->Get("events"); 
  for (int i = 1; i <= maxPvN; ++i) {
    h->Reset();
    cout << gFile->GetName() << " npv = " << i << " cuts: " << Form("%s&&pvn==%d", cuts, i) << endl;
    t->Draw(Form("%s>>isoMeanDist", var), Form("%s&&pvn==%d", cuts, i));
    h->Draw();
    cout << "mean: "  << h->GetMean() << endl;
    c0->Modified();
    c0->Update();
    if (h->GetSumOfWeights() > 20) {
      h2->SetBinContent(i, h->GetMean()); 
      h2->SetBinError(i, h->GetRMS()); 
    } else {
      //      break;
    }
  }


//   // -- norm MC
//   fpMc[fNoMc]->cd(); 
//   h = (TH1D*)gFile->Get("isoMeanDist");  
//   if (0 == h) {
//     h = new TH1D("isoMeanDist", "", 101, 0., 1.01);
//   } else{
//     h->Reset();
//   }

//   t = (TTree*)gFile->Get("events"); 
//   for (int i = 1; i <= maxPvN; ++i) {
//     h->Reset();
//     cout << gFile->GetName() << " npv = " << i << " cuts: " << Form("%s&&pvn==%d", cuts, i) << endl;
//     t->Draw(Form("%s>>isoMeanDist", var), Form("%s&&pvn==%d", cuts, i));
//     h->Draw();
//     cout << "mean: "  << h->GetMean() << endl;
//     c0->Modified();
//     c0->Update();
//     if (h->GetSumOfWeights() > 20) {
//       h3->SetBinContent(i, h->GetMean()); 
//       h3->SetBinError(i, h->GetRMS()); 
//     } else {
//       //      break;
//     }
//   } 


//   // -- norm Data
//   fpData[fNoData]->cd(); 
//   h = (TH1D*)gFile->Get("isoMeanDist");  
//   if (0 == h) {
//     h = new TH1D("isoMeanDist", "", 101, 0., 1.01);
//   } else{
//     h->Reset();
//   }

//   t = (TTree*)gFile->Get("events"); 
//   for (int i = 1; i <= maxPvN; ++i) {
//     h->Reset();
//     cout << gFile->GetName() << " npv = " << i << " cuts: " << Form("%s&&m>5.2&&m<5.35&&pvn==%d", cuts, i) << endl;
//     t->Draw(Form("%s>>isoMeanDist", var), Form("%s&&m>5.2&&m<5.35&&pvn==%d", cuts, i));
//     h->Draw();
//     cout << "mean: "  << h->GetMean() << endl;
//     c0->Modified();
//     c0->Update();

//     if (h->GetSumOfWeights() > 20) {
//       h4->SetBinContent(i, h->GetMean()); 
//       h4->SetBinError(i, h->GetRMS()); 
//     } else {
//       //      break;
//     } 
//   }

  setTitles(h1, "N_{PV}", "<I>");
  setHist(h1, kBlue, 24); h1->Draw();
  setHist(h2, kBlue, 20); h2->Draw("same");

  newLegend(0.65, 0.2, 0.80, 0.35); 
  legg->AddEntry(h1, "Signal MC", "p"); 
  legg->AddEntry(h2, "Signal Data", "p"); 
  legg->Draw();

  if (fDoPrint) c0->SaveAs(Form("%s-npv.pdf", var)); 

//   setHist(h3, kRed, 25); h3->Draw("same");
//   setHist(h4, kRed, 21); h4->Draw("same");
}


// ----------------------------------------------------------------------
void anaBmm::varProcess(const char *var, const char *cuts, double lo, double hi, const char *at, int loop) {

  // -- signal MC
  fpMc[fSgMc]->cd(); 
  TH1D *h40 = (TH1D*)gFile->Get("iso40");  
  if (h40)  delete h40;   
  h40 = new TH1D("iso40", "", 40, lo, hi);
  TH1D *h41 = (TH1D*)gFile->Get("iso41");  
  if (h41)  delete h41;   
  h41 = new TH1D("iso41", "", 40, lo, hi);
  TH1D *h42 = (TH1D*)gFile->Get("iso42");  
  if (h42)  delete h42;   
  h42 = new TH1D("iso42", "", 40, lo, hi);

  double ncnt;
  TTree *t; 
  t = (TTree*)gFile->Get("events"); 
  t = (TTree*)gFile->Get("effTree"); 
  cout << "*********************" << endl;
  cout << gFile->GetName() << " -> t = " << t << endl;
  ncnt = t->Draw(Form("%s>>iso40", var), Form("%s&&procid==%d", cuts, 40));
  cout << "40: " << ncnt << endl;
  ncnt = t->Draw(Form("%s>>iso41", var), Form("%s&&procid==%d", cuts, 41));
  cout << "41: " << ncnt << endl;
  ncnt = t->Draw(Form("%s>>iso42", var), Form("%s&&procid==%d", cuts, 42));
  cout << "42: " << ncnt << endl;
  
  h41->Scale(h40->GetSumOfWeights()/h41->GetSumOfWeights()); 
  h42->Scale(h40->GetSumOfWeights()/h42->GetSumOfWeights()); 

  setTitles(h40, at, "a.u.", 0.06, 0.9, 1.3);
  setHist(h40, kBlue);  
  setHist(h41, kRed);     //h41->SetLineStyle(kDashed); 
  setHist(h42, kBlack);   //h42->SetLineStyle(kDotted); 

  shrinkPad(0.12, 0.15); 
  double maxi(0.); 
  maxi = h40->GetMaximum(); 
  if (h41->GetMaximum() > maxi) maxi = h41->GetMaximum(); 
  if (h42->GetMaximum() > maxi) maxi = h42->GetMaximum(); 
  h40->SetMaximum(1.1*maxi); 
  h40->Draw();
  h41->Draw("same");
  h42->Draw("same");


  newLegend(0.55, 0.7, 0.8, 0.85); 
  legg->AddEntry(h40, Form("<GGF> = %4.3f", h40->GetMean()), "l"); 
  legg->AddEntry(h41, Form("<FEX> = %4.3f", h41->GetMean()), "l"); 
  legg->AddEntry(h42, Form("<GSP> = %4.3f", h42->GetMean()), "l"); 
  legg->Draw();


  string pdfname = Form("%s/%s_process-%s.pdf", fDirectory.c_str(), fSuffix.c_str(), var);
  if (fDoPrint){
    if (c0) c0->SaveAs(pdfname.c_str());
  }

  if (0 == loop) return;


  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << "% --- procid numbers" << endl;  
  for (int i = 40; i <= 42; ++i) {

    loopTree(0, i); 
    loopTree(10, i); 

    for (int ichan = 0; ichan < 2; ++ichan) {
      fTEX <<  Form("\\vdef{%s:procid%iChan%i:effTotSg}   {\\ensuremath{{%6.5f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersBs[ichan]->effTot)
	   << endl;
      fTEX <<  Form("\\vdef{%s:procid%iChan%i:effTotNo}   {\\ensuremath{{%6.5f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersNorm[ichan]->effTot)
	   << endl;
      fTEX <<  Form("\\vdef{%s:procid%iChan%i:effTotRatio}   {\\ensuremath{{%4.3f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersBs[ichan]->effTot/fNumbersNorm[ichan]->effTot)
	   << endl;

      fTEX <<  Form("\\vdef{%s:procid%iChan%i:effAnaSg}   {\\ensuremath{{%6.5f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersBs[ichan]->effAna)
	   << endl;
      fTEX <<  Form("\\vdef{%s:procid%iChan%i:effAnaNo}   {\\ensuremath{{%6.5f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersNorm[ichan]->effAna)
	   << endl;
      fTEX <<  Form("\\vdef{%s:procid%iChan%i:effAnaRatio}   {\\ensuremath{{%4.3f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersBs[ichan]->effAna/fNumbersNorm[ichan]->effAna)
	   << endl;

      fTEX <<  Form("\\vdef{%s:procid%iChan%i:accSg}   {\\ensuremath{{%6.5f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersBs[ichan]->acc)
	   << endl;
      fTEX <<  Form("\\vdef{%s:procid%iChan%i:accNo}   {\\ensuremath{{%6.5f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersNorm[ichan]->acc)
	   << endl;
      fTEX <<  Form("\\vdef{%s:procid%iChan%i:accRatio}   {\\ensuremath{{%4.3f } } }", 
		    fSuffix.c_str(), i, ichan, fNumbersBs[ichan]->acc/fNumbersNorm[ichan]->acc)
	   << endl;
    }
  }
  
}


// ----------------------------------------------------------------------
void anaBmm::procAcc(int mode, int chan) {

  string cuts; 

  if (0 == mode) {
    if (0 == chan) {
      cuts = "g1pt>1&&g2pt>1&&abs(g1eta)<1.4&&abs(g2eta)<1.4&&m1pt>1&&m2pt>1&&abs(m1eta)<1.4&&abs(m2eta)<1.4&&m1gt&&m2gt";
    } else {
      cuts = "g1pt>1&&g2pt>1&&(abs(g1eta)>1.4||abs(g2eta)>1.4)&&m1pt>1&&m2pt>1&&(abs(m1eta)>1.4||abs(m2eta)>1.4)&&m1gt&&m2gt";
    }
  } else if (10 == mode) {
    if (0 == chan) {
      cuts = "g1pt>1&&g2pt>1&&g3pt>0.4&&abs(g1eta)<1.4&&abs(g2eta)<1.4&&abs(g3eta)<2.5&&m1pt>1&&m2pt>1&&k1pt>0.5&&abs(m1eta)<1.4&&abs(m2eta)<1.4&&abs(k1eta)<2.4&&m1gt&&m2gt&&k1gt";
    } else {
      cuts = "g1pt>1&&g2pt>1&&g3pt>0.4&&(abs(g1eta)>1.4||abs(g2eta)>1.4)&&abs(g1eta)<2.5&&abs(g2eta)<2.5&&abs(g3eta)<2.5&&m1pt>1&&m2pt>1&&k1pt>0.5&&(abs(m1eta)>1.4||abs(m2eta)>1.4)&&abs(m1eta)<2.4&&abs(m2eta)<2.4&&abs(k1eta)<2.4&&m1gt&&m2gt&&k1gt";
    }
  }

  TTree *t; 
  t = (TTree*)gFile->Get("effTree"); 
  double c40 = t->Draw("g1pt", Form("%s&&procid==40", cuts.c_str()));
  double n40 = t->Draw("g1pt", "procid==40");
  
  double c41 = t->Draw("g1pt", Form("%s&&procid==41", cuts.c_str()));
  double n41 = t->Draw("g1pt", "procid==41");
  
  double c42 = t->Draw("g1pt", Form("%s&&procid==42", cuts.c_str()));
  double n42 = t->Draw("g1pt", "procid==42");
  
  cout << Form("GGF: %4.3f +/- %4.3f", c40/n40, dEff(static_cast<int>(c40), static_cast<int>(n40))) << endl;
  cout << Form("FEX: %4.3f +/- %4.3f", c41/n41, dEff(static_cast<int>(c41), static_cast<int>(n41))) << endl;
  cout << Form("GSP: %4.3f +/- %4.3f", c42/n42, dEff(static_cast<int>(c42), static_cast<int>(n42))) << endl;

}

// ----------------------------------------------------------------------
void anaBmm::allInvertedIso() {
  histInvertedIso("fls3d>", 10, 10., 20.); 
  histInvertedIso("chi2/dof<", 10, 1., 2.); 
  histInvertedIso("alpha<", 10, 0.01, 0.06); 
  histInvertedIso("m1pt>", 10, 3., 5.0); 
  histInvertedIso("m2pt>", 10, 3., 5.0); 
  histInvertedIso("pt>", 10, 4.0, 8.0); 
  histInvertedIso("docatrk>", 10, 0.00, 0.02); 
}


// ----------------------------------------------------------------------  
void anaBmm::histInvertedIso(const char *var, int n, double lo, double hi) {

  gStyle->SetOptStat(0);
  
  vector<string> cuts;
  vector<double> cutv;
  vector<string> cutt;
  vector<string> cutf;
  string filename; 

  TH1D *BE = new TH1D("BE", "", n, lo, hi); BE->Sumw2(); 
  TH1D *BO = new TH1D("BO", "", n, lo, hi); BO->Sumw2(); 
  TH1D *EE = new TH1D("EE", "", n, lo, hi); EE->Sumw2(); 
  TH1D *EO = new TH1D("EO", "", n, lo, hi); EO->Sumw2(); 

  for (int ichan = 0; ichan < 2; ++ichan) {
    cuts.clear(); 
    cutv.clear(); 
    
    cuts.push_back("fls3d>");    cutv.push_back(fCuts[ichan]->fls3d);   
    cutt.push_back("l_{3d}/#sigma(l_{3d}) >");
    cutf.push_back("fls3d");

    cuts.push_back("chi2/dof<"); cutv.push_back(fCuts[ichan]->chi2dof);
    cutt.push_back("#chi^{2}/dof < ");
    cutf.push_back("chi2dof");

    cuts.push_back("alpha<");    cutv.push_back(fCuts[ichan]->alpha);
    cutt.push_back("#alpha < ");
    cutf.push_back("alpha");

    cuts.push_back("m1pt>");     cutv.push_back(fCuts[ichan]->m1pt);
    cutt.push_back("p_{T,#mu1} > ");
    cutf.push_back("m1pt");

    cuts.push_back("m2pt>");     cutv.push_back(fCuts[ichan]->m2pt);
    cutt.push_back("p_{T,#mu2} > ");
    cutf.push_back("m2pt");

    cuts.push_back("docatrk>");  cutv.push_back(fCuts[ichan]->docatrk);
    cutt.push_back("d^{0}_{trk} > ");
    cutf.push_back("docatrk");

    cuts.push_back("pt>");  cutv.push_back(fCuts[ichan]->pt);
    cutt.push_back("p_{T,B} > ");
    cutf.push_back("pt");
    
    string ocuts, tcuts, rcuts;
    for (unsigned int i = 0; i < cuts.size(); ++i) {
      if (!strcmp(cuts[i].c_str(), var)) {
	setTitles(BE, cutt[i].c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 
	setTitles(EE, cutt[i].c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 
	filename = cutf[i];
	continue;
      }
      tcuts += Form("%s %f && ", cuts[i].c_str(), cutv[i]);
    }
    
    string::size_type n1 = tcuts.find_last_of("&&"); 
    ocuts = tcuts.substr(0, n1-1); 
    
    double steps = (hi-lo)/n;
    double cut; 
    for (int i = 0; i < n; ++i) {
      cut = lo + i*steps; 
      rcuts = ocuts + Form(" && %s %f", var, cut); 
      //      cout << "-->" << rcuts << endl;
      
      invertedIso(ichan, rcuts.c_str()); 
      if (0 == ichan) {
	BE->SetBinContent(i+1, fBlExp); BE->SetBinError(i+1, fBlExpE); 
	BO->SetBinContent(i+1, fBlObs); BO->SetBinError(i+1, fBlObsE); 
      } else {
	EE->SetBinContent(i+1, fBlExp); EE->SetBinError(i+1, fBlExpE); 
	EO->SetBinContent(i+1, fBlObs); EO->SetBinError(i+1, fBlObsE); 
      }	
    }
  }
  
  c0->Clear();
  c0->Divide(1,2);
  
  c0->cd(1);
  shrinkPad(0.2, 0.15); 
  double maxi = BE->GetMaximum(); 
  if (BO->GetMaximum() > maxi) {
    maxi = 1.2*BO->GetMaximum();
  } else {
    maxi = 1.2*BE->GetMaximum();
  }
  setFilledHist(BE);
  BE->SetMaximum(maxi); 
  BE->SetMinimum(0.); 
  BE->Draw("hist");
  BO->Draw("esame");
  tl->DrawLatex(0.15, 0.92, "Barrel");

  maxi = EE->GetMaximum(); 
  if (EO->GetMaximum() > maxi) {
    maxi = 1.2*EO->GetMaximum();  
  } else {
    maxi = 1.2*EE->GetMaximum();  
  }
  c0->cd(2);
  shrinkPad(0.2, 0.15); 
  setFilledHist(EE);
  EE->SetMaximum(maxi); 
  EE->SetMinimum(0.); 
  EE->Draw("hist");
  EO->Draw("esame");
  tl->DrawLatex(0.15, 0.92, "Endcap");

  c0->SaveAs(Form("invertedIso-%s.pdf", filename.c_str())); 

}



// ----------------------------------------------------------------------
void anaBmm::invertedIso(int chan, const char *cuts) {
  string baseCuts = "gmuid&&hlt&&gmupt&&gmueta&&iso5<0.7";
  string chanCut; 
  if (0 == chan) {
    chanCut = "(abs(m1eta)<1.4&&abs(m2eta)<1.4)";
  } else {
    chanCut = "(abs(m1eta)>1.4||abs(m2eta)>1.4)";
  }
  
  string allCuts  = baseCuts + "&&" + chanCut + "&&" + cuts; 
  cout << allCuts << endl;
  
  fF["SgData"]->cd();
  TH1D *h1(0); 
  if (h1) delete h1; 
  h1 = new TH1D("h1", "", 1200, 4.8, 6.0); 
  TTree *t = (TTree*)gFile->Get("events"); 
  t->Draw("m>>h1", allCuts.c_str(), "goff");
  h1->Draw();

  int blo = h1->FindBin(5.2); 
  int bhi = h1->FindBin(5.45); 
  //  cout << "blo: " << blo << " bhi: " << bhi << endl;

  double lo = h1->Integral(1, blo);
  double bl = h1->Integral(blo, bhi);
  double blE= TMath::Sqrt(bl); 
  double hi = h1->Integral(bhi, h1->GetNbinsX());
  double ex = (lo+hi)*0.25/(1.2-0.25);
  double exE= TMath::Sqrt(ex); 
  double df = bl - ex; 
  double dfE=TMath::Sqrt(blE*blE + exE*exE);
  
  //  cout << "bl: " << bl << " lo: " << lo << " hi: " << hi << endl;
  cout << Form("seen: %4.0f+/-%3.1f", bl, blE) 
       << Form(" expected: %4.1f+/-%3.1f", ex, exE) 
       << Form(" difference: %4.1f+/-%3.1f", df, dfE) 
       << endl;

  fBlExp  = ex; 
  fBlExpE = exE; 
  fBlObs  = (bl>0.?bl:0.01);
  fBlObsE = (bl>0.5?blE:1.);
  
}


// ----------------------------------------------------------------------
void anaBmm::rareBg() {

  cout << "normalize to " << fDataLumi[fSgData] << endl;

  c0->Clear();
  gStyle->SetOptStat(0);
  
  THStack *hRareBg0 = new THStack("hRareBg0","");
  THStack *hRareBg1 = new THStack("hRareBg1","");

  std::map<string, int> colors, hatches;
  std::map<string, double> mscale;  // pions: 0.002, kaons: 0.005, protons: 0.001
  std::map<string, double> err;  
  double epsPi = 0.003; 
  double epsKa = 0.003; 
  double epsPr = 0.0005;
  colors.insert(make_pair("bg60", 46)); hatches.insert(make_pair("bg60", 3004)); mscale.insert(make_pair("bg60", epsPi*epsPr)); 
  err.insert(make_pair("bg60", 0.3)); 
  colors.insert(make_pair("bg61", 49)); hatches.insert(make_pair("bg61", 3005)); mscale.insert(make_pair("bg61", epsKa*epsPr)); 
  err.insert(make_pair("bg61", 0.31)); 

  colors.insert(make_pair("bg82", 30)); hatches.insert(make_pair("bg82", 3004)); mscale.insert(make_pair("bg82", epsKa*epsKa)); 
  err.insert(make_pair("bg82", 0.15)); 
  colors.insert(make_pair("bg83", 32)); hatches.insert(make_pair("bg83", 3005)); mscale.insert(make_pair("bg83", epsPi*epsKa)); 
  err.insert(make_pair("bg83", 0.22)); 
  colors.insert(make_pair("bg84", 33)); hatches.insert(make_pair("bg84", 3007)); mscale.insert(make_pair("bg84", epsPi*epsPi)); 
  err.insert(make_pair("bg84", 1.0)); 

  colors.insert(make_pair("bg93", 40)); hatches.insert(make_pair("bg93", 3004)); mscale.insert(make_pair("bg93", epsKa*epsKa)); 
  err.insert(make_pair("bg93", 0.73)); 
  colors.insert(make_pair("bg92", 41)); hatches.insert(make_pair("bg92", 3005)); mscale.insert(make_pair("bg92", epsKa*epsPi)); 
  err.insert(make_pair("bg92", 0.05)); 
  colors.insert(make_pair("bg91", 42)); hatches.insert(make_pair("bg91", 3007)); mscale.insert(make_pair("bg91", epsPi*epsPi)); 
  err.insert(make_pair("bg91", 0.04)); 

  newLegend(0.65, 0.4, 0.80, 0.85); 

  double teff0(0.72), teff1(0.85); 
  double bsRare0(0.), bdRare0(0.), bsRare1(0.), bdRare1(0.); 

  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    if (string::npos == imap->first.find("bg")) {
      continue;
    }

    //     if (string::npos == imap->first.find("bg82")) {
    //       continue;
    //     }

    if (0 == fF[imap->first]) {
      continue; 
    }

    double lscale = fDataLumi[fSgData]/fLumi[imap->first];
    double misid = mscale[imap->first];
    cout << "==>" << imap->first << " lumi: " << fLumi[imap->first] << " -> lscale: " << lscale << " -> misid: " << misid << endl;
 
    fF[imap->first]->cd();
    loopTree(99); 
    TH1D *h1Rare0 = (TH1D*)(fhMassWithAllCuts[0]->Clone("h1Rare0"));  h1Rare0->SetLineColor(kBlack); 
    TH1D *h1Rare1 = (TH1D*)(fhMassWithAllCuts[1]->Clone("h1Rare1"));  h1Rare1->SetLineColor(kBlack);

    double bd0    = fhMassWithAllCutsManyBins[0]->Integral(fhMassWithAllCutsManyBins[0]->FindBin(fCuts[0]->mBdLo)+1, 
							fhMassWithAllCutsManyBins[0]->FindBin(fCuts[0]->mBdHi)-1);
    double bs0    = fhMassWithAllCutsManyBins[0]->Integral(fhMassWithAllCutsManyBins[0]->FindBin(fCuts[0]->mBsLo)+1, 
							fhMassWithAllCutsManyBins[0]->FindBin(fCuts[0]->mBsHi)-1);

    double bd1    = fhMassWithAllCutsManyBins[1]->Integral(fhMassWithAllCutsManyBins[1]->FindBin(fCuts[1]->mBdLo)+1, 
							fhMassWithAllCutsManyBins[1]->FindBin(fCuts[1]->mBdHi)-1);
    double bs1    = fhMassWithAllCutsManyBins[1]->Integral(fhMassWithAllCutsManyBins[1]->FindBin(fCuts[1]->mBsLo)+1, 
							fhMassWithAllCutsManyBins[1]->FindBin(fCuts[1]->mBsHi)-1);

    bsRare0 += bs0*lscale*misid*teff0; 
    bsRare1 += bs1*lscale*misid*teff1; 

    bdRare0 += bd0*lscale*misid*teff0; 
    bdRare1 += bd1*lscale*misid*teff1; 

    cout << imap->first << ": bsRare0 increment by " << bs0*lscale*misid*teff0 << "  " << bs0 << endl;
    cout << imap->first << ": bsRare1 increment by " << bs1*lscale*misid*teff1 << "  " << bs1 << endl;
    cout << imap->first << ": bdRare0 increment by " << bd0*lscale*misid*teff0 << "  " << bd0 << endl;
    cout << imap->first << ": bdRare1 increment by " << bd1*lscale*misid*teff1 << "  " << bd1 << endl;

    
    h1Rare0->SetFillColor(colors[imap->first]);
    h1Rare1->SetFillColor(colors[imap->first]);
    h1Rare0->SetFillStyle(1000);
    h1Rare1->SetFillStyle(1000);
    
    h1Rare0->Scale(lscale*misid);
    h1Rare1->Scale(lscale*misid);

    legg->AddEntry(h1Rare0, fName[imap->first].c_str(), "f"); 
    h1Rare0->Draw();
    c0->Modified();
    c0->Update();

    hRareBg0->Add(h1Rare0); 
    hRareBg1->Add(h1Rare1); 
    
  }
							 
  gStyle->SetOptStat(0);

  shrinkPad(0.12, 0.18); 
  hRareBg0->Draw();
  TH1D *hhRareBg0 = (TH1D*)hRareBg0->GetHistogram(); 
  setTitles(hhRareBg0, "m_{h h} [GeV]", "Cands/bin", 0.06, 0.9, 1.5);
  legg->Draw(); 
  hhRareBg0->Draw("same");
  string pdfname = Form("%s/%s_rare0.pdf", fDirectory.c_str(), fSuffix.c_str());
  if (fDoPrint) {
    if (c0) c0->SaveAs(pdfname.c_str());
  }

  hRareBg1->Draw();
  TH1D *hhRareBg1 = (TH1D*)hRareBg1->GetHistogram(); 
  setTitles(hhRareBg1, "m_{h h} [GeV]", "Cands/bin", 0.06, 0.9, 1.5);
  legg->Draw(); 
  hhRareBg1->Draw("same");
  pdfname = Form("%s/%s_rare1.pdf", fDirectory.c_str(), fSuffix.c_str());
  if (fDoPrint){
    if (c0) c0->SaveAs(pdfname.c_str());
  }

  cout << "0 mBsLo: " << fCuts[0]->mBsLo << " mBsHi: " <<  fCuts[0]->mBsHi << endl;
  cout << "1 mBsLo: " << fCuts[1]->mBsLo << " mBsHi: " <<  fCuts[1]->mBsHi << endl;
  cout << "0 mBdLo: " << fCuts[0]->mBdLo << " mBdHi: " <<  fCuts[0]->mBdHi << endl;
  cout << "1 mBdLo: " << fCuts[1]->mBdLo << " mBdHi: " <<  fCuts[1]->mBdHi << endl;

  cout << "bsRare0: " << bsRare0 <<  " bsRare1: " << bsRare1 << endl;
  cout << "bdRare0: " << bdRare0 <<  " bdRare1: " << bdRare1 << endl;

  numbers *aa = fNumbersBs[0];
  aa->bsRare = bsRare0; 
  aa->bdRare = bdRare0; 

  aa = fNumbersBs[1]; 
  aa->bsRare = bsRare1; 
  aa->bdRare = bdRare1; 

  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << " % --- rare Background numbers" << endl;  
  fTEX <<  Form("\\vdef{%s:bsRare0}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), bsRare0) << endl;
  fTEX <<  Form("\\vdef{%s:bsRare1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), bsRare1) << endl;
  fTEX <<  Form("\\vdef{%s:bdRare0}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), bdRare0) << endl;
  fTEX <<  Form("\\vdef{%s:bdRare1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), bdRare1) << endl;

}


// ----------------------------------------------------------------------
void anaBmm::sbsDistributionOverlaySameFile(std::string file1, std::string r1, std::string r2, std::string sel, std::string L1, std::string L2) {

  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 

  char loption1[100], loption2[100]; 
  int color1(1), marker1(20); 

  sprintf(loption1, "f");
  sprintf(loption2, "p");

  if (string::npos != file1.find("No")) {
    color1 = kBlue; 
  }
  if (string::npos != file1.find("Cs")) {
    color1 = kRed; 
  }

  vector<string> leftList; 
  leftList.push_back("iso"); 


  fF[file1]->cd(); 
  cout << "==> File1: "; gFile->pwd(); 
  TH1D *h = (TH1D*)gFile->Get("analysisDistributions"); 
  if (0 == h) {
    cout << "no histogram analysisDistributions found, returning" << endl;
    return;
  }

  TH1D *h1, *h2;
  AnalysisDistribution a("A_allevents"); 

  string pdfname = Form("%s/%s_%s_%s_%s_sbs_%s.pdf", 
			fDirectory.c_str(), fSuffix.c_str(), file1.c_str(), r1.c_str(), r2.c_str(), sel.c_str());
  
  fF[file1]->cd(); 
  cout << "==> File1: "; gFile->pwd(); 
  cout << "==> pdf: " << pdfname << endl;

  bool mm(false); 

  if (mm) {
    h1 = (TH1D*)gDirectory->Get(Form("%s%s1", r1.c_str(), sel.c_str()));
  } else {
    cout << "r1: "<< r1 << endl;
    h1 = a.sbsDistribution(r1.c_str(), sel.c_str());
    c0->cd();
    h1->Draw();
  }

  if (mm) {
    h2 = (TH1D*)gDirectory->Get(Form("%s%s0", r2.c_str(), sel.c_str()));
  } else {
    cout << "r2: "<< r2 << endl;
    h2 = a.sbsDistribution(r2.c_str(), sel.c_str()) ;
    c0->cd();
    h2->Draw();
  }

  if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  
  c0->cd(); 
  c0->Clear(); 
  h1->SetMinimum(0);
  double dmax = (h1->GetMaximum() > h2->GetMaximum()? 1.1*h1->GetMaximum(): 1.1*h2->GetMaximum()); 
  h1->SetMaximum(dmax); 
  setHist(h1, color1, marker1, 1.5); 
  h1->SetTitle("");
  h1->Draw("hist");
  
  setHist(h2, color1, marker1, 1.5); 
  h2->Draw("samee");

  cout << "r1: " << r1 << endl;
  cout << "h1: " << h1->GetName() << endl;
  cout << "L1: " << L1 << endl;

//  newLegend(0.25, 0.7, 0.50, 0.85); 
//  if (leftList.end() != find(leftList.begin(), leftList.end(), search)) {
//   } else {
//   }
  newLegend(0.45, 0.7, 0.70, 0.85); 
  legg->AddEntry(h1, L1.c_str(), loption1); 
  legg->AddEntry(h2, L2.c_str(), loption2); 
  legg->Draw(); 

  if (fDoPrint) c0->SaveAs(pdfname.c_str()); 
  
}


// ----------------------------------------------------------------------
void anaBmm::sbsDistributionOverlay(std::string file1, std::string file2, const char *selection, const char *region) {

  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 

  char option1[100], option2[100],  loption1[100], loption2[100]; 
  int color1(1), fill1(1000), color2(1), fill2(1000), marker1(20), marker2(25); 
  int mm(0); 

  if (string::npos != file1.find("No")) {
    color1 = kBlue; 
    fill1  = 3004; 
  }
  if (string::npos != file2.find("No")) {
    color2 = kBlue; 
    fill2  = 3004; 
  }
  if (string::npos != file1.find("Cs")) {
    color1 = kRed; 
    fill1  = 3005; 
  }
  if (string::npos != file2.find("Cs")) {
    color2 = kRed; 
    fill2  = 3005; 
  }

  if (string::npos != file2.find("Sg")) {
    // -- data is first, then signal MC
    mm = 1; 
    color1 = kBlack; 
    fill1  = 0; 
    sprintf(option1, "e"); 
    sprintf(loption1, "p");
    color2 = kBlue; 
    fill2  = 3005; 
    sprintf(option2, "hist"); 
    sprintf(loption2, "f");
  }
  
  // -- now fix overlays of two data files
  if (0 == mm ) {
    sprintf(option1, "e");
    sprintf(loption1, "p");
    
    sprintf(option2,  "hist");
    sprintf(loption2, "f");
    
    if (string::npos != file2.find("Data")) {
      sprintf(option1, "hist"); 
      sprintf(loption1, "f");
      sprintf(option2, "e"); 
      sprintf(loption2, "p");
      fill2   = 0; 
      color2  = kBlack; 
      marker2 = 21; 
    }
    
    if (string::npos != file1.find("Mc")) {
      sprintf(option1, "hist"); 
      sprintf(loption1, "f");
      sprintf(loption2, "p");
      fill2   = 0; 
      color2  = kBlack; 
      marker2 = 25; 
    }
  }
 
  //  ofstream OUT("testUL.txt", ios::app);

  fF[file1]->cd(); 
  cout << "==> File1: "; gFile->pwd(); 
  TH1D *h = (TH1D*)gFile->Get("analysisDistributions"); 
  if (0 == h) {
    cout << "no histogram analysisDistributions found, returning" << endl;
    return;
  }

  TH1D *h1, *h2;
  AnalysisDistribution a("A_allevents"); 

  vector<string> skipList; 
  skipList.push_back("mpsi"); 
  skipList.push_back("q"); 
  skipList.push_back("allevents"); 
  skipList.push_back("muonsid"); 
  skipList.push_back("tracksqual"); 
  skipList.push_back("hlt"); 
  skipList.push_back("allevents"); 

  skipList.push_back("pv1"); 
  skipList.push_back("pv2"); 
  skipList.push_back("pv3"); 
  skipList.push_back("pv4"); 
  skipList.push_back("pv5"); 
  skipList.push_back("pv6"); 

  vector<string> nolegendList; 
  nolegendList.push_back("muonseta"); 
  nolegendList.push_back("muon1eta"); 
  nolegendList.push_back("muon2eta"); 
  nolegendList.push_back("trackseta"); 
  nolegendList.push_back("pvz"); 
  nolegendList.push_back("eta"); 
  nolegendList.push_back("ip1"); 
  nolegendList.push_back("ip2"); 
  nolegendList.push_back("pvz"); 

  vector<string> leftList; 
  leftList.push_back("iso"); 
  leftList.push_back("iso1"); 
  leftList.push_back("iso2"); 
  leftList.push_back("iso3"); 
  leftList.push_back("iso4"); 
  leftList.push_back("cosa"); 
  leftList.push_back("cosa0"); 

  TCanvas *c1;
  string skipcut, cut, pdfname; 
  bool skip; 
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i)); 
    skipcut = cut.substr(cut.find("_")+1); 
    skip = false;
    //    cout << "skipcut: " << skipcut << endl;
    for (unsigned int j = 0; j < skipList.size(); ++j) {
      if (string::npos != skipcut.find(skipList[j])) {
	//	cout << "found " << skipcut << " in skiplist" << endl;
	skip = true; 
	break;
      }
    }
    if (skip) continue;
    // -- skip  if a specific region is requested
    if (strcmp(region, "") && string::npos == cut.find(region)) {
      //      cout << "region " << region << " not part of name, skipping" << endl;
      continue;
    }

    if (cut == string("")) {
      //      cout << "empty string found, break ..." << endl;
      //      OUT.close();
      break;
    }
    pdfname = Form("%s/%s_%s-%s_sbs_%s_%s.pdf", fDirectory.c_str(), fSuffix.c_str(),
		   file1.c_str(), file2.c_str(), cut.c_str(), selection);
  
    fF[file1]->cd(); 
    cout << "==> File1: "; gFile->pwd(); 
    cout << "==> pdf: " << pdfname << endl;
    if (mm) {
      // -- normal
      h1 = (TH1D*)gDirectory->Get(Form("%s%s1", cut.c_str(), selection));
      // -- validation of signal MC!
      //       cout << "WARNING MC VALIDATION" << endl;
      //       h1 = (TH1D*)gDirectory->Get(Form("%s%s2", cut.c_str(), selection));
    } else {
      h1 = a.sbsDistribution(cut.c_str(), selection);
    }

    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint){
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file1.c_str(), cut.c_str(), selection));
    }

    fF[file2]->cd(); 
    cout << "==> File2: "; gFile->pwd(); 
    if (mm) {
      // -- normal
      h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection));
      // -- validation of signal MC!
      //       cout << "WARNING MC VALIDATION" << endl;
      //       h2 = (TH1D*)gDirectory->Get(Form("%s%s2", cut.c_str(), selection));
    } else {
      h2 = a.sbsDistribution(cut.c_str(), selection) ;
    }
    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint) {
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file2.c_str(), cut.c_str(), selection));
    }

    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    c0->cd(); 
    c0->Clear(); 
    h1->SetMinimum(0);
    double dmax = (h1->GetMaximum() > h2->GetMaximum()? 1.1*h1->GetMaximum(): 1.1*h2->GetMaximum()); 
    h1->SetMaximum(dmax); 
    setHist(h1, color1, marker1, 1.5); 
    setFilledHist(h1, color1, color1, fill1); 
    h1->SetTitle("");
    h1->Draw(option1);

    setHist(h2, color2, marker2, 1.5); 
    setFilledHist(h2, color2, color2, fill2); 
    h2->Draw(Form("same%s", option2));

    if (nolegendList.end() != find(nolegendList.begin(), nolegendList.end(), cut)) {
      cout << "++++++++++++++++++++++" << endl;
    } else {
      if (leftList.end() != find(leftList.begin(), leftList.end(), cut)) {
	newLegend(0.25, 0.7, 0.50, 0.85); 
      } else {
	newLegend(0.45, 0.7, 0.70, 0.85); 
      }
      legg->AddEntry(h1, fName[file1].c_str(), loption1); 
      legg->AddEntry(h2, fName[file2].c_str(), loption2); 
      legg->Draw(); 
    }
    if (fDoPrint) c0->SaveAs(pdfname.c_str()); 
  }
} 


// ----------------------------------------------------------------------
void anaBmm::breco(TH1D *h) {
  fPeak = 0.;
  fWidth = 0.;
}


// ----------------------------------------------------------------------
void anaBmm::optimizeCut(const char *cut, double lo, double hi, const char *otherCuts) {

  int NSTEPS(10);

  double delta = hi - lo;
  double cutVal; 

  for (int j = 0; j < NSTEPS; ++j) {
    cutVal = lo + j*delta/NSTEPS;
    string cutline(""); 
    cutline += Form("%s %5.4f && %s", cut, cutVal, otherCuts); 
    cout << Form("%3d %s", j, cutline.c_str()) << endl;
    testSimpleUL(cutline.c_str()); 
    c0->Modified(); 
    c0->Update();
  }
} 


// ----------------------------------------------------------------------
void anaBmm::plotVar(const char* plotstring, const char *cuts, const char *options) {

  TTree *td = (TTree*)(fpData[fSgData]->Get("events"));
  td->Draw(plotstring, Form("hlt&&gmuid &&m > 4.8 && m < 6.0 && %s", cuts), options);

}




// ----------------------------------------------------------------------
void anaBmm::playUL(int nruns) {

  int NCUTS(6);
 
  string cuts[]   = {"m2pt>", "cosa>", "iso1>", "chi2<", "fls3d>", "docatrk>"}; 
  double loCuts[] = {2.5,    0.99,    0.0,     1.,       5,        0.0};
  double hiCuts[] = {3.5,    0.9999,  0.99,    5.,       15,       0.05};
  
  string cutline; 
  double cut; 
  for (int j = 0; j < nruns; ++j) {
    cutline = "hlt&&gmuid";
    for (int i = 0; i < NCUTS; ++i) {
      cut = gRandom->Rndm()*(hiCuts[i]-loCuts[i]) + loCuts[i];
      cutline += Form(" && %s %5.4f", cuts[i].c_str(), cut); 
    }

    cout << cutline << endl;
    testSimpleUL(cutline.c_str()); 
  }
}


// ----------------------------------------------------------------------
void anaBmm::filterEff(int mode) {

  double pTCut(0.), pCut(0.), etaCut(0.); 
  if (0 == mode) {
    pTCut = -1.;
    pCut = -1.;
    etaCut = 2.5;
  }

  if (10 == mode) {
    pTCut = 2.5;
    pCut = -1.;
    etaCut = 2.5;
  }

  if (20 == mode) {
    pTCut = -1.;
    pCut = 2.5;
    etaCut = 2.5;
  }

  
  int   bprocid; 
  float bg1pt, bg2pt, bg1eta, bg2eta, bm;

  TTree *t  = (TTree*)(gFile->Get("effTree"));
  t->SetBranchAddress("procid",&bprocid);

  t->SetBranchAddress("m",&bm);
  t->SetBranchAddress("g1pt",&bg1pt);
  t->SetBranchAddress("g2pt",&bg2pt);
  t->SetBranchAddress("g1eta",&bg1eta);
  t->SetBranchAddress("g2eta",&bg2eta);

  int nentries = Int_t(t->GetEntries());
  int nb(0); 
  int ngen(0), nfilter(0);
  TVector3 pm1, pm2;
  cout << "pTCut:  " << pTCut << endl;
  cout << "pCut:   " << pCut << endl;
  cout << "etaCut: " << etaCut << endl;
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    //    if (bm < 0) continue;
    ++ngen;
    pm1.SetPtEtaPhi(bg1pt, bg1eta, 0.);
    pm2.SetPtEtaPhi(bg2pt, bg2eta, 0.);
    // cout << bg1pt << " " << bg2pt << " " << bg1eta << " " << bg2eta << " -> " << pm1.Mag() << " " << pm2.Mag() << endl;
    if (TMath::Abs(bg1eta) > etaCut) continue;
    if (TMath::Abs(bg2eta) > etaCut) continue;
    if (bg1pt < pTCut) continue;
    if (bg2pt < pTCut) continue;
    if (pm1.Mag() < pCut) continue;
    if (pm2.Mag() < pCut) continue;
    ++nfilter;
  }
  
  cout << "nfilter/ngen = " << nfilter << "/" << ngen << " = " << static_cast<double>(nfilter)/ngen << endl;
}


// ----------------------------------------------------------------------
void anaBmm::accEffFromEffTree(string fname, numbers &a, cuts &b, int proc) {

  TFile *f = fF[fname];
  if (0 == f) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): no file " << fname << " found " << endl;
    return;
  }
  TTree *t  = (TTree*)(f->Get("effTree"));
  double effFilter(1.); 
  if (!t) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): no tree `effTree' found " << endl;
    f->pwd(); 
    return;
  } else {
    effFilter = fFilterEff[fname]; 
    cout << "anaBmm::accEffFromEffTree(" << a.name << ")" << endl
	 << " get acceptance from file " << f->GetName() 
	 << " with filterEff = " << effFilter
	 << endl;
  }

  bool sg(false), no(false), cs(false); 

  int   bprocid, bidx; 
  bool  bhlt;
  float bg1pt, bg2pt, bg1eta, bg2eta;
  float bm1pt, bm1eta, bm2pt, bm2eta;
  bool  bm1gt, bm2gt; 
  bool  bm1id, bm2id; 
  float bm;

  float bg3pt, bg4pt, bg3eta, bg4eta; 
  float bk1pt, bk2pt, bk1eta, bk2eta; 
  bool  bk1gt, bk2gt; 

  t->SetBranchAddress("hlt",&bhlt);
  t->SetBranchAddress("procid",&bprocid);
  t->SetBranchAddress("bidx",&bidx);

  t->SetBranchAddress("g1pt",&bg1pt);
  t->SetBranchAddress("g2pt",&bg2pt);
  t->SetBranchAddress("g1eta",&bg1eta);
  t->SetBranchAddress("g2eta",&bg2eta);

  t->SetBranchAddress("m1pt",&bm1pt);
  t->SetBranchAddress("m2pt",&bm2pt);
  t->SetBranchAddress("m1eta",&bm1eta);
  t->SetBranchAddress("m2eta",&bm2eta);

  t->SetBranchAddress("m1gt",&bm1gt);
  t->SetBranchAddress("m2gt",&bm2gt);
  t->SetBranchAddress("m1id",&bm1id);
  t->SetBranchAddress("m2id",&bm2id);


  if (string::npos != a.name.find("signal")) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): SIGNAL " << endl;
    sg = true; 
  }

  if (string::npos != a.name.find("normalization")) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): NORMALIZATION " << endl;
    no = true; 
    t->SetBranchAddress("g3pt", &bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
    t->SetBranchAddress("k1pt", &bk1pt);
    t->SetBranchAddress("k1eta",&bk1eta);
    t->SetBranchAddress("k1gt", &bk1gt);
  }

  if (string::npos != a.name.find("control sample")) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): CONTROL SAMPLE " << endl;
    cs = true; 
    t->SetBranchAddress("g3pt", &bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
    t->SetBranchAddress("k1pt", &bk1pt);
    t->SetBranchAddress("k1eta",&bk1eta);
    t->SetBranchAddress("k1gt", &bk1gt);

    t->SetBranchAddress("g4pt", &bg4pt);
    t->SetBranchAddress("g4eta",&bg4eta);
    t->SetBranchAddress("k2pt", &bk2pt);
    t->SetBranchAddress("k2eta",&bk2eta);
    t->SetBranchAddress("k2gt", &bk2gt);
  }

  t->SetBranchAddress("m",&bm);


  int nentries = Int_t(t->GetEntries());
  int nb(0); 
  int ngen(0), nchangen(0), nreco(0), nchan(0), nmuid(0), nhlt(0), ncand(0); 
  int chan(-1); 
  cout << "channel = " << a.index << endl;
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    if (bidx < 0) continue;
    ++ngen;
    if (proc > 0 && bprocid != proc) continue;
    if (sg) {
      // -- Signal
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1.) {
	    if (bm1pt > 1. && bm2pt > 1.
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4
		&& bm1gt && bm2gt
		) {
	      chan = detChan(bm1eta, bm2eta); 
	      if (chan == a.index) {
		++nreco;
		if (bm1pt > b.m1pt && bm2pt > b.m2pt 
		    ) {
		  ++nchan; 
		  if (bm1id && bm2id) {
		    ++nmuid;
		    if (bhlt) {
		      ++nhlt;
		      if (bm > 0) {
			++ncand;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } else if (no) {
      // -- Normalization
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5 && TMath::Abs(bg3eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1. && bg3pt > 0.4) {
	    if (bm1pt > 1. && bm2pt > 1. && bk1pt > 0.5
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4 && TMath::Abs(bk1eta) < 2.4
		&& bm1gt && bm2gt && bk1gt
		) {
	      chan = detChan(bm1eta, bm2eta); 
	      if (chan == a.index) {
		++nreco;
		if (bm1pt > b.m1pt && bm2pt > b.m2pt
		    ) {
		  ++nchan; 
		  if (bm1id && bm2id) {
		    ++nmuid;
		    if (bhlt) {
		      ++nhlt;
		      if (bm > 0) {
			++ncand;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } else if (cs) {
      // -- Control sample
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5 && TMath::Abs(bg3eta) < 2.5 && TMath::Abs(bg4eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1. && bg3pt > 0.4 && bg4pt > 0.4) {
	    if (bm1pt > 1. && bm2pt > 1. && bk1pt > 0.5 && bk2pt > 0.5
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4 && TMath::Abs(bk1eta) < 2.4 && TMath::Abs(bk2eta) < 2.4
		&& bm1gt && bm2gt && bk1gt && bk2gt
		) {
	      chan = detChan(bm1eta, bm2eta); 
	      if (chan == a.index) {
		++nreco;
		if (bm1pt > b.m1pt && bm2pt > b.m2pt) {
		  ++nchan; 
		  if (bm1id && bm2id) {
		    ++nmuid;
		    if (bhlt) {
		      ++nhlt;
		      if (bm > 0) {
			++ncand;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  a.effGenFilter  = effFilter; 
  a.genFileYield  = ngen;
  a.genYield      = a.genFileYield/effFilter;
  a.genChanYield  = nchangen; 
  a.recoYield     = nreco; // reco'ed in chan, basic global reconstruction cuts 
  a.chanYield     = nchan; // reco'ed in chan, with channel-dependent (pT) cuts
  a.muidYield     = nmuid;
  a.trigYield     = nhlt;
  a.candYield     = ncand;
  
  if (a.genYield > 0) {
    a.cFrac  = a.genChanYield/a.genYield;
    a.cFracE = dEff(static_cast<int>(a.genChanYield), static_cast<int>(a.genYield));
  }
  if (a.genYield > 0) {
    a.acc = a.recoYield/a.genYield;
    a.accE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genYield));
  }
  
  if (a.genChanYield > 0) {
    a.accChan = a.recoYield/a.genChanYield;
    a.accChanE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genChanYield));
  }
  
  if (a.recoYield > 0) {
    a.effChan = a.chanYield/a.recoYield;
    a.effChanE = dEff(static_cast<int>(a.chanYield), static_cast<int>(a.recoYield));
  }


  if (a.trigYield > 0) {
    a.effCand = a.candYield/a.trigYield;
    a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.trigYield));
  } 

  //   if (0) {
  //     if (a.chanYield > 0) {
  //       a.accMuidMC = a.muidYield/a.chanYield;
  //       a.accMuidMCE = dEff(static_cast<int>(a.muidYield), static_cast<int>(a.chanYield));
  //     } 
  //     if (a.muidYield > 0) {
  //       a.accTrigMC = a.trigYield/a.muidYield;
  //       a.accTrigMCE = dEff(static_cast<int>(a.trigYield), static_cast<int>(a.muidYield));
  //     } 
  //   }
}


// ----------------------------------------------------------------------
void anaBmm::computeNormUL() {
  cout << "--> rareBg" << endl;
  rareBg();

  cout << "--> loopTree: signal MC" << endl;
  loopTree(0);  // signal eff
  c0->Modified(); c0->Update();
  loopTree(1);  // Bd2MuMu eff
  c0->Modified(); c0->Update();
  cout << "--> loopTree: signal data" << endl;
  loopTree(5);  // data signal
  c0->Modified(); c0->Update();
  cout << "--> loopTree: norm MC" << endl;
  loopTree(10); // normalization eff
  c0->Modified(); c0->Update();
  cout << "--> loopTree: norm data" << endl;
  loopTree(11); // data normalization 
  c0->Modified(); c0->Update();

  printUlcalcNumbers();

  fNobs = static_cast<int>(fBgExp + 0.5);
  
  double ulbarlow = barlow(fNobs, fBgExp, fBgExpE, 0.2);
  double alpha = 1.; 

  if (fBgExp < 0.001) fBgExp = 0.1;
  double ulBayes90 =  blimit(0.9, fNobs, 1.0, 0.2, fBgExp, fBgExpE, alpha);
  double ulBayes95 =  blimit(0.95, fNobs, 1.0, 0.2, fBgExp, fBgExpE, alpha);
  
  cout << "Bayes: " << ulBayes90 << "(95%CL: " << ulBayes95  << ") vs barlow: " << ulbarlow << endl;

  fNul = ulBayes90;

  cout << "printing fNumbersBs[0]" << endl;
  printNumbers(*fNumbersBs[0], cout); 
  printNumbers(*fNumbersBs[0], fOUT); 

  cout << "printing fNumbersBd[0]" << endl;
  printNumbers(*fNumbersBd[0], cout); 
  printNumbers(*fNumbersBd[0], fOUT); 

  cout << "printing fNumbersNorm[0]" << endl;
  printNumbers(*fNumbersNorm[0], cout); 
  printNumbers(*fNumbersNorm[0], fOUT); 

  fUL = (fNul/fNormSig)
    *(fu/fs)
    *(fNumbersNorm[0]->acc/fNumbersBs[0]->acc)
    *(fNumbersNorm[0]->effCand/fNumbersBs[0]->effCand)     
    *(fNumbersNorm[0]->effMuidMC/fNumbersBs[0]->effMuidMC)
    *(fNumbersNorm[0]->effTrigMC/fNumbersBs[0]->effTrigMC)
    *(fNumbersNorm[0]->effAna/fNumbersBs[0]->effAna)
    * fBF;

  cout << "prod(eff) expected UL: " << fUL << endl;

  fUL = (fNul/fNormSig)
    *(fu/fs)
    *(fNumbersNorm[0]->effTot/fNumbersBs[0]->effTot)
    * fBF;

  cout << "effTot expected UL:    " << fUL << endl;

  //  system(Form("../ulcalc/bin/ulcalc %s", fUlcalcFileName.c_str())); 

}


// ----------------------------------------------------------------------
void anaBmm::computeCsBF() {
  cout << "--> loopTree: CS MC" << endl;
  loopTree(20);  // CS signal eff
  c0->Modified(); c0->Update();
  cout << "--> loopTree: signal data" << endl;
  loopTree(21);  // control sample data 
  c0->Modified(); c0->Update();
  cout << "--> loopTree: norm MC" << endl;
  loopTree(10); // normalization eff
  c0->Modified(); c0->Update();
  cout << "--> loopTree: norm data" << endl;
  loopTree(11); // data normalization 
  c0->Modified(); c0->Update();

  fNobs = static_cast<int>(fBgExp + 0.5);

  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;    
  fTEX << "% -- Control sample branching fraction" << endl;
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;    
  
  double result, resultE; 
  for (int i = 0; i < 2; ++i) {
    result = (fNumbersCS[i]->fitYield/fNumbersNorm[i]->fitYield)
      *(fu/fs)
      *(fNumbersNorm[i]->acc/fNumbersCS[i]->acc)
      *(fNumbersNorm[i]->effCand/fNumbersCS[i]->effCand)     
      *(fNumbersNorm[i]->effAna/fNumbersCS[i]->effAna)
      *(fNumbersNorm[i]->effMuidPid/fNumbersCS[i]->effMuidPid)
      *(fNumbersNorm[i]->effTrigPid/fNumbersCS[i]->effTrigPid)
      * fBF;

    
    resultE = dRatio(fNumbersCS[i]->fitYield, fNumbersCS[i]->fitYieldE, fNumbersNorm[i]->fitYield, fNumbersNorm[i]->fitYieldE)
      *(fu/fs)
      *(fNumbersNorm[i]->acc/fNumbersCS[i]->acc)
      *(fNumbersNorm[i]->effCand/fNumbersCS[i]->effCand)     
      *(fNumbersNorm[i]->effMuidPid/fNumbersCS[i]->effMuidPid)
      *(fNumbersNorm[i]->effTrigPid/fNumbersCS[i]->effTrigPid)
      *(fNumbersNorm[i]->effAna/fNumbersCS[i]->effAna)
      * fBF;

    cout << "chan " << i << ": PID fact branching fraction: " << result << "+/-" << resultE << endl;
    fTEX << formatTex(result, Form("%s:N-CSBF-Pid-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(resultE, Form("%s:N-CSBF-Pid-BS%i:err", fSuffix.c_str(), i), 6) << endl;

    result = (fNumbersCS[i]->fitYield/fNumbersNorm[i]->fitYield)
      *(fu/fs)
      *(fNumbersNorm[i]->acc/fNumbersCS[i]->acc)
      *(fNumbersNorm[i]->effCand/fNumbersCS[i]->effCand)     
      *(fNumbersNorm[i]->effAna/fNumbersCS[i]->effAna)
      *(fNumbersNorm[i]->effMuidMC/fNumbersCS[i]->effMuidMC)
      *(fNumbersNorm[i]->effTrigMC/fNumbersCS[i]->effTrigMC)
      * fBF;

    
    resultE = dRatio(fNumbersCS[i]->fitYield, fNumbersCS[i]->fitYieldE, fNumbersNorm[i]->fitYield, fNumbersNorm[i]->fitYieldE)
      *(fu/fs)
      *(fNumbersNorm[i]->acc/fNumbersCS[i]->acc)
      *(fNumbersNorm[i]->effCand/fNumbersCS[i]->effCand)     
      *(fNumbersNorm[i]->effMuidMC/fNumbersCS[i]->effMuidMC)
      *(fNumbersNorm[i]->effTrigMC/fNumbersCS[i]->effTrigMC)
      *(fNumbersNorm[i]->effAna/fNumbersCS[i]->effAna)
      * fBF;

    cout << "chan " << i << ": MC fact branching fraction: " << result << "+/-" << resultE << endl;
    fTEX << formatTex(result, Form("%s:N-CSBF-MC-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(resultE, Form("%s:N-CSBF-MC-BS%i:err", fSuffix.c_str(), i), 6) << endl;
    
    result = (fNumbersCS[i]->fitYield/fNumbersNorm[i]->fitYield)
      *(fu/fs)
      *(fNumbersNorm[i]->effTot/fNumbersCS[i]->effTot)
      * fBF;

    cout << "chan " << i << ": branching fraction: " << result << "+/-" << resultE << endl;

    fTEX << formatTex(result, Form("%s:N-CSBF-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(resultE, Form("%s:N-CSBF-BS%i:err", fSuffix.c_str(), i), 6) << endl;

  }

  printCsBFNumbers();
}

// ----------------------------------------------------------------------
void anaBmm::effTree(int mode) {

}


// ----------------------------------------------------------------------
void anaBmm::puEff(const char *var, double cut, const char *ylabel, const char *file1, const char *selection) {

  TH1D *h[6];
  
  fF[file1]->cd(); 
  AnalysisDistribution a("A_allevents"); 
  cout << "==> File1: "; gFile->pwd(); 
  c0->Clear();
  c0->Divide(2,3);
  for (int i = 0; i < 6; ++i) {
    cout << "    ==> " << Form("%spv%d", var, i+1) << endl;
    h[i] = a.sbsDistribution(Form("%spv%d", var, i+1), selection);
    c0->cd(i+1);
    // -- check that the bins are not negative. If they are, reset to zero
    for (int ix = 1; ix <= h[i]->GetNbinsX(); ++ix) {
      if (h[i]->GetBinContent(ix) < 0) {
	h[i]->SetBinContent(ix, 0.); 
	h[i]->SetBinError(ix, 0.); 
      }	
    }

    h[i]->Draw();
  }

  TH1D *heff = new TH1D("heff", "", 6, 0., 12.);
  for (int i = 0; i < 6; ++i) {
    double ntot = h[i]->GetSumOfWeights();
    double ncut = h[i]->Integral(h[i]->FindBin(cut)+1, h[i]->GetNbinsX());
    double eff  = ncut/ntot; 
    double stat = dEff(static_cast<int>(ncut), static_cast<int>(ntot));
    double effE = TMath::Sqrt(stat*stat + 0.05*eff*0.05*eff);
    // -- can go above 1, if single bins with too large negative entries are present
    if (eff > 1) {
      eff = 0.; 
      effE = 0.;
    }
    cout << "h[i]->FindBin(cut) = " << h[i]->FindBin(cut) << endl;
    cout << "h[i]->GetNbinsX()  = " << h[i]->GetNbinsX() << endl;
    cout << "ntot = " << ntot << endl;
    cout << "ncut = " << ncut << endl;
    cout << "eff = " << eff << endl;
    cout << "effE = " << effE << endl;
    heff->SetBinContent(i+1, eff); 
    heff->SetBinError(i+1, effE); 
  }

  c0->Clear();
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(1); 
  setTitles(heff, "N_{PV}", Form("%s", ylabel)); 
  heff->SetMaximum(1.2);
  if (heff->GetMinimum() > 0.5) {
    heff->SetMinimum(0.5);
  } else {
    heff->SetMinimum(0.);
  }
  heff->SetMarkerSize(2); 
  heff->Fit("pol1", "FM");
  heff->GetFunction("pol1")->SetLineWidth(2);
  
  if (fDoPrint)  c0->SaveAs(Form("%s/pu-eff-%s-%s-0_%d.pdf", fDirectory.c_str(), var, file1, static_cast<int>(100.*cut)));

}


// ----------------------------------------------------------------------
void anaBmm::plotWithCut(const char *var, const char *cuts, double cut, const char *title, double hmin, double hmax) {

  fpMc[fSgMc]->cd(); 
  TH1D *h1 = new TH1D("plotWithCut1", "", 40, hmin, hmax); setTitles(h1, title, "Entries/Bin"); 
  TTree *t = (TTree*)gFile->Get("events"); 
  t->Draw(Form("%s>>plotWithCut1", var), cuts, "goff"); 

  fpData[fSgData]->cd(); 
  TH1D *h2 = new TH1D("plotWithCut2", "", 40, hmin, hmax); setTitles(h2, title, "Entries/Bin"); 
  t = (TTree*)gFile->Get("events"); 
  t->Draw(Form("%s>>plotWithCut2", var), cuts, "goff"); 

  double imax = 1.2*h1->GetMaximum(); 
  cout << "imax = " << imax << endl;
  if (h2->GetMaximum() > imax) imax = 1.2*h2->GetMaximum(); 
  cout << "imax = " << imax << endl;
  h1->SetMaximum(imax); 

  gStyle->SetOptStat(0); 
  h1->SetLineWidth(2);
  h1->Draw();
  gPad->SetLogy(1); 
  h2->Draw("samee");

  pa->DrawArrow(cut, 0.2*imax, cut, h1->GetMinimum()); 

  newLegend(0.65, 0.75, 0.80, 0.85); 
  legg->SetTextSize(0.035);  
  legg->AddEntry(h1, "Signal MC", "l"); 
  legg->AddEntry(h2, "Data", "p"); 
  legg->Draw();
  
  if (fDoPrint)  c0->SaveAs(Form("%s/plotWithCut-%s.pdf", fDirectory.c_str(), var));
  c0->SetLogy(0); 
}


// ----------------------------------------------------------------------
void anaBmm::compareTrigger(double m1pt, double m2pt) {
  
  for (unsigned int i = 0; i < fNchan; ++i) {
    fCuts[i]->m1pt = m1pt; 
    fCuts[i]->m2pt = m2pt; 
  }
  
  loopTree(0); 
  loopTree(10); 

  cout << Form("pT > %3.1f muid: ", m2pt) 
       << Form("%3.2f %3.2f %3.2f %3.2f", 
	       fNumbersBs[0]->effMuidMC, fNumbersBs[0]->effMuidPid, 
	       fNumbersBs[1]->effMuidMC, fNumbersBs[1]->effMuidPid)
       << " " 
       << Form("%3.2f %3.2f %3.2f %3.2f", 
	       fNumbersNorm[0]->effMuidMC, fNumbersNorm[0]->effMuidPid, 
	       fNumbersNorm[1]->effMuidMC, fNumbersNorm[1]->effMuidPid)
       << " -> " 
       << Form(" rho_B = %3.2f rho_E = %3.2f", 
	       fNumbersNorm[0]->effMuidMC/fNumbersNorm[0]->effMuidPid, 
	       fNumbersNorm[1]->effMuidMC/fNumbersNorm[1]->effMuidPid)
       << Form(" r_MC = %3.2f %3.2f r_Pid = %3.2f %3.2f", 
	       fNumbersBs[0]->effMuidMC/fNumbersNorm[0]->effMuidMC, 
	       fNumbersBs[1]->effMuidMC/fNumbersNorm[1]->effMuidMC, 
	       fNumbersBs[0]->effMuidPid/fNumbersNorm[0]->effMuidPid,
	       fNumbersBs[1]->effMuidPid/fNumbersNorm[1]->effMuidPid
	       )
       << endl;

  cout << Form("pT > %3.1f trig: ", m2pt) 
       << Form("%3.2f %3.2f %3.2f %3.2f", 
	       fNumbersBs[0]->effTrigMC, fNumbersBs[0]->effTrigPid, 
	       fNumbersBs[1]->effTrigMC, fNumbersBs[1]->effTrigPid)
       << " " 
       << Form("%3.2f %3.2f %3.2f %3.2f", 
	       fNumbersNorm[0]->effTrigMC, fNumbersNorm[0]->effTrigPid, 
	       fNumbersNorm[1]->effTrigMC, fNumbersNorm[1]->effTrigPid)
       << " -> " 
       << Form(" rho_B = %3.2f rho_E = %3.2f", 
	       fNumbersNorm[0]->effTrigMC/fNumbersNorm[0]->effTrigPid, 
	       fNumbersNorm[1]->effTrigMC/fNumbersNorm[1]->effTrigPid)
       << Form(" r_MC = %3.2f %3.2f r_Pid = %3.2f %3.2f", 
	       fNumbersBs[0]->effTrigMC/fNumbersNorm[0]->effTrigMC, 
	       fNumbersBs[1]->effTrigMC/fNumbersNorm[1]->effTrigMC, 
	       fNumbersBs[0]->effTrigPid/fNumbersNorm[0]->effTrigPid,
	       fNumbersBs[1]->effTrigPid/fNumbersNorm[1]->effTrigPid
	       )
       << endl;
}


// ----------------------------------------------------------------------
TH1* anaBmm::loopTree(int mode, int proc) {
  // the loopTree
  // -- mode definition
  // 0  Bs2MuMu MC
  // 1  Bd2MuMu MC
  // 5  Bs2MuMu data
  // 10 Bp2JpsiKp MC
  // 11 Bp2JpsiKp data
  // 20 Bs2JpsiPhi MC
  // 21 Bs2JpsiPhi data

  PidTable *ptT1;
  PidTable *ptT2;
  PidTable *ptM; 
  
  bool bp2jpsikp(false), bs2jpsiphi(false), isMC(false); 

  cout << "--> loopTree with mode " << mode << " proc = " << proc << endl;

  string fAcc;

  numbers *aa(0);
  if (0 == mode) {
    isMC = true; 
    fpMc[fSgMc]->cd(); 
    fAcc = "SgMcAcc";
  }  else if (1 == mode) {
    isMC = true; 
    fpMc[fBdMc]->cd(); 
    fAcc = "BdMcAcc";
  } else if (5 == mode) {
    isMC = false;     
    fpData[fSgData]->cd(); 
  } else if (10 == mode) {
    isMC = true; 
    bp2jpsikp = true; 
    fpMc[fNoMc]->cd(); 
    fAcc = "NoMcAcc";
  } else if (11 == mode) {
    isMC = false;     
    bp2jpsikp = true; 
    fpData[fNoData]->cd(); 
  } else if (20 == mode) {
    isMC = true; 
    bs2jpsiphi = true; 
    fpMc[fCsMc]->cd(); 
    fAcc = "CsMcAcc";
  } else if (21 == mode) {
    isMC = false;     
    bs2jpsiphi = true; 
    fpData[fCsData]->cd(); 
  } else {
    cout << "mode 99" << endl;
  }

  if (isMC) {
    ptT1 = new PidTable("pidtables/110625/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 	
    ptT2 = new PidTable("pidtables/110625/H2D_L3Efficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 	
    ptM  = new PidTable("pidtables/110625/H2D_MuonIDEfficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 
//     ptT1 = new PidTable("pidtables/110625/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_mc_MCTRUTH.dat"); 	
//     ptT2 = new PidTable("pidtables/110625/H2D_L3Efficiency_GlbTM_ProbeTrackMatched_mc_MCTRUTH.dat"); 	
//     ptM  = new PidTable("pidtables/110625/H2D_MuonIDEfficiency_GlbTM_ProbeTrackMatched_mc_MCTRUTH.dat"); 
  } else {
    ptT1 = new PidTable("pidtables/110606/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 	
    ptT2 = new PidTable("pidtables/110606/H2D_L3Efficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 	
    ptM  = new PidTable("pidtables/110606/H2D_MuonIDEfficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 
  }

  gFile->pwd();

  // -- reset all histograms
  for (unsigned int i = 0; i < fNchan; ++i) {
    fhMuId[i]->Reset();
    fhMuTr[i]->Reset();

    fh0PidTrigger[i]->Reset();
    fh1PidTrigger[i]->Reset();

    fh0PidMuID[i]->Reset();
    fh1PidMuID[i]->Reset();

    fh0MCTrigger[i]->Reset();
    fh1MCTrigger[i]->Reset();

    fh0MCMuID[i]->Reset();
    fh1MCMuID[i]->Reset();

    fhMassAbsNoCuts[i]->Reset();
    fhMassNoCuts[i]->Reset();

    fhMassWithAnaCuts[i]->Reset();
    fhMassWithAnaCutsManyBins[i]->Reset();

    fhMassWithMuonCuts[i]->Reset();
    fhMassWithMuonCutsManyBins[i]->Reset();

    fhMassWithTriggerCuts[i]->Reset();
    fhMassWithTriggerCutsManyBins[i]->Reset();

    fhMassWithAllCuts[i]->Reset();
    fhMassWithAllCutsManyBins[i]->Reset();

    fhMassWithMassCutsManyBins[i]->Reset();
    fhMassWithMassCuts[i]->Reset();
  }

  // -- set up tree
  double mass(0.); 
  TTree *t;
  t = (TTree*)gFile->Get("events");
  int brun, bevt, bls, btm, bq1, bq2, bprocid; 
  double bg1pt, bg2pt, bg1eta, bg2eta;
  double bm, bcm, bpt, beta, bcosa, biso1, biso4, bchi2, bdof, bdocatrk, bfls3d, bfl3dE, bfl3d;
  double bm1pt, bm1eta, bm2pt, bm2eta, bm1phi, bm2phi;
  double bk1pt, bk1eta, bk2pt, bk2eta; 
  double bg3pt, bg3eta, bg4pt, bg4eta; 
  double bmpsi, bmkk, bdr;
  double bw8mu, bw8tr;
  bool bhlt, bgmuid, bgtqual, bjson;
  double tr1w8(0.), tr2w8(0.), trw8(0.), m1w8(0.), m2w8(0.), mw8(0.0);
  t->SetBranchAddress("run",&brun);
  t->SetBranchAddress("evt",&bevt);
  t->SetBranchAddress("hlt",&bhlt);
  t->SetBranchAddress("ls",&bls);
  t->SetBranchAddress("json",&bjson);
  t->SetBranchAddress("gmuid",&bgmuid);
  t->SetBranchAddress("gtqual",&bgtqual);
  t->SetBranchAddress("w8mu",&bw8mu);
  t->SetBranchAddress("w8tr",&bw8tr);
  t->SetBranchAddress("tm",&btm);
  t->SetBranchAddress("procid",&bprocid);
  t->SetBranchAddress("m",&bm);
  t->SetBranchAddress("cm",&bcm);
  t->SetBranchAddress("pt",&bpt);
  t->SetBranchAddress("eta",&beta);
  t->SetBranchAddress("cosa",&bcosa);
  t->SetBranchAddress("iso1",&biso1);
  t->SetBranchAddress("iso5",&biso4);
  t->SetBranchAddress("chi2",&bchi2);
  t->SetBranchAddress("dof",&bdof);
  t->SetBranchAddress("fls3d",&bfls3d);
  t->SetBranchAddress("fl3d",&bfl3d);
  t->SetBranchAddress("fl3dE",&bfl3dE);
  t->SetBranchAddress("m1pt",&bm1pt);
  t->SetBranchAddress("m1eta",&bm1eta);
  t->SetBranchAddress("m1phi",&bm1phi);
  t->SetBranchAddress("m1q",&bq1);
  t->SetBranchAddress("m2pt",&bm2pt);
  t->SetBranchAddress("m2eta",&bm2eta);
  t->SetBranchAddress("m2phi",&bm2phi);
  t->SetBranchAddress("m2q",&bq2);
  t->SetBranchAddress("docatrk",&bdocatrk);

  t->SetBranchAddress("g1pt",&bg1pt);
  t->SetBranchAddress("g2pt",&bg2pt);
  t->SetBranchAddress("g1eta",&bg1eta);
  t->SetBranchAddress("g2eta",&bg2eta);
  if (bp2jpsikp) {
    if (isMC) {
      t->SetBranchAddress("g3pt",&bg3pt);
      t->SetBranchAddress("g3eta",&bg3eta);
    }
    t->SetBranchAddress("kpt",&bk1pt);
    t->SetBranchAddress("keta",&bk1eta);
    t->SetBranchAddress("mpsi",&bmpsi);
  }

  if (bs2jpsiphi) {
    if (isMC) {
      t->SetBranchAddress("g3pt",&bg3pt);
      t->SetBranchAddress("g3eta",&bg3eta);
      t->SetBranchAddress("g4pt",&bg4pt);
      t->SetBranchAddress("g4eta",&bg4eta);
    }
    t->SetBranchAddress("mpsi",&bmpsi);
    t->SetBranchAddress("mkk",&bmkk);
    t->SetBranchAddress("dr",&bdr);
    t->SetBranchAddress("k1pt",&bk1pt);
    t->SetBranchAddress("k1eta",&bk1eta);
    t->SetBranchAddress("k2pt",&bk2pt);
    t->SetBranchAddress("k2eta",&bk2eta);
  } else {
    bmkk = 999.;
    bdr = 999.;
  }

  int nentries = Int_t(t->GetEntries());
  int nb(0); 
  cuts *pCuts(0); 
  bool cowboy(false); 
  double dphi; 
  TLorentzVector vm1, vm2, vpsi; 
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    // -- require truth matching when on MC
    if (0 == mode && 0 == btm) continue;
    if (1 == mode && 0 == btm) continue;
    if (10 == mode && 0 == btm) continue;
    if (20 == mode && 0 == btm) continue;

    if (isMC && proc > 0) {
      if (bprocid != proc) continue;
    }

    // -- channel index
    fChan = -1; 
    pCuts = 0; 
    
    fChan = detChan(bm1eta, bm2eta); 
    if (fChan > -1) {
      pCuts = fCuts[fChan]; 
    } else {
      //       cout << "event " << jentry << ", fChan = " << fOver << " for eta(B) = " << beta 
      // 	   << " m1eta = " << bm1eta << " m2eta = " << bm2eta
      // 	   << endl;
      continue;
    }

    mass = bm; 
    //     if (10 == mode || 11 == mode) {
    //       mass = bcm; 
    //     } else if (20 == mode || 21 == mode) {
    //       mass = bcm; 
    //     }
    

    fhMassAbsNoCuts[fChan]->Fill(mass);
    // -- require wide mass window
    if (mass < fMassLo) continue;
    if (fMassHi < mass) continue;

    // -- gen-level acceptance cuts
    if (isMC) {
      if (TMath::Abs(bg1eta) > 2.5) continue;
      if (TMath::Abs(bg2eta) > 2.5) continue;
      if (bg1pt < 1.0) continue;
      if (bg2pt < 1.0) continue;
      if (bp2jpsikp) {
	if (TMath::Abs(bg3eta) > 2.5) continue;
	if (bg3pt < 0.4) continue;
      }
      
      if (bs2jpsiphi) {
	if (TMath::Abs(bg3eta) > 2.5) continue;
	if (TMath::Abs(bg4eta) > 2.5) continue;
	if (bg3pt < 0.4) continue;
	if (bg4pt < 0.4) continue;
      }
    } else {
      if (!bjson) {
	//	if (5 == mode) cout << "skipping run: " << brun << " LS: " << bls << endl;
	continue;
      }
    }

    // -- immutable cuts: require basic muon and trackQual cuts
    if (false == bgtqual) continue;
    if (TMath::Abs(bm1eta) > 2.4) continue;
    if (TMath::Abs(bm2eta) > 2.4) continue;

    if (bm1pt < 1.0) continue;
    if (bm2pt < 1.0) continue;

    if (bp2jpsikp) {
      if (TMath::Abs(bk1eta) > 2.4) continue;
      if (bk1pt < 0.5) continue;
    }
    if (bs2jpsiphi) {
      if (TMath::Abs(bk1eta) > 2.4) continue;
      if (TMath::Abs(bk2eta) > 2.4) continue;
      if (bk1pt < 0.5) continue;
      if (bk2pt < 0.5) continue;
    }

    // -- this is the base, after the raw acceptance cuts
    fhMassNoCuts[fChan]->Fill(mass);
    
    // -- analysis cuts
    if (bq1*bq2 > 0) continue;
    if (bm1pt < pCuts->m1pt) continue; 
    if (bm2pt < pCuts->m2pt) continue; 

    if (bfl3d > 2) continue;

    if (bpt < pCuts->pt) continue; 
    if (biso4 < pCuts->iso1) continue; 
    if (bchi2/bdof > pCuts->chi2dof) continue;
    if (TMath::IsNaN(bfls3d)) continue;
    if (bfls3d < pCuts->fls3d) continue;
    if (TMath::ACos(bcosa) > pCuts->alpha) continue;
    if (bdocatrk < pCuts->docatrk) continue;

    if (bs2jpsiphi && bdr >0.3) continue;
    if (bs2jpsiphi && bmkk < 0.995) continue;
    if (bs2jpsiphi && bmkk > 1.045) continue;

    if (bs2jpsiphi || bp2jpsikp) {
      cowboy = false; 
      vm1.SetPtEtaPhiM(bm1pt, bm1eta, bm1phi, MMUON);
      vm2.SetPtEtaPhiM(bm2pt, bm2eta, bm2phi, MMUON);
      vpsi = vm1 + vm2; 
      if (bmpsi > 3.2) continue;
      if (bmpsi < 3.0) continue;
      dphi = vm1.DeltaPhi(vm2); 
      cowboy = bq1*dphi > 0; 
      // -- cowboy veto 
      if (cowboy) {
	//cout << vpsi.Perp() << " and " << bq1*dphi << " killed by cowboy veto" << endl;
 	continue;
      }
      if (vpsi.Perp() < 7) {
	//cout << vpsi.Perp() << " and " << bq1*dphi << " killed by psi pT cut" << endl;
	continue;
      }
    }

    //     if (bs2jpsiphi && cowboy) {
    //       cout << "q1 = " << bq1 << " phi1: " << bm1phi << " phi2: " << bm2phi << " -> " << dphi << "==> " << cowboy << endl;
    //     }

    fhMassWithAnaCuts[fChan]->Fill(mass); 
    fhMassWithAnaCutsManyBins[fChan]->Fill(mass); 

    // -- now check for muon ID and trigger
    if (false == bgmuid) continue;

    m1w8 = ptM->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    m2w8 = ptM->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    mw8  = m1w8*m2w8; 
    
    //    cout << "mw8 = " << mw8 << ": " << m1w8 << "*" << m2w8 
    //         << " for muons = " << bm1pt << " " << bm1eta << " " << bm2pt << " " << bm2eta << endl;

    if (mw8 > 0.) {
      fhMuId[fChan]->Fill(mw8, 1./mw8); 
    } else {
      cout << "mw8 = " << mw8 << " for muons = " << bm1pt << " " << bm1eta << " " << bm2pt << " " << bm2eta << endl;
    }

    fh0PidMuID[fChan]->Fill(bpt, mw8); 
    fh1PidMuID[fChan]->Fill(bpt); 

    fhMassWithMuonCuts[fChan]->Fill(mass); 
    fhMassWithMuonCutsManyBins[fChan]->Fill(mass); 

    fh1MCTrigger[fChan]->Fill(bpt); 

    if (false == bhlt) continue;
    fhMassWithTriggerCuts[fChan]->Fill(mass); 
    fhMassWithTriggerCutsManyBins[fChan]->Fill(mass); 

    // -- weights for trigger
    tr1w8 = ptT1->effD(bm1pt, TMath::Abs(bm1eta), 0.)*ptT2->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    tr2w8 = ptT1->effD(bm2pt, TMath::Abs(bm2eta), 0.)*ptT2->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    trw8  = tr1w8*tr2w8; 
    if (trw8 > 0.) {
      fhMuTr[fChan]->Fill(trw8, 1./trw8); 
    }

    fh0PidTrigger[fChan]->Fill(bpt, trw8); 
    fh1PidTrigger[fChan]->Fill(bpt); 

    fh0MCTrigger[fChan]->Fill(bpt); 

    fhMassWithAllCuts[fChan]->Fill(mass); 
    fhMassWithAllCutsManyBins[fChan]->Fill(mass); 
    
    if (0 == mode && mass < pCuts->mBsLo) continue;
    if (0 == mode && mass > pCuts->mBsHi) continue;
    if (1 == mode && mass < pCuts->mBdLo) continue;
    if (1 == mode && mass > pCuts->mBdHi) continue;
    if (10 == mode && mass < fNormLo) continue;
    if (10 == mode && mass > fNormHi) continue;
    if (20 == mode && mass < fCsLo) continue;
    if (20 == mode && mass > fCsHi) continue;

    fhMassWithMassCuts[fChan]->Fill(mass);
    fhMassWithMassCutsManyBins[fChan]->Fill(mass); 

    if (fDoPrint && 5 == mode && mass > 4.5 && mass < 6.0) {
      cout << Form("m = %4.3f pT = %4.3f eta = %4.3f", mass, bpt, beta)
	   <<	" r = " << brun << "/" << bevt
	   << " chan = " << fChan 
	   << Form(" mpt = %4.3f,%4.3f", bm1pt, bm2pt)
	   << Form(" meta = %4.3f,%4.3f", TMath::Abs(bm1eta), TMath::Abs(bm2eta))
	   << Form(" a = %4.3f iso = %4.3f chi2 = %4.3f fls3d = %4.3f, fl/E=%4.3f/%4.3f", 
		   TMath::ACos(bcosa), biso4, bchi2/bdof, bfls3d, bfl3d, bfl3dE)
	   << endl;
      fOUT << Form("m = %4.3f pT = %4.3f eta = %4.3f", mass, bpt, beta)
	//	   <<	" run = " << brun << " event = " << bevt
	   << " chan = " << fChan 
	   << Form(" mpt = %4.3f,%4.3f", bm1pt, bm2pt)
	   << Form(" meta = %4.3f,%4.3f", TMath::Abs(bm1eta), TMath::Abs(bm2eta))
	   << Form(" a = %4.3f iso = %4.3f chi2 = %4.3f fls3d = %4.3f, fl/E=%4.3f/%4.3f", 
		   TMath::ACos(bcosa), biso4, bchi2/bdof, bfls3d, bfl3d, bfl3dE)
	   << endl;
    }
    
  }

  if (99 == mode) return fhMassWithAllCuts[0];

  for (unsigned int i = 0; i < fNchan; ++i) {
    pCuts = fCuts[i]; 
    aa = 0; 
    if (0 == mode || 5 == mode) {
      aa = fNumbersBs[i];
      aa->mBdHi = pCuts->mBdHi;
      aa->mBdLo = pCuts->mBdLo;
      aa->mBsHi = pCuts->mBsHi;
      aa->mBsLo = pCuts->mBsLo;
    }
    // -- Bd2MuMu only for the MC numbers!
    if (1 == mode) {
      aa = fNumbersBd[i];
      aa->mBdHi = pCuts->mBdHi;
      aa->mBdLo = pCuts->mBdLo;
      aa->mBsHi = pCuts->mBsHi;
      aa->mBsLo = pCuts->mBsLo;
    }

    if (10 == mode || 11 == mode) {
      aa = fNumbersNorm[i];
    }

    if (20 == mode || 21 == mode) {
      aa = fNumbersCS[i];
    }


    if (0 == aa) { 
      cout << "++++++++++++????????????++++++++++ aa unset?????????????" << endl;
      cout << "mode = " << mode << endl;
      cout << "i    = " << i << endl;
      cout << "++++++++++++????????????++++++++++ aa unset?????????????" << endl;

    }


    // -- signal cross feed 
    double tot   = fhMassWithAllCutsManyBins[i]->GetSumOfWeights();
    double bd    = fhMassWithAllCutsManyBins[i]->Integral(fhMassWithAllCutsManyBins[i]->FindBin(pCuts->mBdLo), 
						       fhMassWithAllCutsManyBins[i]->FindBin(pCuts->mBdHi));
    double bs    = fhMassWithAllCutsManyBins[i]->Integral(fhMassWithAllCutsManyBins[i]->FindBin(pCuts->mBsLo), 
						       fhMassWithAllCutsManyBins[i]->FindBin(pCuts->mBsHi));
    
    fhMuId[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/muid-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMuTr[i]->Draw();
    if (fDoPrint)   c0->SaveAs(Form("%s/mutr-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassAbsNoCuts[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/anc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassNoCuts[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/noc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithAllCutsManyBins[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/wac-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithMuonCutsManyBins[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/muc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithTriggerCutsManyBins[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/trc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithMassCutsManyBins[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/wmc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    // -- Efficiency and acceptance
    if (isMC) {
      accEffFromEffTree(fAcc, *aa, *fCuts[i], proc);
      //      accEffFromEffTree(gFile, *aa, *fCuts[i], 0, 1, proc);
      //      double a = fhMassAbsNoCuts[i]->GetSumOfWeights(); 
      double a = fhMassNoCuts[i]->GetSumOfWeights(); 
      double b = fhMassWithAnaCuts[i]->GetSumOfWeights();
      double c = fhMassWithMuonCuts[i]->GetSumOfWeights();
      double d = fhMassWithTriggerCuts[i]->GetSumOfWeights();
      double e = fhMassWithAllCuts[i]->GetSumOfWeights();
      double f = fhMassWithMassCuts[i]->GetSumOfWeights();
      //WRONG if effTree is from different file than events
      //aa->effCand          = a/aa->recoYield;
      //aa->effCandE         = dEff(static_cast<int>(a), static_cast<int>(aa->recoYield));
      aa->ana0Yield        = a;
      aa->ana0YieldE       = TMath::Sqrt(a);
      aa->anaYield         = b; 
      aa->anaYieldE        = TMath::Sqrt(b); 
      aa->anaMuonYield     = c; 
      aa->anaMuonYieldE    = TMath::Sqrt(c); 
      aa->anaTriggerYield  = d; 
      aa->anaTriggerYieldE = TMath::Sqrt(d); 
      aa->anaWmcYield      = f; 
      aa->anaWmcYieldE     = TMath::Sqrt(f);
      aa->effAna           = b/a;
      aa->effAnaE          = dEff(static_cast<int>(b), static_cast<int>(a));
      aa->effMuidMC        = c/b;
      aa->effMuidMCE       = dEff(static_cast<int>(c), static_cast<int>(b));
      aa->effMuidPid       = fhMuId[i]->GetMean();
      aa->effMuidPidE      = fhMuId[i]->GetMeanError();
      aa->effTrigMC        = d/c;
      aa->effTrigMCE       = dEff(static_cast<int>(d), static_cast<int>(c));
      aa->effTrigPid       = fhMuTr[i]->GetMean();
      aa->effTrigPidE      = fhMuTr[i]->GetMeanError();
      aa->effTot           = e/(aa->genYield);
      aa->effTotE          = dEff(static_cast<int>(e), static_cast<int>(aa->genYield));
      aa->effTotChan       = e/(aa->genChanYield);
      aa->effTotChanE      = dEff(static_cast<int>(e), static_cast<int>(aa->genChanYield));
      aa->effProdMC        = aa->effCand * aa->effAna * aa->effMuidMC * aa->effTrigMC;
      aa->effProdMCE       = 0.;
      aa->effProdPid       = aa->effCand * aa->effAna * aa->effMuidPid * aa->effTrigPid;
      aa->effProdPidE      = 0.;
      aa->aEffProdMC       = aa->effProdMC * aa->accChan * aa->cFrac;
      aa->aEffProdMCE      = 0.;

      aa->combGenYield     = e/(aa->acc * aa->effProdMC);
      aa->chanGenYield     = e/(aa->accChan * aa->effProdMC); 
      aa->prodGenYield     = e/(aa->effTot); 
    }

    if (fDoPrint) {
      cout << "bs:  " << bs << endl;
      cout << "bd:  " << bd << endl;
      cout << "tot: " << tot << endl;
      cout << "wmc: " << aa->anaWmcYield << endl;
    }

    // -- compute PSD, PDS, etc
    if (0 == mode) {
      aa->pss = bs/tot;
      aa->pssE = dEff(static_cast<int>(bs), static_cast<int>(tot)); 
      aa->pds = bd/tot;
      aa->pdsE = dEff(static_cast<int>(bd), static_cast<int>(tot)); 
    } 
    
    if (1 == mode) {
      aa->pdd = bd/tot;
      aa->pddE = dEff(static_cast<int>(bd), static_cast<int>(tot)); 
      aa->psd = bs/tot;
      aa->psdE = dEff(static_cast<int>(bs), static_cast<int>(tot)); 
    } 

    if (0 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC SIGNAL, channel " << i << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      if (fDoPrint)      c0->SaveAs(Form("%s/sig-mc-chan%d.pdf", fDirectory.c_str(), i));
      cout << "----> "; gFile->pwd(); 
    } else if (5 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA SIGNAL, channel " << i  << endl;
      bgBlind(fhMassWithAllCuts[i], 1, fBgLo, fBgHi);
      cout << "fBgExp = " << fBgExp << "+/-" << fBgExpE << endl;
      cout << "fBgHist = " << fBgHist << "+/-" << fBgHistE << endl;
      aa->bgObs = fBgHist;
      double scaleBs = (aa->mBsHi-aa->mBsLo)/(fBgHi-fBgLo);
      double scaleBd = (aa->mBdHi-aa->mBdLo)/(fBgHi-fBgLo);
      aa->bgBsExp  = scaleBs*aa->bgObs;
      aa->bgBsExpE = scaleBs*TMath::Sqrt(aa->bgObs);
      aa->bgBdExp  = scaleBd*aa->bgObs;
      aa->bgBdExpE = scaleBd*TMath::Sqrt(aa->bgObs);
      TH1D *h = fhMassWithAllCuts[i]; 
      setHist(h, kBlack, 20, 1.); 
      setTitles(h, "m_{#mu#mu} [GeV]", "Candidates/Bin", 0.05, 1.2, 1.3); 
      h->SetMinimum(0.01); 
      gStyle->SetOptStat(0); 
      gStyle->SetOptTitle(0); 
      h->Draw();
      if (fDoPrint)  c0->SaveAs(Form("%s/sig-data-chan%d.pdf", fDirectory.c_str(), i));
    } else if (10 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC NORMALIZATION, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      if (fDoPrint) c0->SaveAs(Form("%s/norm-mc-chan%d.pdf", fDirectory.c_str(), i));
    } else if (11 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA NORMALIZATION, channel " << i  << endl;
      //      TH1D *h = fhMassWithAllCutsManyBins[i]; 
      TH1D *h = fhMassWithAllCuts[i]; 
      normYield(h, mode, 5.10, 5.5);
      aa->fitYield  = fNormSig; 
      aa->fitYieldE = fNormSigE; 
      setHist(h, kBlack, 20, 1.); 
      setTitles(h, "m_{#mu#muK} [GeV]", "Candidates/Bin", 0.05, 1.2, 1.6); 
      h->SetMinimum(0.01); 
      gStyle->SetOptStat(0); 
      gStyle->SetOptTitle(0); 
      h->Draw();
      if (fDoPrint) c0->SaveAs(Form("%s/norm-data-chan%d.pdf", fDirectory.c_str(), i));
    } else if (20 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC CONTROL SAMPLE, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      if (fDoPrint) c0->SaveAs(Form("%s/cs-mc-chan%d.pdf", fDirectory.c_str(), i));
    } else if (21 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA CONTROL SAMPLE, channel " << i  << endl;
      TH1D *h = fhMassWithAllCuts[i]; 
      csYield(h, mode, 5.20, 5.6);
      aa->fitYield  = fCsSig; 
      aa->fitYieldE = fCsSigE; 
      setHist(h, kBlack, 20, 1.); 
      setTitles(h, "m_{#mu#muKK} [GeV]", "Candidates/Bin", 0.05, 1.2, 1.6); 
      h->SetMinimum(0.01); 
      gStyle->SetOptStat(0); 
      gStyle->SetOptTitle(0); 
      h->Draw();
      if (fDoPrint) c0->SaveAs(Form("%s/cs-data-chan%d.pdf", fDirectory.c_str(), i));
    } 

    printNumbers(*aa, cout); 
    printNumbers(*aa, fOUT); 
    
    // -- Cache the pwd...
    TDirectory *pD = gFile; 
    fHistFile->cd();
    cout << "==> " << i << "  " << mode << " hMassWithMassCuts[i] = " << fhMassWithMassCuts[i] << endl;
    fhMassWithMassCutsManyBins[i]->SetName(Form("hMassWithMassCutsManyBins%d_chan%d", mode, i)); fhMassWithMassCutsManyBins[i]->Write();
    fhMassWithAllCutsManyBins[i]->SetName(Form("hMassWithAllCutsManyBins%d_chan%d", mode, i)); fhMassWithAllCutsManyBins[i]->Write();
    fhMassWithMassCuts[i]->SetName(Form("hMassWithMassCuts%d_chan%d", mode, i)); fhMassWithMassCuts[i]->Write();
    fhMassWithAllCuts[i]->SetName(Form("hMassWithAllCuts%d_chan%d", mode, i)); fhMassWithAllCuts[i]->Write();
    fhMassNoCuts[i]->SetName(Form("hMassNoCuts%d_chan%d", mode, i)); fhMassNoCuts[i]->Write();
    fhMassAbsNoCuts[i]->SetName(Form("hMassAbsNoCuts%d_chan%d", mode, i)); fhMassAbsNoCuts[i]->Write();
    fhMuTr[i]->SetName(Form("hMuTr%d_chan%d", mode, i)); fhMuTr[i]->Write();
    fhMuId[i]->SetName(Form("hMuId%d_chan%d", mode, i)); fhMuId[i]->Write();
    // -- and get back to it
    pD->cd();
  }

  delete ptT1;
  delete ptT2;
  delete ptM;

  return fhMassWithAllCuts[0];
}


// ----------------------------------------------------------------------
void anaBmm::triggerSignal(const char *cuts) {

  static int version(0); 

  string defaultCuts = "4.8<m&&m<6&&gmuid&&gtqual&&m1pt>3&&m2pt>3&&m1q*m2q<0&&abs(eta)<2.4&&abs(m2eta)<2.4&&pt>5&&abs(m1ip)<2&&abs(m2ip)<2";

  string allCuts = defaultCuts + "&&" + cuts;
  string hltCuts = defaultCuts + "&&" + cuts + "&&hlt";

  int NBINS(15); 
  double HLO(0.), HHI(45.); 
  TH1D *M0 = new TH1D("M0", "", 40, 4.8, 6.0); 
  TH1D *H0 = new TH1D("H0", "", NBINS, HLO, HHI); 
  TH1D *H1 = new TH1D("H1", "", NBINS, HLO, HHI); 
  TH1D *EF = new TH1D("EF", "", NBINS, HLO, HHI); EF->Sumw2(); 
  
  const int NFILE(3);
  string names[4]; 
  TFile *f1[NFILE]; 
  f1[0] = TFile::Open("/shome/ursl/scratch/bmm/stuff/cv02-1313-bmt-ht.root");       names[0] = "HT";
  f1[1] = TFile::Open("/shome/ursl/scratch/bmm/stuff/cv02-1313-bmt-photon.root");   names[1] = "Photon";
  f1[2] = TFile::Open("/shome/ursl/scratch/bmm/stuff/cv02-1313-bmt-jet.root");      names[2] = "Jet"; 
  //  f1[3] = TFile::Open("/shome/ursl/scratch/bmm/stuff/cv02-1313-bmt-multiJet.root"); names[3] = "MultiJet"; 

  TTree *t; 

  loopTree(0); 

  c0->Clear(); 
  c0->Divide(2,2); 
  //  c0->Divide(2,2); 
  gStyle->SetOptStat(0); 

  double vEff[NFILE], vEffE[NFILE]; 
  
  for (int i = 0; i < NFILE; ++i) {
    f1[i]->cd();
    TH1D *m0 = new TH1D("m0", "", 40, 4.8, 6.0); 
    TH1D *h0 = new TH1D("h0", "", NBINS, HLO, HHI); 
    TH1D *h1 = new TH1D("h1", "", NBINS, HLO, HHI); 
    t = (TTree*)f1[i]->Get("events"); 
    
    t->Draw("m>>m0", allCuts.c_str(), "goff");
    double n0   = t->Draw("pt>>h0", allCuts.c_str(), "goff");
    double n1   = t->Draw("pt>>h1", hltCuts.c_str(), "goff");
    double eff  = n1/n0; 
    double deff = dEff(static_cast<int>(n1), static_cast<int>(n0)); 

    M0->Add(m0); 
    H0->Add(h0); 
    H1->Add(h1); 
    
    vEff[i]  = eff;
    vEffE[i] = deff;

    c0->cd(i+1); 
    setFilledHist(h0, kBlack, kRed, 3004); 
    setTitles(h0, "p_{T} [GeV]", "Entries/bin"); 
    h0->Draw();

    setFilledHist(h1, kBlack, kBlue, 3005); 
    h1->Draw("same");

    tl->SetTextSize(0.04); tl->DrawLatex(0.2, 0.92, names[i].c_str()); 
    tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.8, Form("eff = %3.0f/%3.0f", n1, n0)); 
    tl->SetTextSize(0.04); tl->DrawLatex(0.55,0.72, Form("= %4.3f #pm %4.3f", eff, deff)); 
    
  }

  double ave(0.), aveE(0.); 
  average(ave, aveE, 3, vEff, vEffE); 

  if (fDoPrint) c0->SaveAs(Form("%s/triggerSignal-pt-%d.pdf", fDirectory.c_str(), version));     

  c0->Clear();
  c0->Divide(1); 

  EF->Divide(H1, H0, 1., 1., "b"); 
  EF->SetMinimum(0.); EF->SetMaximum(1.2);
  setTitles(EF, "p_{T} [GeV]", "Efficiency"); 
  EF->Draw(); 
  tl->SetTextSize(0.03);  tl->DrawLatex(0.2, 0.92, cuts); 
  tl->SetTextSize(0.04); tl->DrawLatex(0.5,0.85, Form("<#epsilon> = %4.3f #pm %4.3f", ave, aveE)); 
  if (fDoPrint) c0->SaveAs(Form("%s/triggerSignal-eff-pt-%d.pdf", fDirectory.c_str(), version));     


  TH1D *h1 = (TH1D*)fh0PidTrigger[0]->Clone();
  h1->SetName("heff"); 
  h1->Divide(fh0PidTrigger[0], fh1PidTrigger[0], 1., 1., "b");
  setHist(h1, kRed, 25, 1.); 
  h1->Draw("samee");

  TH1D *h2 = (TH1D*)fh0MCTrigger[0]->Clone();
  h2->SetName("mceff"); 
  h2->Divide(fh0MCTrigger[0], fh1MCTrigger[0], 1., 1., "b");
  setHist(h2, kRed, 25, 1.); 
  h2->Draw("samehist");

  cout << Form("version: %3d eff = %4.3f+/-%4.3f", version, ave, aveE) << endl;

  EF->Draw();


  h1->Draw("samee");
  h2->Draw("samehist");

  newLegend(0.22, 0.65, 0.30, 0.89); 
  legg->SetTextSize(0.045);  
  legg->AddEntry(EF, "non-muon PD", "p"); 
  legg->AddEntry(h1, "TNP", "p"); 
  legg->AddEntry(h2, "MC", "l"); 
  legg->Draw();


  if (fDoPrint) c0->SaveAs(Form("%s/triggerSignal-eff-pt-overlay-%d.pdf", fDirectory.c_str(), version));     
 
  ++version;

}


// ----------------------------------------------------------------------
void anaBmm::triggerNorm(const char *cuts) {

  static int version(0); 

  string defaultCuts = "4.8<m&&m<6&&gmuid&&gtqual&&m1pt>3&&m2pt>3&&m1q*m2q<0&&abs(eta)<2.4&&abs(m2eta)<2.4&&pt>5&&abs(m1ip)<2&&abs(m2ip)<2";

  string allCuts = defaultCuts + "&&" + cuts;
  string hltCuts = defaultCuts + "&&" + cuts + "&&hlt";

  int NBINS(15); 
  double HLO(0.), HHI(45.); 
  TH1D *M0 = new TH1D("M0", "", 40, 4.8, 6.0); 
  TH1D *H0 = new TH1D("H0", "", NBINS, HLO, HHI); 
  TH1D *H1 = new TH1D("H1", "", NBINS, HLO, HHI); 
  TH1D *EF = new TH1D("EF", "", NBINS, HLO, HHI); EF->Sumw2(); 
  
  const int NFILE(4); 
  string names[4]; 
  TFile *f1[NFILE]; 
  f1[0] = TFile::Open("/shome/ursl/scratch/bmm/stuff/cv01-300521-bmt-ht.root");       names[0] = "HT";	   
  f1[1] = TFile::Open("/shome/ursl/scratch/bmm/stuff/cv01-300521-bmt-photon.root");   names[1] = "Photon";   
  f1[2] = TFile::Open("/shome/ursl/scratch/bmm/stuff/cv01-300521-bmt-jet.root");      names[2] = "Jet"; 	   
  //  f1[3] = TFile::Open("/shome/ursl/scratch/bmm/stuff/cv01-300521-bmt-multiJet.root"); names[3] = "MultiJet"; 

  TTree *t; 

  loopTree(10); 

  c0->Clear(); 
  c0->Divide(2,2); 
  gStyle->SetOptStat(0); 

  double vEff[NFILE], vEffE[NFILE]; 
  
  for (int i = 0; i < NFILE; ++i) {
    f1[i]->cd();
    TH1D *m0 = new TH1D("m0", "", 40, 4.8, 6.0); 
    TH1D *h0 = new TH1D("h0", "", NBINS, HLO, HHI); 
    TH1D *h1 = new TH1D("h1", "", NBINS, HLO, HHI); 
    t = (TTree*)f1[i]->Get("events"); 
    
    t->Draw("m>>m0", allCuts.c_str(), "goff");
    double n0   = t->Draw("pt>>h0", allCuts.c_str(), "goff");
    double n1   = t->Draw("pt>>h1", hltCuts.c_str(), "goff");
    double eff  = n1/n0; 
    double deff = dEff(static_cast<int>(n1), static_cast<int>(n0)); 

    M0->Add(m0); 
    H0->Add(h0); 
    H1->Add(h1); 
    
    vEff[i]  = eff;
    vEffE[i] = deff;

    c0->cd(i+1); 
    setFilledHist(h0, kBlack, kRed, 3004); 
    setTitles(h0, "p_{T} [GeV]", "Entries/bin"); 
    h0->Draw();

    setFilledHist(h1, kBlack, kBlue, 3005); 
    h1->Draw("same");

    tl->SetTextSize(0.03); tl->DrawLatex(0.2, 0.92, gFile->GetName()); 
    tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.8, Form("eff = %3.0f/%3.0f", n1, n0)); 
    tl->SetTextSize(0.04); tl->DrawLatex(0.55,0.72, Form("= %4.3f #pm %4.3f", eff, deff)); 
    
  }

  double ave(0.), aveE(0.); 
  average(ave, aveE, NFILE, vEff, vEffE); 


  if (fDoPrint) c0->SaveAs(Form("%s/triggerNorm-pt-%d.pdf", fDirectory.c_str(), version));     

  c0->Clear();
  c0->Divide(1); 
  M0->Draw();  
  if (fDoPrint) c0->SaveAs(Form("%s/triggerNorm-mass-%d.pdf", fDirectory.c_str(), version));     

  c0->cd(NFILE+2); 
  EF->Divide(H1, H0, 1., 1., "b"); 
  EF->SetMinimum(0.); EF->SetMaximum(1.2);
  setTitles(EF, "p_{T} [GeV]", "Efficiency"); 
  EF->Draw(); 
  tl->SetTextSize(0.03);  tl->DrawLatex(0.2, 0.92, cuts); 
  tl->SetTextSize(0.04); tl->DrawLatex(0.5,0.85, Form("<#epsilon> = %4.3f #pm %4.3f", ave, aveE)); 

  cout << Form("version: %3d eff = %4.3f+/-%4.3f", version, ave, aveE) << endl;
  if (fDoPrint) c0->SaveAs(Form("%s/triggerNorm-eff-pt-%d.pdf", fDirectory.c_str(), version));     


  TH1D *h1 = (TH1D*)fh0PidTrigger[0]->Clone();
  h1->SetName("heff"); 
  h1->Divide(fh0PidTrigger[0], fh1PidTrigger[0], 1., 1., "b");
  setHist(h1, kRed, 25, 1.); 
  h1->Draw("samee");


  TH1D *h2 = (TH1D*)fh0MCTrigger[0]->Clone();
  h2->SetName("mceff"); 
  h2->Divide(fh0MCTrigger[0], fh1MCTrigger[0], 1., 1., "b");
  setHist(h2, kRed, 25, 1.); 
  h2->Draw("samehist");



  EF->Draw();
  h1->Draw("samee");
  h2->Draw("samehist");

  newLegend(0.22, 0.65, 0.30, 0.89); 
  legg->SetTextSize(0.045);  
  legg->AddEntry(EF, "non-muon PD", "p"); 
  legg->AddEntry(h1, "TNP", "p"); 
  legg->AddEntry(h2, "MC", "l"); 
  legg->Draw();
  

  if (fDoPrint)  c0->SaveAs(Form("%s/triggerNorm-eff-pt-overlay-%d.pdf", fDirectory.c_str(), version));     

  ++version;

}


// ----------------------------------------------------------------------
void anaBmm::testUL(int ichan) {
  static int version(0); 
  ++version; 
  readCuts("anaBmm.opt.cuts"); 
  ofstream OUT(Form("testUL-v%d", version)); 
  printCuts(OUT); 

  cout << "--> loopTree: signal MC" << endl;
  loopTree(0);  // signal eff
  c0->Modified(); c0->Update();
  loopTree(1);  // Bd2MuMu eff
  c0->Modified(); c0->Update();
  cout << "--> loopTree: signal data" << endl;
  loopTree(5);  // data signal
  c0->Modified(); c0->Update();

  fhMassWithAllCuts[ichan]->Draw();

  // -- single channel version
  double scale = fDataLumi[fSgData]/39.4;
  double nbs = 2.0e9*(1.0-0.12)*scale;
  double siglo = fNumbersBs[ichan]->mBsLo; 
  double sighi = fNumbersBs[ichan]->mBsHi; 
  double mFactor = (sighi-siglo)/(fBgHi-fBgLo-(5.45-5.20)); 
  fBgHistExp  = fNumbersBs[ichan]->bgObs*mFactor;
  
  fBgExp = fBgHistExp;
  fBgExpE = 0.2*fBgExp;
  fNobs = static_cast<int>(fBgExp + 0.5);
  double nulbayes  = blimit(0.9, fNobs, 1.0, 0.2, fBgExp, fBgExpE, 1);
  double effTot = fNumbersBs[ichan]->effTot*fNumbersBs[0]->pss;
  double ulbayes  = nulbayes/(effTot*nbs);
  cout << "==> Chan:      " << ichan << endl;
  cout << "==> effTot:    " << effTot << endl;
  cout << "==> nObs:      " << fNumbersBs[ichan]->bgObs << endl;
  cout << "==> nExp:      " << fBgHistExp << endl;
  cout << "==> ul(Bayes): " << ulbayes << endl;

  OUT << "==> chan   = " << ichan << endl;
  OUT << "==> effTot = " << effTot << endl;
  OUT << "==> nObs   = " << fNumbersBs[ichan]->bgObs << endl;
  OUT << "==> nExp   = " << fBgHistExp << endl;
  OUT << "==> UL     = " << ulbayes << endl;


//   fUlcalcFileName = Form("testUL-v%d.ulc", version);
//   system(Form("/bin/rm -f %s", fUlcalcFileName.c_str()));
//   printUlcalcNumbers();
//   system("/bin/rm -f bla");
//   system(Form("../ulcalc/bin/ulcalc %s > bla", fUlcalcFileName.c_str())); 
//   system(Form("tail -1 bla >> %s", Form("testUL-v%d", version)));

}


// ----------------------------------------------------------------------
void anaBmm::optimizeULs(int nruns, int seed) {
  int version(-1); 
  fDoPrint = false; 
  ofstream OUT(Form("optimizeUL-%d.txt", seed)); 

  int NCUTS(9);
  //string cuts[] = {"m2pt", "m1pt","pt", "alpha", "chi2", "fls3d", "docatrk", "iso", "mwindow"}; 
  double loCuts[] = {3.0,    3.0,   4.0,  0.02,    0.8,     8,       0.0,       0.7,  0.020 };
  double hiCuts[] = {4.0,    7.0,   10.,  0.07,    2.2,     20,      0.1,       0.9,  0.100};

  if (seed > 0) {
    cout << "Setting random number seed " << seed << endl;
    OUT  << "Setting random number seed " << seed << endl;
    gRandom->SetSeed(seed);
  }

  
  string cutline; 
  double cut; 
  for (int j = 0; j < nruns; ++j) {
    ++version; 
    for (int i = 0; i < NCUTS; ++i) {
      cut = gRandom->Rndm()*(hiCuts[i]-loCuts[i]) + loCuts[i];
      for (int ic = 0; ic < 2; ++ic) {
	if (0 == i) fCuts[ic]->m2pt = cut; 
	if (1 == i) fCuts[ic]->m1pt = cut; 
	if (2 == i) fCuts[ic]->pt = cut; 
	if (3 == i) fCuts[ic]->alpha = cut; 
	if (4 == i) fCuts[ic]->chi2dof = cut; 
	if (5 == i) fCuts[ic]->fls3d = cut; 
	if (6 == i) fCuts[ic]->docatrk = cut; 
	if (7 == i) fCuts[ic]->iso1 = cut;
	if (8 == i) fCuts[ic]->mBsLo = 5.370 - cut; 
	if (8 == i) fCuts[ic]->mBsHi = 5.370 + cut; 
      }
    }
    printCuts(OUT); 

    cout << "--> loopTree: signal MC" << endl;
    loopTree(0);  // signal eff
    c0->Modified(); c0->Update();
    loopTree(1);  // Bd2MuMu eff
    c0->Modified(); c0->Update();
    cout << "--> loopTree: signal data" << endl;
    loopTree(5);  // data signal
    c0->Modified(); c0->Update();
    
    
    // -- simple blimits in two channels
    for (int ichan = 0; ichan < 2; ++ichan) {
	fhMassWithAllCuts[ichan]->Draw();
	double scale   = fDataLumi[fSgData]/39.4;
	double nbs     = 2.0e9*(1.0-0.12)*scale;
	double siglo   = fNumbersBs[ichan]->mBsLo; 
	double sighi   = fNumbersBs[ichan]->mBsHi; 
	double mFactor = (sighi-siglo)/(fBgHi-fBgLo-(5.45-5.20)); 
	double lscale  = fDataLumi[fSgData]/fLumi["SgMc"];
    
	double signal  = fNumbersBs[ichan]->anaWmcYield*lscale;
	fBgHistExp     = fNumbersBs[ichan]->bgObs*mFactor;
	
	fBgExp = fBgHistExp;
	fBgExpE = 0.2*fBgExp;
	fNobs = static_cast<int>(fBgExp + 0.5);
	double nulbayes  = blimit(0.9, fNobs, 1.0, 0.2, fBgExp, fBgExpE, 1);
	double effTot = fNumbersBs[ichan]->effTot*fNumbersBs[0]->pss;
	double ulbayes  = nulbayes/(effTot*nbs);
	cout << "==> effTot:     " << effTot << " chan=" << ichan << " version=" << version << endl;
	cout << "==> nObs:       " << fNumbersBs[ichan]->bgObs << " chan=" << ichan << " version=" << version  << endl;
	cout << "==> nExp:       " << fBgHistExp << " chan=" << ichan << " version=" << version << endl;
	cout << "==> ul(Bayes):  " << ulbayes << " chan=" << ichan << " version=" << version << endl;
	cout << "==> Signal:     " << signal << " chan=" << ichan << " version=" << version << endl;
	cout << "==> SSB:        " << signal/TMath::Sqrt(signal+fBgHistExp) << " chan=" << ichan << " version=" << version << endl;
	
	OUT << "==> effTot = " << effTot << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> nObs   = " << fNumbersBs[ichan]->bgObs << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> nExp   = " << fBgHistExp << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> UL     = " << ulbayes << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> Signal = " << signal << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> SSB    = " << signal/TMath::Sqrt(signal+fBgHistExp) << " chan=" << ichan << " version=" << version << endl;
      }
  }

}


struct bla{
  double ul, nobs, nexp, eff; 
  double mlo, mhi; 
  double m1pt, m2pt, pt; 
  double chi2dof, iso1, alpha, fls3d, docatrk; 
};

// ----------------------------------------------------------------------
void anaBmm::bestUL(const char *fname, int ichan, int nsettings) {

  TFile *f = TFile::Open(fname); 
  
  TTree *t = (TTree*)f->Get("t");

  int chan, file, run; 
  float mlo, mhi, pt, m1pt, m2pt, iso1, chi2dof, alpha, fls3d, docatrk; 
  float ul, nobs, nexp, eff; 

  t->SetBranchAddress("chan", &chan);
  t->SetBranchAddress("file", &file);
  t->SetBranchAddress("run", &run);
  t->SetBranchAddress("ul", &ul);
  t->SetBranchAddress("nobs", &nobs);
  t->SetBranchAddress("nexp", &nexp);
  t->SetBranchAddress("eff", &eff);

  t->SetBranchAddress("mlo", &mlo);
  t->SetBranchAddress("mhi", &mhi);
  t->SetBranchAddress("pt", &pt);
  t->SetBranchAddress("m1pt", &m1pt);
  t->SetBranchAddress("m2pt", &m2pt);
  t->SetBranchAddress("iso1", &iso1);
  t->SetBranchAddress("chi2dof", &chi2dof);
  t->SetBranchAddress("alpha", &alpha);
  t->SetBranchAddress("fls3d", &fls3d);
  t->SetBranchAddress("docatrk", &docatrk);

  bla ini = {1., 0., 0., 0., 
	     0., 0.,
	     0., 0., 0., 
	     0., 0., 0., 0., 0.};

  list<bla> bestList(1, ini) ;

  vector<bla> bestVector;
  for (int i = 0 ; i < nsettings; ++i) {
    bestVector.push_back(ini); 
  }

  bestVector.reserve(200000); 

  int nb(0); 
  int nentries = Int_t(t->GetEntries());
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    if (chan != ichan) continue;
    ini.ul = ul; 
    ini.nobs = nobs; 
    ini.nexp = nexp;
    ini.eff = eff; 
    ini.mlo = mlo; 
    ini.mhi = mhi; 
    ini.pt = pt; 
    ini.m1pt = m1pt; 
    ini.m2pt = m2pt; 
    ini.iso1 = iso1; 
    ini.chi2dof = chi2dof; 
    ini.alpha = alpha; 
    ini.fls3d = fls3d; 
    ini.docatrk = docatrk; 
    for (list<bla>::iterator i = bestList.begin(); i != bestList.end(); ++i) {
      if (ini.ul < i->ul) { 
	bestList.insert(i, ini); 
	break;
      }
    }

//     for (int i = 0; i < bestVector.size(); ++i) {
//       if (ini.ul < bestVector[i].ul) { 
// 	bestVector.insert(bestVector.begin() + i, 1, ini); 
// 	break;
//       }
//     }
  }

  int icnt(0); 
  for (list<bla>::iterator i = bestList.begin(); i != bestList.end(); ++i) {
    ++icnt;
    if (icnt > nsettings) break;
    ini = *i;
    cout << Form("ul=%2.1e nobs=%2.0f nexp=%4.3f e=%5.4f ", ini.ul,   ini.nobs,  ini.nexp, ini.eff)
	 << Form("%4.3f<m<%4.3f pT=%3.2f ",  ini.mlo,  ini.mhi,  ini.pt)
	 << Form("pt1=%3.2f pt2=%3.2f I=%3.2f c/d=%3.2f a=%4.3f f=%4.3f d=%4.3f", 
		 ini.m1pt, ini.m2pt, ini.iso1, ini.chi2dof, ini.alpha, ini.fls3d, ini.docatrk)
	 << endl;
  }

//   for (int i = 0; i < nsettings; ++i) {
//     ini = bestVector[i];
//     cout << Form("ul=%3.2e nobs=%2.0f nexp=%4.3f ", ini.ul,   ini.nobs,  ini.nexp)
// 	 << Form("mlo=%4.3f mhi=%4.3f pT=%4.3f ",  ini.mlo,  ini.mhi,  ini.pt)
// 	 << Form("m1pt=%4.3f m2pt=%4.3f iso=%4.3f c/d=%4.3f a=%4.3f fls3d=%4.3f", ini.m1pt, ini.m2pt, ini.iso1, ini.chi2dof, ini.alpha, ini.fls3d)
// 	 << endl;
//   }

}



// ----------------------------------------------------------------------
void anaBmm::testSimpleUL(const char *cuts) {
  //  ofstream OUT("testUL.txt", ios::app);

  fpMc[fSgMc]->cd();

  TH1D *h = (TH1D*)gFile->Get("uSG");  
  if (0 == h) {
    h = new TH1D("uSG", "", 40, fMassLo, fMassHi);
  } else{
    h->Reset();
  }

  const char *defCuts = "hlt&&gmuid&&gtqual&&m1q*m2q<0"; 

  // -- Signal: event counts *including* muon and trigger cuts!
  TTree *ts = (TTree*)(fpMc[fSgMc]->Get("events"));
  // -- all cuts
  string cutString = Form("%s&&%s&&m>%f&&m<%f", defCuts, cuts, fSigLo, fSigHi);
  cout << cutString << endl;
  ts->Draw("m>>uSG", cutString.c_str(), "goff");
  double afterCuts = h->GetSumOfWeights();

  fNumbersBs[0]->effTot  = afterCuts/(fNumbersBs[0]->genYield);
  fNumbersBs[0]->effTotE = dEff(static_cast<int>(afterCuts), static_cast<int>(fNumbersBs[0]->genYield));

  // -- Background expectation
  TTree *td = (TTree*)(fpData[fSgData]->Get("events"));
  td->Draw("m>>uSG", Form("%s&&%s", defCuts, cuts), "goff");
  h->Draw();
  c0->Modified();
  c0->Update();

  bgBlind(h, 1, 4.7, 6.0);
  
  fNobs = static_cast<int>(fBgExp + 0.5);
  
  fBgExp += 0.01; // to protect against blimit assert

  if (fBgExp < 0) {
    fBgExp = fBgHistExp; 
    fBgExpE = fBgHistExpE; 
  }
    
  cout << "eff:   " << afterCuts << "/" << fNumbersBs[0]->genYield << " = " << fNumbersBs[0]->effTot << "+/-" << fNumbersBs[0]->effTotE << endl;
  cout << "BG:    " << fBgExp << "+/-" << fBgExpE << " histogram counts: " 
       << fBgHistExp << "+/-" << fBgHistExpE
       << endl;
  cout << "Nobs : " << fNobs << endl;
  if (fNobs > 10) return; 

  double nbs = 2.0e9*(1.0-0.12);
  
  double nulbarlow = barlow(fNobs, fBgExp, fBgExpE, 0.2);
  double nulbayes  = blimit(0.9, fNobs, 1.0, 0.2, fBgExp, fBgExpE, 1);

  double ulbarlow = nulbarlow/(fNumbersBs[0]->effTot*nbs);
  double ulbayes  = nulbayes/(fNumbersBs[0]->effTot*nbs);
  cout << Form("Nul: barlow = %3.2f, bayes = %4.2e", nulbarlow, nulbayes) << endl;
  fOUT << "==> " << fUL 
      << ". obs = " << fNobs << " bg = " << fBgExp << "+/-" << fBgExpE << "(fBgHistExp = " << fBgHistExp << ")" << endl
      << " eff = " << afterCuts << "/" << fNumbersBs[0]->genYield << " = " << fNumbersBs[0]->effTot << endl 
      << endl;

  cout << "UL: barlow = " << ulbarlow << " ulbayes " << ulbayes << " BG: " << fBgExp << " hist: " << h->GetSumOfWeights() << endl;
  h->Draw();
  fOUT.flush();
}


// ----------------------------------------------------------------------
double anaBmm::barlow(int nobs, double bg, double bgE, double sE) {
  double ul(-99.); 

  if (nobs > 10)  {
    cout << "barlow: dunno" << endl;
    return -99;
  }

  // -- "unified confidence limits" from PDG as intervals to probe 
  double uclLo[] = {0.00, 0.11, 0.53, 1.10, 1.47, 1.84, 2.21, 3.56, 3.96, 4.36, 5.50};
  double uclHi[] = {2.44, 4.36, 5.91, 7.42, 8.60, 9.99,11.47,12.53,13.99,15.30,16.50};

  double cl(0.1); // this is 90% UL!

  // -- Set up Gaussian pdf for expected background 
  double bgSig= (bgE > 0. ? bgE : bg/3.); 

  // -- Set up Gaussian pdf for sensitivity
  double sSig= (sE > 0. ? sE : 0.1); 

  // -- Set up histogram for counting
  TH1D *h = (TH1D*)gROOT->FindObject("hbarlow"); 
  if (h) {
    h->Reset();
  } else {
    //    h = new TH1D("hbarlow", "", 10000, 0., 100.); 
    h = new TH1D("hbarlow", "", 50, 0., 50.); 
  }

  int bin = h->FindBin(nobs); 

  // -- Run toy
  double tSignal(0.), tBackground(0.), tSensitivity(1.), tCounts(0.);
  double tFraction(0.); 
  double minDeviation(99.), bestAlpha(99.); 

  int NSTEPS(100); 
  double step = (uclHi[nobs] - uclLo[nobs])/NSTEPS; 

  if (nobs <= 2) {
    step   = 0.02; 
  } else {
    step   = 0.1; 
  }
  NSTEPS = int((uclHi[nobs] - uclLo[nobs])/step + 1);

  // -- This loop is quite rough ...?
  for (int k = 0; k < NSTEPS; ++k) {
    h->Reset();
    tSignal = uclLo[nobs] + k*step; 
    
    for (int i = 0; i < 100000; ++i) {
      tBackground  = gRandom->Gaus(bg, bgSig); 
      tSensitivity = gRandom->Gaus(1., sSig); 
    
      tCounts = gRandom->Poisson((tSignal+tBackground)*tSensitivity);
      h->Fill(tCounts); 
    }
    tFraction= h->Integral(1, bin) / h->Integral();
    if (TMath::Abs(tFraction - cl) < minDeviation) {
      minDeviation = TMath::Abs(tFraction - cl); 
      bestAlpha    = tFraction; 
      ul           = tSignal;
    }

    if (tFraction < (cl - minDeviation)) {
      //      cout << "break " << tFraction << " bestAlpha: " << bestAlpha << " (should be " << cl << ")" << endl;
      break;
    }
  }

  return ul; 
}

// ----------------------------------------------------------------------
void anaBmm::rolkeM3() {

  // -- Old interface to TRolke in ROOT 5.22
  int    x   = fNobs;     // observed number of events FIXME
  double bm  = fBgExp;    // measured bg
  double sdb = fBgExpE;   // error of bg
  double tau = 1.0;       // scale factor of bg to sg region
  int     mid = 3;        // model ID
  int m = 0;              // number of MC runs
  int z = 0;              // simulated events
  int y = 0;              // simulated events
  double e = 1.0;         // true efficiency if known
  double em = 1.0;        // measured efficiency
  double sde= fEffE/fEff; // error on measured efficiency
  double b = fBgExp;      //??

  //  fpRolke->SetGaussBkgGaussEff(x,bm,em,sde,sdb);
  //  fNul = fpRolke->GetUpperLimit();
  fNul = fpRolke->CalculateInterval(x,y,z,bm,em,e,mid,sde,sdb,tau,b,m); 
 
}


// ----------------------------------------------------------------------
void anaBmm::rolkeM3(int x, double bm, double em, double sde, double sdb) {
//   fpRolke->SetGaussBkgGaussEff(x,bm,em,sde,sdb);
//   fNul = fpRolke->GetUpperLimit();

  // -- Old interface to TRolke in ROOT 5.22
  int    xl   = x;         // observed number of events FIXME
  double tau = 1.0;       // scale factor of bg to sg region
  int     mid = 3;        // model ID
  int m = 0;              // number of MC runs
  int z = 0;              // simulated events
  int y = 0;              // simulated events
  double e = 1.0;         // true efficiency if known
  double b = fBgExp;      //??

  if (bm < 0.2) {
    xl = 0; 
  } else if (bm < 1.0) {
    xl = 1;
  } else {
    xl = fBgExp;
  }
    

  cout << "xl = " << xl << endl;
  cout << "bm = " << bm << endl;
  cout << "sde = " << sde << endl;

  fNul = fpRolke->CalculateInterval(xl,y,z,bm,em,e,mid,sde,sdb,tau,b,m); 
  cout << "Nul = " << fNul << endl;
}


// ----------------------------------------------------------------------
void anaBmm::bgBlind(TH1 *h, int mode, double lo, double hi) {
  
  if (0 == h) { 
    cout << "anaBmm::bgBlind(...): No histogram passed! mode = " << mode << endl;
    return;
  }
  
  TF1 *lF1(0);

  double histCount = h->Integral(h->FindBin(fBgLo), h->FindBin(fBgHi)-1); 
  cout << "bgBlind: histCount = " << histCount << " starting at " << h->FindBin(fBgLo) << " to " << h->FindBin(fBgHi)-1 << endl;
  fBgHist  = histCount; 
  fBgHistE  = TMath::Sqrt(histCount); 
  fBgHistExp  = histCount*(fSigHi-fSigLo)/(fBgHi-fBgLo-0.25); // FIXME fixed limits
  if (histCount > 0) {
    fBgHistE = TMath::Sqrt(histCount)/histCount*fBgHist;
  } else {
    fBgHistE = 0.2; // FIXME?!
    fBgExp = 0.;
    fBgExpE = 0.2; 
    return;
  }

  if (0 == mode) {
    fBgExp = fBgHist;
    fBgExpE = fBgHistE;
    return;
  }
  else if (1 == mode) { 
    //    lF1 = fpFunc->pol0(h); 
    lF1 = fpFunc->pol0BsBlind(h); 
  } else if (2 == mode) {
    lF1 = fpFunc->pol1BsBlind(h); 
  }
  
  //  setErrors(h);
  h->Fit(lF1, "lr", "", lo, hi); 

  double c  = h->GetFunction("f1")->GetParameter(0); 
  double cE = h->GetFunction("f1")->GetParError(0); 
  fBgExp  = c * (fSigHi - fSigLo)/h->GetBinWidth(1);
  fBgExpE = cE/c*fBgExp; 
  cout << "bgBlind: c = " << c << " sig width = " << (fSigHi - fSigLo) << endl;
}


// ----------------------------------------------------------------------
void anaBmm::normYield(TH1 *h, int mode, double lo, double hi) {

  TF1 *lF1(0);
  
  fpFunc->fLo = lo;
  fpFunc->fHi = hi;
  lF1 = fpFunc->pol1Gauss(h); 
  h->Fit(lF1, "lr", "", lo, hi); 

  double c  = h->GetFunction("f1")->GetParameter(0); 
  double cE = h->GetFunction("f1")->GetParError(0); 

  fNormSig = c/h->GetBinWidth(1);
  fNormSigE = cE/c*fNormSig;

  cout << "N(Sig) = " << fNormSig << " +/- " << fNormSigE << endl;
  
  delete lF1; 

}


// ----------------------------------------------------------------------
void anaBmm::csYield(TH1 *h, int mode, double lo, double hi) {

  TF1 *lF1(0);
  
  fpFunc->fLo = lo;
  fpFunc->fHi = hi;
  lF1 = fpFunc->pol1Gauss(h); 
  h->Fit(lF1, "lr", "", lo, hi); 

  double c  = h->GetFunction("f1")->GetParameter(0); 
  double cE = h->GetFunction("f1")->GetParError(0); 

  fCsSig = c/h->GetBinWidth(1);
  fCsSigE = cE/c*fCsSig;

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  
  delete lF1; 

}


// ----------------------------------------------------------------------
string anaBmm::formatTex(double n, std::string name, int digits) {

  char line[200]; 
  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str());
  } else if (0 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.0f } } }", name.c_str(), n);
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.1f } } }", name.c_str(), n);
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.2f } } }", name.c_str(), n);
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.3f } } }", name.c_str(), n);
  } else if (4 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.4f } } }", name.c_str(), n);
  } else if (5 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%6.5f } } }", name.c_str(), n);
  } else if (6 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%7.6f } } }", name.c_str(), n);
  } else {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%f } } }", name.c_str(), n);
  }

  return string(line); 
}

// ----------------------------------------------------------------------
string anaBmm::formatTex(double n, string name, string tag) {
  
  char line[200]; 
  
  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str(), tag.c_str());
    //   } else if ( n > 1.e10) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm2(n)).Data());
    //   } else if ( n > 1.e4) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm(n)).Data());
  } else if ( n > 100. ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%6.0f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1. ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.2f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1.e-1) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1.e-3) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name.c_str(), tag.c_str(), n);
    //   } else if ( n > 1.e-9 ) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm(n)).Data());
    //   } else if ( n > 1.e-19 ){
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm2(n)).Data());
  } else {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{0.0 } } }", name.c_str(), tag.c_str());
  }
  
  string result(line); 
  return result;
}



// ----------------------------------------------------------------------
void anaBmm::makeCanvas(int i) {
  if (i & 16) { 
    c5 = new TCanvas("c5", "c5", 210,   0, 800, 1000);
    c5->ToggleEventStatus();
  }
  if (i & 8) { 
    c4 = new TCanvas("c4", "c4", 210,   0, 800, 600);
    c4->ToggleEventStatus();
  }
  if (i & 4) {
    c3 = new TCanvas("c3", "c3", 200,  20, 800, 800);
    c3->ToggleEventStatus();
  }
  if (i & 1) {
    //    c1 = new TCanvas("c1", "c1", 20,  60, 1200, 400);
    c1 = new TCanvas("c1", "c1", 20,  60, 1000, 400);
    c1->ToggleEventStatus();
  }
  if (i & 2) { 
    c2 = new TCanvas("c2", "c2", 300, 200, 400, 800);
    c2->ToggleEventStatus();
  }
}


// ----------------------------------------------------------------------
int anaBmm::wait() {
  cout << " Continue [<RET>|q]?  "; 
  char x;
  x = getchar();
  if ((x == 'q') || (x == 'Q')) return 1;
  return 0;
}


// ----------------------------------------------------------------------
void anaBmm::replaceAll(std::string &s, std::string a, std::string b) {
  
  TString ts(s.c_str()); 
  ts.ReplaceAll(a.c_str(), b.c_str()); 
  s = ts.Data(); 

}


// ----------------------------------------------------------------------
void anaBmm::setErrors(TH1D *h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinLowEdge(i) > 5.1 && h->GetBinLowEdge(i+1) <= 5.5) continue;
    if (h->GetBinContent(i) < 1) {
      h->SetBinError(i, 1.); 
    }
  }
}


// ----------------------------------------------------------------------
void anaBmm::texUlcalcNumbers(const char *filename, const char *output) {
  cout << "==> anaBmm: TeXing " << filename << " for UL numbers" << endl;
  vector<string> cutLines; 
  char  buffer[200];
  ifstream is(filename);
  while (is.getline(buffer, 200, '\n')) {
    cutLines.push_back(string(buffer));
  }

  ofstream TEX(output);

  
  string sname; 
  char name[100]; 
  int digits(-1); 
  int channel; 
  float val, err; 
  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    if (buffer[0] == '#') continue;
    val = -1.; 
    err = 0.;
    sscanf(buffer, "%s\t%d\t%f\t%f", name, &channel, &val, &err);
      //    sscanf(buffer, "%s\t%d\t%f", name, &channel, &val);
    sname = name; 
    replaceAll(sname, "_", "-"); 
    cout << sname << "[" << channel << "] : " << val << " +/- " << err << endl;
    if (string::npos != sname.find("OBS")) digits = 0; 
    if (string::npos != sname.find("EFF")) digits = 4; 
    if (string::npos != sname.find("ACC")) digits = 4; 
    if (string::npos != sname.find("LOW")) digits = 3; 
    if (string::npos != sname.find("HIGH")) digits = 3; 
    TEX << formatTex(val, Form("%s:%s%d:val", fSuffix.c_str(), sname.c_str(), channel), digits) << endl;
    TEX << formatTex(err, Form("%s:%s%d:err", fSuffix.c_str(), sname.c_str(), channel), digits) << endl;
  }
  
  TEX.close();
}


// ----------------------------------------------------------------------
void anaBmm::printCsBFNumbers() {

  fTEX << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  for (unsigned int i = 0; i < fNchan; ++i) {
    fTEX << "% -- CONTROL SAMPLE " << i << endl;
    fTEX << formatTex(fNumbersCS[i]->effTot, Form("%s:N-EFF-TOT-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(fNumbersCS[i]->effTotE, Form("%s:N-EFF-TOT-BS%i:err", fSuffix.c_str(), i), 6) << endl;

    fTEX << formatTex(fNumbersCS[i]->acc, Form("%s:N-ACC-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCS[i]->accE, Form("%s:N-ACC-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCS[i]->effMuidPid, Form("%s:N-EFF-MU-PID-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCS[i]->effMuidPidE, Form("%s:N-EFF-MU-PID-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCS[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCS[i]->effMuidMCE, Form("%s:N-EFF-MU-MC-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCS[i]->effTrigPid, Form("%s:N-EFF-TRIG-PID-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCS[i]->effTrigPidE, Form("%s:N-EFF-TRIG-PID-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCS[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCS[i]->effTrigMCE, Form("%s:N-EFF-TRIG-MC-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCS[i]->effCand, Form("%s:N-EFF-CAND-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCS[i]->effCandE, Form("%s:N-EFF-CAND-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCS[i]->effAna, Form("%s:N-EFF-ANA-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCS[i]->effAnaE, Form("%s:N-EFF-ANA-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCS[i]->fitYield, Form("%s:N-OBS-BS%i:val", fSuffix.c_str(), i), 0) << endl;
    fTEX << formatTex(fNumbersCS[i]->fitYieldE, Form("%s:N-OBS-BS%i:err", fSuffix.c_str(), i), 0) << endl;
  }
}

// ----------------------------------------------------------------------
void anaBmm::printUlcalcNumbers() {
  ofstream OUT(fUlcalcFileName.c_str());

  OUT << "######################################################################" << endl;
  fTEX << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  //  if (fNumbersNorm[0]->fitYield > 0) {
  if (1) {
    for (unsigned int i = 0; i < fNchan; ++i) {
      OUT << "# -- NORMALIZATION " << i << endl;
      fTEX << "% -- NORMALIZATION " << i << endl;

      OUT << "#EFF_TOT_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effTot << endl;
      fTEX << formatTex(fNumbersNorm[i]->effTot, Form("%s:N-EFF-TOT-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNorm[i]->effTotE, Form("%s:N-EFF-TOT-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;

      OUT << "ACC_BPLUS\t" << i << "\t" << fNumbersNorm[i]->acc << endl;
      fTEX << formatTex(fNumbersNorm[i]->acc, Form("%s:N-ACC-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNorm[i]->accE, Form("%s:N-ACC-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;

      //    OUT << "EFF_MU_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effMuidPid << endl;
      fTEX << formatTex(fNumbersNorm[i]->effMuidPid, Form("%s:N-EFF-MU-PID-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNorm[i]->effMuidPidE, Form("%s:N-EFF-MU-PID-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;

      OUT << "EFF_MU_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effMuidMC << endl;
      fTEX << formatTex(fNumbersNorm[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNorm[i]->effMuidMCE, Form("%s:N-EFF-MU-MC-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;

      //    OUT << "EFF_TRIG_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effTrigPid << endl;
      fTEX << formatTex(fNumbersNorm[i]->effTrigPid, Form("%s:N-EFF-TRIG-PID-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNorm[i]->effTrigPidE, Form("%s:N-EFF-TRIG-PID-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;

      OUT << "EFF_TRIG_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effTrigMC << endl;
      fTEX << formatTex(fNumbersNorm[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNorm[i]->effTrigMCE, Form("%s:N-EFF-TRIG-MC-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;

      OUT << "EFF_CAND_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effCand << endl;
      fTEX << formatTex(fNumbersNorm[i]->effCand, Form("%s:N-EFF-CAND-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNorm[i]->effCandE, Form("%s:N-EFF-CAND-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;

      //    OUT << "EFF_ANA_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effAna << endl;
      OUT << "EFF_ANA_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effAna << endl;
      fTEX << formatTex(fNumbersNorm[i]->effAna, Form("%s:N-EFF-ANA-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNorm[i]->effAnaE, Form("%s:N-EFF-ANA-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;

      OUT << "OBS_BPLUS\t" << i << "\t" << fNumbersNorm[i]->fitYield << endl;
      fTEX << formatTex(fNumbersNorm[i]->fitYield, Form("%s:N-OBS-BPLUS%i:val", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(fNumbersNorm[i]->fitYieldE, Form("%s:N-OBS-BPLUS%i:err", fSuffix.c_str(), i), 0) << endl;
    } 
  } else {
    OUT << "TOT_BPLUS\t" << "0\t" << (fDataLumi[fSgData]/39.4)*440000 << endl;
    OUT << "TOT_BPLUS\t" << "1\t" << (fDataLumi[fSgData]/39.4)*383000 << endl;
  }

  OUT << "######################################################################" << endl;
  fTEX << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  double scaledSig, scaledSigE; 
  for (unsigned int i = 0; i < fNchan; ++i) {
    OUT << "# -- SIGNAL " << i << endl;
    fTEX << "% -- SIGNAL " << i << endl;

    scaledSig  = fNumbersBs[i]->anaWmcYield*fDataLumi[fSgData]/fLumi["SgMc"];
    scaledSigE = scaledSig*TMath::Sqrt(fNumbersBs[i]->anaWmcYield)/fNumbersBs[i]->anaWmcYield;
    fTEX << formatTex(scaledSig, Form("%s:N-EXP-SIG-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(scaledSigE, Form("%s:N-EXP-SIG-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    scaledSig  = fNumbersBd[i]->anaWmcYield*fDataLumi[fSgData]/fLumi["BdMc"];
    scaledSigE = scaledSig*TMath::Sqrt(fNumbersBd[i]->anaWmcYield)/fNumbersBd[i]->anaWmcYield;
    fTEX << formatTex(scaledSig, Form("%s:N-EXP-SIG-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(scaledSigE, Form("%s:N-EXP-SIG-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "OBS_BKG\t" << i << "\t" << fNumbersBs[i]->bgObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bgObs, Form("%s:N-OBS-BKG%d:val", fSuffix.c_str(), i), 0) << endl;
    fTEX << formatTex(fNumbersBs[i]->bgBsExp, Form("%s:N-EXP-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->bgBsExpE, Form("%s:N-EXP-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->bgBdExp, Form("%s:N-EXP-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->bgBdExpE, Form("%s:N-EXP-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;

    OUT << "LOW_BD\t" << i << "\t" << fNumbersBs[i]->mBdLo << endl;
    fTEX << formatTex(fNumbersBs[i]->mBdLo, Form("%s:N-LOW-BD%d:val", fSuffix.c_str(), i), 3) << endl;

    OUT << "HIGH_BD\t" << i << "\t" << fNumbersBs[i]->mBdHi << endl;
    fTEX << formatTex(fNumbersBs[i]->mBdHi, Form("%s:N-HIGH-BD%d:val", fSuffix.c_str(), i), 3) << endl;

    OUT << "LOW_BS\t" << i << "\t" << fNumbersBs[i]->mBsLo << endl;
    fTEX << formatTex(fNumbersBs[i]->mBsLo, Form("%s:N-LOW-BS%d:val", fSuffix.c_str(), i), 3) << endl;

    OUT << "HIGH_BS\t" << i << "\t" << fNumbersBs[i]->mBsHi << endl;
    fTEX << formatTex(fNumbersBs[i]->mBsHi, Form("%s:N-HIGH-BS%d:val", fSuffix.c_str(), i), 3) << endl;

    OUT << "PSS\t" << i << "\t" << fNumbersBs[i]->pss << endl;
    fTEX << formatTex(fNumbersBs[i]->pss, Form("%s:N-PSS%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->pssE, Form("%s:N-PSS%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "PSD\t" << i << "\t" << fNumbersBd[i]->psd << endl;
    fTEX << formatTex(fNumbersBs[i]->psd, Form("%s:N-PSD%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->psdE, Form("%s:N-PSD%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "PDS\t" << i << "\t" << fNumbersBs[i]->pds << endl;
    fTEX << formatTex(fNumbersBs[i]->pds, Form("%s:N-PDS%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->pdsE, Form("%s:N-PDS%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "PDD\t" << i << "\t" << fNumbersBd[i]->pdd << endl;
    fTEX << formatTex(fNumbersBs[i]->pdd, Form("%s:N-PDD%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->pddE, Form("%s:N-PDD%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "#EFF_TOT_BSMM\t" << i << "\t" << fNumbersBs[i]->effTot << endl;
    fTEX << formatTex(fNumbersBs[i]->effTot, Form("%s:N-EFF-TOT-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTotE, Form("%s:N-EFF-TOT-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "ACC_BSMM\t" << i << "\t" << fNumbersBs[i]->acc << endl;
    fTEX << formatTex(fNumbersBs[i]->acc, Form("%s:N-ACC-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->accE, Form("%s:N-ACC-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    //    OUT << "EFF_MU_BSMM\t" << i << "\t" << fNumbersBs[i]->effMuidPid << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidPid, Form("%s:N-EFF-MU-PID-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidPidE, Form("%s:N-EFF-MU-PID-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_MU_BSMM\t" << i << "\t" << fNumbersBs[i]->effMuidMC << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidMCE, Form("%s:N-EFF-MU-MC-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    //    OUT << "EFF_TRIG_BSMM\t" << i << "\t" << fNumbersBs[i]->effTrigPid << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigPid, Form("%s:N-EFF-TRIG-PID-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigPidE, Form("%s:N-EFF-TRIG-PID-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_TRIG_BSMM\t" << i << "\t" << fNumbersBs[i]->effTrigMC << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigMCE, Form("%s:N-EFF-TRIG-MC-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_CAND_BSMM\t" << i << "\t" << fNumbersBs[i]->effCand << endl;
    fTEX << formatTex(fNumbersBs[i]->effCand, Form("%s:N-EFF-CAND-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effCandE, Form("%s:N-EFF-CAND-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_ANA_BSMM\t" << i << "\t" << fNumbersBs[i]->effAna << endl;
    fTEX << formatTex(fNumbersBs[i]->effAna, Form("%s:N-EFF-ANA-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effAnaE, Form("%s:N-EFF-ANA-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "#EFF_TOT_BDMM\t" << i << "\t" << fNumbersBd[i]->effTot << endl;
    fTEX << formatTex(fNumbersBd[i]->effTot, Form("%s:N-EFF-TOT-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTotE, Form("%s:N-EFF-TOT-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "ACC_BDMM\t" << i << "\t" << fNumbersBd[i]->acc << endl;
    fTEX << formatTex(fNumbersBd[i]->acc, Form("%s:N-ACC-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->accE, Form("%s:N-ACC-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    //    OUT << "EFF_MU_BDMM\t" << i << "\t" << fNumbersBd[i]->effMuidPid << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidPid, Form("%s:N-EFF-MU-PID-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidPidE, Form("%s:N-EFF-MU-PID-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_MU_BDMM\t" << i << "\t" << fNumbersBd[i]->effMuidMC << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidMCE, Form("%s:N-EFF-MU-MC-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    //    OUT << "EFF_TRIG_BDMM\t" << i << "\t" << fNumbersBd[i]->effTrigPid << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigPid, Form("%s:N-EFF-TRIG-PID-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigPidE, Form("%s:N-EFF-TRIG-PID-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_TRIG_BDMM\t" << i << "\t" << fNumbersBd[i]->effTrigMC << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigMCE, Form("%s:N-EFF-TRIG-MC-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_CAND_BDMM\t" << i << "\t" << fNumbersBd[i]->effCand << endl;
    fTEX << formatTex(fNumbersBd[i]->effCand, Form("%s:N-EFF-CAND-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effCandE, Form("%s:N-EFF-CAND-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_ANA_BDMM\t" << i << "\t" << fNumbersBd[i]->effAna << endl;
    fTEX << formatTex(fNumbersBd[i]->effAna, Form("%s:N-EFF-ANA-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effAnaE, Form("%s:N-EFF-ANA-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;

    OUT << "# Observed in signal boxes" << endl;
    OUT << "OBS_BSMM\t" << i << "\t" << fNumbersBs[i]->bsObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bsObs, Form("%s:N-OBS-BSMM%d:val", fSuffix.c_str(), i), 0) << endl;

    OUT << "OBS_BDMM\t" << i << "\t" << fNumbersBs[i]->bdObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bdObs, Form("%s:N-OBS-BDMM%d:val", fSuffix.c_str(), i), 0) << endl;
    
  }

  OUT.close();

}


// ----------------------------------------------------------------------
void anaBmm::printNumbers(numbers &a, ostream &OUT) {
  OUT << "======================================================================" << endl;
  OUT << "numbers for \"" << a.name.c_str() << "\"" << endl;
  OUT << "fitYield     = " << a.fitYield << "+/-" << a.fitYieldE << endl;
  OUT << "genFileYield = " << a.genFileYield << endl;
  OUT << "genYield     = " << a.genYield << endl;
  OUT << "genChanYield = " << a.genChanYield << endl;
  OUT << "recoYield    = " << a.recoYield << endl;
  OUT << "chanYield    = " << a.chanYield << endl;
  OUT << "muidYield    = " << a.muidYield << endl;
  OUT << "trigYield    = " << a.trigYield << endl;
  OUT << "candYield    = " << a.candYield << endl;
  OUT << "ana0Yield    = " << a.ana0Yield << endl;
  OUT << "anaYield     = " << a.anaYield << endl;
  OUT << "anaMuYield   = " << a.anaMuonYield << endl;
  OUT << "anaTrigYield = " << a.anaTriggerYield << endl;
  OUT << "anaWmcYield  = " << a.anaWmcYield << endl;
  OUT << "mBsLo        = " << a.mBsLo << endl;
  OUT << "mBsHi        = " << a.mBsHi << endl;
  OUT << "mBdLo        = " << a.mBdLo << endl;
  OUT << "mBdHi        = " << a.mBdHi << endl;
  OUT << "PSS          = " << a.pss << endl;
  OUT << "PDS          = " << a.pds << endl;
  OUT << "PSD          = " << a.psd << endl;
  OUT << "PDD          = " << a.pdd << endl;
  OUT << "bsRare       = " << a.bsRare << endl;
  OUT << "bdRare       = " << a.bdRare << endl;
  OUT << "gen filter   = " << a.effGenFilter << endl;
  OUT << "acceptance   = " << a.acc << "+/-" << a.accE << endl;
  OUT << "accChan      = " << a.accChan << "+/-" << a.accChanE << endl; 
  OUT << "cFrac        = " << a.cFrac << "+/-" << a.cFracE << endl; 
  OUT << "effChan      = " << a.effChan << "+/-" << a.effChanE << endl; 
  OUT << "effCand      = " << a.effCand << "+/-" << a.effCandE << endl;
  OUT << "effAna       = " << a.effAna << "+/-" << a.effAnaE << endl; 
  OUT << "effMuidMC    = " << a.effMuidMC << "+/-" << a.effMuidMCE << endl;
  OUT << "effMuidPid   = " << a.effMuidPid << "+/-" << a.effMuidPidE << endl;
  OUT << "accMuidMC    = " << a.accMuidMC << "+/-" << a.accMuidMCE << endl;
  OUT << "effTrigMC    = " << a.effTrigMC << "+/-" << a.effTrigMCE << endl;
  OUT << "effTrigPid   = " << a.effTrigPid << "+/-" << a.effTrigPidE << endl;
  OUT << "accTrigMC    = " << a.accTrigMC << "+/-" << a.accTrigMCE << endl;
  OUT << "effProd(MC)  = " << a.effProdMC << endl;
  OUT << "effProd(Pid) = " << a.effProdPid << endl;
  OUT << "effProd(MC)A = " << a.effProdMC*a.acc << endl;
  OUT << "effTot       = " << a.effTot << "+/-" << a.effTotE << endl; 
  OUT << "effTotChan   = " << a.effTotChan << "+/-" << a.effTotChanE << endl; 
  OUT << "combGenYield     = " << a.combGenYield << endl; 
  OUT << "prodGenYield     = " << a.prodGenYield << endl; 
  OUT << "prodChanGenYield = " << a.chanGenYield << endl; 
  OUT.flush();

  // -- dump into fHistFile
  int bin(0); 
  // -- Cache the pwd...
  TDirectory *pD = gFile; 
  fHistFile->cd();
  TH1D *hn = new TH1D(Form("numbers: %s", a.name.c_str()), Form("numbers: %s", a.name.c_str()), 100, 0., 100.); hn->Sumw2();
  bin = 1; 
  hn->GetXaxis()->SetBinLabel(bin, "fitYield"); hn->SetBinContent(bin, a.fitYield);  hn->SetBinError(bin, a.fitYieldE); 

  bin = 2; 
  hn->GetXaxis()->SetBinLabel(bin, "genFileYield"); hn->SetBinContent(bin, a.genFileYield);  

  bin = 3; 
  hn->GetXaxis()->SetBinLabel(bin, "genYield"); hn->SetBinContent(bin, a.genYield);  

  bin = 4; 
  hn->GetXaxis()->SetBinLabel(bin, "recoYield"); hn->SetBinContent(bin, a.recoYield);  

  bin = 5; 
  hn->GetXaxis()->SetBinLabel(bin, "muidYield"); hn->SetBinContent(bin, a.muidYield);  

  bin = 6; 
  hn->GetXaxis()->SetBinLabel(bin, "trigYield"); hn->SetBinContent(bin, a.trigYield);  

  bin = 7; 
  hn->GetXaxis()->SetBinLabel(bin, "candYield"); hn->SetBinContent(bin, a.candYield);  

  bin = 8; 
  hn->GetXaxis()->SetBinLabel(bin, "ana0Yield"); hn->SetBinContent(bin, a.ana0Yield);  

  bin = 9; 
  hn->GetXaxis()->SetBinLabel(bin, "anaWmcYield"); hn->SetBinContent(bin, a.anaWmcYield);  

  bin = 10; 
  hn->GetXaxis()->SetBinLabel(bin, "anaYield"); hn->SetBinContent(bin, a.anaYield);  

  bin = 11; 
  hn->GetXaxis()->SetBinLabel(bin, "effGenFilter"); hn->SetBinContent(bin, a.effGenFilter);  

  bin = 12; 
  hn->GetXaxis()->SetBinLabel(bin, "acc"); hn->SetBinContent(bin, a.acc);  

  bin = 13; 
  hn->GetXaxis()->SetBinLabel(bin, "effMuidMC"); hn->SetBinContent(bin, a.effMuidMC);  

  bin = 14; 
  hn->GetXaxis()->SetBinLabel(bin, "effTrigMC"); hn->SetBinContent(bin, a.effTrigMC);  

  bin = 15; 
  hn->GetXaxis()->SetBinLabel(bin, "effMuidPid"); hn->SetBinContent(bin, a.effMuidPid);  

  bin = 16; 
  hn->GetXaxis()->SetBinLabel(bin, "effTrigPid"); hn->SetBinContent(bin, a.effTrigPid);  

  bin = 17; 
  hn->GetXaxis()->SetBinLabel(bin, "effCand"); hn->SetBinContent(bin, a.effCand);  

  bin = 18; 
  hn->GetXaxis()->SetBinLabel(bin, "effAna"); hn->SetBinContent(bin, a.effAna);  

  bin = 19; 
  hn->GetXaxis()->SetBinLabel(bin, "effTot"); hn->SetBinContent(bin, a.effTot);  

  bin = 30; 
  hn->GetXaxis()->SetBinLabel(bin, "pss"); hn->SetBinContent(bin, a.pss);  

  bin = 31; 
  hn->GetXaxis()->SetBinLabel(bin, "pdd"); hn->SetBinContent(bin, a.pdd);  

  bin = 32; 
  hn->GetXaxis()->SetBinLabel(bin, "psd"); hn->SetBinContent(bin, a.psd);  

  bin = 33; 
  hn->GetXaxis()->SetBinLabel(bin, "pds"); hn->SetBinContent(bin, a.pds);  

  // -- and get back to it. 
  pD->cd();
}


// ----------------------------------------------------------------------
void anaBmm::newLegend(double x1, double y1, double x2, double y2) {
  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2);
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 
}


// ----------------------------------------------------------------------
void anaBmm::readCuts(const char *filename) {
  cout << "==> anaBmm: Reading " << filename << " for cut settings" << endl;
  vector<string> cutLines; 
  char  buffer[200];
  ifstream is(filename);
  while (is.getline(buffer, 200, '\n')) {
    cutLines.push_back(string(buffer));
  }

  char CutName[100];
  float CutValue;
  int dump(0), ok(0);

  cuts *a = 0;
  
  fCuts.clear();

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "index")) {
      ok = 1;
      if (dump) cout << "index:            " << CutValue << endl;
      if (a) fCuts.push_back(a); 
      a = new cuts; 
      a->index = static_cast<int>(CutValue); 
    }
    
    if (!strcmp(CutName, "mBdLo")) {
      a->mBdLo = CutValue; ok = 1;
      if (dump) cout << "mBdLo:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBdHi")) {
      a->mBdHi = CutValue; ok = 1;
      if (dump) cout << "mBdHi:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBsLo")) {
      a->mBsLo = CutValue; ok = 1;
      if (dump) cout << "mBsLo:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBsHi")) {
      a->mBsHi = CutValue; ok = 1;
      if (dump) cout << "mBsHi:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMin")) {
      a->etaMin = CutValue; ok = 1;
      if (dump) cout << "etaMin:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMax")) {
      a->etaMax = CutValue; ok = 1;
      if (dump) cout << "etaMax:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "pt")) {
      a->pt = CutValue; ok = 1;
      if (dump) cout << "pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1pt")) {
      a->m1pt = CutValue; ok = 1;
      if (dump) cout << "m1pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2pt")) {
      a->m2pt = CutValue; ok = 1;
      if (dump) cout << "m2pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1eta")) {
      a->m1eta = CutValue; ok = 1;
      if (dump) cout << "m1eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2eta")) {
      a->m2eta = CutValue; ok = 1;
      if (dump) cout << "m2eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "iso1")) {
      a->iso1 = CutValue; ok = 1;
      if (dump) cout << "iso1:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "chi2dof")) {
      a->chi2dof = CutValue; ok = 1;
      if (dump) cout << "chi2dof:             " << CutValue << endl;
    }

    if (!strcmp(CutName, "alpha")) {
      a->alpha = CutValue; ok = 1;
      if (dump) cout << "alpha:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "fls3d")) {
      a->fls3d = CutValue; ok = 1;
      if (dump) cout << "fls3d:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "docatrk")) {
      a->docatrk = CutValue; ok = 1;
      if (dump) cout << "docatrk:               " << CutValue << endl;
    }

  }

  if (a) fCuts.push_back(a); 

  if (!ok) cout << "==> what about " << CutName << endl;
  

}


// ----------------------------------------------------------------------
void anaBmm::printCuts(ostream &OUT) {

  OUT << "----------------------------------------------------------------------" << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    cuts *a = fCuts[i]; 
    OUT << "# -- channel " << a->index << endl;
    OUT << "index   " << a->index << endl;
    fTEX << "% ----------------------------------------------------------------------" << endl;
    fTEX << "% -- Cuts for channel " << a->index << endl;

    OUT << "mBdLo   " << Form("%4.3f", a->mBdLo) << endl;
    OUT << "mBdHi   " << Form("%4.3f", a->mBdHi) << endl;
    fTEX <<  Form("\\vdef{%s:mBdLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBdLo) << endl;
    fTEX <<  Form("\\vdef{%s:mBdHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBdHi) << endl;

    OUT << "mBsLo   " << Form("%4.3f", a->mBsLo) << endl;
    OUT << "mBsHi   " << Form("%4.3f", a->mBsHi) << endl;
    fTEX <<  Form("\\vdef{%s:mBsLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBsLo) << endl;
    fTEX <<  Form("\\vdef{%s:mBsHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBsHi) << endl;

    OUT << "etaMin  " << Form("%3.1f", a->etaMin) << endl;
    OUT << "etaMax  " << Form("%3.1f", a->etaMax) << endl;
    fTEX <<  Form("\\vdef{%s:etaMin:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->etaMin) << endl;
    fTEX <<  Form("\\vdef{%s:etaMax:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->etaMax) << endl;

    OUT << "pt      " << Form("%3.1f", a->pt) << endl;
    fTEX <<  Form("\\vdef{%s:pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->pt) << endl;
    OUT << "m1pt    " << a->m1pt << endl;
    fTEX <<  Form("\\vdef{%s:m1pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m1pt) << endl;
    OUT << "m2pt    " << a->m2pt << endl;
    fTEX <<  Form("\\vdef{%s:m2pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m2pt) << endl;
    OUT << "m1eta   " << a->m1eta << endl;
    fTEX <<  Form("\\vdef{%s:m1eta:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m1eta) << endl;
    OUT << "m2eta   " << a->m2eta << endl;
    fTEX <<  Form("\\vdef{%s:m2eta:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m2eta) << endl;

    OUT << "iso1    " << a->iso1 << endl;
    fTEX <<  Form("\\vdef{%s:iso1:%d}   {\\ensuremath{{%3.2f } } }", fSuffix.c_str(), a->index, a->iso1) << endl;
    OUT << "chi2dof " << a->chi2dof << endl;
    fTEX <<  Form("\\vdef{%s:chi2dof:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->chi2dof) << endl;
    OUT << "alpha   " << a->alpha << endl;
    fTEX <<  Form("\\vdef{%s:alpha:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->alpha) << endl;
    OUT << "fls3d   " << a->fls3d << endl;
    fTEX <<  Form("\\vdef{%s:fls3d:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->fls3d) << endl;
    OUT << "docatrk   " << a->docatrk << endl;
    fTEX <<  Form("\\vdef{%s:docatrk:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->docatrk) << endl;
  }
  OUT.flush();
}


// ----------------------------------------------------------------------
void anaBmm::initNumbers(numbers *a) {

  a->name = "";
  a->effGenFilter = a->effGenFilterE = 1.;
  a->fitYield = a->fitYieldE = 0.;
  a->genFileYield = a->genYield = a->recoYield = a->muidYield = a->trigYield = a->candYield = a->ana0Yield = a->anaYield = a->anaWmcYield = 0; 
  a->acc = a->accE = 0; 
  a->effMuidMC =  a->effMuidMCE = a->effTrigMC = a->effTrigMCE = 0; 
  a->effMuidPid = a->effMuidPidE = a->effTrigPid = a->effTrigPidE = 0; 
  a->effCand = a->effCandE = 0; 
  a->effAna = a->effAnaE = 0; 
  // -- this is only relevant for the signal(s)
  a->pss   = a->pdd = 1.;
  a->psd   = a->pds = 1.;
  a->bgObs = a->bgBsExp = a->bgBsExpE = a->bgBdExp = a->bgBdExpE =0; 
  a->bsObs = a->bdObs = 0; 
  a->mBdLo = a->mBdHi = a->mBsLo = a->mBsHi = 0.;
}


// ----------------------------------------------------------------------
int anaBmm::detChan(double m1eta, double m2eta) {
  // -- simple two channel analysis: channel 0 if both muons in barrel, channel 1 else
  if (TMath::Abs(m1eta) < fCuts[0]->etaMax && TMath::Abs(m2eta) < fCuts[0]->etaMax) return 0; 
  if (TMath::Abs(m1eta) < 2.4 && TMath::Abs(m2eta) < 2.4) return 1; 
  return -1; 
}
