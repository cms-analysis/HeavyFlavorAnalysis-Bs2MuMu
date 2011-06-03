#include "anaBmm.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
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
#include "TMath.h"
#include "TVirtualPad.h"  // access to gPad
#include "TMinuit.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"

#include <iomanip>
#include <string>

// Usage: root[0] anaBmm(

using namespace std; 
using std::string; 

ClassImp(anaBmm)

// ----------------------------------------------------------------------
anaBmm::anaBmm(const char *files, const char *dir, int mode) { 

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

  init(files, dir, mode);
}

// ----------------------------------------------------------------------
anaBmm::~anaBmm() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void anaBmm::init(const char *files, const char *dir, int mode) {

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
  cout << "Reading cuts from " << Form("anaBmm.%s.cuts", dir) << endl;
  readCuts(Form("anaBmm.%s.cuts", dir)); 

  printCuts(cout); 

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
    a->effGenFilter  = 0.63;
    //    a->effGenFilter  = 1.0;
    a-> effGenFilterE = 0.03;
    fNumbersBs.push_back(a); 

    // -- signal Bd2MuMu
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("signal Bd2MuMu %i", i); 
    a->effGenFilter  = 0.63;
    //    a->effGenFilter  = 1.0;
    a-> effGenFilterE = 0.03;
    fNumbersBd.push_back(a); 

    // --  normalization
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name = Form("normalization %i", i); 
    a->effGenFilter  = 0.24;
    //    a->effGenFilter  = 1.0;
    a->effGenFilterE = 0.013;
    fNumbersNorm.push_back(a); 

    // --  control sample
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("control sample %i", i); 
    a->effGenFilter  = 0.22;
    //    a->effGenFilter  = 1.0;
    a-> effGenFilterE = 0.015;
    fNumbersCS.push_back(a); 

    h = new TH1D(Form("hMassWithMassCuts%d", i), Form("hMassWithMassCuts%d", i), NBINS, fMassLo, fMassHi);
    //    h->SetLineColor(kBlue); 
    fhMassWithMassCuts.push_back(h); 

    h = new TH1D(Form("fhMassWithMassCutsManyBins%d", i), Form("fhMassWithMassCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    //    h->SetLineColor(kBlue); 
    fhMassWithMassCutsManyBins.push_back(h); 
  
    h = new TH1D(Form("hMassWithCuts%d", i), Form("hMassWithCuts%d", i), NBINS, fMassLo, fMassHi);
    //    h->SetLineColor(kBlue); 
    fhMassWithCuts.push_back(h); 

    h = new TH1D(Form("hMassWithCutsManyBins%d", i), Form("hMassWithCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    //    h->SetLineColor(kBlue); 
    fhMassWithCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassNoCuts%d", i), Form("hMassNoCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassNoCuts.push_back(h); 

    h = new TH1D(Form("hMassChan%d", i), Form("hMassChan%d", i), 100, 0, 10);
    fhMassChan.push_back(h); 

    h = new TH1D(Form("hMassAcc%d", i), Form("hMassAcc%d", i), 100, 0, 10);
    fhMassAcc.push_back(h); 

    h = new TH1D(Form("hMassAbsNoCuts%d", i), Form("hMassAbsNoCuts%d", i), 100, 0, 10);
    fhMassAbsNoCuts.push_back(h); 

    h = new TH1D(Form("hMuId%d", i), Form("hMuId%d", i), 100, 0., 1.);
    fhMuId.push_back(h); 

    h = new TH1D(Form("hMuTr%d", i), Form("hMuTr%d", i), 100, 0., 1.);
    fhMuTr.push_back(h); 
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
  fSuffix    = dir; 
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
      fpMc[fNMc] = loadFile(sfile, stype); 
      fMcLumi[fNMc] = atof(slumi.c_str()); 
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	fSgMc = fNMc;
	sname = "SgMc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
      }
      if (string::npos != stype.find("sg,2011")) {
	sname = "SgMc11";
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (Spring11)")); 
      }
      if (string::npos != stype.find("sg,2010")) {
	sname = "SgMc10"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (Fall10)")); 
      }
      if (string::npos != stype.find("default") && string::npos != stype.find("bd")) {
	fBdMc = fNMc;
	sname = "BdMc"; 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (Spring11)")); 
      }
      if (string::npos != stype.find("bd,2011")) {
	sname = "BdMc11"; 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (Spring11)")); 
      }
      if (string::npos != stype.find("bd,2010")) {
	sname = "BdMc10"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (Fall10)")); 
      }
      if (string::npos != stype.find("default") && string::npos != stype.find("no")) {
	fNoMc = fNMc;
	sname = "NoMc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #mu^{+}#mu^{-}K^{+} (MC)")); 
      }
      if (string::npos != stype.find("2010") && string::npos != stype.find("no")) {
	sname = "NoMc10"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #mu^{+}#mu^{-}K^{+} (Fall10)")); 
      }
      if (string::npos != stype.find("2011") && string::npos != stype.find("no")) {
	sname = "NoMc11"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #mu^{+}#mu^{-}K^{+} (Spring11)")); 
      }
      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	fCsMc = fNMc;
	sname = "CsMc"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-}K^{+}K^{-} (MC)")); 
      }	
      if (string::npos != stype.find("2010") && string::npos != stype.find("cs")) {
	sname = "CsMc10"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-}K^{+}K^{-} (Fall10)")); 
      }	
      if (string::npos != stype.find("2011") && string::npos != stype.find("cs")) {
	sname = "CsMc11"; 
	fF.insert(make_pair(sname, fpMc[fNMc])); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-}K^{+}K^{-} (Spring11)")); 
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

  dumpSamples();
  fF["SgMc"]->cd();
  dumpCutNames();

  // -- some special plots
  if ((0x1<<0) & channel) {
    plotWithCut("fl3d", "", 2.0, "l_{3d} [cm]", 0., 5);
    plotWithCut("fl3dE", "", -1.0, "#sigma(l_{3d}) [cm]", 0., 0.5);
    
    isoEff("iso4", 0.70, "NoData", "Ao");
    isoEff("iso4", 0.75, "NoData", "Ao");
    isoEff("iso4", 0.80, "NoData", "Ao");
      
    isoEff("iso4", 0.70, "CsData", "Ao");
    isoEff("iso4", 0.75, "CsData", "Ao");
    isoEff("iso4", 0.80, "CsData", "Ao");
      
    isoEff("iso", 0.70, "CsData", "Ao");
    isoEff("iso", 0.75, "CsData", "Ao");
    isoEff("iso", 0.80, "CsData", "Ao");
  }

  if (0) {
    
    cout << "1: " << endl;
    fF["SgMc"]->cd();
    histAcceptanceAndPreselection(*fNumbersBs[0]);
    printNumbers(*fNumbersBs[0], cout); 
    printNumbers(*fNumbersBs[0], fOUT); 
    cout << "2: " << endl;
    fF["SgMc"]->cd();
    accEffFromEffTree(*fNumbersBs[0], *fCuts[0]);
    printNumbers(*fNumbersBs[0], cout); 
    printNumbers(*fNumbersBs[0], fOUT); 
    cout << "3: " << endl;
    fF["BdMc10"]->cd();
    accEffFromEffTree(*fNumbersBd[0], *fCuts[0]);
    printNumbers(*fNumbersBs[0], cout); 
    printNumbers(*fNumbersBs[0], fOUT); 
  }

  // -- 2
  if ((0x1<<1) & channel) {
    sbsDistributionOverlay("SgData", "SgMc", "Ao");
    sbsDistributionOverlay("NoData", "NoMc", "Ao");
    sbsDistributionOverlay("CsData", "CsMc", "Ao");

    sbsDistributionOverlay("NoData", "NoMc", "Presel");
    sbsDistributionOverlay("CsData", "CsMc", "Presel");

    sbsDistributionOverlay("NoData", "NoMc", "Nm");
    sbsDistributionOverlay("CsData", "CsMc", "Nm");

    //     sbsDistributionOverlay("NoData", "NoMc", "HLT");
    //     sbsDistributionOverlay("CsData", "CsMc", "HLT");

    //??  sbsDistributionOverlay("NoData", "NoData11", "Ao");
    //??  sbsDistributionOverlay("CsData", "CsData2011", "Ao");
    sbsDistributionOverlay("NoData", "CsData", "Ao");
    sbsDistributionOverlay("NoMc",   "CsMc", "Ao");
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

  ngen.insert(make_pair("SgMc10", -1)); 
  ngen.insert(make_pair("SgMc11", 100.e6)); 

  ngen.insert(make_pair("BdMc11", 40.e6)); 
  ngen.insert(make_pair("BdMc10", -1)); 

  ngen.insert(make_pair("NoMc11", 160.e6)); 
  ngen.insert(make_pair("CsMc11", 40.e6)); 
  
  //  ofstream OUT(fNumbersFileName.c_str(), ios::app);

  fTEX << "% ----------------------------------------------------------------------" << endl;
  string name; 
  double lumi, n, f; 
  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    cout << "===> " << imap->first << " -> " << imap->second << endl;
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
    fTEX <<  Form("\\vdef{%s:%s:cutLine}   {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
    fTEX <<  Form("\\vdef{%s:%s:cutValue}  {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
  }
  fTEX.flush();
}


// ----------------------------------------------------------------------
void anaBmm::allEffTables() {
  // -- BMM
  effTable("SgMc");
  effTable("SgData");

  // -- B+
  effTable("NoData");
  effTable("NoData10");
  effTable("NoMc");

  // -- Bs2JpsiPhi
  effTable("CsData");
  effTable("CsData10");
  effTable("CsMc");

  fTEX.flush();
  // -- Compute difference
  double deltaIso     = computeDelta(Form("{%s:iso1eff:NoData}", fSuffix.c_str()), Form("{%s:iso1eff:NoMc}", fSuffix.c_str()));
  double deltaAlpha   = computeDelta(Form("{%s:alphaeff:NoData}", fSuffix.c_str()), Form("{%s:alphaeff:NoMc}", fSuffix.c_str()));
  double deltaFls     = computeDelta(Form("{%s:fls3deff:NoData}", fSuffix.c_str()), Form("{%s:fls3deff:NoMc}", fSuffix.c_str()));
  double deltaChi2dof = computeDelta(Form("{%s:chi2dofeff:NoData}", fSuffix.c_str()), Form("{%s:chi2dofeff:NoMc}", fSuffix.c_str()));
  double deltaTot     = TMath::Sqrt(deltaIso*deltaIso + deltaAlpha*deltaAlpha + deltaFls*deltaFls + deltaChi2dof*deltaChi2dof); 
  double deltaTotNoIso= TMath::Sqrt(deltaAlpha*deltaAlpha + deltaFls*deltaFls + deltaChi2dof*deltaChi2dof); 

  fTEX << Form("%s", (formatTex(deltaIso, fSuffix+":deltaIso", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaAlpha, fSuffix+":deltaAlpha", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaFls, fSuffix+":deltaFls", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaChi2dof, fSuffix+":deltaChi2dof", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(100*deltaTot, fSuffix+":deltaTot", "No")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(100*deltaTotNoIso, fSuffix+":deltaTotNoIso", "No")).c_str()) << endl;


  deltaIso     = computeDelta(Form("{%s:iso1eff:CsData}", fSuffix.c_str()), Form("{%s:iso1eff:CsMc}", fSuffix.c_str()));
  deltaAlpha   = computeDelta(Form("{%s:alphaeff:CsData}", fSuffix.c_str()), Form("{%s:alphaeff:CsMc}", fSuffix.c_str()));
  deltaFls     = computeDelta(Form("{%s:fls3deff:CsData}", fSuffix.c_str()), Form("{%s:fls3deff:CsMc}", fSuffix.c_str()));
  deltaChi2dof = computeDelta(Form("{%s:chi2dofeff:CsData}", fSuffix.c_str()), Form("{%s:chi2dofeff:CsMc}", fSuffix.c_str()));
  deltaTot     = TMath::Sqrt(deltaIso*deltaIso + deltaAlpha*deltaAlpha + deltaFls*deltaFls + deltaChi2dof*deltaChi2dof); 
  deltaTotNoIso= TMath::Sqrt(deltaAlpha*deltaAlpha + deltaFls*deltaFls + deltaChi2dof*deltaChi2dof); 
  
  fTEX << Form("%s", (formatTex(deltaIso, fSuffix+":deltaIso", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaAlpha, fSuffix+":deltaAlpha", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaFls, fSuffix+":deltaFls", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(deltaChi2dof, fSuffix+":deltaChi2dof", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(100*deltaTot, fSuffix+":deltaTot", "Cs")).c_str()) << endl;
  fTEX << Form("%s", (formatTex(100*deltaTotNoIso, fSuffix+":deltaTotNoIso", "Cs")).c_str()) << endl;
  
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
void anaBmm::effTable(string smode) {
  gStyle->SetOptTitle(1); 
  gStyle->SetOptFit(1112); 
  tl->SetTextSize(0.05); 
  fF[smode]->cd();
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
    massSigma = 0.030;
  }

  if (string::npos != smode.find("Cs")) {
    massLo = 5.15; 
    massHi = 5.6;
    massPeak = 5.37;
    massSigma = 0.030;
  }
  
  //  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << " % --- " << smode << endl;
  int mode(0); 
  if (string::npos != smode.find("SgMc")) {
    cout << "==> This is signal MC, counting events in histogram" << endl;
    mode = 0; 
  }
  if (string::npos != smode.find("SgData")) {
    cout << "==> This is signal Data, counting events in blind signal box from sidebands" << endl;
    mode = 1; 
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
  AnalysisDistribution *an = new AnalysisDistribution("docatrk");
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
  c0->SaveAs(pdfname.c_str(), "Portrait");

  delete an;
  
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i)); 
    
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
    tl->DrawLatex(0.22, 0.75, Form("%4.2f+/-%4.2f", n, nE)); 
    pdfname = Form("%s/%s_%s_%s_hMassAo.pdf", fDirectory.c_str(), fSuffix.c_str(), smode.c_str(), cut.c_str());
    cout << "AD for " << cut << " results in " << pdfname << endl;
    c0->SaveAs(pdfname.c_str(), "Portrait");

    delete a; 

    if ((string::npos != cut.find("tracks")) || (string::npos != cut.find("muons"))) {
      n *=0.5; 
      nE *=0.5; 
    }

    if (nE < normE) nE = normE;
    eff    = norm/n;
    if (eff > 1.0) eff = 1.0; 
    effE   = dEff(norm, normE, n, nE); 

    //     relEff = n/nprev;
    //     relEffE= dEff(n, nE, nprev, nprevE);
    //     cumEff = n/norm;
    //     cumEffE= dEff(n, nE, norm, normE);
    
    cout << cut << " n = " << n << "+/-" << nE 
	 << " eff = " << eff << "+/-" << effE
      // 	 << " rel eff = " << relEff << "+/-" << relEffE
      // 	 << " cum eff = " << cumEff << "+/-" << cumEffE
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

  c0->SaveAs(Form("%s-npv.pdf", var)); 

//   setHist(h3, kRed, 25); h3->Draw("same");
//   setHist(h4, kRed, 21); h4->Draw("same");
}


// ----------------------------------------------------------------------
void anaBmm::isoProcess(const char *var, const char *cuts) {

  // -- signal MC
  fpMc[fSgMc]->cd(); 
  TH1D *h40 = (TH1D*)gFile->Get("iso40");  
  if (0 == h40) {h40 = new TH1D("iso40", "", 40, 0., 1.01); } else  h40->Reset();

  TH1D *h41 = (TH1D*)gFile->Get("iso41");  
  if (0 == h41) {h41 = new TH1D("iso41", "", 40, 0., 1.01); } else  h41->Reset();

  TH1D *h42 = (TH1D*)gFile->Get("iso42");  
  if (0 == h42) {h42 = new TH1D("iso42", "", 40, 0., 1.01); } else  h42->Reset();


  TTree *t; 
  t = (TTree*)gFile->Get("events"); 
  cout << "*********************" << endl;
  cout << gFile->GetName() << " -> t = " << t << endl;
  t->Draw(Form("%s>>iso40", var), Form("%s&&procid==%d", cuts, 40));
  cout << "2" << endl;
  t->Draw(Form("%s>>iso41", var), Form("%s&&procid==%d", cuts, 41));
  t->Draw(Form("%s>>iso42", var), Form("%s&&procid==%d", cuts, 42));
  
  h41->Scale(h40->GetSumOfWeights()/h41->GetSumOfWeights()); 
  h42->Scale(h40->GetSumOfWeights()/h42->GetSumOfWeights()); 

  setTitles(h40, "I", "a.u.", 0.06, 0.9, 1.3);
  setHist(h40, kBlue);  
  setHist(h41, kRed);     //h41->SetLineStyle(kDashed); 
  setHist(h42, kBlack);   //h42->SetLineStyle(kDotted); 

  shrinkPad(0.1, 0.15); 
  double maxi = h40->GetMaximum(); 
  if (h41->GetMaximum() > maxi) maxi = h41->GetMaximum(); 
  if (h42->GetMaximum() > maxi) maxi = h42->GetMaximum(); 
  h40->SetMaximum(1.1*maxi); 
  h40->Draw();
  h41->Draw("same");
  h42->Draw("same");


  newLegend(0.25, 0.6, 0.4, 0.75); 
  legg->AddEntry(h40, Form("GGF, mean = %4.3f", h40->GetMean()), "l"); 
  legg->AddEntry(h41, Form("FEX, mean = %4.3f", h41->GetMean()), "l"); 
  legg->AddEntry(h42, Form("GSP, mean = %4.3f", h42->GetMean()), "l"); 
  legg->Draw();


  string pdfname = Form("%s/%s_isoProcess.pdf", fDirectory.c_str(), fSuffix.c_str());
  if (c0) c0->SaveAs(pdfname.c_str());

  loopTree(0, 40); 
  double ggf0 = fNumbersBs[0]->effAna; 
  double ggf1 = fNumbersBs[1]->effAna;
  loopTree(0, 41); 
  double fex0 = fNumbersBs[0]->effAna;
  double fex1 = fNumbersBs[1]->effAna;
  loopTree(0, 42); 
  double gsp0 = fNumbersBs[0]->effAna;
  double gsp1 = fNumbersBs[1]->effAna;

  loopTree(0); 
  double all0 = fNumbersBs[0]->effAna;
  double all1 = fNumbersBs[1]->effAna;

  cout << "GGF: " << ggf0  << " " << ggf1 << endl;
  cout << "FEX: " << fex0  << " " << fex1 << endl;
  cout << "GSP: " << gsp0  << " " << gsp1 << endl;
  cout << "ALL: " << all0  << " " << all1 << endl;
  
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
  double epsPi = 0.003; 
  double epsKa = 0.004; 
  double epsPr = 0.001;
  colors.insert(make_pair("bg60", 46)); hatches.insert(make_pair("bg60", 3004)); mscale.insert(make_pair("bg60", epsPi*epsPr)); 
  colors.insert(make_pair("bg61", 49)); hatches.insert(make_pair("bg61", 3005)); mscale.insert(make_pair("bg61", epsKa*epsPr)); 

  colors.insert(make_pair("bg82", 30)); hatches.insert(make_pair("bg82", 3004)); mscale.insert(make_pair("bg82", epsKa*epsKa)); 
  colors.insert(make_pair("bg83", 32)); hatches.insert(make_pair("bg83", 3005)); mscale.insert(make_pair("bg83", epsPi*epsKa)); 
  colors.insert(make_pair("bg84", 33)); hatches.insert(make_pair("bg84", 3007)); mscale.insert(make_pair("bg84", epsPi*epsPi)); 

  colors.insert(make_pair("bg93", 40)); hatches.insert(make_pair("bg93", 3004)); mscale.insert(make_pair("bg93", epsKa*epsKa)); 
  colors.insert(make_pair("bg92", 41)); hatches.insert(make_pair("bg92", 3005)); mscale.insert(make_pair("bg92", epsKa*epsPi)); 
  colors.insert(make_pair("bg91", 42)); hatches.insert(make_pair("bg91", 3007)); mscale.insert(make_pair("bg91", epsPi*epsPi)); 

  newLegend(0.65, 0.4, 0.80, 0.85); 

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
    TH1D *h1Rare0 = (TH1D*)(fhMassWithCuts[0]->Clone("h1Rare0"));  h1Rare0->SetLineColor(kBlack); 
    TH1D *h1Rare1 = (TH1D*)(fhMassWithCuts[1]->Clone("h1Rare1"));  h1Rare1->SetLineColor(kBlack);

    double bd0    = fhMassWithCutsManyBins[0]->Integral(fhMassWithCutsManyBins[0]->FindBin(fCuts[0]->mBdLo)+1, 
							fhMassWithCutsManyBins[0]->FindBin(fCuts[0]->mBdHi)-1);
    double bs0    = fhMassWithCutsManyBins[0]->Integral(fhMassWithCutsManyBins[0]->FindBin(fCuts[0]->mBsLo)+1, 
							fhMassWithCutsManyBins[0]->FindBin(fCuts[0]->mBsHi)-1);

    double bd1    = fhMassWithCutsManyBins[1]->Integral(fhMassWithCutsManyBins[1]->FindBin(fCuts[1]->mBdLo)+1, 
							fhMassWithCutsManyBins[1]->FindBin(fCuts[1]->mBdHi)-1);
    double bs1    = fhMassWithCutsManyBins[1]->Integral(fhMassWithCutsManyBins[1]->FindBin(fCuts[1]->mBsLo)+1, 
							fhMassWithCutsManyBins[1]->FindBin(fCuts[1]->mBsHi)-1);

    bsRare0 += bs0*lscale*misid; 
    bsRare1 += bs1*lscale*misid; 

    bdRare0 += bd0*lscale*misid; 
    bdRare1 += bd1*lscale*misid; 

    cout << imap->first << ": bsRare0 increment by " << bs0*lscale*misid << "  " << bs0 << endl;
    cout << imap->first << ": bsRare1 increment by " << bs1*lscale*misid << "  " << bs1 << endl;
    cout << imap->first << ": bdRare0 increment by " << bd0*lscale*misid << "  " << bd0 << endl;
    cout << imap->first << ": bdRare1 increment by " << bd1*lscale*misid << "  " << bd1 << endl;

    
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
  if (c0) c0->SaveAs(pdfname.c_str());
  
  hRareBg1->Draw();
  TH1D *hhRareBg1 = (TH1D*)hRareBg1->GetHistogram(); 
  setTitles(hhRareBg1, "m_{h h} [GeV]", "Cands/bin", 0.06, 0.9, 1.5);
  legg->Draw(); 
  hhRareBg1->Draw("same");
  pdfname = Form("%s/%s_rare1.pdf", fDirectory.c_str(), fSuffix.c_str());
  if (c0) c0->SaveAs(pdfname.c_str());

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
void anaBmm::sbsDistributionOverlay(std::string file1, std::string file2, const char *selection) {

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
  AnalysisDistribution a("allevents"); 

  vector<string> skipList; 
  skipList.push_back("mpsi"); 
  skipList.push_back("q"); 
  skipList.push_back("allevents"); 
  skipList.push_back("muonsid"); 
  skipList.push_back("tracksqual"); 
  skipList.push_back("hlt"); 
  skipList.push_back("allevents"); 


  vector<string> nolegendList; 
  nolegendList.push_back("muonseta"); 
  nolegendList.push_back("muon1eta"); 
  nolegendList.push_back("muon2eta"); 
  nolegendList.push_back("trackseta"); 
  nolegendList.push_back("zpv"); 
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
  string cut, pdfname; 
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i)); 
    
    if (skipList.end() != find(skipList.begin(), skipList.end(), cut)) continue;
    
    if (cut == string("")) {
      cout << "empty string found, break ..." << endl;
      //      OUT.close();
      break;
    }
    pdfname = Form("%s/%s_%s-%s_sbs_%s_%s.pdf", fDirectory.c_str(), fSuffix.c_str(),
		   file1.c_str(), file2.c_str(), cut.c_str(), selection);
  
    fF[file1]->cd(); 
    cout << "==> File1: "; gFile->pwd(); 
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
    if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file1.c_str(), cut.c_str(), selection));
    
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
    if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file2.c_str(), cut.c_str(), selection));

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
    c0->SaveAs(pdfname.c_str()); 
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
void anaBmm::histAcceptanceAndPreselection(numbers &a) {

  // -- Acceptance and MC efficiency numbers based on prefilled histogram
  //    (could be changed to use effTree!)
  TH1D *hacc  = (TH1D*)(gFile->Get("efficiency"));
  if (!hacc) {
    cout << "anaBmm::histAcceptanceAndPreselection(" << a.name << "): no histogram `efficiency' found " << endl;
    return;
  }
  a.genFileYield  = hacc->GetBinContent(2);
  a.genYield  = a.genFileYield/a.effGenFilter;
  a.recoYield = hacc->GetBinContent(5);
  a.muidYield = hacc->GetBinContent(6);
  a.trigYield = hacc->GetBinContent(7);
  a.candYield = hacc->GetBinContent(8);
  
  if (a.genYield > 0) {
    a.acc = a.recoYield/a.genYield;
    a.accE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genYield));
  }

  if (a.recoYield > 0) {
    a.effMuidMC = a.muidYield/a.recoYield;
    a.effMuidMCE = dEff(static_cast<int>(a.muidYield), static_cast<int>(a.recoYield));
  } 

  if (a.muidYield > 0) {
    a.effTrigMC = a.trigYield/a.muidYield;
    a.effTrigMCE = dEff(static_cast<int>(a.trigYield), static_cast<int>(a.muidYield));
  } 

  if (a.trigYield > 0) {
    a.effCand = a.candYield/a.trigYield;
    a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.trigYield));
  } 

}


// ----------------------------------------------------------------------
void anaBmm::accEffFromEffTree(numbers &a, cuts &b) {

  TTree *t  = (TTree*)(gFile->Get("effTree"));
  if (!t) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): no tree `effTree' found " << endl;
    gFile->pwd(); 
    return;
  }

  bool sg(false), no(false), cs(false); 

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
    ++ngen;
    if (sg) {
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1.) {
	    if (bm1pt > 1. && bm2pt > 1.
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4
		&& bm1gt && bm2gt
		) {
	      ++nreco;
	      chan = detChan(bm1eta, bm2eta); 
	      if (bm1pt > b.m1pt && bm2pt > b.m2pt 
		  && chan == a.index
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
    } else if (no) {
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5 && TMath::Abs(bg3eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1. && bg3pt > 0.4) {
	    if (bm1pt > 1. && bm2pt > 1. && bk1pt > 0.4
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4 && TMath::Abs(bk1eta) < 2.4
		&& bm1gt && bm2gt && bk1gt
		) {
	      ++nreco;
	      chan = detChan(bm1eta, bm2eta); 
	      if (bm1pt > b.m1pt && bm2pt > b.m2pt
		  && chan == a.index
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
    } else if (cs) {
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1. && bg3pt > 1. && bg4pt > 1.) {
	    if (bm1pt > 1. && bm2pt > 1. && bk1pt > 0.4 && bk2pt > 0.4
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4 && TMath::Abs(bk1eta) < 2.4 && TMath::Abs(bk2eta) < 2.4
		&& bm1gt && bm2gt && bk1gt && bk2gt
		) {
	      ++nreco;
	      if (bm1pt > b.m1pt && bm2pt > b.m2pt
		  && chan == a.index
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
  }

  a.genFileYield  = ngen;
  a.genYield      = a.genFileYield/a.effGenFilter;
  a.genChanYield  = nchangen; 
  a.recoYield     = nreco;
  a.chanYield     = nchan;
  a.muidYield     = nmuid;
  a.trigYield     = nhlt;
  a.candYield     = ncand;
  
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

  if (a.chanYield > 0) {
    a.effMuidMC = a.muidYield/a.chanYield;
    a.effMuidMCE = dEff(static_cast<int>(a.muidYield), static_cast<int>(a.chanYield));
  } 

  if (a.muidYield > 0) {
    a.effTrigMC = a.trigYield/a.muidYield;
    a.effTrigMCE = dEff(static_cast<int>(a.trigYield), static_cast<int>(a.muidYield));
  } 

  if (a.trigYield > 0) {
    a.effCand = a.candYield/a.trigYield;
    a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.trigYield));
  } 
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

  system(Form("../ulcalc/bin/ulcalc %s", fUlcalcFileName.c_str())); 

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
  
  double result = (fCsSig/fNormSig)
    *(fu/fs)
    *(fNumbersNorm[0]->acc/fNumbersCS[0]->acc)
    *(fNumbersNorm[0]->effCand/fNumbersCS[0]->effCand)     
    *(fNumbersNorm[0]->effMuidMC/fNumbersCS[0]->effMuidMC)
    *(fNumbersNorm[0]->effTrigMC/fNumbersCS[0]->effTrigMC)
    *(fNumbersNorm[0]->effAna/fNumbersCS[0]->effAna)
    * fBF;

  cout << "fact branching fraction: " << result << endl;

  result = (fCsSig/fNormSig)
    *(fu/fs)
    *(fNumbersNorm[0]->effTot/fNumbersCS[0]->effTot)
    * fBF;

  cout << "branching fraction: " << result << endl;

}

// ----------------------------------------------------------------------
void anaBmm::effTree(int mode) {

}


// ----------------------------------------------------------------------
void anaBmm::isoEff(const char *var, double cut, const char *file1, const char *selection) {

  TH1D *h[6];
  
  fF[file1]->cd(); 
  AnalysisDistribution a("allevents"); 
  cout << "==> File1: "; gFile->pwd(); 
  c0->Clear();
  c0->Divide(2,3);
  for (int i = 0; i < 6; ++i) {
    h[i] = a.sbsDistribution(Form("%spv%d", var, i+1), selection);
    c0->cd(i+1);
    h[i]->Draw();
  }

  TH1D *heff = new TH1D("heff", "", 6, 0., 12.);
  for (int i = 0; i < 6; ++i) {
    double ntot = h[i]->GetSumOfWeights();
    double ncut = h[i]->Integral(h[i]->FindBin(cut)+1, h[i]->GetNbinsX());
    double eff  = ncut/ntot; 
    double effE = dEff(static_cast<int>(ncut), static_cast<int>(ntot));
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
  setTitles(heff, "N_{PV}", Form("#epsilon(I>%3.2f)", cut)); 
  heff->SetMaximum(1.2);
  if (heff->GetMinimum() > 0.5) {
    heff->SetMinimum(0.5);
  } else {
    heff->SetMinimum(0.);
  }
  heff->SetMarkerSize(2); 
  heff->Fit("pol1");
  heff->GetFunction("pol1")->SetLineWidth(2);
  
  c0->SaveAs(Form("%s/isoeff-%s-%s-0_%d.pdf", fDirectory.c_str(), var, file1, static_cast<int>(100.*cut)));

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
  
  c0->SaveAs(Form("%s/plotWithCut-%s.pdf", fDirectory.c_str(), var));
  c0->SetLogy(0); 
}


// ----------------------------------------------------------------------
TH1* anaBmm::loopTree(int mode, int proc) {
  // -- mode definition
  // 0  Bs2MuMu MC
  // 1  Bd2MuMu MC
  // 5  Bs2MuMu data
  // 10 Bp2JpsiKp MC
  // 11 Bp2JpsiKp data
  // 20 Bs2JpsiPhi MC
  // 21 Bs2JpsiPhi data

  bool bp2jpsikp(false), bs2jpsiphi(false), isMC(false); 

  cout << "--> loopTree with mode " << mode << endl;

  numbers *aa(0);
  if (0 == mode) {
    isMC = true; 
    fpMc[fSgMc]->cd(); 
  }  else if (1 == mode) {
    isMC = true; 
    fpMc[fBdMc]->cd(); 
  } else if (5 == mode) {
    isMC = false;     
    fpData[fSgData]->cd(); 
  } else if (10 == mode) {
    isMC = true; 
    bp2jpsikp = true; 
    fpMc[fNoMc]->cd(); 
  } else if (11 == mode) {
    isMC = false;     
    fpData[fNoData]->cd(); 
  } else if (20 == mode) {
    isMC = true; 
    bs2jpsiphi = true; 
    fpMc[fCsMc]->cd(); 
    //    aa = fNumbersCS[0];
  } else if (21 == mode) {
    isMC = false;     
    bs2jpsiphi = true; 
    fpData[fCsData]->cd(); 
  } else {
    cout << "mode 99" << endl;
  }

  gFile->pwd();

  // -- reset all histograms
  for (unsigned int i = 0; i < fNchan; ++i) {
    fhMassAbsNoCuts[i]->Reset();
    fhMuId[i]->Reset();
    fhMassNoCuts[i]->Reset();
    fhMuTr[i]->Reset();
    fhMassWithCuts[i]->Reset();
    fhMassWithCutsManyBins[i]->Reset();
    fhMassWithMassCutsManyBins[i]->Reset();
    fhMassWithMassCuts[i]->Reset();
  }

  // -- set up tree
  TTree *t;
  t = (TTree*)gFile->Get("events");
  int brun, bevt, btm, bq1, bq2, bprocid; 
  double bg1pt, bg2pt, bg1eta, bg2eta;
  double bm, bpt, beta, bcosa, biso1, bchi2, bdof, bdocatrk, bfls3d, bfl3dE, bfl3d, bm1pt, bm1eta, bm2pt, bm2eta;
  double bg3pt, bg3eta, bg4pt, bg4eta; 
  double bmkk, bdr;
  double bw8mu, bw8tr;
  bool bhlt, bgmuid, bgtqual;
  t->SetBranchAddress("run",&brun);
  t->SetBranchAddress("evt",&bevt);
  t->SetBranchAddress("hlt",&bhlt);
  t->SetBranchAddress("gmuid",&bgmuid);
  t->SetBranchAddress("gtqual",&bgtqual);
  t->SetBranchAddress("w8mu",&bw8mu);
  t->SetBranchAddress("w8tr",&bw8tr);
  t->SetBranchAddress("tm",&btm);
  t->SetBranchAddress("procid",&bprocid);
  t->SetBranchAddress("m",&bm);
  t->SetBranchAddress("pt",&bpt);
  t->SetBranchAddress("eta",&beta);
  t->SetBranchAddress("cosa",&bcosa);
  t->SetBranchAddress("iso1",&biso1);
  t->SetBranchAddress("chi2",&bchi2);
  t->SetBranchAddress("dof",&bdof);
  t->SetBranchAddress("fls3d",&bfls3d);
  t->SetBranchAddress("fl3d",&bfl3d);
  t->SetBranchAddress("fl3dE",&bfl3dE);
  t->SetBranchAddress("m1pt",&bm1pt);
  t->SetBranchAddress("m1eta",&bm1eta);
  t->SetBranchAddress("m1q",&bq1);
  t->SetBranchAddress("m2pt",&bm2pt);
  t->SetBranchAddress("m2eta",&bm2eta);
  t->SetBranchAddress("m2q",&bq2);
  t->SetBranchAddress("docatrk",&bdocatrk);

  t->SetBranchAddress("g1pt",&bg1pt);
  t->SetBranchAddress("g2pt",&bg2pt);
  t->SetBranchAddress("g1eta",&bg1eta);
  t->SetBranchAddress("g2eta",&bg2eta);
  if (bp2jpsikp) {
    t->SetBranchAddress("g3pt",&bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
  }


  if (bs2jpsiphi) {
    t->SetBranchAddress("g3pt",&bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
    t->SetBranchAddress("g4pt",&bg4pt);
    t->SetBranchAddress("g4eta",&bg4eta);
    t->SetBranchAddress("mkk",&bmkk);
    t->SetBranchAddress("dr",&bdr);
  } else {
    bmkk = 999.;
    bdr = 999.;
  }

  int nentries = Int_t(t->GetEntries());
  int nb(0); 
  cuts *pCuts(0); 
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    // -- require truth matching when on MC
    if (0 == mode && 0 == btm) continue;
    if (1 == mode && 0 == btm) continue;
    if (10 == mode && 0 == btm) continue;
    if (20 == mode && 0 == btm) continue;

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

    fhMassAbsNoCuts[fChan]->Fill(bm);
    // -- require wide mass window
    if (bm < fMassLo) continue;
    if (fMassHi < bm) continue;

    // -- gen-level acceptance cuts
    if (isMC) {
      if (proc > 0 && bprocid != proc) continue;
      if (TMath::Abs(bg1eta) > 2.5) continue;
      if (TMath::Abs(bg2eta) > 2.5) continue;
      if (TMath::Abs(bg1pt) < 1.0) continue;
      if (TMath::Abs(bg2pt) < 1.0) continue;
      if (bp2jpsikp) {
	if (TMath::Abs(bg3eta) > 2.5) continue;
	if (TMath::Abs(bg3pt) < 0.4) continue;
      }
      
      if (bs2jpsiphi) {
	if (TMath::Abs(bg3eta) > 2.5) continue;
	if (TMath::Abs(bg3pt) < 0.4) continue;
	if (TMath::Abs(bg4eta) > 2.5) continue;
	if (TMath::Abs(bg4pt) < 0.4) continue;
      }
    }

    fhMassAcc[fChan]->Fill(bm);

    // -- immutable cuts: require basic muon and trackQual cuts
    if (TMath::Abs(bm1eta) > 2.4) continue;
    if (TMath::Abs(bm2eta) > 2.4) continue;
    if (false == bgtqual) continue;
    if (bq1*bq2 > 0) continue;

    if (bm1pt < pCuts->m1pt) continue; 
    if (bm2pt < pCuts->m2pt) continue; 

    if (bfl3d > 2) continue;

    // -- Channel 
    fhMassChan[fChan]->Fill(bm);

    
    // -- must fill this BEFORE the trigger requirement!
    if (bw8mu > 0.) {
      fhMuId[fChan]->Fill(bw8mu, 1./bw8mu); 
    }
    
    // -- now check for muon ID and trigger
    if (false == bhlt) continue;
    if (false == bgmuid) continue;

    fhMassNoCuts[fChan]->Fill(bm);

    // -- weights for trigger
    if (bw8tr > 0.) {
      fhMuTr[fChan]->Fill(bw8tr, 1./bw8tr); 
    }

    // -- apply analysis cand selection 
    if (bpt < pCuts->pt) continue; 
    if (biso1 < pCuts->iso1) continue; 
    if (bchi2/bdof > pCuts->chi2dof) continue;
    if (TMath::IsNaN(bfls3d)) continue;
    if (bfls3d < pCuts->fls3d) continue;
    if (TMath::ACos(bcosa) > pCuts->alpha) continue;

    if (bs2jpsiphi && bdr >0.3) continue;
    if (bs2jpsiphi && bmkk < 0.995) continue;
    if (bs2jpsiphi && bmkk > 1.045) continue;

    fhMassWithCuts[fChan]->Fill(bm); 
    fhMassWithCutsManyBins[fChan]->Fill(bm); 
    
    if (0 == mode && bm < pCuts->mBsLo) continue;
    if (0 == mode && bm > pCuts->mBsHi) continue;
    if (1 == mode && bm < pCuts->mBdLo) continue;
    if (1 == mode && bm > pCuts->mBdHi) continue;
    if (10 == mode && bm < fNormLo) continue;
    if (10 == mode && bm > fNormHi) continue;
    if (20 == mode && bm < fCsLo) continue;
    if (20 == mode && bm > fCsHi) continue;

    fhMassWithMassCuts[fChan]->Fill(bm);
    fhMassWithMassCutsManyBins[fChan]->Fill(bm); 

    if (5 == mode && bm > 4.8 && bm < 6.0) {
      cout << Form("m = %4.3f pT = %4.3f eta = %4.3f", bm, bpt, beta)
	//	   <<	" run = " << brun << " event = " << bevt
	   << " chan = " << fChan 
	   << Form(" mpt = %4.3f,%4.3f", bm1pt, bm2pt)
	   << Form(" meta = %4.3f,%4.3f", TMath::Abs(bm1eta), TMath::Abs(bm2eta))
	   << Form(" a = %4.3f iso = %4.3f chi2 = %4.3f fls3d = %4.3f, fl/E=%4.3f/%4.3f", 
		   TMath::ACos(bcosa), biso1, bchi2/bdof, bfls3d, bfl3d, bfl3dE)
	   << endl;
    }
    
  }

  if (99 == mode) return fhMassWithCuts[0];


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
    double tot   = fhMassWithCutsManyBins[i]->GetSumOfWeights();
    double bd    = fhMassWithCutsManyBins[i]->Integral(fhMassWithCutsManyBins[i]->FindBin(pCuts->mBdLo), 
						       fhMassWithCutsManyBins[i]->FindBin(pCuts->mBdHi));
    double bs    = fhMassWithCutsManyBins[i]->Integral(fhMassWithCutsManyBins[i]->FindBin(pCuts->mBsLo), 
						       fhMassWithCutsManyBins[i]->FindBin(pCuts->mBsHi));
    
    fhMuId[i]->Draw();
    c0->SaveAs(Form("%s/muid-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMuTr[i]->Draw();
    c0->SaveAs(Form("%s/mutr-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithCutsManyBins[i]->Draw();
    c0->SaveAs(Form("%s/nmc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithMassCutsManyBins[i]->Draw();
    c0->SaveAs(Form("%s/wmc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    // -- Efficiency and acceptance
    if (isMC) {
      accEffFromEffTree(*aa, *fCuts[i]);
      double a = fhMassNoCuts[i]->GetSumOfWeights(); 
      double b = fhMassWithCuts[i]->GetSumOfWeights();
      aa->ana0Yield   = a;
      aa->anaYield    = b; 
      aa->anaWmcYield = fhMassWithMassCuts[i]->GetSumOfWeights(); // "with mass cut"
      aa->effAna      = b/a;
      aa->effAnaE     = dEff(static_cast<int>(b), static_cast<int>(a));
      aa->effMuidPid  = fhMuId[i]->GetMean();
      aa->effMuidPidE = fhMuId[i]->GetMeanError();
      aa->effTrigPid  = fhMuTr[i]->GetMean();
      aa->effTrigPidE = fhMuTr[i]->GetMeanError();
      aa->effTot      = b/(aa->genYield);
      aa->effTotE     = dEff(static_cast<int>(b), static_cast<int>(aa->genYield));
      aa->effTotChan  = b/(aa->genChanYield);
      aa->effTotChanE = dEff(static_cast<int>(b), static_cast<int>(aa->genChanYield));
      aa->combGenYield= b/(aa->effTot);
      aa->prodGenYield= b/(aa->acc * aa->effChan * aa->effMuidMC * aa->effTrigMC * aa->effCand * aa->effAna); 
      aa->chanGenYield= b/(aa->accChan * aa->effChan * aa->effMuidMC * aa->effTrigMC * aa->effCand * aa->effAna); 
    }


    cout << "bs:  " << bs << endl;
    cout << "bd:  " << bd << endl;
    cout << "tot: " << tot << endl;
    cout << "wmc: " << aa->anaWmcYield << endl;

    // -- compute PSD, PDS, etc
    if (0 == mode) {
      aa->pss = bs/tot;
      aa->pds = bd/tot;
    } 
    
    if (1 == mode) {
      aa->pdd = bd/tot;
      aa->psd = bs/tot;
    } 

    if (0 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC SIGNAL, channel " << i << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithCuts[i]->Draw("same");
      c0->SaveAs(Form("%s/sig-mc-chan%d.pdf", fDirectory.c_str(), i));
      cout << "----> "; gFile->pwd(); 
    } else if (5 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA SIGNAL, channel " << i  << endl;
      bgBlind(fhMassWithCuts[i], 1, 4.7, 6.0);
      cout << "fBgExp = " << fBgExp << "+/-" << fBgExpE << endl;
      cout << "fBgHist = " << fBgHist << "+/-" << fBgHistE << endl;
      aa->bgObs = fBgHist;
      fhMassWithCuts[i]->Draw();
      c0->SaveAs(Form("%s/sig-data-chan%d.pdf", fDirectory.c_str(), i));
    } else if (10 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC NORMALIZATION, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithCuts[i]->Draw("same");
      c0->SaveAs(Form("%s/norm-mc-chan%d.pdf", fDirectory.c_str(), i));
    } else if (11 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA NORMALIZATION, channel " << i  << endl;
      normYield(fhMassWithCuts[i], mode, 5.0, 5.5);
      aa->fitYield  = fNormSig; 
      aa->fitYieldE = fNormSigE; 
      fhMassWithCuts[i]->Draw();
      c0->SaveAs(Form("%s/norm-data-chan%d.pdf", fDirectory.c_str(), i));
    } else if (20 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC CONTROL SAMPLE, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithCuts[i]->Draw("same");
      c0->SaveAs(Form("%s/cs-mc-chan%d.pdf", fDirectory.c_str(), i));
    } else if (21 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA CONTROL SAMPLE, channel " << i  << endl;
      csYield(fhMassWithCuts[i], mode, 5.0, 5.6);
      aa->fitYield  = fCsSig; 
      aa->fitYieldE = fCsSigE; 
      fhMassWithCuts[i]->Draw();
      c0->SaveAs(Form("%s/cs-data-chan%d.pdf", fDirectory.c_str(), i));
    } 

    printNumbers(*aa, cout); 
    printNumbers(*aa, fOUT); 
    
    // -- Cache the pwd...
    TDirectory *pD = gFile; 
    fHistFile->cd();
    cout << "==> " << i << "  " << mode << " hMassWithMassCuts[i] = " << fhMassWithMassCuts[i] << endl;
    fhMassWithMassCutsManyBins[i]->SetName(Form("hMassWithMassCutsManyBins%d_chan%d", mode, i)); fhMassWithMassCutsManyBins[i]->Write();
    fhMassWithCutsManyBins[i]->SetName(Form("hMassWithCutsManyBins%d_chan%d", mode, i)); fhMassWithCutsManyBins[i]->Write();
    fhMassWithMassCuts[i]->SetName(Form("hMassWithMassCuts%d_chan%d", mode, i)); fhMassWithMassCuts[i]->Write();
    fhMassWithCuts[i]->SetName(Form("hMassWithCuts%d_chan%d", mode, i)); fhMassWithCuts[i]->Write();
    fhMassNoCuts[i]->SetName(Form("hMassNoCuts%d_chan%d", mode, i)); fhMassNoCuts[i]->Write();
    fhMassAbsNoCuts[i]->SetName(Form("hMassAbsNoCuts%d_chan%d", mode, i)); fhMassAbsNoCuts[i]->Write();
    fhMuTr[i]->SetName(Form("hMuTr%d_chan%d", mode, i)); fhMuTr[i]->Write();
    fhMuId[i]->SetName(Form("hMuId%d_chan%d", mode, i)); fhMuId[i]->Write();
    // -- and get back to it
    pD->cd();
  }

  return fhMassWithCuts[0];
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

  fhMassWithCuts[ichan]->Draw();

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
  ofstream OUT(Form("optimizeUL-%d.txt", seed)); 

  int NCUTS(7);
//string cuts[]   = {"m2pt", "m1pt","pt", "alpha", "chi2", "fls3d", "mwindow"}; 
  double loCuts[] = {2.5,    3.0,   4,    0.02,    1.,     8,        0.020 };
  double hiCuts[] = {3.5,    7.0,   10,   0.06,    3.,     20,       0.100};

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
	fCuts[ic]->iso1 = 0.75;
	if (0 == i) fCuts[ic]->m2pt = cut; 
	if (1 == i) fCuts[ic]->m1pt = cut; 
	if (2 == i) fCuts[ic]->pt = cut; 
	if (3 == i) fCuts[ic]->alpha = cut; 
	if (4 == i) fCuts[ic]->chi2dof = cut; 
	if (5 == i) fCuts[ic]->fls3d = cut; 
	if (6 == i) fCuts[ic]->mBsLo = 5.370 - cut; 
	if (6 == i) fCuts[ic]->mBsHi = 5.370 + cut; 
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
	fhMassWithCuts[ichan]->Draw();
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
	cout << "==> effTot:    " << effTot << " chan=" << ichan << " version=" << version << endl;
	cout << "==> nObs:      " << fNumbersBs[ichan]->bgObs << " chan=" << ichan << " version=" << version  << endl;
	cout << "==> nExp:      " << fBgHistExp << " chan=" << ichan << " version=" << version << endl;
	cout << "==> ul(Bayes): " << ulbayes << " chan=" << ichan << " version=" << version << endl;
	
	OUT << "==> effTot = " << effTot << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> nObs   = " << fNumbersBs[ichan]->bgObs << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> nExp   = " << fBgHistExp << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> UL     = " << ulbayes << " chan=" << ichan << " version=" << version << endl;
      }
  }

}



// ----------------------------------------------------------------------
void anaBmm::testSimpleUL(const char *cuts) {
  //  ofstream OUT("testUL.txt", ios::app);

  fpMc[fSgMc]->cd();
  histAcceptanceAndPreselection(*fNumbersBs[0]); 

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
  fCsSigE = cE/c*fNormSig;

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  
  delete lF1; 

}


// ----------------------------------------------------------------------
string anaBmm::formatTex(double n, std::string name, int digits) {

  char line[200]; 
  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str());
  } else if (0 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%.0f } } }", name.c_str(), n);
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%.1f } } }", name.c_str(), n);
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%.2f } } }", name.c_str(), n);
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%.3f } } }", name.c_str(), n);
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
      fTEX << formatTex(fNumbersNorm[i]->effTot, Form("%s:EFF_TOT_BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;

      OUT << "ACC_BPLUS\t" << i << "\t" << fNumbersNorm[i]->acc << endl;
      fTEX << formatTex(fNumbersNorm[i]->acc, Form("%s:ACC_BPLUS:val", fSuffix.c_str()), 4) << endl;

      //    OUT << "EFF_MU_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effMuidPid << endl;
      fTEX << formatTex(fNumbersNorm[i]->effMuidPid, Form("%s:EFF_MU_PID_BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;

      OUT << "EFF_MU_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effMuidMC << endl;
      fTEX << formatTex(fNumbersNorm[i]->effMuidMC, Form("%s:EFF_MU_MC_BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;

      //    OUT << "EFF_TRIG_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effTrigPid << endl;
      fTEX << formatTex(fNumbersNorm[i]->effTrigPid, Form("%s:EFF_TRIG_PID_BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;

      OUT << "EFF_TRIG_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effTrigMC << endl;
      fTEX << formatTex(fNumbersNorm[i]->effTrigMC, Form("%s:EFF_TRIG_MC_BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;

      OUT << "EFF_CAND_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effCand << endl;
      fTEX << formatTex(fNumbersNorm[i]->effCand, Form("%s:EFF_CAND_BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;

      //    OUT << "EFF_ANA_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effAna << endl;
      OUT << "EFF_ANA_BPLUS\t" << i << "\t" << fNumbersNorm[i]->effAna*fNumbersNorm[i]->effChan << endl;
      fTEX << formatTex(fNumbersNorm[i]->effAna*fNumbersNorm[i]->effChan, Form("%s:EFF_ANA_BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;

      OUT << "OBS_BPLUS\t" << i << "\t" << fNumbersNorm[i]->fitYield << endl;
      fTEX << formatTex(fNumbersNorm[i]->fitYield, Form("%s:OBS_BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
    } 
  } else {
    OUT << "TOT_BPLUS\t" << "0\t" << (fDataLumi[fSgData]/39.4)*440000 << endl;
    OUT << "TOT_BPLUS\t" << "1\t" << (fDataLumi[fSgData]/39.4)*383000 << endl;
  }

  OUT << "######################################################################" << endl;
  fTEX << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  for (unsigned int i = 0; i < fNchan; ++i) {
    OUT << "# -- SIGNAL " << i << endl;
    fTEX << "% -- SIGNAL " << i << endl;

    OUT << "OBS_BKG\t" << i << "\t" << fNumbersBs[i]->bgObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bgObs, Form("%s:OBS_BKG%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "LOW_BD\t" << i << "\t" << fNumbersBs[i]->mBdLo << endl;
    fTEX << formatTex(fNumbersBs[i]->mBdLo, Form("%s:LOW_BD%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "HIGH_BD\t" << i << "\t" << fNumbersBs[i]->mBdHi << endl;
    fTEX << formatTex(fNumbersBs[i]->mBdHi, Form("%s:HIGH_BD%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "LOW_BS\t" << i << "\t" << fNumbersBs[i]->mBsLo << endl;
    fTEX << formatTex(fNumbersBs[i]->mBsLo, Form("%s:LOW_BS%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "HIGH_BS\t" << i << "\t" << fNumbersBs[i]->mBsHi << endl;
    fTEX << formatTex(fNumbersBs[i]->mBsHi, Form("%s:HIGH_BS%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "PSS\t" << i << "\t" << fNumbersBs[i]->pss << endl;
    fTEX << formatTex(fNumbersBs[i]->pss, Form("%s:PSS%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "PSD\t" << i << "\t" << fNumbersBd[i]->psd << endl;
    fTEX << formatTex(fNumbersBs[i]->psd, Form("%s:PSD%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "PDS\t" << i << "\t" << fNumbersBs[i]->pds << endl;
    fTEX << formatTex(fNumbersBs[i]->pds, Form("%s:PDS%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "PDD\t" << i << "\t" << fNumbersBd[i]->pdd << endl;
    fTEX << formatTex(fNumbersBs[i]->pdd, Form("%s:PDD%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "#EFF_TOT_BSMM\t" << i << "\t" << fNumbersBs[i]->effTot << endl;
    fTEX << formatTex(fNumbersBs[i]->effTot, Form("%s:EFF_TOT_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "ACC_BSMM\t" << i << "\t" << fNumbersBs[i]->acc << endl;
    fTEX << formatTex(fNumbersBs[i]->acc, Form("%s:ACC_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    //    OUT << "EFF_MU_BSMM\t" << i << "\t" << fNumbersBs[i]->effMuidPid << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidPid, Form("%s:EFF_MU_PID_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_MU_BSMM\t" << i << "\t" << fNumbersBs[i]->effMuidMC << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidMC, Form("%s:EFF_MU_MC_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    //    OUT << "EFF_TRIG_BSMM\t" << i << "\t" << fNumbersBs[i]->effTrigPid << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigPid, Form("%s:EFF_TRIG_PID_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_TRIG_BSMM\t" << i << "\t" << fNumbersBs[i]->effTrigMC << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigMC, Form("%s:EFF_TRIG_MC_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_CAND_BSMM\t" << i << "\t" << fNumbersBs[i]->effCand << endl;
    fTEX << formatTex(fNumbersBs[i]->effCand, Form("%s:EFF_CAND_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_ANA_BSMM\t" << i << "\t" << fNumbersBs[i]->effAna*fNumbersBs[i]->effChan << endl;
    fTEX << formatTex(fNumbersBs[i]->effAna*fNumbersBs[i]->effChan, Form("%s:EFF_ANA_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "#EFF_TOT_BDMM\t" << i << "\t" << fNumbersBd[i]->effTot << endl;
    fTEX << formatTex(fNumbersBd[i]->effCand, Form("%s:EFF_TOT_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "ACC_BDMM\t" << i << "\t" << fNumbersBd[i]->acc << endl;
    fTEX << formatTex(fNumbersBd[i]->acc, Form("%s:ACC_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;

    //    OUT << "EFF_MU_BDMM\t" << i << "\t" << fNumbersBd[i]->effMuidPid << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidPid, Form("%s:EFF_MU_PID_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_MU_BDMM\t" << i << "\t" << fNumbersBd[i]->effMuidMC << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidMC, Form("%s:EFF_MU_MC_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;

    //    OUT << "EFF_TRIG_BDMM\t" << i << "\t" << fNumbersBd[i]->effTrigPid << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigPid, Form("%s:EFF_TRIG_PID_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_TRIG_BDMM\t" << i << "\t" << fNumbersBd[i]->effTrigMC << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigMC, Form("%s:EFF_TRIG_MC_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_CAND_BDMM\t" << i << "\t" << fNumbersBd[i]->effCand << endl;
    fTEX << formatTex(fNumbersBd[i]->effCand, Form("%s:EFF_CAND_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "EFF_ANA_BDMM\t" << i << "\t" << fNumbersBd[i]->effAna*fNumbersBd[i]->effChan << endl;
    fTEX << formatTex(fNumbersBd[i]->effAna*fNumbersBd[i]->effChan, Form("%s:EFF_ANA_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "# Observed in signal boxes" << endl;
    OUT << "OBS_BSMM\t" << i << "\t" << fNumbersBs[i]->bsObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bsObs, Form("%s:OBS_BSMM%d:val", fSuffix.c_str(), i), 4) << endl;

    OUT << "OBS_BDMM\t" << i << "\t" << fNumbersBs[i]->bdObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bdObs, Form("%s:OBS_BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    
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
  OUT << "anaWmcYield  = " << a.anaWmcYield << endl;
  OUT << "anaYield     = " << a.anaYield << endl;
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
  OUT << "effChan      = " << a.effChan << "+/-" << a.effChanE << endl; 
  OUT << "effMuidMC    = " << a.effMuidMC << "+/-" << a.effMuidMCE << endl;
  OUT << "effMuidPid   = " << a.effMuidPid << "+/-" << a.effMuidPidE << endl;
  OUT << "effTrigMC    = " << a.effTrigMC << "+/-" << a.effTrigMCE << endl;
  OUT << "effTrigPid   = " << a.effTrigPid << "+/-" << a.effTrigPidE << endl;
  OUT << "effCand      = " << a.effCand << "+/-" << a.effCandE << endl;
  OUT << "effAna       = " << a.effAna << "+/-" << a.effAnaE << endl; 
  OUT << "prod(eff)    = " << a.acc*a.effChan*a.effMuidMC*a.effTrigMC*a.effCand*a.effAna << endl;
  OUT << "prod(effPid) = " << a.acc*a.effChan*a.effMuidPid*a.effTrigPid*a.effCand*a.effAna << endl;
  OUT << "effTot       = " << a.effTot << "+/-" << a.effTotE << endl; 
  OUT << "effTotChan   = " << a.effTotChan << "+/-" << a.effTotChanE << endl; 
  OUT << "combGenYield = " << a.combGenYield << endl; 
  OUT << "prodGenYield         = " << a.prodGenYield << endl; 
  OUT << "chanGenYield(prod)   = " << a.chanGenYield << endl; 
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
    fTEX <<  Form("\\vdef{%s:iso1:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->iso1) << endl;
    OUT << "chi2dof " << a->chi2dof << endl;
    fTEX <<  Form("\\vdef{%s:chi2dof:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->chi2dof) << endl;
    OUT << "alpha   " << a->alpha << endl;
    fTEX <<  Form("\\vdef{%s:alpha:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->alpha) << endl;
    OUT << "fls3d   " << a->fls3d << endl;
    fTEX <<  Form("\\vdef{%s:fls3d:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->fls3d) << endl;
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
  a->bgObs = a->bgExp = a->bgExpE = 0; 
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
