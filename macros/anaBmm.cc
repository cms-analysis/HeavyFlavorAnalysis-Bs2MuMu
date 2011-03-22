#include "anaBmm.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "AnalysisDistribution.hh"
#include "initFunc.hh"
#include "mclimit_csm.hh"
#include "bayesianlimit.hh"

#include "TF1.h"
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

  // -- Christoph's cuts
  M2PT  = 3.;
  ISO1  = 0.55;
  CHI2  = 1.7;
  FLS3D = 10.;
  ALPHA = 0.04; 

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

  fBgLo = 4.7;
  fBgHi = 6.0;

  fNData = fNMc = 0; 
  fSgData = fSgMc = fNoData = fNoMc = -1; 
  fCsData = fCsMc = -1; 

  // -- initialize "numbers"
  numbers a;
  a.name = "";
  a.effGenFilter = a.effGenFilterE = 0.;
  a.fitYield = a.fitYieldE = 0.;
  a.genFileYield = a.genYield = a.recoYield = a.muidYield = a.trigYield = a.candYield = a.ana0Yield = a.anaYield = a.anaNmcYield = 0; 
  a.acc = a.accE = 0; 
  a.effMuidMC =  a.effMuidMCE = a.effTrigMC = a.effTrigMCE = 0; 
  a.effMuidPid = a.effMuidPidE = a.effTrigPid = a.effTrigPidE = 0; 
  a.effCand = a.effCandE = 0; 
  a.effAna = a.effAnaE = 0; 
  
  fNumbersSig = fNumbersNorm = fNumbersCS = a;

  fNumbersSig.name  = "signal"; 
  fNumbersSig.effGenFilter  = 0.63;
  fNumbersSig.effGenFilterE = 0.03;

  fNumbersNorm.name = "normalization"; 
  fNumbersNorm.effGenFilter  = 0.24;
  fNumbersNorm.effGenFilterE = 0.013;

  fNumbersCS.name   = "control sample"; 
  fNumbersCS.effGenFilter  = 0.22;
  fNumbersCS.effGenFilterE = 0.015;

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
  if (fMode > 0)  fSuffix += Form("%d", fMode); 
  cout << "--> Dumping output into " << fDirectory << endl;
  fNumbersFileName = fDirectory + "/anaBmm." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  loadFiles(files);
  string hfname  = fDirectory + "/anaBmm." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");
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
    
    if (string::npos != sdset.find("data")) {
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) fSgData = fNData;
      if (string::npos != stype.find("default") && string::npos != stype.find("no")) fNoData = fNData;
      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) fCsData = fNData;
      fpData[fNData] = loadFile(sfile, stype); 
      fDataLumi[fNData] = atof(slumi.c_str()); 
      cout << "open data " << sfile << " as " << stype << " with lumi = " << fDataLumi[fNData] << endl;
      ++fNData;
    } else {
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) fSgMc = fNMc;
      if (string::npos != stype.find("default") && string::npos != stype.find("no")) fNoMc = fNMc;
      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) fCsMc = fNMc;
      fpMc[fNMc] = loadFile(sfile, stype); 
      fMcLumi[fNMc] = atof(slumi.c_str()); 
      cout << "open MC " << sfile << " as " << stype << " with lumi = " << fMcLumi[fNMc] << endl;
      ++fNMc;
    }
  }
}


// ----------------------------------------------------------------------
TFile* anaBmm::loadFile(string file, string type) {
  TFile *f = TFile::Open(file.c_str());
  return f; 
}


// ----------------------------------------------------------------------
void anaBmm::makeAll(int channel) {
  fpMc[fSgMc]->cd();
  dumpCutNames();

  allEffTables();


  //  computeNormUL();
}


// ----------------------------------------------------------------------
void anaBmm::dumpCutNames() {

  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  
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
    OUT <<  Form("\\vdef{%s:%s:cutLine}   {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
    OUT <<  Form("\\vdef{%s:%s:cutValue}  {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
  }
}


// ----------------------------------------------------------------------
void anaBmm::allEffTables() {
  // -- BMM
  fpMc[fSgMc]->cd(); 
  effTable("SgMc");
  fpData[fSgData]->cd(); 
  effTable("SgData");

  // -- B+
  fpData[fNoData]->cd(); 
  effTable("NoData");
  fpMc[fNoMc]->cd(); 
  effTable("NoMc");

  // -- Bs2JpsiPhi
  fpData[fCsData]->cd(); 
  effTable("CsData");
  fpMc[fCsMc]->cd(); 
  effTable("CsMc");
}

// ----------------------------------------------------------------------
void anaBmm::effTable(string smode) {
  TH1D *h = (TH1D*)gFile->Get("analysisDistributions"); 
  if (0 == h) {
    cout << "no histogram analysisDistributions found, returning" << endl;
    return;
  }

  double massLo(5.1), massHi(5.6);
  if (string::npos != smode.find("No")) {
    massLo = 5.0; 
    massHi = 5.5;
  }

  if (string::npos != smode.find("Cs")) {
    massLo = 5.15; 
    massHi = 5.6;
  }
  
  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << " % --- " << smode << endl;
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
  //  double nprev(0.), nprevE(0.), relEff(0.), relEffE(0.), cumEff(0.), cumEffE(0.); 
  double n(0.), nE(0.); 
  double norm(0.),  normE(0.), eff(0.), effE(0.);

  // -- normalization
  AnalysisDistribution *an = new AnalysisDistribution("docatrk");
  an->fMassLo = massLo; 
  an->fMassHi = massHi; 
  an->hMassCu->SetMinimum(0.);
  norm = an->fitMass(an->hMassCu, normE, mode); 
  OUT << Form("%s", (formatTex(norm, fSuffix+":"+cut+"Norm", smode)).c_str()) << endl;
  OUT << Form("%s", (formatTex(normE, fSuffix+":"+cut+"NormE", smode)).c_str()) << endl;
  delete an;
  
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i)); 
    
    if (cut == string("")) {
      cout << "empty string found, break ..." << endl;
      OUT.close();
      break;
    }
    pdfname = Form("%s/%s_%s.pdf", fDirectory.c_str(), smode.c_str(), cut.c_str());
    cout << "AD for " << cut << " results in " << pdfname << endl;
    AnalysisDistribution *a = new AnalysisDistribution(cut.c_str());
    a->fMassLo = massLo; 
    a->fMassHi = massHi; 
    //    a.hMassCu->SetMinimum(0.);
    //    n = a.fitMass(a.hMassCu, nE, mode); 
    a->hMassAo->SetMinimum(0.);
    n = a->fitMass(a->hMassAo, nE, mode); 
    delete a; 

    if ((string::npos != cut.find("tracks")) || (string::npos != cut.find("muons"))) {
      n *=0.5; 
      nE *=0.5; 
    }

    //     if ("hlt" == cut) {
    //       cout << "initialize normalization numbers" << endl;
    //       norm  = n; 
    //       normE = nE; 
    //       nprev = norm; 
    //       nprevE= normE; 
    //     }

    if (nE < normE) nE = normE;
    eff    = norm/n;
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
    //    pdfname = Form("%s/%s_%s_%s_hMassCu.pdf", fDirectory.c_str(), fSuffix.c_str(), smode.c_str(), cut.c_str());
    pdfname = Form("%s/%s_%s_%s_hMassAo.pdf", fDirectory.c_str(), fSuffix.c_str(), smode.c_str(), cut.c_str());
    c0->SaveAs(pdfname.c_str(), "Portrait");

    OUT << Form("%s", (formatTex(n, fSuffix+":"+cut+"N", smode)).c_str()) << endl;
    OUT << Form("%s", (formatTex(nE, fSuffix+":"+cut+"NE", smode)).c_str()) << endl;
    OUT << Form("%s", (formatTex(eff, fSuffix+":"+cut+"eff", smode)).c_str()) << endl;
    OUT << Form("%s", (formatTex(effE, fSuffix+":"+cut+"effE", smode)).c_str()) << endl;
    //     OUT << Form("%s", (formatTex(relEff, fSuffix+":"+cut+"eRel", smode)).c_str()) << endl;
    //     OUT << Form("%s", (formatTex(relEffE, fSuffix+":"+cut+"eRelE", smode)).c_str()) << endl;
    //     OUT << Form("%s", (formatTex(cumEff, fSuffix+":"+cut+"eCum", smode)).c_str()) << endl;
    //     OUT << Form("%s", (formatTex(cumEffE, fSuffix+":"+cut+"eCumE", smode)).c_str()) << endl;

    //     nprev = n; 
    //     nprevE= nE; 
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
    testUL(cutline.c_str()); 
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
void anaBmm::optimizeUL(int nruns) {

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
    testUL(cutline.c_str()); 
  }
}

// ----------------------------------------------------------------------
void anaBmm::acceptanceAndPreselection(numbers &a) {

  // -- Acceptance and MC efficiency numbers based on prefilled histogram
  //    (could be changed to use effTree!)
  TH1D *hacc  = (TH1D*)(gFile->Get("efficiency"));
  if (!hacc) return;
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
void anaBmm::computeNormUL() {
  cout << "--> loopTree: signal MC" << endl;
  loopTree(0);  // signal eff
  c0->Modified(); c0->Update();
  cout << "--> loopTree: signal data" << endl;
  loopTree(1);  // data signal
  c0->Modified(); c0->Update();
  cout << "--> loopTree: norm MC" << endl;
  loopTree(10); // normalization eff
  c0->Modified(); c0->Update();
  cout << "--> loopTree: norm data" << endl;
  loopTree(11); // data normalization 
  c0->Modified(); c0->Update();

  fNobs = static_cast<int>(fBgExp + 0.5);
  
  double ulbarlow = barlow(fNobs, fBgExp, fBgExpE, 0.2);
  double alpha = 1.; 

  if (fBgExp < 0.001) fBgExp = 0.1;
  double ulBayes90 =  blimit(0.9, fNobs, 1.0, 0.2, fBgExp, fBgExpE, alpha);
  double ulBayes95 =  blimit(0.95, fNobs, 1.0, 0.2, fBgExp, fBgExpE, alpha);
  
  cout << "Bayes: " << ulBayes90 << "(95%CL: " << ulBayes95  << ") vs barlow: " << ulbarlow << endl;

  fNul = ulBayes90;

  fUL = (fNul/fNormSig)
    *(fu/fs)
    *(fNumbersNorm.acc/fNumbersSig.acc)
    *(fNumbersNorm.effCand/fNumbersSig.effCand)     
    *(fNumbersNorm.effMuidMC/fNumbersSig.effMuidMC)
    *(fNumbersNorm.effTrigMC/fNumbersSig.effTrigMC)
    *(fNumbersNorm.effAna/fNumbersSig.effAna)
    * fBF;

  cout << "prod(eff) expected UL: " << fUL << endl;

  fUL = (fNul/fNormSig)
    *(fu/fs)
    *(fNumbersNorm.effTot/fNumbersSig.effTot)
    * fBF;

  cout << "effTot expected UL:    " << fUL << endl;

  TH1D *hn = (TH1D*)fHistFile->Get(Form("numbers: %s", fNumbersSig.name.c_str()));
  int bin = 50; 
  hn->GetXaxis()->SetBinLabel(bin, "fNul"); hn->SetBinContent(bin, fNul); 

  bin = 51; 
  hn->GetXaxis()->SetBinLabel(bin, "fNormSig"); hn->SetBinContent(bin, fNormSig);   hn->SetBinError(bin, fNormSigE);  

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
    *(fNumbersNorm.acc/fNumbersCS.acc)
    *(fNumbersNorm.effCand/fNumbersCS.effCand)     
    *(fNumbersNorm.effMuidMC/fNumbersCS.effMuidMC)
    *(fNumbersNorm.effTrigMC/fNumbersCS.effTrigMC)
    *(fNumbersNorm.effAna/fNumbersCS.effAna)
    * fBF;

  cout << "fact branching fraction: " << result << endl;

  result = (fCsSig/fNormSig)
    *(fu/fs)
    *(fNumbersNorm.effTot/fNumbersCS.effTot)
    * fBF;

  cout << "branching fraction: " << result << endl;

}

// ----------------------------------------------------------------------
void anaBmm::effTree(int mode) {

}

// ----------------------------------------------------------------------
TH1* anaBmm::loopTree(int mode) {
  // -- mode definition
  // 0  Bs2MuMu MC
  // 1  Bs2MuMu data
  // 10 Bp2JpsiKp MC
  // 11 Bp2JpsiKp data
  // 20 Bs2JpsiPhi MC
  // 21 Bs2JpsiPhi data

  bool bp2jpsikp(false), bs2jpsiphi(false), isMC(false); 

  numbers *aa(0);
  if (0 == mode) {
    isMC = true; 
    fpMc[fSgMc]->cd(); 
    aa = &fNumbersSig;
  }
  if (1 == mode) {
    fpData[fSgData]->cd(); 
  }

  if (10 == mode) {
    isMC = true; 
    bp2jpsikp = true; 
    fpMc[fNoMc]->cd(); 
    aa = &fNumbersNorm;
  }

  if (11 == mode) {
    fpData[fNoData]->cd(); 
  }

  if (20 == mode) {
    isMC = true; 
    bs2jpsiphi = true; 
    fpMc[fCsMc]->cd(); 
    aa = &fNumbersCS;
  }

  if (21 == mode) {
    bs2jpsiphi = true; 
    fpData[fCsData]->cd(); 
  }

  // -- set up histograms
  int NBINS = (fMassHi - fMassLo)/0.025+1;
  TH1D *hMassWithMassCuts = (TH1D*)gFile->Get("hMassWithMassCuts");  
  if (0 == hMassWithMassCuts) {
    hMassWithMassCuts = new TH1D("hMassWithMassCuts", "", NBINS, fMassLo, fMassHi);
    hMassWithMassCuts->SetLineColor(kBlue); 
  } else {
    hMassWithMassCuts->Reset(); 
  }
  TH1D *hMassWithCuts = (TH1D*)gFile->Get("hMassWithCuts");  
  if (0 == hMassWithCuts) {
    hMassWithCuts = new TH1D("hMassWithCuts", "", NBINS, fMassLo, fMassHi);
    hMassWithCuts->SetLineColor(kBlue); 
  } else {
    hMassWithCuts->Reset(); 
  }
  TH1D *hMassNoCuts = (TH1D*)gFile->Get("hMassNoCuts");  
  if (0 == hMassNoCuts) {
    hMassNoCuts = new TH1D("hMassNoCuts", "", NBINS, fMassLo, fMassHi);
  } else {
    hMassNoCuts->Reset();
  }
  TH1D *hMassAbsNoCuts = (TH1D*)gFile->Get("hMassAbsNoCuts");  
  if (0 == hMassAbsNoCuts) {
    hMassAbsNoCuts = new TH1D("hMassAbsNoCuts", "", 100, 0, 10);
  } else {
    hMassAbsNoCuts->Reset();
  }
  TH1D *hMuId = (TH1D*)gFile->Get("hMuId");  
  if (0 == hMuId) {
    hMuId = new TH1D("hMuId", "", 100, 0., 1.);
  } else {
    hMuId->Reset();
  }
  TH1D *hMuTr = (TH1D*)gFile->Get("hMuTr");  
  if (0 == hMuTr) {
    hMuTr = new TH1D("hMuTr", "", 100, 0., 1.);
  } else {
    hMuTr->Reset();
  }

  // -- set up tree
  TTree *t;
  t = (TTree*)gFile->Get("events");
  int brun, bevt, btm, bq1, bq2; 
  double bg1pt, bg2pt, bg1eta, bg2eta;
  double bm, bcosa, biso1, bchi2, bdof, bdocatrk, bfls3d, bm1pt, bm1eta, bm2pt, bm2eta;
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
  t->SetBranchAddress("m",&bm);
  t->SetBranchAddress("cosa",&bcosa);
  t->SetBranchAddress("iso1",&biso1);
  t->SetBranchAddress("chi2",&bchi2);
  t->SetBranchAddress("dof",&bdof);
  t->SetBranchAddress("fls3d",&bfls3d);
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
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    // -- require truth matching when on MC
    if (0 == mode && 0 == btm) continue;
    if (10 == mode && 0 == btm) continue;
    if (20 == mode && 0 == btm) continue;

    hMassAbsNoCuts->Fill(bm);
    // -- require wide mass window
    if (bm < fMassLo) continue;
    if (fMassHi < bm) continue;

    // -- gen-level cuts
    if (isMC) {
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

    // -- require basic muon and trackQual cuts
    if (TMath::Abs(bm1eta) > 2.4) continue;
    if (TMath::Abs(bm2eta) > 2.4) continue;
    if (false == bgtqual) continue;
    if (bq1*bq2 > 0) continue;
    if (bm2pt < M2PT) continue; 
    
    // -- must fill this BEFORE the trigger requirement!
    if (bw8mu > 0.) hMuId->Fill(bw8mu, 1./bw8mu); 
    
    // -- now check for muon ID and trigger
    if (false == bhlt) continue;
    if (false == bgmuid) continue;

    hMassNoCuts->Fill(bm);

    // -- weights for trigger
    if (bw8tr > 0.) hMuTr->Fill(bw8tr, 1./bw8tr); 

    // -- apply analysis cand selection 
    if (biso1 < ISO1) continue; 
    if (bchi2/bdof > CHI2) continue;
    if (TMath::IsNaN(bfls3d)) continue;
    if (bfls3d < FLS3D) continue;
    if (TMath::ACos(bcosa) > ALPHA) continue;

    if (bs2jpsiphi && bdr >0.3) continue;
    if (bs2jpsiphi && bmkk < 0.995) continue;
    if (bs2jpsiphi && bmkk > 1.045) continue;

    hMassWithCuts->Fill(bm); 

    if (0 == mode && bm < fSigLo) continue;
    if (0 == mode && bm > fSigHi) continue;
    if (10 == mode && bm < fNormLo) continue;
    if (10 == mode && bm > fNormHi) continue;
    if (20 == mode && bm < fCsLo) continue;
    if (20 == mode && bm > fCsHi) continue;
    hMassWithMassCuts->Fill(bm);

    if (1 == mode && bm > 4.8 && bm < 6.0) {
      cout << "m = " << bm << " run = " << brun << " event = " << bevt 
	   << " mpt = " << bm1pt << "," << bm2pt 
	   << " a = " << TMath::ACos(bcosa) << " iso = " << biso1 << " chi2 = " << bchi2 << " fls3d = " << bfls3d
	   << endl;
    }
    
  }

  // -- Efficiency and acceptance
  if (aa) {
    acceptanceAndPreselection(*aa);
    double a = hMassNoCuts->GetSumOfWeights(); 
    double b = hMassWithMassCuts->GetSumOfWeights();
    aa->anaNmcYield = hMassWithCuts->GetSumOfWeights(); // "no mass cut"
    aa->ana0Yield   = a;
    aa->anaYield    = b; 
    aa->effAna      = b/a;
    aa->effAnaE     = dEff(static_cast<int>(b), static_cast<int>(a));
    aa->effMuidPid  = hMuId->GetMean();
    aa->effMuidPidE = hMuId->GetMeanError();
    aa->effTrigPid  = hMuTr->GetMean();
    aa->effTrigPidE = hMuTr->GetMeanError();
    aa->effTot      = b/(aa->genYield);
    aa->effTotE     = dEff(static_cast<int>(b), static_cast<int>(aa->genYield));
  }

  if (0 == mode) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "==> MC SIGNAL" << endl;
    printNumbers(*aa); 
    hMassNoCuts->Draw();
    hMassWithCuts->Draw("same");
    c0->SaveAs(Form("%s/sig-mc.pdf", fDirectory.c_str()));
  } else if (1 == mode) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "==> Data SIGNAL" << endl;
    bgBlind(hMassWithCuts, 1, 4.7, 6.0);
    cout << "fBgExp = " << fBgExp << "+/-" << fBgExpE << endl;
    hMassWithCuts->Draw();
    c0->SaveAs(Form("%s/sig-data.pdf", fDirectory.c_str()));
  } else if (10 == mode) {
    // -- efficiency filling
    cout << "----------------------------------------------------------------------" << endl;
    cout << "==> MC NORMALIZATION" << endl;
    printNumbers(*aa); 
    hMassNoCuts->Draw();
    hMassWithCuts->Draw("same");
    c0->SaveAs(Form("%s/norm-mc.pdf", fDirectory.c_str()));
  } else if (11 == mode) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "==> Data NORMALIZATION" << endl;
    normYield(hMassWithCuts, mode, 5.0, 5.5);
    hMassWithCuts->Draw();
    c0->SaveAs(Form("%s/norm-data.pdf", fDirectory.c_str()));
  } else if (20 == mode) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "==> MC control sample" << endl;
    printNumbers(*aa); 
    hMassNoCuts->Draw();
    hMassWithCuts->Draw("same");
    c0->SaveAs(Form("%s/cs-mc.pdf", fDirectory.c_str()));
  } else if (21 == mode) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "==> Data control sample" << endl;
    csYield(hMassWithCuts, mode, 5.0, 5.6);
    hMassWithCuts->Draw();
    c0->SaveAs(Form("%s/cs-data.pdf", fDirectory.c_str()));
  } 

  fHistFile->cd();
  hMassWithMassCuts->SetName(Form("hMassWithMassCuts_%d", mode)); hMassWithMassCuts->Write();
  hMassWithCuts->SetName(Form("hMassWithCuts_%d", mode)); hMassWithCuts->Write();
  hMassNoCuts->SetName(Form("hMassNoCuts_%d", mode)); hMassNoCuts->Write();
  hMassAbsNoCuts->SetName(Form("hMassAbsNoCuts_%d", mode)); hMassAbsNoCuts->Write();
  hMuTr->SetName(Form("hMuTr_%d", mode)); hMuTr->Write();
  hMuId->SetName(Form("hMuId_%d", mode)); hMuId->Write();

  return hMassWithCuts;
}


// ----------------------------------------------------------------------
void anaBmm::testUL(const char *cuts) {
  ofstream OUT("testUL.txt", ios::app);

  fpMc[fSgMc]->cd();
  acceptanceAndPreselection(fNumbersSig); 

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

  fNumbersSig.effTot  = afterCuts/(fNumbersSig.genYield);
  fNumbersSig.effTotE = dEff(static_cast<int>(afterCuts), static_cast<int>(fNumbersSig.genYield));

  // -- Background expectation
  TTree *td = (TTree*)(fpData[fSgData]->Get("events"));
  td->Draw("m>>uSG", Form("%s&&%s", defCuts, cuts), "goff");
  bgBlind(h, 1, 4.7, 6.0);
  
  fNobs = static_cast<int>(fBgExp + 0.5);
  
  fBgExp += 0.01; // to protect against blimit assert

  if (fBgExp < 0) {
    fBgExp = fBgHist; 
    fBgExpE = fBgHistE; 
  }
    
  cout << "eff:   " << afterCuts << "/" << fNumbersSig.genYield << " = " << fNumbersSig.effTot << "+/-" << fNumbersSig.effTotE << endl;
  cout << "BG:    " << fBgExp << "+/-" << fBgExpE << " histogram counts: " 
       << fBgHist << "+/-" << fBgHistE
       << endl;
  cout << "Nobs : " << fNobs << endl;
  if (fNobs > 10) return; 

  double nbs = 2.0e9*(1.0-0.12);
  
  double nulbarlow = barlow(fNobs, fBgExp, fBgExpE, 0.2);
  double nulbayes  = blimit(0.9, fNobs, 1.0, 0.2, fBgExp, fBgExpE, 1);

  double ulbarlow = nulbarlow/(fNumbersSig.effTot*nbs);
  double ulbayes  = nulbayes/(fNumbersSig.effTot*nbs);
  cout << Form("Nul: barlow = %3.2f, bayes = %4.2e", nulbarlow, nulbayes) << endl;
  OUT << "==> " << fUL 
      << ". obs = " << fNobs << " bg = " << fBgExp << "+/-" << fBgExpE << "(fBgHist = " << fBgHist << ")" << endl
      << " eff = " << afterCuts << "/" << fNumbersSig.genYield << " = " << fNumbersSig.effTot << endl 
      << endl;

  cout << "UL: barlow = " << ulbarlow << " ulbayes " << ulbayes << " BG: " << fBgExp << " hist: " << h->GetSumOfWeights() << endl;
  h->Draw();

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

  TF1 *lF1(0);

  double histCount = h->Integral(h->FindBin(fBgLo), h->FindBin(fBgHi)-1); 
  cout << "bgBlind: histCount = " << histCount << " starting at " << h->FindBin(fBgLo) << " to " << h->FindBin(fBgHi)-1 << endl;
  fBgHist  = histCount*(fSigHi-fSigLo)/(fBgHi-fBgLo-0.4);
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
void anaBmm::printNumbers(numbers &a) {
  cout << "numbers for \"" << a.name.c_str() << "\"" << endl;
  cout << "fitYield     = " << a.fitYield << "+/-" << a.fitYieldE << endl;
  cout << "genFileYield = " << a.genFileYield << endl;
  cout << "genYield     = " << a.genYield << endl;
  cout << "recoYield    = " << a.recoYield << endl;
  cout << "muidYield    = " << a.muidYield << endl;
  cout << "trigYield    = " << a.trigYield << endl;
  cout << "candYield    = " << a.candYield << endl;
  cout << "ana0Yield    = " << a.ana0Yield << endl;
  cout << "anaNmcYield  = " << a.anaNmcYield << endl;
  cout << "anaYield     = " << a.anaYield << endl;
  cout << "gen filter   = " << a.effGenFilter << endl;
  cout << "acceptance   = " << a.acc << "+/-" << a.accE << endl;
  cout << "effMuidMC    = " << a.effMuidMC << "+/-" << a.effMuidMCE << endl;
  cout << "effTrigMC    = " << a.effTrigMC << "+/-" << a.effTrigMCE << endl;
  cout << "effMuidPid   = " << a.effMuidPid << "+/-" << a.effMuidPidE << endl;
  cout << "effTrigPid   = " << a.effTrigPid << "+/-" << a.effTrigPidE << endl;
  cout << "effCand      = " << a.effCand << "+/-" << a.effCandE << endl;
  cout << "effAna       = " << a.effAna << "+/-" << a.effAnaE << endl; 
  cout << "prod(eff)    = " << a.acc*a.effMuidMC*a.effTrigMC*a.effCand*a.effAna << endl;
  cout << "prod(eff')   = " << a.acc*a.effMuidPid*a.effTrigPid*a.effCand*a.effAna << endl;
  cout << "effTot       = " << a.effTot << "+/-" << a.effTotE << endl; 


  // -- dump into fHistFile
  int bin(0); 
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
  hn->GetXaxis()->SetBinLabel(bin, "anaNmcYield"); hn->SetBinContent(bin, a.anaNmcYield);  

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
}
