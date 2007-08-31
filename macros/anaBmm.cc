#include "anaBmm.hh"

#include "TF1.h"
#include "THStack.h"
#include "TTree.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TMath.h"

#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include <vector.h>

// -- CDF bayesian limit calculator
#include "bayesianlimit.h"
#include "blimit.c"
#include "incompletebeta.c"
#include "incompletegamma.c"
#include "posterior.c"  
#include "postint.c"

// -- CMS sCp calculator (this is very bad code; 
//    should be replaced with http://cmsdoc.cern.ch/~bityukov/durham/ScP.cc)
#include "scp.cc"


ClassImp(anaBmm)


// ----------------------------------------------------------------------
double f_expo(double *x, double *par) {
  return par[0]*TMath::Exp(-x[0]*par[1]);
}

// ----------------------------------------------------------------------
double f_Gauss(double *x, double *par) {
  // par[0] -> area
  // par[1] -> mean
  // par[2] -> sigma

  double sqrt2pi = 2.506628275;

  if (par[2] > 0.) {
    Double_t arg = (x[0] - par[1]) / par[2];
    Double_t fitval =  (par[0]/(sqrt2pi*par[2])) * TMath::Exp(-0.5*arg*arg);
    return fitval;
  }
  else {
    return -1.;
  }
}


// ----------------------------------------------------------------------
double f_gauss(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma

  if (par[2] > 0.) {
    Double_t arg = (x[0] - par[1]) / par[2];
    Double_t fitval =  par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
  }
  else {
    return -1.;
  }
}


// ----------------------------------------------------------------------
double f_2gauss(double *x, double *par) {
  // par[0] -> area
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> area
  // par[4] -> mean
  // par[5] -> sigma


  Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.); 
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[5] > 0.) {
    arg2 = (x[0] - par[4]) / par[5];
    fitval2 =  par[3]*TMath::Exp(-0.5*arg2*arg2);
  }
  Double_t fitval = fitval1 + fitval2;
  return fitval;
}


// ----------------------------------------------------------------------
double f_2G(double *x, double *par) {
  // par[0] -> area
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> fraction in second gaussian
  // par[4] -> mean
  // par[5] -> sigma

  double sqrt2pi = 2.506628275;

  Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.); 
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  (par[0]/(sqrt2pi*par[2]))*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[5] > 0.) {
    arg2 = (x[0] - par[4]) / par[5];
    fitval2 =  (par[3]*par[0]/(sqrt2pi*par[2]))*TMath::Exp(-0.5*arg2*arg2);
  }
  Double_t fitval = fitval1 + fitval2;
  return fitval;
}


// ----------------------------------------------------------------------
double f_2g(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> fraction in second gaussian
  // par[4] -> mean
  // par[5] -> sigma
  Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.); 
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[5] > 0.) {
    arg2 = (x[0] - par[4]) / par[5];
    fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
  }
  Double_t fitval = fitval1 + fitval2;
  return fitval;
}


// ----------------------------------------------------------------------
// pol0 and double Gauss
double f_p0ag(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] = par 0 of pol0

  return  (par[3] + f_Gauss(x, &par[0]));
}

// ----------------------------------------------------------------------
anaBmm::anaBmm(const char *files) { 
  init(files);
}


// ----------------------------------------------------------------------
void anaBmm::init(const char *files) {
  fFont = 132; 

  sprintf(inDir, "bmmroot");
  sprintf(outDir, "anabmm");

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (c0 == 0) {
    cout << "TCanvas c0 not found. Creating my own version." << endl;
    c0 = new TCanvas("c0","--c0--",356,0,656,700);
  }
  tl  = new TLatex();
  pl  = new TLine();
  pa  = new TArrow();
  box = new TBox();

  // -- functions
  f0 = new TF1("f0", f_gauss,  5.0,  6.0, 3);
  f1 = new TF1("f1", f_2gauss, 5.0,  6.0, 6);
  f2 = new TF1("f2", f_2G,     5.0,  6.0, 6);
  f3 = new TF1("f3", f_2g,     5.0,  6.0, 6);
  f4 = new TF1("f4", f_p0ag,   4.9,  5.9, 4);

  f10= new TF1("f10", f_expo,  4.0, 10.0, 2);

  
  // -- setup all files
  nSg = nMc = nDa = 0;
  for (int i = 0; i < 10; ++i) {
    fNevtS[i] = 0.;
    fLumiS[i] = 0.;
    fS[i] = 0; 
    fD[i] = 0; 
  }

  for (int i = 0; i < 20; ++i) {
    fNevtM[i] = 0.;
    fLumiM[i] = 0.;
    fM[i] = 0; 
  }

  cout << "--> Loading rootfiles" << endl;
  fLumiD[0] = 10; 
  if (!strcmp(files, "nada")) {
    loadFiles();
  } else {
    loadFiles(files);
  }
}


// ----------------------------------------------------------------------
void anaBmm::loadFiles() {

  // ----------------------------------------------------------------------
  // vXs = visible cross-section 
  //     = \sigma * (NSEL/NGEN) * BF-corrections * kinematic-cuts
  //       +-----------------------------------+   +------------+
  //         -> from PYTHIA logfiles                -> possibly add.
  //                                                cuts in bmmTree
  //                                                (to unify kinematic cuts)
  // ----------------------------------------------------------------------

  //   TH1D *h; 

  // -- private MC
  //   fS[0]     = new TFile(Form("%s/csg.default-001.root",inDir));  
  //   fS[0]     = new TFile(Form("%st/csg-001.hard2.root", inDir));  
  //   fS[0]     = new TFile(Form("%/csg-001.root", inDir));  
  //   cout << "Loaded " << fS[0]->GetName() << endl;
  //   h         = (TH1D*)fS[0]->Get("ER1");
  //   fvXsS[0]  = 5.522E+01*1.E12*(5000./24.8e6)*3.5e-9;
  //   fNevtS[0] = h->GetBinContent(h->FindBin(1.1));
  //   fLumiS[0] = fNevtS[0]/fvXsS[0] ;

  //   // -- official MC
  //   fS[1]     = new TFile("/data/ursl/060316/sg-001.root");  
  //   h         = (TH1D*)fS[1]->Get("ER1");
  //   fvXsS[1]  = 5.522E+01*1.E12*(5000./24.8e6)*3.5e-9;
  //   fNevtS[1] =  h->GetBinContent(h->FindBin(1.1));
  //   fLumiS[1] = fNevtS[1]/fvXsS[1] ;

  //   // -- default
  //   fM[0]     = new TFile(Form("%s/cbg.default-001.root",inDir));  
  //   //  fM[0]     = new TFile(Form("%s/cbg-001.hard2.root", inDir));  
  //   //  fM[0]     = new TFile(Form("%s/cbg-001.root", inDir));  
  //   cout << "Loaded " << fM[0]->GetName() << endl;
  //   h         = (TH1D*)fM[0]->Get("ER1");
  //   fvXsM[0]  = 1.1E+01*1.E12*1.584e-6;
  //   fNevtM[0] =  h->GetBinContent(h->FindBin(1.1));
  //   fLumiM[0] = fNevtM[0]/fvXsM[0] ;

  //   fNumbersFileName = TString(Form("%s/anaBmm.tex", outDir));
  //   sprintf(line, "rm -f %s", fNumbersFileName.Data());
  //   system(line);
  //   dumpFiles();
  
}


// ----------------------------------------------------------------------
void anaBmm::loadFiles(const char *filename) {
  char buffer[200];
  char type[100];
  char file[100];
  char signature[100];
  float visXsection;
  ifstream is(filename);

  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}

    visXsection = -1.;

    sscanf(buffer, "type=%s file=%s vxs=%e signature=%s", type, file, &visXsection, signature);  

    TString cn(file);
    cn.ReplaceAll("bmmroot", "");
    cn.ReplaceAll("default", "");
    cn.ReplaceAll("treebmm", "");
    cn.ReplaceAll("root", "");
    cn.ReplaceAll(".", "");
    cn.ReplaceAll("/", "");

    if (!strcmp(type, "mc")) {
      sprintf(signature, "bbbar");
      sprintf(line, "fM[%d] = ", nMc);
      cout << "anaBmm::loadFiles> Loading MC file   " << line << file << " and vis x-section: " << visXsection << endl;
      loadMc(file, visXsection, signature);
    }
    else if (!strcmp(type, "hmc")) {
      sprintf(signature, "%s", cn.Data());
      sprintf(line, "fM[%d] = ", nMc);
      cout << "anaBmm::loadFiles> Loading MC file   " << line << file << " and vis x-section: " << visXsection << endl;
      loadMc(file, visXsection, signature);
    }
    else if (!strcmp(type, "sg")) {
      sprintf(line, "fS[%d] = ", nSg);
      cout << "anaBmm::loadFiles> Loading signal file" << line << file << " and vis x-section: " << visXsection << endl;
      loadSg(file, visXsection, signature);
    }
    else if (!strcmp(type, "da")) {
      sprintf(signature, "Data");
      sprintf(line, "fD[%d] = ", nDa);
      cout << "anaBmm::loadFiles> Loading data file  " << line << file << " and vis x-section: " << visXsection  << endl;
      loadDa(file, visXsection, signature);
    }
    else {
      //
    }
  }

  fMu = 1.;
  fPi = 5E-3;
  fKa = 1E-2;
  fProt = 1E-3;

  TString fn(filename);
  fn.ReplaceAll("bmm", "");
  fn.ReplaceAll("files", "");
  fn.ReplaceAll(".", "");

  fNumbersFileName = TString(Form("%s/anaBmm.%s.tex", outDir, fn.Data()));
  sprintf(line, "rm -f %s", fNumbersFileName.Data());
  system(line);
  dumpFiles();
}



// ----------------------------------------------------------------------
void anaBmm::loadSg(const char *name, double lumi, const char *sign) {
  if (nSg > 10) {
    cout << "Too many open Signal files. Increase nSg. " << endl;
    return;
  } 
  fS[nSg] = new TFile(name);
  fLumiS[nSg] = lumi;

  cout << "Loaded " << fS[nSg]->GetName() << " with vis x-section " << lumi << " and signature " << sign << endl;
  TH1 *h      = (TH1D*)fS[nSg]->Get("ER1");
  fvXsS[nSg]  = lumi;
  fNevtS[nSg] = h->GetBinContent(h->FindBin(1.1));
  fLumiS[nSg] = fNevtS[nSg]/fvXsS[nSg] ;
  fSignS[nSg] = TString(sign);
  getSignature(fSignS[nSg], fSignTitleS[nSg]);

  ++nSg; 
}

// ----------------------------------------------------------------------
void anaBmm::loadMc(const char *name, double lumi, const char *sign) {
  if (nMc > 20) {
    cout << "Too many open MC files. Increase nMc. " << endl;
    return;
  } 
  fM[nMc] = new TFile(name);
  fLumiM[nMc] = lumi;
  
  cout << "Loaded " << fM[nMc]->GetName() << " with vis x-section " << lumi << " and signature " << sign << endl;
  TH1 *h      = (TH1D*)fM[nMc]->Get("ER1");
  fvXsM[nMc]  = lumi;
  fNevtM[nMc] = h->GetBinContent(h->FindBin(1.1));
  fLumiM[nMc] = fNevtM[nMc]/fvXsM[nMc] ;
  fSignM[nMc] = TString(sign);
  getSignature(fSignM[nMc], fSignTitleM[nMc]);

  ++nMc; 
}

// ----------------------------------------------------------------------
void anaBmm::loadDa(const char *name, double lumi, const char *sign) {
  if (nDa > 10) {
    cout << "Too many open DATA files. Increase nDa. " << endl;
    return;
  } 
  fD[nDa] = new TFile(name);
  fLumiD[nDa] = lumi;

  cout << "Loaded " << fD[nDa]->GetName() << " with vis x-section " << lumi << " and signature " << sign << endl;
  TH1 *h      = (TH1D*)fD[nDa]->Get("ER1");
  fvXsD[nDa]  = lumi;
  fNevtD[nDa] = h->GetBinContent(h->FindBin(1.1));
  fLumiD[nDa] = fNevtD[nDa]/fvXsD[nDa] ;
  fSignD[nDa] = TString(sign);
  getSignature(fSignD[nDa], fSignTitleD[nDa]);

  ++nDa; 
}


// ----------------------------------------------------------------------
void anaBmm::dumpFiles() {

  cout << "Assuming data luminosity " << Form("%2.1f", fLumiD[0]) << " /fb" << endl;
  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% -- X-sections, NEVT, and lumi" << endl;
  OUT << Form("\\vdef{Lumi:d0}    {\\ensuremath{ {%2.1f } } }", fLumiD[0]) << endl;
  for (int i = 0; i < 10; ++i) {

    if (fS[i]) {

      fNexpS[i] = fLumiD[0]*fvXsS[i];

       cout << Form("Signal MC[%d]     visible cross-section: %3.2e fb, Nevt: %8.0f -> Lumi: %3.2e /fb   Signature: %s Title: %s", 
 		   i, fvXsS[i], fNevtS[i], fLumiS[i], fSignS[i].Data(), fSignTitleS[i].Data()) 
 	   << endl;


      OUT << Form("\\vdef{channel:s%i} {\\ensuremath{ {%s } } }", i, fSignTitleS[i].Data()) << endl;
      OUT << Form("\\vdef{vNevt:s%i}   {\\ensuremath{ {%5.0f } } }", i, fNevtS[i]) << endl;
      OUT << Form("\\vdef{vNexp:s%i}   {\\ensuremath{ {%5.0f } } }", i, fNexpS[i]) << endl;
      OUT << Form("\\vdef{vLumi:s%i}   {\\ensuremath{ {%5.1f } } }", i, fLumiS[i]) << endl;
    }
  }
 
  for (int i = 0; i < 20; ++i) {

    if (fM[i]) {

      fMisIdM[i] = getMisID(fSignM[i]);
      fNexpM[i] = fLumiD[0]*fvXsM[i]*fMisIdM[i];


      cout << Form("Background MC[%d] visible cross-section: %3.2e fb, Nevt: %8.0f -> lumi: %3.2e /fb   Signature: %s Title: %s", 
		   i, fvXsM[i], fNevtM[i], fLumiM[i], fSignM[i].Data(), fSignTitleM[i].Data()) 
	   << endl;

      OUT << Form("\\vdef{channel:m%i} {\\ensuremath{ {%s } } }", i, fSignTitleM[i].Data()) << endl;

    //   if ( fvXsM[i] < 1.e4) {
// 	OUT << Form("\\vdef{vXs:m%i}     {\\ensuremath{ {%4.1f } } }", i, fvXsM[i]) << endl;
//       } else if ( fvXsM[i] < 1.e9) {
// 	OUT << Form("\\vdef{vXs:m%i}     {\\ensuremath{ {%s } } }", i, (texForm(fvXsM[i])).Data()) << endl;
//       } else {
// 	OUT << Form("\\vdef{vXs:m%i}     {\\ensuremath{ {%s } } }", i, (texForm2(fvXsM[i])).Data()) << endl;
//       }

      OUT << Form("\\vdef{vXs:m%i}     {\\ensuremath{ \\tt {%3.2E } } }", i, fvXsM[i]) << endl;

//       if ( fLumiM[i] < 1.e4 ) {
// 	OUT << Form("\\vdef{vLumi:m%i}   {\\ensuremath{ {%3.2E } } }", i, fLumiM[i]) << endl;
//       } else if ( fLumiM[i] < 1.e9 ) {
// 	OUT << Form("\\vdef{vLumi:m%i}   {\\ensuremath{ {%s } } }", i, (texForm(fLumiM[i]).Data())) << endl;
//       } else {
// 	OUT << Form("\\vdef{vLumi:m%i}   {\\ensuremath{ {%s } } }", i, (texForm2(fLumiM[i]).Data())) << endl;
//       }

      OUT << Form("\\vdef{vLumi:m%i}   {\\ensuremath{ \\tt  {%3.2E } } }", i, fLumiM[i]) << endl;
      OUT << Form("\\vdef{vNevt:m%i}   {\\ensuremath{ {%5.0f } } }", i, fNevtM[i]) << endl;

      if ( fNexpM[i] < 1.e4 ) {
	OUT << Form("\\vdef{vNexp:m%i}   {\\ensuremath{ {%5.0f } } }", i, fNexpM[i]+0.5) << endl;
      } else if ( fNexpM[i] < 1.e9 ) {	
	OUT << Form("\\vdef{vNexp:m%i}   {\\ensuremath{ {%s } } }", i, (texForm(fNexpM[i]).Data())) << endl;
      } else {
	OUT << Form("\\vdef{vNexp:m%i}   {\\ensuremath{ {%s } } }", i, (texForm2(fNexpM[i]).Data())) << endl;
      }

      //      OUT << Form("\\vdef{vNexp:m%i}   {\\ensuremath{ \\tt  {%5.0f } } }", i, fNexpM[i]) << endl;
    }

  }
}


// ----------------------------------------------------------------------
void anaBmm::makeAllPlots() {
   

  mcValidation();
  breco(0);
  breco(2);
  breco(3);
  breco(4);

  jreco(3);
  jreco(4);
   
  muonMisId();

//   effTables();
//   calculateUpperLimit();
//   dumpCuts();
   
//   plotRareDecays();

//   showDistributions(2, 0);
   
//   showProcesses(1);
//   showProcesses(0);
  
}

// -----------------------------------------------------
void anaBmm::plotRareDecays() {

  singleHBG("c030");  // plots for every rare background together with signal
  singleHBG("c033");  // plots for every rare background together with signal

  overlayHBG("c030", 6); // overlay of different rare background with sg + bg
  overlayHBG("c430", 4); // overlay of different rare background with sg + bg
  overlayHBG("c530", 5); // overlay of different rare background with sg + bg

  //  stackHBG(hist);   // stack of different rare background with sg + bg
  plotAll("c030");    // plot sg + (bg & hbg) normalized to lumi 

}


// ----------------------------------------------------------------------
void anaBmm::calculateUpperLimit() {

  double nBs  = fLumiD[0] * 500.*1.e9 * 0.107 * 2.;
  double ekin = 1000. / (4.8e6 * (500./55.e3) * 0.107 * 2.);

  // -- BF < N_UL(Nobs) / scaleUL  =  N_UL(Nobs) / (epsilon * N_Bs)
  double atlasUL = expUL(7, 20);
  double nUL = expUL(fNsg, fNbg);
  double scaleUL = ekin * fEsg * nBs; // this is: eff_kin * eff_total * N_Bs 

  double expectedUL = nUL/scaleUL;


  cout << "Nsg: " << fNsg << " +/- " << fNsgE
       << " Nbg: " << fNbg << " +/- " << fNbgE
       << endl;

  double Scp = scp(fNsg, fNbg, fNbgE, 0.);

  for (int i = 2; i < 50.; i += 2) {
    cout << "i = " << i << "  "; 
    scp(i*fNsg, i*fNbg, i*fNbgE, 0.);
  }


  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT  << Form("\\vdef{nBs}   {\\ensuremath{{%4.1e} } }", nBs) << endl;
  OUT  << Form("\\vdef{ekin}  {\\ensuremath{{%4.3f} } }", ekin) << endl;
  OUT  << Form("\\vdef{ExpectedNobs}        {\\ensuremath{{%4.1f} } }", fNsg + fNbg) << endl;
  OUT  << Form("\\vdef{ExpectedNul}         {\\ensuremath{{%4.1f} } }", nUL) << endl;
  OUT  << Form("\\vdef{ExpectedUpperLimit}  {\\ensuremath{{%s} } }", (texForm31(expectedUL)).Data()) << endl;
  OUT  << Form("\\vdef{Scp}                 {\\ensuremath{{%4.1f} } }", Scp) << endl;
  OUT  << Form("\\vdef{AtlasUpperLimit}     {\\ensuremath{{%s} } }", (texForm31(atlasUL/scaleUL)).Data()) << endl;
  OUT.close();
}
  


// ----------------------------------------------------------------------
double anaBmm::expUL(double s0, double b0) {

  double beta  = 0.9;
  double e0    = 1.; 
  double esig  = 0.25; 
  double bsig  = 1.6*b0;  
  double alpha = 1.; 

  int nbins = int(2*(b0+s0)); 
  TH1D *h1 = new TH1D("h1", "Events observed", nbins, 0., nbins);

  double ulmean =  blimit(beta, s0+b0, e0, esig, b0, bsig, alpha);

  int nobs(0); 
  double ul(0.), ulSum(0.);
  int ulSample(1000);
  // -- This loop results in an UL that is very close to the mean above...
  for (int i = 0; i < ulSample; ++i) {
    nobs  = gRandom->Poisson(e0*s0 + b0);
    h1->Fill(nobs);
    ul    = blimit(beta, nobs, e0, esig, b0, bsig, alpha);
    ulSum += ul;
  }
  ul = ulSum/ulSample;


  cout << "----------------------------------------------------------------------" << endl;
  cout << Form("Mean UL  for n = s(%4.3f) + b(%4.3f): %5.4f", s0, b0, ulmean)
       << endl;

  cout << Form("Expected UL for s = %4.3f, b = %4.3f: %5.4f", s0, b0, ul)
       << endl;
  cout << "----------------------------------------------------------------------" << endl;


  h1->Draw();
  
  return ul;
}


// ----------------------------------------------------------------------
void anaBmm::effTables() {

  effTable(fS[0], "s0");
  effTable(fM[0], "m0");
  effTable(fM[0], "2mId");
  effTable(fM[0], "mIdMu+");
  effTable(fM[0], "2mu+");
  effTable(fM[0], "qcd");
  effTable(fM[0], "r0");
  effTable(fM[0], "c0");
  

}


// ----------------------------------------------------------------------
void anaBmm::dumpCuts() {
  TH1D *hS = (TH1D*)fS[0]->Get("hcuts");
  TH1D *hM = (TH1D*)fM[0]->Get("hcuts");

  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% -- Cuts" << endl;

  char label[100];
  sprintf(label, "%s", hS->GetXaxis()->GetBinLabel(1));

  double sVal, mVal;
  for (int i = 1; i < hS->GetNbinsX(); ++i) {
    if (strcmp(hS->GetXaxis()->GetBinLabel(i), "")) {
      cout << hS->GetXaxis()->GetBinLabel(i) << endl;
      sVal = hS->GetBinContent(i);
      mVal = hM->GetBinContent(i);

      if (sVal != mVal) {
	cout << "====> Error: Signal and BG MC run with different cut values!" << endl;
      } else {

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "p_{T}(B_{s}) [GeV]")) {
	  OUT  << Form("\\vdef{cut:%s:ptbs}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "p_{T}^{min}(l) [GeV]")) {
	  OUT  << Form("\\vdef{cut:%s:ptlo}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "p_{T}^{max}(l) [GeV]")) {
	  OUT  << Form("\\vdef{cut:%s:pthi}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "R_{#mu#mu}^{min}")) {
	  OUT  << Form("\\vdef{cut:%s:rmmlo}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "R_{#mu#mu}^{max}")) {
	  OUT  << Form("\\vdef{cut:%s:rmmhi}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "#eta_{T}^{min}(l)")) {
	  OUT  << Form("\\vdef{cut:%s:etalo}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "#eta_{T}^{max}(l)")) {
	  OUT  << Form("\\vdef{cut:%s:etahi}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "TIP(l) [cm]")) {
	  OUT  << Form("\\vdef{cut:%s:tip}      {\\ensuremath{%5.3f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "#chi^2")) {
	  OUT  << Form("\\vdef{cut:%s:chi2}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "l_{3d} [cm]")) {
	  OUT  << Form("\\vdef{cut:%s:l3d}    {\\ensuremath{%5.3f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "cos(#alpha)")) {
	  OUT  << Form("\\vdef{cut:%s:cosalpha}    {\\ensuremath{%6.5f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "l_{xy} [cm]")) {
	  OUT  << Form("\\vdef{cut:%s:lxy}    {\\ensuremath{%5.3f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "l_{xy}/#sigma_{xy}")) {
	  OUT  << Form("\\vdef{cut:%s:lxy/sxy}    {\\ensuremath{%3.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "I_{veto}")) {
	  OUT  << Form("\\vdef{cut:%s:isoveto}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "R_{I}")) {
	  OUT  << Form("\\vdef{cut:%s:isocone}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "I")) {
	  OUT  << Form("\\vdef{cut:%s:isolation}    {\\ensuremath{%5.3f } } ", label, sVal) << endl;
	}


	
      }
    }
  }
}

// ----------------------------------------------------------------------
void anaBmm::effTable(TFile *f, const char *tag) {

  TH1D *h = (TH1D*)f->Get("ER1");
  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% -- Efficiency" << endl;

  double n;

  // -- Fix lumi normalisation scaling factor
  double SF(0.), comb(0.);
  if (!strcmp(tag, "s0")) {
    SF = fLumiD[0]/fLumiS[0];
  }

  else if (!strcmp(tag, "m0")) {
    SF = fLumiD[0]/fLumiM[0];
  }

  else {
    SF = 1; comb = 1;
    
    TH1 *sumER1 = sumHistMC("ER1", 2, tag);  
    h = (TH1D*)sumER1;
  }

  //  SF = 1.;


  double nchain = h->GetBinContent(h->FindBin(0.1));
  OUT  << Form("\\vdef{nchain:%s}    {\\ensuremath{{%5.1f } } }", tag, nchain) << endl;

  // -- kinematic cuts
  n = h->GetBinContent(h->FindBin(1.1));
  OUT  << Form("\\vdef{cutKin:%s}    {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  OUT  << Form("\\vdef{effKin:%s}    {\\ensuremath{{%4.3f } } }", tag, n/nchain) << endl;
  OUT  << Form("\\vdef{effKinE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nchain))) << endl;
  double nkin = n; 

  // -- L1
  n = h->GetBinContent(h->FindBin(101.1));
  OUT  << Form("\\vdef{cutL1:%s}    {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  OUT  << Form("\\vdef{effL1:%s}    {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{effL1E:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;
  double nl1 = n;

  // -- PV and two signal tracks
  n = h->GetBinContent(h->FindBin(203.1));
  OUT  << Form("\\vdef{cutJunk:%s}    {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  OUT  << Form("\\vdef{effJunk:%s}    {\\ensuremath{{%4.3f } } }", tag, n/nl1) << endl;
  OUT  << Form("\\vdef{effJunkE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nl1))) << endl;
  OUT  << Form("\\vdef{effcJunk:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{effcJunkE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nkin))) << endl;

  // -- HLT
  n = h->GetBinContent(h->FindBin(305.1));
  OUT  << Form("\\vdef{cutHLT:%s}    {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  OUT  << Form("\\vdef{effHLT:%s}    {\\ensuremath{{%4.3f } } }", tag, n/nl1) << endl;
  OUT  << Form("\\vdef{effHLTE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nl1))) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effcHLT:%s} {\\ensuremath{{%4.3E } } }", tag, n/nkin) << endl;
    OUT  << Form("\\vdef{effcHLTE:%s}{\\ensuremath{{%4.3E } } }", tag, dEff(int(n),int(nkin))) << endl;
  } else {
    OUT  << Form("\\vdef{effcHLT:%s} {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
    OUT  << Form("\\vdef{effcHLTE:%s}{\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nkin))) << endl;
  }
  double nhlt = n;


  // -- Tracking cuts
  double nt = h->GetBinContent(h->FindBin(420.1));
  n  = h->GetBinContent(h->FindBin(421.1));
  OUT  << Form("\\vdef{cutL0PT:%s}     {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effL0PT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL0PTE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  } else {
    OUT  << Form("\\vdef{effL0PT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL0PTE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  }

  n  = h->GetBinContent(h->FindBin(422.1));
  OUT  << Form("\\vdef{cutL0ETA:%s}     {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effL0ETA:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL0ETAE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  } else {
    OUT  << Form("\\vdef{effL0ETA:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL0ETAE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  }

  n  = h->GetBinContent(h->FindBin(423.1));
  OUT  << Form("\\vdef{cutL0TIP:%s}     {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effL0TIP:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL0TIPE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  } else {
    OUT  << Form("\\vdef{effL0TIP:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL0TIPE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  }

  n  = h->GetBinContent(h->FindBin(431.1));
  OUT  << Form("\\vdef{cutL1PT:%s}     {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effL1PT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL1PTE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  } else {
    OUT  << Form("\\vdef{effL1PT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL1PTE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  }

  n  = h->GetBinContent(h->FindBin(432.1));
  OUT  << Form("\\vdef{cutL1ETA:%s}     {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effL1ETA:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL1ETAE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  } else {
    OUT  << Form("\\vdef{effL1ETA:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL1ETAE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  }

  n  = h->GetBinContent(h->FindBin(433.1));
  OUT  << Form("\\vdef{cutL1TIP:%s}     {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effL1TIP:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL1TIPE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  } else {
    OUT  << Form("\\vdef{effL1TIP:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effL1TIPE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
  }

  n  = h->GetBinContent(h->FindBin(407.1));
  OUT  << Form("\\vdef{cutTRK:%s}     {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effTRK:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nhlt) << endl;
    OUT  << Form("\\vdef{effTRKE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nhlt))) << endl;
    OUT  << Form("\\vdef{effcTRK:%s}  {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
    OUT  << Form("\\vdef{effcTRKE:%s} {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nkin))) << endl;
  } else {
    OUT  << Form("\\vdef{effTRK:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nhlt) << endl;
    OUT  << Form("\\vdef{effTRKE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nhlt))) << endl;
    OUT  << Form("\\vdef{effcTRK:%s}  {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
    OUT  << Form("\\vdef{effcTRKE:%s} {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nkin))) << endl;
  }
  nt = n; // nt now contains nevt with goodtt

  // -- PID cuts
  n  = h->GetBinContent(h->FindBin(509.1));
  OUT  << Form("\\vdef{cutPID:%s}      {\\ensuremath{{%5.1f } } }", tag, SF*n) << endl;
  if (strstr(tag, "m")) {
    OUT  << Form("\\vdef{effPID:%s}    {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effPIDE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
    OUT  << Form("\\vdef{effcPID:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
    OUT  << Form("\\vdef{effcPIDE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nkin))) << endl;
  } else {
    OUT  << Form("\\vdef{effPID:%s}    {\\ensuremath{{%4.3f } } }", tag, n/nt) << endl;
    OUT  << Form("\\vdef{effPIDE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nt))) << endl;
    OUT  << Form("\\vdef{effcPID:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
    OUT  << Form("\\vdef{effcPIDE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n),int(nkin))) << endl;
  }

  // -- HLT selection
  double normc   = h->GetBinContent(h->FindBin(330.1));
  double norm1   = h->GetBinContent(h->FindBin(330.1));
  OUT  << Form("\\vdef{nH0:%s}      {\\ensuremath{{%5.1f } } }", tag, normc) << endl;
  n  = h->GetBinContent(h->FindBin(331.1));
  OUT  << Form("\\vdef{eHpt:%s}     {\\ensuremath{{%4.3f } } }", tag, n/norm1) << endl;
  OUT  << Form("\\vdef{eHptE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm1))) << endl;
  OUT  << Form("\\vdef{cHpt:%s}     {\\ensuremath{{%4.3f } } }", tag, n/normc) << endl;
  OUT  << Form("\\vdef{cHptE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(normc))) << endl;

  norm1 = n;
  n  = h->GetBinContent(h->FindBin(332.1));
  OUT  << Form("\\vdef{eHlrap:%s}     {\\ensuremath{{%4.3f } } }", tag, n/norm1) << endl;
  OUT  << Form("\\vdef{eHlrapE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm1))) << endl;
  OUT  << Form("\\vdef{cHlrap:%s}     {\\ensuremath{{%4.3f } } }", tag, n/normc) << endl;
  OUT  << Form("\\vdef{cHlrapE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(normc))) << endl;

  norm1 = n;
  n  = h->GetBinContent(h->FindBin(333.1));
  OUT  << Form("\\vdef{eHqq:%s}     {\\ensuremath{{%4.3f } } }", tag, n/norm1) << endl;
  OUT  << Form("\\vdef{eHqqE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm1))) << endl;
  OUT  << Form("\\vdef{cHqq:%s}     {\\ensuremath{{%4.3f } } }", tag, n/normc) << endl;
  OUT  << Form("\\vdef{cHqqE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(normc))) << endl;

  norm1 = n;
  n  = h->GetBinContent(h->FindBin(334.1));
  OUT  << Form("\\vdef{eHtip:%s}     {\\ensuremath{{%4.3f } } }", tag, n/norm1) << endl;
  OUT  << Form("\\vdef{eHtipE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm1))) << endl;
  OUT  << Form("\\vdef{cHtip:%s}     {\\ensuremath{{%4.3f } } }", tag, n/normc) << endl;
  OUT  << Form("\\vdef{cHtipE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(normc))) << endl;

  norm1 = n;
  n  = h->GetBinContent(h->FindBin(335.1));
  OUT  << Form("\\vdef{eHchi2:%s}     {\\ensuremath{{%4.3f } } }", tag, n/norm1) << endl;
  OUT  << Form("\\vdef{eHchi2E:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm1))) << endl;
  OUT  << Form("\\vdef{cHchi2:%s}     {\\ensuremath{{%4.3f } } }", tag, n/normc) << endl;
  OUT  << Form("\\vdef{cHchi2E:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(normc))) << endl;

  norm1 = n;
  n  = h->GetBinContent(h->FindBin(336.1));
  OUT  << Form("\\vdef{eHl3d:%s}     {\\ensuremath{{%4.3f } } }", tag, n/norm1) << endl;
  OUT  << Form("\\vdef{eHl3dE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm1))) << endl;
  OUT  << Form("\\vdef{cHl3d:%s}     {\\ensuremath{{%4.3f } } }", tag, n/normc) << endl;
  OUT  << Form("\\vdef{cHl3dE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(normc))) << endl;

  norm1 = n;
  n  = h->GetBinContent(h->FindBin(337.1));
  OUT  << Form("\\vdef{eHm:%s}     {\\ensuremath{{%4.3f } } }", tag, n/norm1) << endl;
  OUT  << Form("\\vdef{eHmE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm1))) << endl;
  OUT  << Form("\\vdef{cHm:%s}     {\\ensuremath{{%4.3f } } }", tag, n/normc) << endl;
  OUT  << Form("\\vdef{cHmE:%s}    {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(normc))) << endl;


  // -- Analysis
  h = (TH1D*)f->Get("AR1");

  if ( strcmp(tag, "s0") && strcmp(tag, "m0") ) {
    
    TH1 *sumAR1 = sumHistMC("AR1", 2, tag);
    h = (TH1D*)sumAR1; 
  }
  
  double norm  = h->GetBinContent(h->FindBin(1.1)) + h->GetBinContent(h->FindBin(2.1));
  OUT  << Form("\\vdef{nA0:%s}      {\\ensuremath{{%5.1f } } }", tag, SF*norm) << endl;
  n  = h->GetBinContent(h->FindBin(1.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nAkin:%s}    {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10){
    OUT  << Form("\\vdef{nAkin:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nAkin:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eAkin:%s}    {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eAkinE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cAkin:%s}    {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cAkinE:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(11.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nAL1:%s}    {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nAL1:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nAL1:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eAL1:%s}    {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eAL1E:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cAL1:%s}    {\\ensuremath{{%4.3f } } }", tag, SF*n/nkin) << endl;
  OUT  << Form("\\vdef{cAL1E:%s}   {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(21.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nAHLT:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nAHLT:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nAHLT:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eAHLT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eAHLTE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cAHLT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cAHLTE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(100.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nALPT:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nALPT:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nALPT:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eALPT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eALPTE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cALPT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cALPTE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(101.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nARmm:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nARmm:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nARmm:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eARmm:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eARmmE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cARmmT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cARmmE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(110.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nABX:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nABX:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nABX:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eABX:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eABXE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cABX:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cABXE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(120.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nABpT:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nABpT:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nABpT:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eABpT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eABpTE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cABpT:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cABpTE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(121.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nABeta:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nABeta:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nABeta:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eABeta:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eABetaE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cABeta:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cABetaE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(122.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nABcos:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nABcos:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nABcos:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eABcos:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eABcosE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cABcos:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cABcosE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(123.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nABlxy:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nABlxy:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nABlxy:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  if (n/norm < 1.e-3 && n/norm) {
    OUT  << Form("\\vdef{eABlxy:%s}   {\\ensuremath{{%s } } }", tag, (texForm(n/norm)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{eABlxy:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  }
  OUT  << Form("\\vdef{eABlxyE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cABlxy:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cABlxyE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;
  double eff   = n/nkin;
  double fnorm = n;

  n  = h->GetBinContent(h->FindBin(124.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nABiso:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nABiso:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nABiso:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eABiso:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eABisoE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cABiso:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cABisoE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  n  = h->GetBinContent(h->FindBin(125.1));
  if (SF*n < 1.e4) {
    OUT  << Form("\\vdef{nABchi2:%s}   {\\ensuremath{{%4.1f } } }", tag, SF*n) << endl;
  } else if (SF*n < 1.e10) {
    OUT  << Form("\\vdef{nABchi2:%s}   {\\ensuremath{{%s } } }", tag, (texForm(SF*n)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{nABchi2:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(SF*n)).Data()) << endl;
  }
  OUT  << Form("\\vdef{eABchi2:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
  OUT  << Form("\\vdef{eABchi2E:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
  OUT  << Form("\\vdef{cABchi2:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
  OUT  << Form("\\vdef{cABchi2E:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;

  if (!strcmp(tag, "s0")) {
    fEsg = n/nkin;
    fEsgE= dEff(int(n), int(nkin));
    OUT  << Form("\\vdef{eAllCuts:%s}   {\\ensuremath{{%4.3f } } }", tag, fEsg) << endl;
    OUT  << Form("\\vdef{eAllCutsE:%s}  {\\ensuremath{{%4.3f } } }", tag, fEsgE) << endl;
  }

  // -- factorizing cuts: Iso & Vertex efficiency
  norm = h->GetBinContent(h->FindBin(224.1)) + h->GetBinContent(h->FindBin(324.1));
  n  = h->GetBinContent(h->FindBin(224.1));
  double eff1 = n/norm;
  if ( norm ) {
    OUT  << Form("\\vdef{eAfBiso:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
    OUT  << Form("\\vdef{eAfBisoE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
    OUT  << Form("\\vdef{cAfBiso:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
    OUT  << Form("\\vdef{cAfBisoE:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;
  }
  else {   // QCD events do not survive tightened preselection
    OUT  << Form("\\vdef{eAfBiso:%s}   {\\ensuremath{{- } } }", tag) << endl;
    OUT  << Form("\\vdef{eAfBisoE:%s}  {\\ensuremath{{- } } }", tag) << endl;
    OUT  << Form("\\vdef{cAfBiso:%s}   {\\ensuremath{{- } } }", tag) << endl;
    OUT  << Form("\\vdef{cAfBisoE:%s}  {\\ensuremath{{- } } }", tag) << endl;
  }


  norm = h->GetBinContent(h->FindBin(226.1)) + h->GetBinContent(h->FindBin(326.1));
  n  = h->GetBinContent(h->FindBin(226.1));
  double eff2 = n/norm;
  if ( norm ) {
    OUT  << Form("\\vdef{eAfBchi2:%s}   {\\ensuremath{{%4.3f } } }", tag, n/norm) << endl;
    OUT  << Form("\\vdef{eAfBchi2E:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(norm))) << endl;
    OUT  << Form("\\vdef{cAfBchi2:%s}   {\\ensuremath{{%4.3f } } }", tag, n/nkin) << endl;
    OUT  << Form("\\vdef{cAfBchi2E:%s}  {\\ensuremath{{%4.3f } } }", tag, dEff(int(n), int(nkin))) << endl;
  }
  else {   // QCD events do not survive tightened preselection
    OUT  << Form("\\vdef{eAfBchi2:%s}   {\\ensuremath{{- } } }", tag) << endl;
    OUT  << Form("\\vdef{eAfBchi2E:%s}  {\\ensuremath{{- } } }", tag) << endl;
    OUT  << Form("\\vdef{cAfBchi2:%s}   {\\ensuremath{{- } } }", tag) << endl;
    OUT  << Form("\\vdef{cAfBchi2E:%s}  {\\ensuremath{{- } } }", tag) << endl;
  }


  // -- factorizing cuts 
  double nIsoFact  = SF*fnorm*eff1;
  double nIsoFactE = SF*TMath::Sqrt(fnorm)*eff1; 

  double nChiFact  = SF*fnorm*eff2; 
  double nChiFactE = SF*TMath::Sqrt(fnorm)*eff2;

  double nAllCutsFact = SF*fnorm*eff1*eff2;          
  double nAllCutsFactE = SF*TMath::Sqrt(fnorm)*eff1*eff2;

  // double eIsoFact = eff*eff1;
  // double eChiFact = eff*eff2;
  // double eAllCutsFact = eff*eff1*eff2;

  //  double eAllCutsFactE = (fEsgE/fEsg)*eff*eff1*eff2;


  if ( comb ) {

    nIsoFact = h->GetBinContent(h->FindBin(424.1));
    nIsoFactE = -9999;  
    nChiFact = h->GetBinContent(h->FindBin(426.1)); 
    nChiFactE = -9999; 
    nAllCutsFact = h->GetBinContent(h->FindBin(428.1));
    nAllCutsFactE = -9999; 
  }

  double eIsoFact = nIsoFact/(SF*nkin);
  double eChiFact = nChiFact/(SF*nkin);
  double eAllCutsFact = nAllCutsFact/(SF*nkin);

  double eIsoFactE = (fEsgE/fEsg)*eIsoFact;
  double eChiFactE = (fEsgE/fEsg)*eChiFact;
  double eAllCutsFactE = (fEsgE/fEsg)*eAllCutsFact;


  // -- fact. Isolation
  OUT  << Form("\\vdef{nAfBiso:%s}   {\\ensuremath{{%4.1f } } }", tag, nIsoFact) << endl;
  OUT  << Form("\\vdef{nAfBisoE:%s}  {\\ensuremath{{%5.1f} } }", tag, nIsoFactE) << endl;

  if (!strcmp(tag, "s0")) {
    OUT  << Form("\\vdef{eTotEffIso:%s}   {\\ensuremath{{%s } } }", tag, (texForm(eIsoFact)).Data()) << endl;
  } else if (eIsoFact > 1.e-9) {
    OUT  << Form("\\vdef{eTotEffIso:%s}   {\\ensuremath{{%s } } }", tag, (texForm(eIsoFact)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{eTotEffIso:%s}   {\\ensuremath{{%s } } }", tag, (texForm2(eIsoFact)).Data()) << endl;
  }

  // -- fact. Vertex  
  OUT  << Form("\\vdef{nAfBchi2:%s}   {\\ensuremath{{%4.1f } } }", tag, nChiFact) << endl;
  OUT  << Form("\\vdef{nAfBchi2E:%s}  {\\ensuremath{{%5.1f} } }", tag, nChiFactE) << endl;

  if (!strcmp(tag, "s0")) {
    OUT  << Form("\\vdef{eTotEffChi2:%s}  {\\ensuremath{{%s } } }", tag, (texForm(eChiFact)).Data()) << endl;
  } else if (eChiFact > 1.e-9) {
    OUT  << Form("\\vdef{eTotEffChi2:%s}  {\\ensuremath{{%s } } }", tag, (texForm(eChiFact)).Data()) << endl;
  } else {
    OUT  << Form("\\vdef{eTotEffChi2:%s}  {\\ensuremath{{%s } } }", tag, (texForm2(eChiFact)).Data()) << endl;
  }

  // -- all Cuts
  OUT  << Form("\\vdef{nAllCutsFact:%s}  {\\ensuremath{{%5.1f} } }", tag, nAllCutsFact) << endl;
  OUT  << Form("\\vdef{nAllCutsFactE:%s}  {\\ensuremath{{%5.1f} } }", tag, nAllCutsFactE) << endl;

  if (!strcmp(tag, "s0")) {
    OUT  << Form("\\vdef{eAllCutsFact:%s}  {\\ensuremath{{%4.3f } } }", tag, eAllCutsFact) << endl;
  } else if (eAllCutsFact > 1.e-9 ) {
    OUT  << Form("\\vdef{eAllCutsFact:%s}  {\\ensuremath{{%s } } }", tag, (texForm(eAllCutsFact)).Data()) << endl;
  } else if (eAllCutsFact ) {
    OUT  << Form("\\vdef{eAllCutsFact:%s}  {\\ensuremath{{%s } } }", tag, (texForm2(eAllCutsFact)).Data()) << endl;
  } else  {
    OUT  << Form("\\vdef{eAllCutsFact:%s}  {\\ensuremath{{%4.3f } } }", tag, eAllCutsFact) << endl;
  }



  if (!strcmp(tag, "s0")) {
    OUT  << Form("\\vdef{eAllCutsFactE:%s} {\\ensuremath{{%4.3f } } }", tag, eAllCutsFactE) << endl;
  } else if (eAllCutsFactE > 1.e-9) {
    OUT  << Form("\\vdef{eAllCutsFactE:%s}  {\\ensuremath{{%s } } }", tag, (texForm(eAllCutsFactE).Data())) << endl;
  } else {
    OUT  << Form("\\vdef{eAllCutsFactE:%s}  {\\ensuremath{{%s } } }", tag, (texForm2(eAllCutsFactE).Data())) << endl;
  }



  // -- mass reduction 
  double massRed0  = massReduction("c030", tag);
  double massRed1  = massReduction("c130", tag);
  double massRed2  = massReduction("c230", tag);
  double massRed5  = massReduction("c530", tag);
  double massRed   = massReduction("c330", tag);

  

  double nMassAllCuts = massRed*SF*fnorm*eff1*eff2;
printf("  ****************nMassAllCuts = %4.4f massRed * %4.4f SF * %4.4f fnorm * %4.4f eff1 * %4.4f eff2\n", massRed, SF, fnorm, eff1, eff2);
  double nMassAllCutsE = 1.6*massRed*SF*TMath::Sqrt(fnorm)*eff1*eff2; // XXXX FIXME XXXX

  if ( comb ) {

    massRed = massRed5;

    if ( !strcmp(tag, "c0") ) {

      nMassAllCuts = fNbg + fNrbg;
      nMassAllCutsE = -9999;
    }
    else {

      nMassAllCuts = massRed*nAllCutsFact;
      nMassAllCutsE = -9999;
    }
    
    printf("  *****NEW********nMassAllCuts = %4.4f massRed * %4.4f = %4.4f\n", massRed, nAllCutsFact, nMassAllCuts);
  }

  OUT  << Form("\\vdef{massReduction:%s}  {\\ensuremath{{%3.2f} } }", tag, massRed) << endl;
  OUT  << Form("\\vdef{nMassAllCutsFact:%s}  {\\ensuremath{{%5.1f} } }", tag, nMassAllCuts) << endl;
  OUT  << Form("\\vdef{nMassAllCutsFactE:%s} {\\ensuremath{{%5.1f} } }", tag, nMassAllCutsE) << endl;

  if (!strcmp(tag, "s0")) {
    fNsg  = nAllCutsFact;
    fNsgE = nAllCutsFactE;
  }

  if (!strcmp(tag, "m0")) {
    fNbg = nMassAllCuts;
    fNbgE = nMassAllCutsE;
  }

  if (!strcmp(tag, "r0")) {
    fNrbg = nMassAllCuts;
    fNrbgE = nMassAllCutsE;
  }

}


// ----------------------------------------------------------------------
void anaBmm::singleHBG(const char *hist) {

  int EColor[]     = {   2,    4,    108,  6,  107,   93,   1, 
			 2,    4,    108,  6,  107,   93,   1,  
			 2,    4,    108,  6,  107,   93,   1   }; 
 
 
  cout << endl << endl;
    
  TH1D *hs = (TH1D*)fS[0]->Get(hist)->Clone();
  setFilledHist(hs, kBlack, kYellow, 1001, 2);
  //  hs->Scale(1./hs->GetSumOfWeights());
  //  hs->GetYaxis()->SetRangeUser(0.0001,2);

  setTitles(hs, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin", 0.06, 1.1, 1.3);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetTopMargin(0.2);
  shrinkPad(0.2, 0.2);

  gPad->SetLogy(1);

  legg = new TLegend(0.45,0.79,0.8,0.89);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.035);  legg->SetFillColor(0); 
  
  const int nhist = nMc;
  TH1D *hbg[nhist];

  for (int i = 1; i < nhist; i++) {
 
    printf("anaBmm::singleHBG> Getting histogram %s for %s\n", hist, fSignM[i].Data());
    hbg[i] = (TH1D*)fM[i]->Get(hist)->Clone();
    //    hbg[i]->Scale(1./hbg[i]->GetSumOfWeights());

    printf("anaBmm::singleHBG> Scaling %4.2f to luminosity of %4.2e * %4.1f / %4.6f\n"
	   ,  hbg[i]->GetSumOfWeights(), fMisIdM[i], fLumiD[0], fLumiM[i]);

    hbg[i]->Scale(fMisIdM[i]*fLumiD[0]/fLumiM[i]);

    printf(" = %4.2f\n",  hbg[i]->GetSumOfWeights());

    setHist(hbg[i],  EColor[i]);
     
    if ( hbg[i]->GetSumOfWeights() > 0  ) {

      hs->Scale(hbg[i]->GetSumOfWeights()/hs->GetSumOfWeights());

      if ( hbg[i]->GetMaximum() > 0 ) {

	hs->GetYaxis()->SetRangeUser(0.1*hbg[i]->GetMinimum(1.e-20), 100*hbg[i]->GetMaximum());
	hs->DrawCopy("hist");
	hbg[i]->DrawCopy("same");
      }
    }
    else {

      hs->DrawCopy("hist");
    }

    legg->Clear();
    
    legge = legg->AddEntry(hs, "Signal", "f"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hbg[i], "Background", "p");

//     if ( hbg[i]->GetSumOfWeights() < 1000 ) { 

//       legge = legg->AddEntry(hbg[i], Form("Background (%4.1f)", hbg[i]->GetSumOfWeights()), "p");
//     }
//     else {

//       legge = legg->AddEntry(hbg[i], Form("Background (%2.1e)", hbg[i]->GetSumOfWeights()), "p");
//     }
 
    legge->SetTextColor(EColor[i]);
    legg->Draw();
    
    tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
    tl->SetTextColor(kBlack);    
    tl->DrawLatex(0.22, 0.92, Form("%s",  fSignTitleM[i].Data()));

    c0->SaveAs(Form("%s/mass-%s-%s.eps", outDir, hist, fSignM[i].Data()));
    
    cout << endl;
  }

  cout << endl << endl;

}  
// ----------------------------------------------------------------------
void anaBmm::overlayHBG(const char *hist, const int npers) {

  int EColor[]     = {   13,    2,    4,    108,  6,  107,   93,   1, 
			 13,    2,    4,    108,  6,  107,   93,   1   }; 
  
  double minMID(9999.), maxMID(0.), maxSoW(0.);
  double fnorm_Ch, eff1_Ch, eff2_Ch;  

  // -- signal
  printf("anaBmm::overlayHBG> Getting histogram %s for %s\n", hist, fSignS[0].Data());
  TH1D *hs = (TH1D*)fS[0]->Get(hist)->Clone();
  printf("anaBmm::overlayHBG> Scaling %4.2f to luminosity of  %f / %f\n",  hs->GetSumOfWeights(), fLumiD[0], fLumiS[0]);
  if ( !strcmp(hist,"c530") ) {

    channelEff(fS[0], fnorm_Ch, eff1_Ch, eff2_Ch);
    hs->Scale((fLumiD[0]/fLumiS[0])*eff1_Ch*eff2_Ch);
  }
  else {
    
    hs->Scale(fLumiD[0]/fLumiS[0]);
  } 
  printf(" = %4.2f\n\n",  hs->GetSumOfWeights());
  setFilledHist(hs, kBlack, kYellow, 1001, 2);
  setTitles(hs, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin", 0.06, 1.1, 1.3); 

  if (hs->GetMinimum(1.e-20) < minMID ) { minMID = hs->GetMinimum(1.e-20);      }
  if (hs->GetMaximum()       > maxMID ) { maxMID = hs->GetMaximum();            }
 
 

  // -- background
  printf("anaBmm::overlayHBG> Getting histogram %s for %s\n", hist, fSignM[0].Data());
  TH1D *hb = (TH1D*)fM[0]->Get(hist)->Clone();
  printf("anaBmm::overlayHBG> Scaling %4.2f to luminosity of %f / %f\n",  hb->GetSumOfWeights(), fLumiD[0], fLumiM[0]);
  if ( !strcmp(hist,"c530") ) {

    channelEff(fM[0], fnorm_Ch, eff1_Ch, eff2_Ch);
    hb->Scale((fLumiD[0]/fLumiM[0])*eff1_Ch*eff2_Ch);
  }
  else {
    
    hb->Scale(fLumiD[0]/fLumiM[0]);
  }
  printf(" = %4.2f\n\n",  hb->GetSumOfWeights());
  hb->SetLineStyle(kDashed);
  setFilledHist(hb, 13, 13, 3005, 2);

  if (hb->GetMinimum(1.e-20) < minMID ) { minMID = hb->GetMinimum(1.e-20);      }
  if (hb->GetMaximum()       > maxMID ) { maxMID = hb->GetMaximum();            }



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetTopMargin(0.25);
  gPad->SetRightMargin(0.2);

  gPad->SetLogy(1);

  // -- HBG
  const int nhist = nMc;
  const int sets = int(nhist/npers)+1;

  TH1D *hMid[nhist];


  for (int i = 1; i < nhist; ++i) {

    printf("anaBmm::overlayHBG> Getting histogram %s for %s\n", hist, fSignM[i].Data());

    hMid[i] = (TH1D*)fM[i]->Get(hist)->Clone();

    printf("anaBmm::overlayHBG> Scaling %4.2f to luminosity of %4.2e * %4.1f / %4.4f\n"
	   ,hMid[i]->GetSumOfWeights() , fMisIdM[i], fLumiD[0], fLumiM[i]);

    if ( !strcmp(hist,"c530") ) {

      channelEff(fM[i], fnorm_Ch, eff1_Ch, eff2_Ch);
      hMid[i]->Scale((fMisIdM[i]*fLumiD[0]/fLumiM[i])*eff1_Ch*eff2_Ch);
    
    }
    else {

      hMid[i]->Scale(fMisIdM[i]*fLumiD[0]/fLumiM[i]);
    }

    printf(" = %4.2f\n\n",  hMid[i]->GetSumOfWeights());

    if (hMid[i]->GetSumOfWeights()  > maxSoW )  { maxSoW = hMid[i]->GetSumOfWeights();  }
    if (hMid[i]->GetMinimum(1.e-20) < minMID )  { minMID = hMid[i]->GetMinimum(1.e-20); }
    if (hMid[i]->GetMaximum()       > maxMID )  { maxMID = hMid[i]->GetMaximum();       }
  }

  // -- Draw & Legends
  legg = new TLegend(0.49,0.62,0.99,0.99);
  legg->SetFillStyle(1001); legg->SetBorderSize(0); legg->SetTextSize(0.035);  legg->SetFillColor(10); 

  int i(0), cont (0);

  for (int k = 0; (k < sets) && (i < nhist); k++ ) {

    legge = legg->AddEntry(hs, Form("%s", fSignTitleS[0].Data()), "f");
    legge = legg->AddEntry(hb, Form("%s",fSignTitleM[0].Data()), "f");

  //   legge = legg->AddEntry(hs, Form("%s (%4.1f)", fSignTitleS[0].Data(), hs->GetSumOfWeights()), "f");

//     if ( hb->GetSumOfWeights() > 1000 ) {
      
//       legge = legg->AddEntry(hb, Form("%s (%4.1e)",fSignTitleM[0].Data(), hb->GetSumOfWeights()), "f");
//     }
//     else if ( hb->GetSumOfWeights() < 0.1 ) {
      
//       legge = legg->AddEntry(hb, Form("%s (%4.1e)",fSignTitleM[0].Data(), hb->GetSumOfWeights()), "f");
//     }
//     else if ( hb->GetSumOfWeights() ) {
      
//       legge = legg->AddEntry(hb, Form("%s (%4.1f)",fSignTitleM[0].Data(), hb->GetSumOfWeights()), "f");
//     }

    // -- Draw plot
    hs->GetYaxis()->SetRangeUser(0.1*minMID, 100*maxMID);
 
    hs->DrawCopy("hist"); 
    hb->DrawCopy("histsame");    
    
    cont = 1;
    for (i = npers*k+1; (cont < npers+1) && (i < nhist); ++i) {
      

      if ( hMid[i]->GetSumOfWeights() > 0 && hMid[i]->GetMaximum() > 0 ) {
	
	cont++;

	setHist(hMid[i], EColor[cont]);

	legge = legg->AddEntry(hMid[i], Form("%s", fSignTitleM[i].Data()), "p");

// 	if ( hMid[i]->GetSumOfWeights() > 1000 ) {
// 	  legge = legg->AddEntry(hMid[i], Form("%s (%4.1e)", fSignTitleM[i].Data(), hMid[i]->GetSumOfWeights()), "p"); 
// 	}
// 	else if (hMid[i]->GetSumOfWeights() < 0.1 ) {
// 	  legge = legg->AddEntry(hMid[i], Form("%s (%4.1e)", fSignTitleM[i].Data(), hMid[i]->GetSumOfWeights()), "p"); 
// 	}
// 	else {
// 	  legge = legg->AddEntry(hMid[i], Form("%s (%4.1f)", fSignTitleM[i].Data(), hMid[i]->GetSumOfWeights()), "p"); 
// 	}
	
	legge->SetTextColor(kBlack);  legge->SetFillStyle(1001);  legge->SetFillColor(10); 
		
	hMid[i]->DrawCopy("same");

      }
    }
    
    legg->Draw();

    c0->SaveAs(Form("%s/hbgOverlay%i-%s.eps", outDir, k+1, hist));

    legg->Clear();

  }

  gPad->SetLogy(0);

}

// ----------------------------------------------------------------------
void anaBmm::stackHBG(const char *hist) {
 
  int EColor[]     = {   13,    2,    4,    108,  6,  107,   93,   1, 
			 13,    2,    4,    108,  6,  107,   93,   1   }; 

  int EFillStyle[] = { 3004, 3005, 3004, 3005, 3004, 3005, 3004, 3005, 3004, 3005, 3004, 3005, 3004, 3005 }; 

 
  THStack hStackHBG("hStackHBG","stack with rare bg");

  // -- signal
  TH1D *hs = (TH1D*)fS[0]->Get(hist)->Clone();
  hs->Scale(fLumiD[0]/fLumiS[0]);
  setFilledHist(hs, kBlack, kYellow, 1001, 2);
  setTitles(hs, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin", 0.06, 1.1, 1.3); 

  // -- background
  TH1D *hb = (TH1D*)fM[0]->Get(hist)->Clone();
  hb->Scale(fLumiD[0]/fLumiM[0]);
  setFilledHist(hb, 13, 13, 3005, 2);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetTopMargin(0.25);
  gPad->SetRightMargin(0.2);

  gPad->SetLogy(1);
  
  const int nhist = nMc;
  TH1D *hMid[nhist];

  double minMID(9999.), maxMID(0.), maxSoW(0.);

  // -- HBG
  for (int i = 1; i < nhist; ++i) {

    hMid[i] = (TH1D*)fM[i]->Get(hist)->Clone();

    if (hMid[i]->GetSumOfWeights()  > maxSoW )  { maxSoW = hMid[i]->GetSumOfWeights(); }
    if (hMid[i]->GetMinimum(10e-20) < minMID )  { minMID = hMid[i]->GetMinimum(10e-20); }
    if (hMid[i]->GetMaximum()       > maxMID )  { maxMID = hMid[i]->GetMaximum();       }
    
    hMid[i]->Scale(fMisIdM[i]*fLumiD[0]/fLumiM[i]);
    setFilledHist(hMid[i], EColor[i], EColor[i], 0, 2);
  }

  // -- Scaling
  int scaleS = int(TMath::Log(maxSoW/hs->GetSumOfWeights())/TMath::Ln10());
  hs->Scale(TMath::Power(10,scaleS));

  int scaleM = int(TMath::Log(maxSoW/hb->GetSumOfWeights())/TMath::Ln10());
  hb->Scale(TMath::Power(10,scaleM));

  if (hs->GetMinimum(10e-20) < minMID ) { minMID = hs->GetMinimum(10e-20);      }
  if (hs->GetMaximum()       > maxMID ) { maxMID = hs->GetMaximum();            }
  if (hb->GetMinimum(10e-20) < minMID ) { minMID = hb->GetMinimum(10e-20);      }
  if (hb->GetMaximum()       > maxMID ) { maxMID = hb->GetMaximum();            }

  // -- Legends
  legg = new TLegend(0.65,0.6,0.8,0.93);
  legg->SetFillStyle(1001); legg->SetBorderSize(0); legg->SetTextSize(0.035);  legg->SetFillColor(10); 

  TString fact = scaleFactor(scaleS);
  legge = legg->AddEntry(hs, Form("%s %s", fSignTitleS[0].Data(), fact.Data()), "f");
  fact  =  scaleFactor(scaleM);
  legge = legg->AddEntry(hb, Form("%s %s", fSignTitleM[0].Data(), fact.Data()), "f");
 
  for (int i = 1; i < nhist; ++i) {

    legge = legg->AddEntry(hMid[i], Form("%s", fSignTitleM[i].Data()), "f"); 
    legge->SetTextColor(kBlack);  legge->SetFillStyle(1001);  legge->SetFillColor(10); 

  }

  // -- Stack histogams
  for (int i = 1; i < nhist; ++i) {
    
    if (!TMath::IsNaN(hMid[i]->GetSumOfWeights())) { hStackHBG.Add(hMid[i]); }
    
  }
    
  if ( hStackHBG.GetMaximum() > maxMID ) { maxMID = hStackHBG.GetMaximum();          }
  
  // -- Draw plot
  hs->GetYaxis()->SetRangeUser(0.1*minMID, 100*maxMID);

  hs->DrawCopy("Ahist"); 
  hStackHBG.Draw("histsame");
 
  if (!TMath::IsNaN(hb->GetSumOfWeights())) { hb->DrawCopy("histsame"); }


  legg->Draw();

  c0->SaveAs(Form("%s/hbgStack-%s.eps", outDir, hist));

  gPad->SetLogy(0);

}


// ----------------------------------------------------------------------
void anaBmm::plotAll(const char *hist) {
 
  // -- signal
  TH1D *hs = (TH1D*)fS[0]->Get(hist)->Clone();
  hs->Scale(fLumiD[0]/fLumiS[0]);
  setFilledHist(hs, kBlue, kBlue, 3004, 2);
  setTitles(hs, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin", 0.06, 1.1, 1.3); 

  // -- background
  TH1 *hb = sumHistMC(hist, 2);
  setFilledHist(hb, kBlack, kBlack, 3005, 2);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetTopMargin(0.2);
  gPad->SetLogy(1);

  // -- Scaling
  double minMID(9999.), maxMID(0.);
  int scaleS = int(TMath::Log(hb->GetSumOfWeights()/hs->GetSumOfWeights())/TMath::Ln10());
  hs->Scale(TMath::Power(10,scaleS));

  if (hs->GetMinimum(10e-20) < minMID )  { minMID = hs->GetMinimum(10e-20);      }
  if (hs->GetMaximum()       > maxMID )  { maxMID = hs->GetMaximum();            }
  if (hb->GetMinimum(10e-20) < minMID )  { minMID = hb->GetMinimum(10e-20);      }
  if (hb->GetMaximum()       > maxMID )  { maxMID = hb->GetMaximum();            }

  // -- Legends
  legg = new TLegend(0.65,0.8,0.8,0.89);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.035);  legg->SetFillColor(0);   
  if ( scaleS == 0 ) {

    legge = legg->AddEntry(hs, Form("%s", fSignTitleS[0].Data()), "f"); 
  }
  else if ( scaleS == 0 ) {

    legge = legg->AddEntry(hs, Form("%s (x 10)", fSignTitleS[0].Data()), "f"); 
  }
  else { 

    legge = legg->AddEntry(hs, Form("%s (x10^{%i})", fSignTitleS[0].Data(), scaleS), "f");
  }

  legge->SetTextColor(kBlack);  legge->SetFillStyle(1001);
  legge = legg->AddEntry(hb, "Background (sum)", "f"); 
  legge->SetTextColor(kBlack);  legge->SetFillStyle(3004);  


  // -- Draw plots
  hs->GetYaxis()->SetRangeUser(0.1*minMID, 10*maxMID);

  hs->DrawCopy("hist");
  hb->DrawCopy("histsame");
  legg->Draw();

  c0->SaveAs(Form("%s/mass-%s.eps", outDir, hist));

  gPad->SetLogy(0);
}

// ----------------------------------------------------------------------
double anaBmm::massReduction(const char *hist, const char *sel) {
  TTree *t = (TTree*)fM[0]->Get("events");

  if ( !strcmp(sel, "m0") ) {

    TH1D *h1 = (TH1D*)fM[0]->Get(hist)->Clone();
    h1->GetXaxis()->SetRangeUser(5., 5.7);
    h1->Fit("pol1");
    c0->SaveAs(Form("%s/massReduction-%s.eps", outDir, hist));
    
    TF1 *f1 = (TF1*)h1->GetFunction("pol1");
    double total = f1->Integral(h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()))/h1->GetBinWidth(2); 
    double massW = f1->Integral(fMassBs - 0.1, fMassBs + 0.1)/h1->GetBinWidth(2);
    double massRed = massW/total;
    
    double histT =  h1->GetSumOfWeights();
    double histW =  h1->Integral(h1->FindBin(fMassBs - 0.1), h1->FindBin(fMassBs + 0.1));
    double histRed = histW/histT;

    cout << "->" << hist<< " massRed: " << massRed << " hRed = " << histRed 
	 << " massW = "  << massW 
	 << " total: " << total 
	 << " hW: " << histW
	 << " htotal: " << histT
	 << " fMassBs: " << fMassBs
	 << endl;

    h1->GetListOfFunctions()->Clear();
  
    return massRed;

  }  

  else if ( !strcmp(sel, "s0") ) {

    TH1D *h1 = (TH1D*)fS[0]->Get(hist)->Clone();
    h1->GetXaxis()->SetRangeUser(5., 5.7);
    h1->Fit("gaus"); // XXX FIX ME > Use proper fitting function XXX
    c0->SaveAs(Form("%s/massReduction-%s-sg.eps", outDir, hist));
    
    TF1 *f1 = (TF1*)h1->GetFunction("gaus");
    double total = f1->Integral(h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()))/h1->GetBinWidth(2); 
    double massW = f1->Integral(fMassBs - 0.1, fMassBs + 0.1)/h1->GetBinWidth(2);
    double massRed = massW/total;
    
    double histT =  h1->GetSumOfWeights();
    double histW =  h1->Integral(h1->FindBin(fMassBs - 0.1), h1->FindBin(fMassBs + 0.1));
    double histRed = histW/histT;
    
    cout << "->" << hist<< " massRed: " << massRed << " hRed = " << histRed 
	 << " massW = "  << massW 
	 << " total: " << total 
	 << " hW: " << histW
	 << " htotal: " << histT
	 << " fMassBs: " << fMassBs
	 << endl;
    
    h1->GetListOfFunctions()->Clear();
    
    return massRed;

  }
  else {  
    
    double histT(0.), histW(0.), histRed(0.);
    double fnorm_Ch(0.), eff1_Ch(0.), eff2_Ch(0.);
  
    int nhist = nMc;
    TH1D *h1[nhist];

    for (int i = 0; i < nhist; ++i) {

      if (  !strcmp(sel,"r0") ) {

	if ( i == 0 ) {
	
	  continue;
	}
      }

      else if ( strcmp(sel,"c0") ) {
	
	if ( strcmp((getRareType(fSignM[i])).Data(), sel) ) {
	  
	  continue;
	}
      }

      h1[i]  = (TH1D*)fM[i]->Get(hist)->Clone();

      channelEff(fM[i], fnorm_Ch, eff1_Ch, eff2_Ch);
      h1[i]->Scale((fMisIdM[i]*fLumiD[0]/fLumiM[i])*eff1_Ch*eff2_Ch);

      histT +=  h1[i]->GetSumOfWeights();
      histW +=  h1[i]->Integral(h1[i]->FindBin(fMassBs - 0.1), h1[i]->FindBin(fMassBs + 0.1));

    }

    if ( histT ) { histRed = histW/histT; }

    cout <<  "\t*** " << sel << "->" << hist << " hRed = " << histRed 
	 << " hW: " << histW
	 << " htotal: " << histT
	 << " fMassBs: " << fMassBs
	 << endl;
  
    return histRed;
  }
}

// ----------------------------------------------------------------------
void anaBmm::runOptimization(const char *aoCuts, const char *extraVar, int nbin, double min, double max) {
  char line[2000]; 
  double cut, maxFom(-99.), maxCut(0.); 
  for (int i = 0; i < nbin; ++i) {
    cut = min + i*(max-min)/nbin; 
    handOptimization(aoCuts, Form("%s %5.4f", extraVar, cut));
    if (fFom > maxFom) {
      maxFom = fFom; 
      maxCut = cut;
    }
  }

  cout << "==> maximum fom: " << maxFom << " at cut " << extraVar << maxCut << endl;
}



// ----------------------------------------------------------------------
void anaBmm::handOptimization(const char *aoCuts, const char *extraCuts) {
 
  char filename[200];
  sprintf(filename, "%s/handOptimization.txt", outDir); 
  
  //   system(Form("/bin/mv %s old.%s", filename, filename));
  ofstream OUT(filename, ios::app);

  // -- Fix lumi normalisation scaling factor
  double SF = fLumiD[0]/fLumiM[0];
  
  // -- Run on signal MC to determine efficiency 
  fS[0]->cd(); 
  TH1D *hSG = (TH1D*)gROOT->FindObject("hSG"); 
  if (!hSG) hSG = new TH1D("hSG", "", 50, 5., 6.);
  TTree *s = (TTree*)gFile->Get("events");
  s->Draw("mass>>hSG", Form("goodKinematics"), "goff"); 
  double sNorm = hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("%s", aoCuts), "goff"); 
  double s1Norm = hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("%s", extraCuts), "goff"); 
  double s2Norm  =  hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("%s && %s", aoCuts, extraCuts), "goff"); 
  double s12Norm  =  hSG->GetSumOfWeights();

  // -- Run on bg MC
  fM[0]->cd(); 
  TH1D *hBG = (TH1D*)gROOT->FindObject("hBG"); 
  if (!hBG) hBG = new TH1D("hBG", "", 50, 5., 6.);
  TTree *b = (TTree*)gFile->Get("events");
  b->Draw("mass>>hBG", Form("goodKinematics"), "goff"); 
  double bNorm   = hBG->GetSumOfWeights();
  b->Draw("mass>>hBG", Form("%s", aoCuts), "goff"); 
  double b1Norm  = hBG->GetSumOfWeights();
  b->Draw("mass>>hBG", Form("%s", extraCuts), "goff"); 
  double b2Norm  = hBG->GetSumOfWeights();
  double e2      = b2Norm/bNorm;
  b->Draw("mass>>hBG", Form("%s && %s", aoCuts, extraCuts), "goff"); 
  double b12Norm  =  hSG->GetSumOfWeights();
  
  fFom = (s12Norm/sNorm) / (1. + TMath::Sqrt(e2*b1Norm*SF)) ;
  
  cout << Form("%4.2e s:%4.3f/%4.3f b:%4.3e/%4.3e #bg:%5.2f (%3.1e,%3.1e)", 
	       fFom, s12Norm/sNorm, (s1Norm/sNorm)*(s2Norm/sNorm), b12Norm/bNorm, (b1Norm/bNorm)*(b2Norm/bNorm), e2*b1Norm*SF, b1Norm, b2Norm)
       << " cuts: \"" << aoCuts << "\", \"" << extraCuts << "\" " << endl;

  OUT << Form("%4.2e e12/e1e2:%4.3f/%4.3f bg:%5.2f (%3.1e,%3.1e)", 
	      fFom, s12Norm/sNorm, (s1Norm/sNorm)*(s2Norm/sNorm), e2*b1Norm*SF, b1Norm, b2Norm)
      << " cuts: \"" << aoCuts << "\", \"" << extraCuts << "\" " << endl;

}


// ----------------------------------------------------------------------
// -- loop for stupid loop over all  cuts
// -- ptl1, chi2, lxy/sxy, iso, cosa
void anaBmm::loopOptimization(double pt) {

  double pTlo;
  TH1D *h = (TH1D*)fS[0]->Get("hcuts");
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    if (!strcmp(h->GetXaxis()->GetBinLabel(i), "p_{T}^{min}(l) [GeV]")) {
      pTlo = h->GetBinContent(i);
    }
  }

  char filename[200];
  sprintf(filename, "%s/loopOptimization-pt:%3.2f.txt", outDir, pt); 

  system(Form("/bin/rm -f %s", filename));
  ofstream OUT(filename, ios::app);


  char line[2000];
  char precuts[2000]; 


  sprintf(precuts, Form("goodL1 && goodKinematics && ptl1 > %f", pt));
  cout << "Running with precuts: " << precuts << endl;

  double eff, bg, sg;
  double cut1(0.), cut2(0.),   cut3(0.),   cut4(0.),   cut5(0.); 
  char scut1[200], scut2[200], scut3[200], scut4[200], scut5[200];
  double fom(-1);

  int IMAX1(5), IMAX2(5), IMAX3(5), IMAX4(5),  IMAX5(1); 
  TH1D *hOpt = new TH1D("hOpt", "", IMAX1*IMAX2*IMAX3*IMAX4*IMAX5, 0., IMAX1*IMAX2*IMAX3*IMAX4*IMAX5);
  
  // -- Significance of UL: figure of merit should be maximized
  double a(1.64); // a = 2:    95% CL UL
                  // a = 1.64: 90% CL UL

  double maxFom(-9.), i1max, i2max; 

  for (int i1 = 0; i1 < IMAX1; ++i1) {
    cut1 = 10. + i1*3.;
    sprintf(scut1, "lxy/sxy > %f", cut1);
    cout << "i1: " << scut1 << endl;

    for (int i2 = 0; i2 < IMAX2; ++i2) {
      cut2 = 0.95 + i2*0.005;
      sprintf(scut2, "cosa > %f", cut2);
      cout << " i2: " << scut2 << endl;

      for (int i3 = 0; i3 < IMAX3; ++i3) {
	cut3 = 0.7 + i3*0.05;
	sprintf(scut3, "iso > %f", cut3);
	cout << "i3: " << scut3 << endl;

	for (int i4 = 0; i4 < IMAX4; ++i4) {
	  cut4 = 1.0 + i4*0.5;
	  sprintf(scut4, "chi2 < %f", cut4);
	
	  sprintf(line, "%s && %s && %s && %s", scut1, scut2, scut3, scut4); 
	  optimizerNumbers(precuts, line, eff, sg, bg);
	  
	  fom = eff/(a + TMath::Sqrt(bg));
	  hOpt->SetBinContent(i1*IMAX1 + i2, fom);

	  if (fom > maxFom) {
	    maxFom = fom; 
	    i1max = i1;
	    i2max = i2; 
	  }
      
	  OUT << Form("lxy: %4.3f /cosa: %6.5f /iso: %4.3f /chi2: %4.3f /eff: %5.2f /SG: %5.2f /BG: %5.2f /fom: %7.6f", 
		      cut1, cut2, cut3, cut4, eff, sg, bg, fom) 
	      << endl;
	}
      }
    }
  }

  hOpt->Draw();

  cout << "Maximum fom: " << maxFom << " at lxy/sxy > " << 5. + i1max*2 << " and cos > " << 0.95 + i2max*0.005 << endl;

}





// ----------------------------------------------------------------------
void anaBmm::optimizerNumbers(const char *precuts, const char *line, double &eff, double &sg, double &bg) {

  // -- Run on signal MC to determine efficiency 
  fS[0]->cd(); 
  TH1D *hSG = (TH1D*)gROOT->FindObject("hSG"); 
  if (!hSG) hSG = new TH1D("hSG", "", 50, 5., 6.);
  TTree *s = (TTree*)gFile->Get("events");
  s->Draw("mass>>hSG", Form("%s", precuts), "goff"); 
  double norm = hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("%s && %s", precuts, line), "goff"); 
  sg =  hSG->GetSumOfWeights();
  eff = sg/norm;

  // -- And now MC for background expectation
  fM[0]->cd(); 
  TH1D *hBG = (TH1D*)gROOT->FindObject("hBG"); 
  if (!hBG) hBG = new TH1D("hBG", "", 50, 5., 6.);
  TTree *b = (TTree*)gFile->Get("events");
  b->Draw("mass>>hBG", Form("%s && %s", precuts, line), "goff"); 
  bg =  hBG->GetSumOfWeights();

  //   cout << " sg: " << sg << " and norm: " << norm << "bg: " << bg << endl;

}


// ----------------------------------------------------------------------
void anaBmm::breco(int o) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TH1D *h1;
  h1 = (TH1D*)(fS[o]->Get("c030"))->Clone();

  cout << endl << " ENTRIES " << h1->GetSumOfWeights() << endl << endl;
  h1->Scale(fLumiD[0]/fLumiS[o]);
  cout << endl << " ENTRIES " << h1->GetSumOfWeights() << endl << endl;

  if ( o < 3 ) {
    setTitles(h1, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
    h1->SetMaximum(1.2*h1->GetMaximum());
    h1->SetMaximum(45);
  }
  else {
    setTitles(h1, "m_{#mu^{+}#mu^{-}K^{+}} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
    h1->SetMaximum(1.2*h1->GetMaximum());
    h1->SetMaximum(135);
  }

  setFilledHist(h1, kBlack, kYellow, 1000);

  c0->Clear();
  shrinkPad(0.15, 0.15);
  setFilledHist(h1, kBlack, kYellow, 1000);
  h1->DrawCopy("hist");

 
//     f1->SetParameters(0.8*h1->GetMaximum(), h1->GetMean(), h1->GetRMS(), 
// 		      0.2*h1->GetMaximum(), h1->GetMean(), 3.*h1->GetRMS());
//     h1->Fit(f1);
 
  f2->SetParameters(h1->GetMaximum()*h1->GetRMS(), h1->GetMean(), 0.5*h1->GetRMS(), 
		    1.0, h1->GetMean(), 3.*h1->GetRMS());
    
  f1->SetParameters(h1->GetMaximum()*0.8, h1->GetMean(), 0.5*h1->GetRMS(), 
		    h1->GetMaximum()*0.2, h1->GetMean(), 3.*h1->GetRMS());

  f4->SetParameters(h1->GetMaximum()*0.8, h1->GetMean(), 0.5*h1->GetRMS(), 
		    h1->GetMaximum()*0.2, h1->GetMean(), 3.*h1->GetRMS());

  double mean(0.), sigma(0.);
 
  if ( o < 3 ) {
    h1->Fit(f1, "0");
    f1->DrawCopy("same");
    mean = TMath::Sqrt((f1->GetParameter(0)*f1->GetParameter(0)*f1->GetParameter(1)*f1->GetParameter(1)
			     + f1->GetParameter(3)*f1->GetParameter(3)*f1->GetParameter(4)*f1->GetParameter(4))
			    /
			    (f1->GetParameter(0)*f1->GetParameter(0) + f1->GetParameter(3)*f1->GetParameter(3))
			    );
    sigma = TMath::Sqrt((f1->GetParameter(0)*f1->GetParameter(0)*f1->GetParameter(2)*f1->GetParameter(2)
			 + f1->GetParameter(3)*f1->GetParameter(3)*f1->GetParameter(5)*f1->GetParameter(5))
			/
			(f1->GetParameter(0)*f1->GetParameter(0) + f1->GetParameter(3)*f1->GetParameter(3))
			);
  } else {
    h1->Fit(f4, "0");
    f4->DrawCopy("same");
    mean  = f4->GetParameter(1);
    sigma = f4->GetParameter(2);
  }

  tl->SetTextColor(kBlack);  
  tl->SetNDC(kTRUE); tl->SetTextSize(0.06);


  tl->DrawLatex(0.16, 0.85, Form("#mu: %5.3f#pm%5.3f GeV", mean, 0.001));
  tl->DrawLatex(0.16, 0.79, Form("#sigma: %5.3f#pm%5.3f GeV", sigma, 0.001));
  tl->DrawLatex(0.16, 0.73, Form("N: %4.0f", h1->GetSumOfWeights()));
//   tl->DrawLatex(0.16, 0.85, Form("#mu_{1} = %5.3f#pm%5.3f GeV", f2->GetParameter(1), f2->GetParError(1)));
//   tl->DrawLatex(0.16, 0.80, Form("#mu_{2} = %5.3f#pm%5.3f GeV", f2->GetParameter(4), f2->GetParError(4)));
//   tl->DrawLatex(0.16, 0.75, Form("#sigma_{1} = %5.3f#pm%5.3f GeV", f2->GetParameter(2), f2->GetParError(2)));
//   tl->DrawLatex(0.16, 0.70, Form("#sigma_{2} = %5.3f#pm%5.3f GeV", f2->GetParameter(5), f2->GetParError(5)));
//   tl->DrawLatex(0.16, 0.65, Form("g_{2}/g_{1} = %5.3f#pm%5.3f ", f2->GetParameter(3), f2->GetParError(3)));

  tl->SetTextSize(0.06); tl->SetTextColor(kBlue);

  if ( o == 1 ) {
    tl->DrawLatex(0.50, 0.72, "official MC");
  } else if ( o == 2 ) {
    tl->DrawLatex(0.60, 0.70, "CMSSW");
  } else if ( o == 0 ) {
    tl->DrawLatex(0.48, 0.72, "ORCA/OSCAR");
  } else if ( o == 3 ) {
    tl->DrawLatex(0.56, 0.72, "perfect align.");
    tl->DrawLatex(0.60, 0.65, "B #rightarrow J/#Psi K");
  } else if ( o == 4 ) {
    tl->DrawLatex(0.52, 0.72, "short-term align.");
    tl->DrawLatex(0.56, 0.65, "B #rightarrow J/#Psi K");
  } else {  
    tl->SetTextSize(0.1);
    tl->DrawLatex(0.70, 0.82, Form("MC%i",o));
  }

  delete h1;

  fMassBs = mean;

  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% -- BRECO mass peak numbers" << endl;
  OUT << Form("\\vdef{mBgMean:s0} {\\ensuremath {%4.1f } }", 1000.*mean) << endl;
  OUT << Form("\\vdef{mBgSigma:s0}{\\ensuremath {%4.1f } }", 1000.*sigma) << endl;

  OUT << Form("\\vdef{mBg1n:s0} {\\ensuremath {%4.2f } }", f1->GetParameter(0)) << endl;
  OUT << Form("\\vdef{mBg1nE:s0}{\\ensuremath {%4.2f } }", f1->GetParError(0)) << endl;

  OUT << Form("\\vdef{mBg1m:s0}  {\\ensuremath {%4.1f } }", 1000.*f1->GetParameter(1)) << endl;
  OUT << Form("\\vdef{mBg1mE:s0} {\\ensuremath {%4.1f } }", 1000.*f1->GetParError(1)) << endl;

  OUT << Form("\\vdef{mBg1s:s0}  {\\ensuremath {%4.1f } }", 1000.*f1->GetParameter(2)) << endl;
  OUT << Form("\\vdef{mBg1sE:s0} {\\ensuremath {%4.1f } }", 1000.*f1->GetParError(2)) << endl;

  OUT << Form("\\vdef{mBg2n:s0} {\\ensuremath {%4.2f } }", f1->GetParameter(3)) << endl;
  OUT << Form("\\vdef{mBg2nE:s0}{\\ensuremath {%4.2f } }", f1->GetParError(3)) << endl;

  OUT << Form("\\vdef{mBg2m:s0}  {\\ensuremath {%4.1f } }", 1000.*f1->GetParameter(4)) << endl;
  OUT << Form("\\vdef{mBg2mE:s0} {\\ensuremath {%4.1f } }", 1000.*f1->GetParError(4)) << endl;

  OUT << Form("\\vdef{mBg2s:s0}  {\\ensuremath {%4.1f } }", 1000.*f1->GetParameter(5)) << endl;
  OUT << Form("\\vdef{mBg2sE:s0} {\\ensuremath {%4.1f } }", 1000.*f1->GetParError(5)) << endl;


  c0->SaveAs(Form("%s/s%d-breco.eps", outDir, o));

}

// ----------------------------------------------------------------------
void anaBmm::jreco(int o) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TH1D *h1;
  h1 = (TH1D*)(fS[o]->Get("j101"))->Clone();

  cout << endl << " ENTRIES " << h1->GetSumOfWeights() << endl << endl;
  h1->Scale(fLumiD[0]/fLumiS[o]);
  cout << endl << " ENTRIES " << h1->GetSumOfWeights() << endl << endl;

  setTitles(h1, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
  setFilledHist(h1, kBlack, kYellow, 1000);

  c0->Clear();
  shrinkPad(0.15, 0.15);
  setFilledHist(h1, kBlack, kYellow, 1000);
  h1->SetMaximum(1.3*h1->GetMaximum());
  h1->SetMaximum(600);
  h1->DrawCopy("hist");

  f0->SetParameters(0.8*h1->GetMaximum(), h1->GetMean(), h1->GetRMS(), 
 		    0.2*h1->GetMaximum(), h1->GetMean(), 3.*h1->GetRMS());
  h1->Fit(f0, "0");

//   f2->SetParameters(h1->GetMaximum()*h1->GetRMS(), h1->GetMean(), 0.5*h1->GetRMS(), 
// 		    1.0, h1->GetMean(), 3.*h1->GetRMS());

//   f0->SetParameters(h1->GetMaximum()*0.8, h1->GetMean(), 0.5*h1->GetRMS(), 
// 		    h1->GetMaximum()*0.2, h1->GetMean(), 3.*h1->GetRMS());
//   h1->Fit(f0, "0");

  f0->DrawCopy("same");

  tl->SetTextColor(kBlack);  
  tl->SetNDC(kTRUE); tl->SetTextSize(0.06);

  double mean  = f0->GetParameter(1);
  double sigma = f0->GetParameter(2);
  
  tl->DrawLatex(0.16, 0.85, Form("#mu: %5.3f#pm%5.3f GeV", mean, 0.001));
  tl->DrawLatex(0.16, 0.79, Form("#sigma: %5.3f#pm%5.3f GeV", sigma, 0.001));
//   tl->DrawLatex(0.16, 0.85, Form("#mu_{1} = %5.3f#pm%5.3f GeV", f2->GetParameter(1), f2->GetParError(1)));
//   tl->DrawLatex(0.16, 0.80, Form("#mu_{2} = %5.3f#pm%5.3f GeV", f2->GetParameter(4), f2->GetParError(4)));
//   tl->DrawLatex(0.16, 0.75, Form("#sigma_{1} = %5.3f#pm%5.3f GeV", f2->GetParameter(2), f2->GetParError(2)));
//   tl->DrawLatex(0.16, 0.70, Form("#sigma_{2} = %5.3f#pm%5.3f GeV", f2->GetParameter(5), f2->GetParError(5)));
//   tl->DrawLatex(0.16, 0.65, Form("g_{2}/g_{1} = %5.3f#pm%5.3f ", f2->GetParameter(3), f2->GetParError(3)));

  tl->SetTextSize(0.06); tl->SetTextColor(kBlue);

  if ( o == 3 ) {
    tl->DrawLatex(0.19, 0.74, "perfect align.");
    tl->DrawLatex(0.2, 0.67, "J/#Psi #rightarrow #mu #mu");
  } else if ( o == 4 ) {
    tl->DrawLatex(0.19, 0.74, "short-term align.");
    tl->DrawLatex(0.2, 0.67, "J/#Psi #rightarrow #mu #mu");
  } else {  
    tl->SetTextSize(0.1);
    tl->DrawLatex(0.70, 0.82, Form("MC%i",o));
  }

  c0->SaveAs(Form("%s/s%d-jreco.eps", outDir, o));

  delete h1;
}


// ----------------------------------------------------------------------
void anaBmm::mcValidation(int wiat) {

  mcVal("c000"); if (wiat) if (wait()) goto end;
  mcVal("c001"); if (wiat) if (wait()) goto end;
  mcVal("c010"); if (wiat) if (wait()) goto end;
  mcVal("c011"); if (wiat) if (wait()) goto end;
  mcVal("c020"); if (wiat) if (wait()) goto end;
  mcVal("c021"); if (wiat) if (wait()) goto end;
  mcVal("c024"); if (wiat) if (wait()) goto end;
  mcVal("c025"); if (wiat) if (wait()) goto end;
  mcVal("c030"); if (wiat) if (wait()) goto end;
  mcVal("c023"); if (wiat) if (wait()) goto end;
  mcVal("c026"); if (wiat) if (wait()) goto end;
  mcVal("c027"); if (wiat) if (wait()) goto end;
  mcVal("c035"); if (wiat) if (wait()) goto end;
  mcVal("c022"); if (wiat) if (wait()) goto end;

  mcVal2("c000"); if (wiat) if (wait()) goto end;
  mcVal2("c001"); if (wiat) if (wait()) goto end;
  mcVal2("c010"); if (wiat) if (wait()) goto end;
  mcVal2("c011"); if (wiat) if (wait()) goto end;
  mcVal2("c020"); if (wiat) if (wait()) goto end;
  mcVal2("c021"); if (wiat) if (wait()) goto end;
  mcVal2("c024"); if (wiat) if (wait()) goto end;
  mcVal2("c025"); if (wiat) if (wait()) goto end;
  mcVal2("c030"); if (wiat) if (wait()) goto end;
  mcVal2("c023"); if (wiat) if (wait()) goto end;
  mcVal2("c026"); if (wiat) if (wait()) goto end;
  mcVal2("c027"); if (wiat) if (wait()) goto end;
  mcVal2("c035"); if (wiat) if (wait()) goto end;
  mcVal2("c022"); if (wiat) if (wait()) goto end;
 
  mcVal3("c000"); if (wiat) if (wait()) goto end;
  mcVal3("c001"); if (wiat) if (wait()) goto end;
  mcVal3("c010"); if (wiat) if (wait()) goto end;
  mcVal3("c011"); if (wiat) if (wait()) goto end;
  mcVal3("c020"); if (wiat) if (wait()) goto end;
  mcVal3("c021"); if (wiat) if (wait()) goto end;
  mcVal3("c024"); if (wiat) if (wait()) goto end;
  mcVal3("c025"); if (wiat) if (wait()) goto end;
  mcVal3("c030"); if (wiat) if (wait()) goto end;
  mcVal3("c023"); if (wiat) if (wait()) goto end;
  mcVal3("c026"); if (wiat) if (wait()) goto end;
  mcVal3("c027"); if (wiat) if (wait()) goto end;
  mcVal3("c035"); if (wiat) if (wait()) goto end;
  mcVal3("c022"); if (wiat) if (wait()) goto end;
 end:;

}

// ----------------------------------------------------------------------
void anaBmm::mcVal(const char *hname) {

  gStyle->SetOptTitle(0);

  TH1D *h0, *h1; 
  fS[0]->cd();
  h0 = (TH1D*)(fS[0]->Get(hname))->Clone();
  if (0 == h0) {
    cout << "Histogram with name " << hname << " does not exist" << endl;
    return;
  }

  fS[1]->cd();
  h1 = (TH1D*)(fS[1]->Get(hname))->Clone();
  if (0 == h1) {
    cout << "Histogram with name " << hname << " does not exist" << endl;
    return;
  }

  h0->Scale(1./h0->GetSumOfWeights());
  h1->Scale(1./h1->GetSumOfWeights());
  
  c0->Clear();
  shrinkPad(0.15, 0.2);
  setHist(h0, kBlack, 20, 2);
  setTitles(h0, h0->GetXaxis()->GetTitle(), h0->GetYaxis()->GetTitle(), 0.06, 1.1, 1.8);
  gPad->SetLogy(1);
  h0->DrawCopy();

  setHist(h1, kBlue, 24, 2);
  h1->DrawCopy("same");

  if (!strcmp(hname, "c000")) {
    legg = new TLegend(0.5,0.7,0.8,0.89);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05);  legg->SetFillColor(0); 
    legge = legg->AddEntry(h1, Form("sm06_bsmumu"), "p"); legge->SetTextColor(kBlue);
    legge = legg->AddEntry(h0, Form("private prod."), "p"); legge->SetTextColor(kBlack);
    legg->Draw();
  }

  if (!strcmp(hname, "c021")) {
    legg = new TLegend(0.5,0.7,0.8,0.89);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05);  legg->SetFillColor(0); 
    legge = legg->AddEntry(h1, Form("sm06_bsmumu"), "p"); legge->SetTextColor(kBlue);
    legge = legg->AddEntry(h0, Form("private prod."), "p"); legge->SetTextColor(kBlack);
    legg->Draw();
  }

  c0->Draw(); 
  c0->Update();

  c0->SaveAs(Form("%s/mcval-%s.eps", outDir, hname));

  gPad->SetLogy(0);

  delete h0;
  delete h1;

}
// ----------------------------------------------------------------------
void anaBmm::mcVal2(const char *hname) {

  gStyle->SetOptTitle(0);

  TH1D *h0, *h1; 
  fS[0]->cd();
  h0 = (TH1D*)(fS[0]->Get(hname))->Clone();
  if (0 == h0) {
    cout << "Histogram with name " << hname << " does not exist" << endl;
    return;
  }

  fS[2]->cd();
  h1 = (TH1D*)(fS[2]->Get(hname))->Clone();
  if (0 == h1) {
    cout << "Histogram with name " << hname << " does not exist" << endl;
    return;
  }

  h0->Scale(1./h0->GetSumOfWeights());
  h1->Scale(1./h1->GetSumOfWeights());
  
  c0->Clear();
  shrinkPad(0.15, 0.2);
  setHist(h1, kBlue, 24, 2);
  setTitles(h1, h0->GetXaxis()->GetTitle(), h0->GetYaxis()->GetTitle(), 0.06, 1.1, 1.8);
  gPad->SetLogy(1);
  h1->DrawCopy();

  setHist(h0, kBlack, 20, 2);
  h0->DrawCopy("same");

  if (!strcmp(hname, "c000")) {
    legg = new TLegend(0.5,0.7,0.8,0.89);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05);  legg->SetFillColor(0); 
    legge = legg->AddEntry(h1, Form("CMSSW"), "p"); legge->SetTextColor(kBlue);
    legge = legg->AddEntry(h0, Form("ORCA/OSCAR"), "p"); legge->SetTextColor(kBlack);
    legg->Draw();
  }

  if (!strcmp(hname, "c021")) {
    legg = new TLegend(0.5,0.7,0.8,0.89);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05);  legg->SetFillColor(0); 
    legge = legg->AddEntry(h1, Form("CMSSW"), "p"); legge->SetTextColor(kBlue);
    legge = legg->AddEntry(h0, Form("ORCA/OSCAR"), "p"); legge->SetTextColor(kBlack);
    legg->Draw();
  }

  c0->Draw(); 
  c0->Update();

  c0->SaveAs(Form("%s/mcval2-%s.eps", outDir, hname));

  gPad->SetLogy(0);

  delete h0;
  delete h1;

}

// ----------------------------------------------------------------------
void anaBmm::mcVal3(const char *hname) {

  gStyle->SetOptTitle(0);

  TH1D *h0, *h1; 
  fS[3]->cd();
  h0 = (TH1D*)(fS[3]->Get(hname))->Clone();
  if (0 == h0) {
    cout << "Histogram with name " << hname << " does not exist" << endl;
    return;
  }

  fS[4]->cd();
  h1 = (TH1D*)(fS[4]->Get(hname))->Clone();
  if (0 == h1) {
    cout << "Histogram with name " << hname << " does not exist" << endl;
    return;
  }

  h0->Scale(1./h0->GetSumOfWeights());
  h1->Scale(1./h1->GetSumOfWeights());
  
  c0->Clear();
  shrinkPad(0.15, 0.2);
  setHist(h1, kBlue, 24, 2);
  setTitles(h1, h0->GetXaxis()->GetTitle(), h0->GetYaxis()->GetTitle(), 0.06, 1.1, 1.8);
  gPad->SetLogy(1);
  h1->DrawCopy();

  setHist(h0, kBlack, 20, 2);
  h0->DrawCopy("same");

  if (!strcmp(hname, "c000")) {
    legg = new TLegend(0.4,0.7,0.8,0.89);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05);  legg->SetFillColor(0); 
    legge = legg->AddEntry(h1, Form("short-term misal."), "p"); legge->SetTextColor(kBlue);
    legge = legg->AddEntry(h0, Form("perfect alignment"), "p"); legge->SetTextColor(kBlack);
    legg->Draw();
  }

  if (!strcmp(hname, "c021")) {
    legg = new TLegend(0.4,0.7,0.8,0.89);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05);  legg->SetFillColor(0); 
    legge = legg->AddEntry(h1, Form("short-term misal."), "p"); legge->SetTextColor(kBlue);
    legge = legg->AddEntry(h0, Form("perfect alignment"), "p"); legge->SetTextColor(kBlack);
    legg->Draw();
  }

  c0->Draw(); 
  c0->Update();

  c0->SaveAs(Form("%s/mcval3-%s.eps", outDir, hname));

  gPad->SetLogy(0);

  delete h0;
  delete h1;

}

// ----------------------------------------------------------------------
void anaBmm::showDistributions(int offset, int wiat) { 

  int logy(100);   

  showDistribution(Form("c%d00", offset), 2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d01", offset), 2);            if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d10", offset), 2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d11", offset), 2);            if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d12", offset), 2, 0.6, 0.75); if (wiat) if (wait()) goto end;

  showDistribution(Form("c%d20", offset), 2,      0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d21", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d22", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d23", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d26", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d24", offset), 2,      0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d25", offset), 2,      0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d27", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;

  showDistribution(Form("c%d40", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d41", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d42", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d43", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;


 end:;

}


// ----------------------------------------------------------------------
void anaBmm::showDistribution(const char *hname, int mode, double x, double y) { 

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  c0->Clear();
  shrinkPad(0.15, 0.2);
  
  //  TH1 *hm = sumHistMC(hname);
  TH1 *hm = (TH1D*)fM[0]->Get(hname);

  //  setFilledHist(hm, kBlack, kBlack, 3000, 2);
  setHist(hm, kBlack, kBlack);
  setTitles(hm, hm->GetXaxis()->GetTitle(), hm->GetYaxis()->GetTitle(), 0.06, 1.1, 1.5);

  TH1D *hs = new TH1D(*((TH1D*)fS[0]->Get(hname)));
  hs->SetName(Form("SIG:%s", hname)); 
  setFilledHist(hs, kBlue, kBlue, 3004, 2);
  setTitles(hs, hs->GetXaxis()->GetTitle(), hs->GetYaxis()->GetTitle(), 0.06, 1.1, 1.5);

  int logy(0); 
  if (mode > 99) {
    mode -= 100;
    logy = 1;
  }  
  
  c0->SetLogy(logy);
  
  if (mode == 0) {
    cout << hname << " drawn as is: sg, bg" << endl;
    hs->Draw("hist");
    hm->Draw("samehist");
  } else if (mode == 1) {
    cout << hname << " drawn as is: bg, sg" << endl;
    hm->Draw("hist");
    hs->Draw("samehist");
  } else if (mode == 2) {
    cout << hname << " scale to unity" << endl;
    hm->Scale(1./hm->GetSumOfWeights()); 
    hs->Scale(1./hs->GetSumOfWeights()); 
    hs->SetMaximum(1.1*(hm->GetMaximum() > hs->GetMaximum()? hm->GetMaximum(): hs->GetMaximum()));
    hs->Draw("hist");
    hm->Draw("samehist");
  } else if (mode == 3) {
    cout << hname << " scale sg to bg" << endl;
    hs->Scale(hm->GetSumOfWeights()/hs->GetSumOfWeights()); 
    hm->Draw("hist");
    hs->Draw("samehist");
  } else if (mode == 4) {
    cout << hname << " scale to L = " << fLumiD[0] << "/fb" << endl;
    hs->Scale(fLumiD[0]/fLumiS[0]); 
    hm->Scale(fLumiD[0]/fLumiM[0]); 
    hm->Draw("hist");
    hs->Draw("samehist");
  }

  if (x > 0) {
    legg = new TLegend(x, y, x+0.15, y+0.1);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
    //     legge = legg->AddEntry(hs, Form("S: %5.1f", hs->GetSumOfWeights()), "f"); legge->SetTextColor(kBlack);
    //     legge = legg->AddEntry(hm, Form("B: %5.1f", hm->GetSumOfWeights()), "f"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hs, Form("Signal"), "f"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hm, Form("Background"), "f"); legge->SetTextColor(kBlack);
    legg->Draw();
  }
 
  c0->SaveAs(Form("%s/dist-%s.eps", outDir, hname));

  c0->SetLogy(0);

}
  


// ----------------------------------------------------------------------
void anaBmm::showProcesses(int signal) {

  TString label;
  if (1 == signal) {
    fS[0]->cd();
    label = TString("sproc");
  } else {
    fM[0]->cd();
    label = TString("bproc");
  }


  const char *cuts = "pt>5 && ptl1>4";

  plotProcesses(cuts, "pt", "p_{T, B} [GeV]", 10, 0., 50., 1);
  c0->SaveAs(Form("%s/%s-pt.eps", outDir, label.Data()));

  plotProcesses(cuts, "iso", "I", 20, 0., 1., 0);
  c0->SaveAs(Form("%s/%s-isolation.eps", outDir, label.Data()));

  plotProcesses(cuts, "rmm", "#Delta R(#mu,#mu)", 10, 0., 2., 0);
  c0->SaveAs(Form("%s/%s-rmm.eps", outDir, label.Data()));

  plotProcesses(cuts, "mass", "m_{#mu#mu} [GeV]", 20, 5., 6., 0);
  c0->SaveAs(Form("%s/%s-mass.eps", outDir, label.Data()));

}

// ----------------------------------------------------------------------
void anaBmm::plotProcesses(const char *cut, const char *var, const char *axis, int bin, double lo, double hi, int legend) {

  TTree *tS = (TTree*)gFile->Get("events");
  TH1D *hGSP = new TH1D("hGSP", "", bin, lo, hi); //hGSP->Sumw2();
  TH1D *hFEX = new TH1D("hFEX", "", bin, lo, hi); //hFEX->Sumw2();
  TH1D *hGGF = new TH1D("hGGF", "", bin, lo, hi); //hGGF->Sumw2();
  TH1D *hSum = new TH1D("hSum", "", bin, lo, hi); //hSum->Sumw2();

  tS->Draw(Form("%s >> hGSP", var), Form("%s && process == 42", cut)); 
  tS->Draw(Form("%s >> hFEX", var), Form("%s && process == 41", cut)); 
  tS->Draw(Form("%s >> hGGF", var), Form("%s && process == 40", cut)); 

  hSum->Add(hGSP);
  hSum->Add(hFEX);
  hSum->Add(hGGF);

  //   hGSP->Divide(hGSP, hSum, 1., 1.);
  //   hFEX->Divide(hFEX, hSum, 1., 1.);
  //   hGGF->Divide(hGGF, hSum, 1., 1.);
  
  //   hGSP->SetMinimum(0.);
  //   hGSP->SetMaximum(1.);

  double max = hGSP->GetMaximum();
  if (hFEX->GetMaximum() > max) {
    max = hFEX->GetMaximum();
  } 
  if (hGGF->GetMaximum() > max) {
    max = hGGF->GetMaximum();
  } 
  
  gStyle->SetOptStat(0);
  shrinkPad(0.15, 0.15);
  setTitles(hGSP, axis, "arbitrary units", 0.05, 1.1, 1.5);

  hGSP->SetMaximum(1.1*max);
  hGSP->SetMinimum(0.0);
  setHist(hGSP, kBlue, 20, 2); 
  setHist(hFEX, kRed, 20, 2);    hFEX->SetLineStyle(kDashed);
  setHist(hGGF, kBlack, 20, 2);  hGGF->SetLineStyle(kDotted);

  hGSP->Draw("hist");
  hFEX->Draw("histsame");
  hGGF->Draw("histsame");

  if (1 == legend) {
    legg = new TLegend(0.5,0.7,0.8,0.89);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
    legge = legg->AddEntry(hGSP, Form("gluon splitting"), "l"); legge->SetTextColor(kBlue);
    legge = legg->AddEntry(hFEX, Form("flavor excitation"), "l"); legge->SetTextColor(kRed);
    legge = legg->AddEntry(hGGF, Form("gluon-gluon fusion"), "l"); legge->SetTextColor(kBlack);
    legg->Draw();
  }


}



// ----------------------------------------------------------------------
TH1* anaBmm::sumHistMC(const char *hname, int mode, const char *selection) {
 
  // mode: 0 weighted by xsec
  //       1 weigthed to correspond to Lumi(data)
  //       2 weigthed by xsec incl. Muon-MisId Rate

  cout << endl << endl;

  TH1D *h1 = (TH1D*)gDirectory->Get(Form("MC:%s", hname));
  if (h1) {
    cout << "anaBmm::sumHistMC> Deleting " << Form("MC:%s", hname) << endl;
    delete h1; 
  }

  TH1D *h2 = (TH1D*)fM[0]->Get(hname); 
  if (h2) {
    h1 = new TH1D(*h2);
    h1->SetName(Form("MC:%s", hname)); 
    h1->Reset(); 
  }
    
  const int nhist = nMc;

  for (int i = 0; i < nhist; ++i) {

    if ( !strcmp(selection,"r0") ) {  // skip generic background to
				      // get combined rare background

      if ( i == 0 ) {

	continue;
      }
    }

    else if (  strcmp(selection,"c0") ) {
      
      if ( strcmp((getRareType(fSignM[i])).Data(), selection) ) {
	
	continue;
      }
    }

    if (fM[i]) {

      h2 = (TH1D*)fM[i]->Get(hname);  
   
      if (mode == 0) {
	if (TMath::IsNaN(h2->GetSumOfWeights())) {
	  h2->Reset();
	  cout << "anaBmm::sumHistMC> ***** Problems with histogram " << hname << " from file " <<
	    fM[i]->GetName() << endl;
	}
	else {
	  h1->Add(h2, fLumiD[0]/fLumiM[i]);
	}
      }

      if (mode == 2) {

	if (TMath::IsNaN(h2->GetSumOfWeights())) {
	  h2->Reset();
	  cout << "anaBmm::sumHistMC> ***** Problems with histogram " << hname << " from file " <<
	    fM[i]->GetName() << endl;
	}
	else {
	  sumHistMC_Add(i, hname, h1, h2);
	  //h1->Add(h2,  (fMisIdM[i]*fLumiD[0]/fLumiM[i]));
	}
     }     
    }
  }

  cout << endl << endl;
  return h1; 
}
// -----------------------------------------------------
void anaBmm::channelEff(TFile *f, double &fnorm_Ch, double &eff1_Ch, double &eff2_Ch) {
  
  TH1D *tAR1 = (TH1D*)f->Get("AR1")->Clone();
  fnorm_Ch = 0.; eff1_Ch = 0.; eff2_Ch = 0.;

  fnorm_Ch = tAR1->GetBinContent(tAR1->FindBin(123.1));
  
  if ( tAR1->GetBinContent(tAR1->FindBin(224.1)) || tAR1->GetBinContent(tAR1->FindBin(324.1)) ) {
    
    eff1_Ch  = tAR1->GetBinContent(tAR1->FindBin(224.1))/
      (tAR1->GetBinContent(tAR1->FindBin(224.1))+tAR1->GetBinContent(tAR1->FindBin(324.1)));
  }
  
  if ( tAR1->GetBinContent(tAR1->FindBin(226.1)) || tAR1->GetBinContent(tAR1->FindBin(326.1)) ) {
    
    eff2_Ch  = tAR1->GetBinContent(tAR1->FindBin(226.1))/
      (tAR1->GetBinContent(tAR1->FindBin(226.1))+tAR1->GetBinContent(tAR1->FindBin(326.1)));
  }
}
// -----------------------------------------------------
void anaBmm::sumHistMC_Add(int ch, const char *hist, TH1 *hist1, TH1 *hist2) {
  

  double SF_Ch(0.);
  double fnorm_Ch(0.), eff1_Ch(0.), eff2_Ch(0.);

  SF_Ch = fMisIdM[ch]*fLumiD[0]/fLumiM[ch];
  
  hist1->Add(hist2, SF_Ch );

  if ( !strcmp(hist,"AR1") ) {

    channelEff(fM[ch], fnorm_Ch, eff1_Ch, eff2_Ch);
    
    for (int m = 0; m < hist2->GetNbinsX(); m++ ) { 

      if ( m == hist2->FindBin(424.1) ) {

	hist1->AddBinContent(m, SF_Ch*fnorm_Ch*eff1_Ch);

      }
      else if ( m == hist2->FindBin(426.1) ) {

	hist1->AddBinContent(m, SF_Ch*fnorm_Ch*eff2_Ch);
      }
      else if ( m == hist2->FindBin(428.1) ) {

	hist1->AddBinContent(m, SF_Ch*fnorm_Ch*eff1_Ch*eff2_Ch);
      } 
    }
  }
}
// -----------------------------------------------------
void anaBmm::muonMisId() {

  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.06);
  tl->SetTextColor(kBlack);

  const int mhist = nMc - 1;
  const int shist = nSg - 1;
  //  const int nhist = 8;
  double mis(0.), tot(0.);
  
  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT <<  "% -- " << "Muon misidentification"<< endl;
  
  // --- Mis-ID. histograms for Pions, Kaons, Protons & Muons
 
  // -- from signal

  printf("anaBmm::muonMisId> Getting MisID histogram from %s\n", fSignS[0].Data());
 
  TH1D *h = (TH1D*)fS[0]->Get("MisID");

  TH1D *PI_pT = (TH1D*)fS[0]->Get("m100");    
  TH1D *PI_pT_Mu = (TH1D*)fS[0]->Get("m101");      
  TH1D *PI_Eta = (TH1D*)fS[0]->Get("m110");
  TH1D *PI_Eta_Mu = (TH1D*)fS[0]->Get("m111");
   
  TH1D *K_pT = (TH1D*)fS[0]->Get("m200");    
  TH1D *K_pT_Mu = (TH1D*)fS[0]->Get("m201");      
  TH1D *K_Eta = (TH1D*)fS[0]->Get("m210");
  TH1D *K_Eta_Mu = (TH1D*)fS[0]->Get("m211");
  
  TH1D *P_pT = (TH1D*)fS[0]->Get("m300");    
  TH1D *P_pT_Mu = (TH1D*)fS[0]->Get("m301");      
  TH1D *P_Eta = (TH1D*)fS[0]->Get("m310");
  TH1D *P_Eta_Mu = (TH1D*)fS[0]->Get("m311");
  
  TH1D *MU_pT = (TH1D*)fS[0]->Get("m400");    
  TH1D *MU_pT_Mu = (TH1D*)fS[0]->Get("m401");      
  TH1D *MU_Eta = (TH1D*)fS[0]->Get("m410");
  TH1D *MU_Eta_Mu = (TH1D*)fS[0]->Get("m411");

  TH2D *PI = (TH2D*)fS[0]->Get("M100");
  TH2D *PImis = (TH2D*)fS[0]->Get("M101"); 
  TH2D *K = (TH2D*)fS[0]->Get("M200");
  TH2D *Kmis = (TH2D*)fS[0]->Get("M201");
  TH2D *P = (TH2D*)fS[0]->Get("M300");
  TH2D *Pmis = (TH2D*)fS[0]->Get("M301");
  TH2D *MU = (TH2D*)fS[0]->Get("M400");
  TH2D *MUeff = (TH2D*)fS[0]->Get("M401");  

  // -- from background
  for (int i = 2; i < shist; i++) {

    printf("anaBmm::muonMisId> Getting MisID histogram from %s\n", fSignS[i].Data());    

    h->Add((TH1D*)fS[i]->Get("MisID"));

    PI_pT->Add((TH1D*)fS[i]->Get("m100"));    
    PI_pT_Mu->Add((TH1D*)fS[i]->Get("m101"));      
    PI_Eta->Add((TH1D*)fS[i]->Get("m110"));
    PI_Eta_Mu->Add((TH1D*)fS[i]->Get("m111"));
    
    K_pT->Add((TH1D*)fS[i]->Get("m200"));
    K_pT_Mu->Add((TH1D*)fS[i]->Get("m201"));
    K_Eta->Add((TH1D*)fS[i]->Get("m210"));
    K_Eta_Mu->Add((TH1D*)fS[i]->Get("m211"));
    
    P_pT->Add((TH1D*)fS[i]->Get("m300"));    
    P_pT_Mu->Add((TH1D*)fS[i]->Get("m301"));      
    P_Eta->Add((TH1D*)fS[i]->Get("m310"));
    P_Eta_Mu->Add((TH1D*)fS[i]->Get("m311"));
    
    MU_pT->Add((TH1D*)fS[i]->Get("m400"));    
    MU_pT_Mu->Add((TH1D*)fS[i]->Get("m401"));      
    MU_Eta->Add((TH1D*)fS[i]->Get("m410"));
    MU_Eta_Mu->Add((TH1D*)fS[i]->Get("m411"));

    PI->Add((TH2D*)fS[i]->Get("M100"));
    PImis->Add((TH2D*)fS[i]->Get("M101"));
    K->Add((TH2D*)fS[i]->Get("M200"));
    Kmis->Add((TH2D*)fS[i]->Get("M201"));
    P->Add((TH2D*)fS[i]->Get("M300"));
    Pmis->Add((TH2D*)fS[i]->Get("M301"));
    MU->Add((TH2D*)fS[i]->Get("M400"));
    MUeff->Add((TH2D*)fS[i]->Get("M401"));      
  }

  for (int i = 0; i < mhist; i++) {

    printf("anaBmm::muonMisId> Getting MisID histogram from %s\n", fSignM[i].Data());

    h->Add((TH1D*)fM[i]->Get("MisID"));

    PI_pT->Add((TH1D*)fM[i]->Get("m100"));    
    PI_pT_Mu->Add((TH1D*)fM[i]->Get("m101"));      
    PI_Eta->Add((TH1D*)fM[i]->Get("m110"));
    PI_Eta_Mu->Add((TH1D*)fM[i]->Get("m111"));
    
    K_pT->Add((TH1D*)fM[i]->Get("m200"));
    K_pT_Mu->Add((TH1D*)fM[i]->Get("m201"));
    K_Eta->Add((TH1D*)fM[i]->Get("m210"));
    K_Eta_Mu->Add((TH1D*)fM[i]->Get("m211"));
    
    P_pT->Add((TH1D*)fM[i]->Get("m300"));    
    P_pT_Mu->Add((TH1D*)fM[i]->Get("m301"));      
    P_Eta->Add((TH1D*)fM[i]->Get("m310"));
    P_Eta_Mu->Add((TH1D*)fM[i]->Get("m311"));
    
    MU_pT->Add((TH1D*)fM[i]->Get("m400"));    
    MU_pT_Mu->Add((TH1D*)fM[i]->Get("m401"));      
    MU_Eta->Add((TH1D*)fM[i]->Get("m410"));
    MU_Eta_Mu->Add((TH1D*)fM[i]->Get("m411"));

    PI->Add((TH2D*)fM[i]->Get("M100"));
    PImis->Add((TH2D*)fM[i]->Get("M101"));
    K->Add((TH2D*)fM[i]->Get("M200"));
    Kmis->Add((TH2D*)fM[i]->Get("M201"));
    P->Add((TH2D*)fM[i]->Get("M300"));
    Pmis->Add((TH2D*)fM[i]->Get("M301"));
    MU->Add((TH2D*)fM[i]->Get("M400"));
    MUeff->Add((TH2D*)fM[i]->Get("M401"));
      
  }

   // --- RATE OF ...

  // ... mis-identified pions
  tot = 0.;  mis = 0.; double PiMID(0.);
  tot = h->GetBinContent(h->FindBin(0.1));
  mis = h->GetBinContent(h->FindBin(1.1));
  if ( tot != 0 ) { PiMID = mis/tot; }

  // ... mis-identified kaons
  tot = 0.;  mis = 0.; double KaMID(0.);
  tot = h->GetBinContent(h->FindBin(10.1));
  mis = h->GetBinContent(h->FindBin(11.1));
  if ( tot != 0 ) { KaMID = mis/tot; }

  // ... mis-identified protons
  tot = 0.;  mis = 0.; double ProtMID(0.);
  tot = h->GetBinContent(h->FindBin(20.1));
  mis = h->GetBinContent(h->FindBin(21.1));
  if ( tot != 0 ) { ProtMID = mis/tot; }

  // ... not identified _Muons
  tot = 0.;  mis = 0.; double MuEff(0.);
  tot = h->GetBinContent(h->FindBin(30.1));
  mis = h->GetBinContent(h->FindBin(32.1));
  if ( tot != 0 ) { MuEff = mis/tot; }


  // --- Draw combined Mis-ID. histograms

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetLogy(0); 
  shrinkPad(0.15, 0.15);

  // -- 1D histograms
  // -- Pions
  PI_pT->Draw("hist");
  PI_pT_Mu->Draw("psame");
  tl->DrawLatex(0.16, 0.92, "Misidentified Pions");
  c0->SaveAs(Form("%s/MisID_PionsPt.eps", outDir));

  // - pT
  fPiMid = DivideHisto(PI_pT_Mu, PI_pT);
  fPiMid->GetYaxis()->SetRangeUser(0,0.022);
  fPiMid->GetYaxis()->SetNdivisions(6, kTRUE);
  setTitles(fPiMid, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  fPiMid->Draw();
  tl->DrawLatex(0.16, 0.92,  "Misidentified pions");
  //  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (p_{T}): Pions");
  c0->SaveAs(Form("%s/MisID_rate_PionsPt.eps", outDir));

  // - eta
  TH1D* MisID = DivideHisto(PI_Eta_Mu, PI_Eta);
  MisID->GetYaxis()->SetRangeUser(0,0.02);
  MisID->Draw();
  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (#eta): Pions");
  c0->SaveAs(Form("%s/MisID_rate_PionsEta.eps", outDir));


  // -- Kaons
  K_pT->GetYaxis()->SetRangeUser(0,0.015);
  K_pT->Draw("hist");
  K_pT_Mu->Draw("psame");
  tl->DrawLatex(0.16, 0.92, "Misidentified Kaons");
  c0->SaveAs(Form("%s/MisID_KaonsPt.eps", outDir));

  // - pT
  fKaMid = DivideHisto(K_pT_Mu, K_pT);
  fKaMid->GetYaxis()->SetRangeUser(0,0.022); 
  fKaMid->GetYaxis()->SetNdivisions(6, kTRUE);
  setTitles(fKaMid, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  fKaMid->Draw();
  tl->DrawLatex(0.16, 0.92,  "Misidentified kaons");
  //  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (p_{T}): Kaons");
  c0->SaveAs(Form("%s/MisID_rate_KaonsPt.eps", outDir));

  // - eta
  MisID = DivideHisto(K_Eta_Mu, K_Eta);
  MisID->GetYaxis()->SetRangeUser(0,0.02);
  MisID->Draw();
  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (#eta): Kaons");
  c0->SaveAs(Form("%s/MisID_rate_KaonsEta.eps", outDir));


  // -- Protons
  P_pT->GetYaxis()->SetRangeUser(0,0.02);
  P_pT->Draw("hist");
  P_pT_Mu->Draw("psame");
  tl->DrawLatex(0.16, 0.92, "Misidentified Protons");
  c0->SaveAs(Form("%s/MisID_ProtonsPt.eps", outDir));

  // - pT
  fProtMid = DivideHisto(P_pT_Mu, P_pT);
  fProtMid->GetYaxis()->SetRangeUser(0,0.022);
  fProtMid->GetYaxis()->SetNdivisions(6, kTRUE);
  setTitles(fProtMid, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  fProtMid->Draw();
  tl->DrawLatex(0.16, 0.92,  "Misidentified protons");
  //  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (p_{T}): Protons");
  c0->SaveAs(Form("%s/MisID_rate_ProtonsPt.eps", outDir));

  // - eta
  MisID = DivideHisto(P_Eta_Mu, P_Eta);
  MisID->GetYaxis()->SetRangeUser(0,0.022);
  setTitles(MisID, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  MisID->Draw();
  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (#eta): Protons");
  c0->SaveAs(Form("%s/MisID_rate_ProtonsEta.eps", outDir));


  // -- Muons
  MU_pT->Draw("hist");
  MU_pT_Mu->Draw("psame");
  tl->DrawLatex(0.16, 0.92, "Identified Muons");
  c0->SaveAs(Form("%s/MisID_MuonEffPt.eps", outDir));

  // - pT
  fMuEff = DivideHisto(MU_pT_Mu, MU_pT);
  fMuEff->GetYaxis()->SetRangeUser(0,1.1);
  setTitles(fMuEff, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  //  fMuEff->GetYaxis()->SetNdivisions(6, kTRUE);
  fMuEff->Draw();
  tl->DrawLatex(0.16, 0.92, "Muon efficiency");
  c0->SaveAs(Form("%s/MisID_rate_MuonEffPt.eps", outDir));
  // - eta
  MisID = DivideHisto(MU_Eta_Mu, MU_Eta);
  MisID->GetYaxis()->SetRangeUser(0,1.1);
  MisID->Draw();
  tl->DrawLatex(0.16, 0.92, "Muon efficiency  (#eta)");
  c0->SaveAs(Form("%s/MisID_rate_MuonEffEta.eps", outDir));


  // -- 2D histograms

  fPiMid2 = DivideHisto2(PImis, PI);
  setTitles2(fPiMid2, "p_{T}", "#eta", 0.06, 1.1, 1.1);
  fPiMid2->Draw("colz");
  tl->DrawLatex(0.16, 0.92, "Misidentified Pions");
  //  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (#eta vs. p_{T}): Pions");
  c0->SaveAs(Form("%s/MisID_rate_Pions.eps", outDir));

  fKaMid2 = DivideHisto2(Kmis, K);
  setTitles2(fKaMid2, "p_{T}", "#eta", 0.06, 1.1, 1.1);
  fKaMid2->Draw("colz");
  tl->DrawLatex(0.16, 0.92, "Misidentified Kaons");
  //  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (#eta vs. p_{T}): Kaons");
  c0->SaveAs(Form("%s/MisID_rate_Kaons.eps", outDir));

  fProtMid2 = DivideHisto2(Pmis, P);
  setTitles2(fProtMid2, "p_{T}", "#eta", 0.06, 1.1, 1.1);
  fProtMid2->Draw("colz");
  tl->DrawLatex(0.16, 0.92, "Misidentified protons");
  //  tl->DrawLatex(0.16, 0.92, "Mis-Id. spectrum (#eta vs. p_{T}): Protons");
  c0->SaveAs(Form("%s/MisID_rate_Protons.eps", outDir));

  fMuEff2 = DivideHisto2(MUeff, MU);
  setTitles2(fMuEff2, "p_{T}", "#eta", 0.06, 1.1, 1.1); 
  fMuEff2->GetZaxis()->SetRangeUser(0.,1.);
  fMuEff2->GetZaxis()->SetNdivisions(6, kTRUE);
  fMuEff2->GetYaxis()->SetNdivisions(6, kTRUE);
  fMuEff2->Draw("colz");
  tl->DrawLatex(0.16, 0.92, Form("Muon efficiency"));
  c0->SaveAs(Form("%s/MisID_rate_MuonEff.eps", outDir));
  
  OUT << Form("\\vdef{effIdMuonsSim}    {\\ensuremath{ {%4.3E } } }", MuEff)   << endl;
  OUT << Form("\\vdef{misIdPionsSim}    {\\ensuremath{ {%4.3E } } }", PiMID)   << endl;
  OUT << Form("\\vdef{misIdKaonsSim}    {\\ensuremath{ {%4.3E } } }", KaMID)   << endl;
  OUT << Form("\\vdef{misIdProtonsSim}  {\\ensuremath{ {%4.3E } } }", ProtMID) << endl;
 
  OUT << Form("\\vdef{effIdMuons}    {\\ensuremath{ {%4.3E } } }", fMu)   << endl;
  OUT << Form("\\vdef{misIdPions}    {\\ensuremath{ {%4.3E } } }", fPi)   << endl;
  OUT << Form("\\vdef{misIdKaons}    {\\ensuremath{ {%4.3E } } }", fKa)   << endl;
  OUT << Form("\\vdef{misIdProtons}  {\\ensuremath{ {%4.3E } } }", fProt) << endl;



  OUT.close();
}


// ================================================================================
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
void anaBmm::shrinkPad(double b, double l, double r, double t) {
  gPad->SetBottomMargin(b); 
  gPad->SetLeftMargin(l);
  gPad->SetRightMargin(r);
  gPad->SetTopMargin(t);
}


// ----------------------------------------------------------------------
void anaBmm::setTitles(TH1 *h, const char *sx, const char *sy, float size, 
               float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy); 
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}

// ----------------------------------------------------------------------
void anaBmm::setTitles2(TH2 *h, const char *sx, const char *sy, float size, 
               float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy); 
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}


// ----------------------------------------------------------------------
void anaBmm::setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(int(width));
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
}


// ----------------------------------------------------------------------
void anaBmm::setFilledHist(TH1 *h, Int_t color, Int_t fillcolor, Int_t fillstyle, Int_t width) {
  // Note: 3004, 3005 are crosshatches
  // ----- 1000       is solid
  //       kYellow    comes out gray on bw printers
  h->SetLineColor(color);     h->SetLineWidth(width);   
  h->SetFillStyle(fillstyle); h->SetFillColor(fillcolor);
}




// ----------------------------------------------------------------------
double anaBmm::nul(int nobs) {
  if (nobs <= 0) return  2.3;
  if (nobs == 1) return  3.89; 
  if (nobs == 2) return  5.32;
  if (nobs == 3) return  6.68;
  if (nobs == 4) return  7.99;
  if (nobs == 5) return  9.27;
  if (nobs == 6) return 10.53;
  if (nobs == 7) return 11.77;
  if (nobs == 8) return 12.99;
  if (nobs == 9) return 14.21;
  if (nobs ==10) return 15.41;
  return -1.;
}


// ----------------------------------------------------------------------
double anaBmm::barlow(int nobs, double bg, double bgE, double sE) {

  double ul(-99.); 

  if (nobs > 10)  {
    cout << "barlow: dunno" << endl;
    return ul; 
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
      cout << "break " << tFraction << " bestAlpha: " << bestAlpha << " (should be " << cl << ")" << endl;
      break;
    }
  }

  return ul; 
}
// ----------------------------------------------------------------------
TString anaBmm::texForm(double e) {
  TString seff1(Form("%4.2e", e)); 
  seff1.ReplaceAll("e-0", "\\times 10^{-"); 
  seff1.ReplaceAll("e+0", "\\times 10^{"); 
  seff1 += TString("}");
  return seff1;
}
// ----------------------------------------------------------------------
TString anaBmm::texForm2(double e) {
  TString seff1(Form("%4.2e", e)); 
  seff1.ReplaceAll("e-", "\\times 10^{-"); 
  seff1.ReplaceAll("e+", "\\times 10^{"); 
  seff1 += TString("}");
  return seff1;
}

// ----------------------------------------------------------------------
TString anaBmm::texForm31(double e) {
  TString seff1(Form("%3.1e", e)); 
  seff1.ReplaceAll("e-0", "\\times 10^{-"); 
  seff1.ReplaceAll("e+0", "\\times 10^{"); 
  seff1 += TString("}");
  return seff1;
}

// ----------------------------------------------------------------------
TH1D* anaBmm::DivideHisto(TH1D *hn, TH1D *hN) {
 
  TH1D *hn_new = (TH1D*)hn->Rebin(4,"hn_new"); 
  TH1D *hN_new = (TH1D*)hN->Rebin(4,"hN_new"); 

  TH1D *h = new TH1D(*hn_new);
  h->Reset(); 
  
  double rate = 0.;
  double error = 0.;
  double tmp1 = 0.;
  double tmp2 = 0.;

  for (int i = 0; i < hn_new->GetNbinsX(); ++i) {
    
    tmp1 = 0.;
    tmp2 = 0.;

    tmp1 = hn_new->GetBinContent(i);
    tmp2 = hN_new->GetBinContent(i);

    
    if ( tmp2 ) { rate = tmp1/tmp2; error = dBinomial(tmp1, tmp2); }
    else { rate = 0; error = 0; }
    
    h->SetBinContent(i, rate);
    h->SetBinError(i, error);
    
  }

  return h;

}
// ----------------------------------------------------------------------
TH2D* anaBmm::DivideHisto2(TH2D *hn, TH2D *hN) {
 
  TH2D *hn_new = (TH2D*)hn->Rebin2D(4, 5,"hn_new"); 
  TH2D *hN_new = (TH2D*)hN->Rebin2D(4, 5,"hN_new"); 
  
  TH2D *h = new TH2D(*hn_new);
  h->Reset(); 
  
  double rate = 0.;
  double error = 0.;
  double tmp1 = 0.;
  double tmp2 = 0.;

  for (int i = 0; i < hn_new->GetNbinsX(); ++i) {
    for (int j = 0; j < hn_new->GetNbinsY(); ++j) {
    
      tmp1 = 0.;
      tmp2 = 0.;
      
      tmp1 = hn_new->GetBinContent(i, j);
      tmp2 = hN_new->GetBinContent(i, j);
         
      if ( tmp2 ) { rate = tmp1/tmp2; error = dBinomial(tmp1, tmp2); }
      else { rate = 0; error = 0; }
      
      h->SetBinContent(i, j, rate);
      h->SetBinError(i, j, error);
    }
  }

  return h;

}

// ----------------------------------------------------------------------
double anaBmm::dBinomial(double n, double N) {

  double w = n/N;
  return TMath::Sqrt(TMath::Abs(w*(1-w)/N));
}

// ----------------------------------------------------------------------
double anaBmm::dEff(int in, int iN) {

  double n = (double)in;
  double N = (double)iN;
  return TMath::Sqrt(((n+1)*(N-n+1))/((N+3)*(N+2)*(N+2)));
}

// ----------------------------------------------------------------------
TString anaBmm::scaleFactor(int exp) {

  TString factor;

  if ( exp == 0 ) {
    
    factor = TString(""); 
  }
  else if ( exp == 1 ) {
    
    factor = TString("(x 10)"); 
  }
  else if ( exp > -100 ) { 
    
    factor = TString(Form("(10^{%i})", exp));
  }
  else {

    factor = TString("");
  }

  return factor;

}

// ----------------------------------------------------------------------
void anaBmm::getSignature(TString signIn, TString &signOut) {

 signOut = signIn;

 TString titles[]    = {   TString("B_{s} \\rightarrow \\mu^{-} \\mu^{+}")
			 , TString("B^{+} \\rightarrow J \\Psi K^{+}")
			 , TString("b\\bar{b}")

			 , TString("B_{d} \\rightarrow \\pi^{-} \\pi^{+}")
			 , TString("B_{d} \\rightarrow \\pi^{-} K^{+}")
			 , TString("B_{d} \\rightarrow \\pi^{-} \\mu^{+} \\nu")

			 , TString("B_{s} \\rightarrow \\pi^{-} \\pi^{+}")
			 , TString("B_{s} \\rightarrow K^{-} K^{+}")
			 , TString("B_{s} \\rightarrow K^{-} \\pi^{+}")
			 , TString("B_{s} \\rightarrow K^{-} \\mu^{+} \\nu")
			 , TString("B_{s} \\rightarrow \\mu^{+} \\mu^{-} \\gamma")
			 , TString("B_{s} \\rightarrow \\mu^{+} \\mu^{-} \\pi_{0}")

			 , TString("B_{u} \\rightarrow \\mu^{+} \\mu^{-} \\mu^{+} \\nu")

			 , TString("B_{c} \\rightarrow \\mu^{+} \\mu^{-} \\mu^{+} \\nu")
			 , TString("B_{c} \\rightarrow J/\\Psi \\mu^{+} \\nu")

			 , TString("\\Lambda_{b} \\rightarrow p K^{-}")
			 , TString("\\Lambda_{b} \\rightarrow p \\pi^{-}")

			 , TString("M(hh)")
                         };
 
 TString sign[]       = {  TString("bsmumu")
			 , TString("bpjpsikp")
			 , TString("bbbar")

			 , TString("bd2pi")
			 , TString("bdpik")
			 , TString("bdpimunu")

			 , TString("bs2pi")
			 , TString("bskk")
			 , TString("bskpi")
			 , TString("bskmunu")
			 , TString("bsmumug")
			 , TString("bsmumup0")

			 , TString("bu3munu")

			 , TString("bc3munu")
			 , TString("bcjpmunu")

			 , TString("lbpk")
			 , TString("lbppi")

			 , TString("qcd")

                        };


 for (int i = 0; i < 17; i++) {

   if ( !strcmp(signIn.Data(), sign[i].Data()) ) {
   
     signOut = TString(titles[i]);
     
   }
 }

}

// ----------------------------------------------------------------------
TString anaBmm::getRareType(TString signIn) {


  TString rtype;

  TString sign[]      = {  TString("bsmumu")
			 , TString("bpjpsikp")
			 , TString("bbbar")

			 , TString("bd2pi")
			 , TString("bdpik")
			 , TString("bdpimunu")

			 , TString("bs2pi")
			 , TString("bskk")
			 , TString("bskpi")
			 , TString("bskmunu")
			 , TString("bsmumug")
			 , TString("bsmumup0")

			 , TString("bu3munu")

			 , TString("bc3munu")
			 , TString("bcjpmunu")

			 , TString("lbpk")
			 , TString("lbppi")

			 , TString("qcd")

                        };

  TString type[]      = {   TString("2mu")
			  , TString("2mu")
		   
			  , TString("2mId")
			  , TString("2mId")
			  , TString("mIdMu+")
		   
			  , TString("2mId")
			  , TString("2mId")
			  , TString("2mId")
			  , TString("mIdMu+")
			  , TString("2mu+")
			  , TString("2mu+")
	   
			  , TString("2mu+")

			  , TString("2mu+")
			  , TString("2mu+")

			  , TString("2mId")
			  , TString("2mId")

			  , TString("qcd")
                          };
 

 for (int i = 0; i < 18; i++) {

   if ( !strcmp(signIn.Data(), sign[i].Data()) ) {
  
     rtype = type[i];
  
   }
 }

 return rtype;
}

// ----------------------------------------------------------------------
double anaBmm::getMisID(TString signIn) {

  double id(-1.);

  TString sign[]      = {  TString("bsmumu")
			 , TString("bpjpsikp")
			 , TString("bbbar")

			 , TString("bd2pi")
			 , TString("bdpik")
			 , TString("bdpimunu")

			 , TString("bs2pi")
			 , TString("bskk")
			 , TString("bskpi")
			 , TString("bskmunu")
			 , TString("bsmumug")
			 , TString("bsmumup0")

			 , TString("bu3munu")

			 , TString("bc3munu")
			 , TString("bcjpmunu")

			 , TString("lbpk")
			 , TString("lbppi")

			 , TString("qcd")

                        };


  double eff[]        = {   fMu*fMu
			  , fMu*fMu
			  , fMu*fMu
		   
			  , fPi*fPi
			  , fPi*fKa
			  , fPi*fMu
	   
			  , fPi*fPi
			  , fKa*fKa
			  , fPi*fKa
			  , fKa*fMu
			  , fMu*fMu
			  , fMu*fMu

			  , fMu*fMu

			  , fMu*fMu
			  , fMu*fMu

			  , fProt*fKa
			  , fProt*fPi
	   
			  , fPi*fKa
                          };

 
 for (int i = 0; i < 17; i++) {

   if ( !strcmp(signIn.Data(), sign[i].Data()) ) {
    
     id = eff[i];
  
   }
 }
 
 return id;

}

