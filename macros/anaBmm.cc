#include "anaBmm.hh"

#include "TF1.h"
#include "THStack.h"
#include "TTree.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
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

//===========================================================================================
// -- Fit functions
//===========================================================================================

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
double f_p1a2gauss(double *x, double *par) {
  // par[0] -> area
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> area
  // par[4] -> mean
  // par[5] -> sigma
  // par[6] = par 0 of pol1
  // par[7] = par 1 of pol1

    return ( par[6] + par[7]*x[0] + f_2gauss(x, &par[0]) );
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
// pol0 and Gauss
double f_p0ag(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] = par 0 of pol0

  return  (par[3] + f_Gauss(x, &par[0]));
}

// ----------------------------------------------------------------------
// pol1 and Gauss
double f_p1ag(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] = par 0 of pol1
  // par[4] = par 1 of pol1

  return  (par[3] + par[4]*x[0] + f_Gauss(x, &par[0]));
}

// ----------------------------------------------------------------------
// exp and Gauss
double f_eag(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] = par 0 of exp
  // par[4] = par 1 of exp

  return (par[3]*TMath::Exp(-x[0]*par[4]) + f_Gauss(x, &par[0]));
}


// ----------------------------------------------------------------------
double f_expo(double *x, double *par) {
  return par[0]*TMath::Exp(-x[0]*par[1]);
}

//===========================================================================================

anaBmm::anaBmm(const char *files) { 
  init(files);
}


// ----------------------------------------------------------------------
void anaBmm::init(const char *files) {
  fFont = 132; 

  fMassBs = 5.369;
  fMassBp = 5.279;

  fSgReduction = -1.;
  fBgReduction = -1.;

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
  f0 = new TF1("f0", f_gauss,  4.8,  6.0, 3);
  f1 = new TF1("f1", f_2gauss, 4.8,  6.0, 6);
  f2 = new TF1("f2", f_2G,     5.0,  6.0, 6);
  f3 = new TF1("f3", f_2g,     4.8,  6.0, 6);
  f4 = new TF1("f4", f_p1ag,   4.8,  6.0, 5);
  f5 = new TF1("f5", f_eag,    4.8,  6.0, 5);

  f6 = new TF1("f6", f_p1a2gauss,  4.8,  6.0, 8);

  f10= new TF1("f10", f_expo,  4.8,  6.0, 2);
  f11= new TF1("f11", f_Gauss, 4.8,  6.0, 3);

  
  // -- setup all files
  nSg = nMc = nDa = 0;
  for (int i = 0; i < 10; ++i) {
    fNevtS[i] = 0.;
    fLumiS[i] = 0.;
    fS[i] = 0; 
    fD[i] = 0; 
  }

  for (int i = 0; i < 30; ++i) {
    fNevtM[i] = 0.;
    fLumiM[i] = 0.;
    fM[i] = 0; 
  }

  fLumiD[0] = 10.;
  cout << "================================" << endl;
  cout << "--> Setting Lumi to " << fLumiD[0] << " /fb <--" << endl;
  cout << "================================" << endl<< endl;

  cout << "--> Loading rootfiles" << endl;
  loadFiles(files);
}


//===========================================================================================
// -- load files
//===========================================================================================

void anaBmm::loadFiles(const char *filename) {

  char buffer[200];
  char type[100];
  char file[100];
  char name[100];
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
    cn.ReplaceAll("default2D", "");
    cn.ReplaceAll("default3D", "");
    cn.ReplaceAll("default2a", "");
    cn.ReplaceAll("default2b", "");
    cn.ReplaceAll("default3", "");
    cn.ReplaceAll("default4", "");
    cn.ReplaceAll("default", "");
    cn.ReplaceAll("norm", "");
    cn.ReplaceAll("treebmm", "");
    cn.ReplaceAll("BAK", "");
    cn.ReplaceAll("root", "");
    cn.ReplaceAll(".", "");
    cn.ReplaceAll("/", "");

    sprintf(name, "%s", cn.Data());

    if (!strcmp(type, "mysg")) {

      sgIndex = nSg;
      sprintf(type, "sg");
      sprintf(line, "fS[%d] = ", nSg);
      cout << "Loading signal file   " << line << file << " and vis x-section: " << visXsection 
	   << " (" << signature << "). \t\t <<<<<< my signal sample" <<  endl;
      loadSg(file, visXsection, signature, type, name);
      
    } else if (!strcmp(type, "mynsg")) {
      
      normSgIndex = nSg;
      sprintf(type, "nsg");
      sprintf(line, "fS[%d] = ", nSg);
      cout << "Loading signal file   " << line << file << " and vis x-section: " << visXsection 
	   << " (" << signature << "). \t <<<<<< my norm. sample" <<  endl;
      loadSg(file, visXsection, signature, type, name);
    
    } else if (!strcmp(type, "mymc")) {
      
      bgIndex = nMc;
      sprintf(type, "mc");
      sprintf(line, "fM[%d] = ", nMc);
      cout << "Loading MC file   " << line << file << " and vis x-section: " << visXsection 
	   << " (" << signature << "). \t\t <<<<<< my background sample" <<  endl;
      loadMc(file, visXsection, signature, type, name);
      
    } else if (!strcmp(type, "mynmc")) {
      
      normBgIndex = nMc;
      sprintf(type, "nmc");
      sprintf(line, "fM[%d] = ", nMc);
      cout << "Loading MC file   " << line << file << " and vis x-section: " << visXsection 
	   << " (" << signature << "). \t\t <<<<<< my norm. background sample" <<  endl;
      loadMc(file, visXsection, signature, type, name);
    
    } else if (!strcmp(type, "sg") || !strcmp(type, "nsg") ) {
      
      sprintf(line, "fS[%d] = ", nSg);
      cout << "Loading signal file   " << line << file << " and vis x-section: " << visXsection 
	   << " (" << signature << ")." <<  endl;
      loadSg(file, visXsection, signature, type, name);
    
    } else if (!strcmp(type, "mc") || !strcmp(type, "nmc") ) {

      sprintf(line, "fM[%d] = ", nMc);
      cout << "Loading MC file   " << line << file << " and vis x-section: " << visXsection 
	   << " (" << signature << ")." <<  endl;
      loadMc(file, visXsection, signature, type, name);


    } else if (!strcmp(type, "rmc")) {

      sprintf(signature, "%s", cn.Data());
      sprintf(line, "fM[%d] = ", nMc);
      cout << "Loading MC file   " << line << file << " and vis x-section: " << visXsection 
	   << " (" << signature << ")." <<  endl;
      loadMc(file, visXsection, signature, type, name);

    } else if (!strcmp(type, "da")) {

      sprintf(signature, "Data");
      sprintf(line, "fD[%d] = ", nDa);
      cout << "Loading data file  " << line << file << " and vis x-section: " << visXsection  
	   << " (" << signature << ")." <<  endl;
      loadDa(file, visXsection, signature, type, name);
    
    }  else {
      //
    }
  }


  fMu = 1.;
  fPi = 0.3E-2;
  fKa = 0.7E-2;
  fProt = 0.1E-2;

  TString fn(filename);
  fn.ReplaceAll("bmm", "");
  fn.ReplaceAll("files", "");
  fn.ReplaceAll(".", "");

  fNumbersFileName = TString(Form("%s/anaBmm.%s.tex", outDir, fn.Data()));
  sprintf(line, "rm -f %s", fNumbersFileName.Data());
  system(line);

  dumpFiles();

  dumpCuts();
}



// ----------------------------------------------------------------------
void anaBmm::loadSg(const char *name, double lumi, const char *sign, const char *type, const char *filename) {

  if (nSg > 10) {
    cout << " **** !!!! Too many open Signal files. Increase nSg. !!!! **** " << endl;
    return;
  }

  fFileS[nSg] = TString(filename);

  fS[nSg] = new TFile(name);
  fLumiS[nSg] = lumi;

  TH1 *h      = (TH1D*)fS[nSg]->Get("AR1");
  fvXsS[nSg]  = lumi;
  fNevtS[nSg] = h->GetBinContent(h->FindBin(0.1));
  if ( fFileS[nSg].Contains("sg-001") || fFileS[nSg].Contains("sg-002") || fFileS[nSg].Contains("sg-004m") ) { 
    fNevtS[nSg] = h->GetBinContent(h->FindBin(1.1)); 
  } 
  fLumiS[nSg] = fNevtS[nSg]/fvXsS[nSg] ;
  fSignS[nSg] = TString(sign);
  fTypeS[nSg] = TString(type);

  getSignature(fSignS[nSg], fSignTexS[nSg], fSignLeggS[nSg]);

  ++nSg; 
}

// ----------------------------------------------------------------------
void anaBmm::loadMc(const char *name, double lumi, const char *sign, const char *type, const char *filename) {

  if (nMc > 30) {
    cout << " **** !!!! Too many open MC files. Increase nMc. !!!! **** " << endl;
    return;
  } 

  fFileM[nMc] = TString(filename);

  fM[nMc] = new TFile(name);
  fLumiM[nMc] = lumi;

  TH1 *h      = (TH1D*)fM[nMc]->Get("AR1");
  fvXsM[nMc]  = lumi;
  fNevtM[nMc] = h->GetBinContent(h->FindBin(0.1));
  if ( fFileM[nMc].Contains("bg-001") || fFileM[nMc].Contains("bg-002") ) { 
    fNevtM[nMc] = h->GetBinContent(h->FindBin(1.1)); 
  } 
  fLumiM[nMc] = fNevtM[nMc]/fvXsM[nMc] ;
  fSignM[nMc] = TString(sign);
  fTypeM[nMc] = TString(type);

  getSignature(fSignM[nMc], fSignTexM[nMc], fSignLeggM[nMc]);
 
  ++nMc; 
}

// ----------------------------------------------------------------------
void anaBmm::loadDa(const char *name, double lumi, const char *sign, const char *type, const char *filename) {

  if (nDa > 10) {
    cout << " **** !!!! Too many open DATA files. Increase nDa.  !!!! **** " << endl;
    return;
  } 

  fD[nDa] = new TFile(name);
  fLumiD[nDa] = lumi;

  TH1 *h      = (TH1D*)fD[nDa]->Get("AR1");
  fvXsD[nDa]  = lumi;
  fNevtD[nDa] = h->GetBinContent(h->FindBin(0.1));
  //  fNevtD[nDa] = h->GetBinContent(h->FindBin(1.1));
  fLumiD[nDa] = fNevtD[nDa]/fvXsD[nDa] ;
  fSignD[nDa] = TString(sign);
  fTypeD[nDa] = TString(type);
  fFileD[nDa] = TString(filename);

  getSignature(fSignD[nDa], fSignTexD[nDa], fSignLeggD[nDa]);

  ++nDa; 
}


//===========================================================================================
// -- produce plots & number for note
//===========================================================================================

void anaBmm::makeAllPlots() {
   

  validations();
  distributions();

  processes();
  fakeMuons();

  brecos();  

  effTables();

  bgOverlays();

  double ul = calculateUpperLimit();
  normalizedUpperLimit();
  
}

// ----------------------------------------------------------------------

void anaBmm::validations() {

  mcValidation(0);
  mcValidation(3);

}

// -- validation

void anaBmm::mcValidation(int offset, int wiat) {

  int logy(100);   

  mcVal(Form("c%d00", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d01", offset), logy+2, 0.4, 0.15); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d10", offset), 2,      0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d11", offset), 2,      0.4, 0.15); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d12", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;

  mcVal(Form("c%d20", offset), 2,      0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d21", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d22", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d23", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d26", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d24", offset), 2,      0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d25", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d27", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;

  mcVal(Form("c%d30", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d35", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d40", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d41", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d42", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;

  mcVal(Form("c%d14", offset), logy+10, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d15", offset), logy+10, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d16", offset), logy+10, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d37", offset), logy+10, 0.6, 0.75); if (wiat) if (wait()) goto end;
  mcVal(Form("c%d62", offset), logy+10, 0.3, 0.75); if (wiat) if (wait()) goto end;


 end:;

}

// ----------------------------------------------------------------------

void anaBmm::distributions() {

  showDistributions(0);
  showDistributions(1);  

  showDistributions(2);
  showDistributions(3);

  showDistributions(4);
  showDistributions(5);
}

// -- distribution
void anaBmm::showDistributions(int offset, int wiat) { 

  int logy(100);   

  showDistribution(Form("c%d00", offset), 2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d01", offset), 2);            if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d10", offset), 2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d11", offset), 2);            if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d12", offset), 2, 0.6, 0.75); if (wiat) if (wait()) goto end;

  showDistribution(Form("c%d14", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d15", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;

  showDistribution(Form("c%d20", offset), 2,      0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d21", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d22", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d23", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d26", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d16", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d24", offset), 2,      0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d25", offset), 2,      0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d27", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d37", offset), logy+2, 0.6, 0.75); if (wiat) if (wait()) goto end;

  showDistribution(Form("c%d30", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d35", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d40", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d41", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d42", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;
  showDistribution(Form("c%d62", offset), logy+2, 0.3, 0.75); if (wiat) if (wait()) goto end;


 end:;

}

// ----------------------------------------------------------------------
void anaBmm::processes() {

  showProcesses(1);
  showProcesses(0);
}


// ----------------------------------------------------------------------

void anaBmm::brecos() {

  breco(-1);

  char f_n[200];
  int n(0);

  // -- J/Psi-peaks from other samples
  int njf = 2;
  TString jfiles[]    = { TString("csg-004") , TString("csg-004m") };

  // -- B-peaks from signal samples
  int nbf = 4;
  TString bfiles[]    = {  TString("csg-003") , TString("csg-003h") 
			 , TString("csg-004h") , TString("csg-004m") };

  // -- B-peaks from other samples
  int nnf = 2;
  TString nfiles[]    = { TString("csg-004"), TString("csg-004j") };


  for (int i = 0; i < nbf; i++ ) {

    sprintf(f_n,"%s", bfiles[i].Data());
    n = findIndex(f_n);
    
    if ( n < 0 ) {
      cout << " ====> Error: couldn't find index of file " << f_n << endl;
      continue;
    }

    breco(n, "c030");
    breco(n, "c230");
    breco(n, "c330");
    breco(n, "c430");
    breco(n, "c530");
  }

  for (int i = 0; i < nnf; i++ ) {

    sprintf(f_n,"%s", nfiles[i].Data());
    n = findIndex(f_n);
    
    if ( n < 0 ) {
      cout << " ====> Error: couldn't find index of file " << f_n << endl;
      continue;
    }

    nreco(n, "c030");
    nreco(n, "c230");
    nreco(n, "c330");
    nreco(n, "c430");
    nreco(n, "c530"); 

  }

  for (int i = 0; i < njf; i++ ) {

    sprintf(f_n,"%s", jfiles[i].Data());
    n = findIndex(f_n); 
    
    if ( n < 0 ) {
      cout << " ====> Error: couldn't find index of file " << f_n << endl;
      continue;
    }

    jreco(n);
  }
}

// ----------------------------------------------------------------------

void anaBmm::effTables() {

  effTable(fS[sgIndex], "mysg");
  effTable(fM[bgIndex], "mymc");

  effTable(fS[normSgIndex], "mynsg");
  effTable(fM[normBgIndex], "mynmc");
  
  effTable(fM[bgIndex], "2mId");
  effTable(fM[bgIndex], "mIdMu+");
  effTable(fM[bgIndex], "2mu+");
  
  //  effTable(fM[bgIndex], "qcd");
  
  effTable(fM[bgIndex], "r0"); // has to be before c0 !!!
  effTable(fM[bgIndex], "c0");
}

// ----------------------------------------------------------------------

void anaBmm::bgOverlays() {

  bgOverlay("c033", 6); // overlay of different rare background with sg + bg
  bgOverlay("c030", 6); // overlay of different rare background with sg + bg
  bgOverlay("c430", 4); // overlay of different rare background with sg + bg
  bgOverlay("c530", 5); // overlay of different rare background with sg + bg

  nbgOverlay("c033", 5); // overlay of sg + bg (norm)
  nbgOverlay("c030", 5); // overlay of sg + bg (norm)
  nbgOverlay("c430", 5); // overlay of sg + bg (norm)
  nbgOverlay("c530", 5); // overlay of sg + bg (norm)
}

// ----------------------------------------------------------------------

void anaBmm::cuts() {
  
 
  TString histos[] = { TString("PTLO"),  TString("RMM"),  TString("MASSBAND"),  TString("PTBS"),  TString("ETABS"),  TString("COSALPHA"),  TString("LXYSXY"),  TString("L3DS3D"),  TString("ISOLATION"),  TString("VTXCHI"),   TString("MASSWI") } ;
  int logy[]  = {0,0,1,0,0,1,1,1,0,1,0};
  double lo[] = {ptmulo, rmmlo, masslo, ptbs, etalo, coslo, lxylo, l3dlo, isolo, vtxhi, 5.369-masswi };
  double hi[] = {25.,    rmmhi, masshi,  30., etahi,    1.,   50.,   50.,   1.1,    0., 5.369+masswi };

  TString histosF[] = { TString("LXYSXY_PRE"),   TString("L3DS3D_PRE"),   TString("ISOLATION_PRE"),  TString("VTXCHI_PRE"), TString("MASSWI_F"),  TString("ISOLATION_F"),   TString("VTXCHI_F") };
  int logyF[] = {1,1,0,1,0,0,1};
  double loF[] = { 7.,  7., isolo, vtxhi, 5.369-masswi, isolo, vtxhi };
  double hiF[] = {50., 50.,   1.1,    0., 5.369+masswi,   1.1,    0.,};

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  c0->Clear();

  TCanvas *c1 = new TCanvas("c1", "cuts", 1000, 800);
  c1->Clear();
  c1->Divide(4,3);

  TH1D *hs, *hm;
  double x(0.6), y(0.75);
  pl->SetLineColor(2);
  pl->SetLineWidth(2);
  for (int i = 0; i < 11; i++ ) {
    
    hs = (TH1D*)fS[sgIndex]->Get(Form("%s", histos[i].Data()));
    hm = (TH1D*)fM[bgIndex]->Get(Form("%s", histos[i].Data()));

    //    c0->cd((i+1)%6);
    
    setHist(hm, kBlack, kBlack);
    setTitles(hm, hm->GetXaxis()->GetTitle(), hm->GetYaxis()->GetTitle(), 0.06, 1.1, 1.5);
    
    setFilledHist(hs, kBlue, kBlue, 3004, 2);
    setTitles(hs, hs->GetXaxis()->GetTitle(), hs->GetYaxis()->GetTitle(), 0.06, 1.1, 1.5);
    
    hm->Scale(1./hm->GetSumOfWeights()); 
    hs->Scale(1./hs->GetSumOfWeights()); 
    hs->SetMaximum(1.1*(hm->GetMaximum() > hs->GetMaximum()? hm->GetMaximum(): hs->GetMaximum()));

    c0->cd();
    gPad->SetLogy(logy[i]);
    shrinkPad(0.15, 0.2);
    hs->DrawCopy("hist");
    hm->DrawCopy("samehist");

    legg = new TLegend(x, y, x+0.15, y+0.1);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
    legge = legg->AddEntry(hs, Form("Signal"), "f"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hm, Form("Background"), "f"); legge->SetTextColor(kBlack);
    legg->Draw();

    pl->DrawLine(lo[i], 0, lo[i], hs->GetMaximum());
    pl->DrawLine(hi[i], 0, hi[i], hs->GetMaximum());


    c1->cd(i+1);
    gPad->SetLogy(logy[i]);
    shrinkPad(0.15, 0.2);
    hs->DrawCopy("hist");
    hm->DrawCopy("samehist");

    legg = new TLegend(x, y, x+0.15, y+0.1);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
    legge = legg->AddEntry(hs, Form("Signal"), "f"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hm, Form("Background"), "f"); legge->SetTextColor(kBlack);
    legg->Draw();

    pl->DrawLine(lo[i], 0, lo[i], hs->GetMaximum());
    pl->DrawLine(hi[i], 0, hi[i], hs->GetMaximum());

    c0->SaveAs(Form("%s/dist/dist-before-c%i%s.eps", outDir, i, histos[i].Data()));
    
  }

  c1->SaveAs(Form("%s/dist/dist-before-cuts.eps", outDir));
  //  delete c1;

  TCanvas *c2 = new TCanvas("c2", "fact", 1000, 500);
  c2->Clear();
  c2->Divide(4,2);

  for (int i = 0; i < 7; i++ ) {
    
    hs = (TH1D*)fS[sgIndex]->Get(Form("%s", histosF[i].Data()));
    hm = (TH1D*)fM[bgIndex]->Get(Form("%s", histosF[i].Data()));
    
    //   c0->cd(i+1);    
    setHist(hm, kBlack, kBlack);
    setTitles(hm, hm->GetXaxis()->GetTitle(), hm->GetYaxis()->GetTitle(), 0.06, 1.1, 1.5);
    
    setFilledHist(hs, kBlue, kBlue, 3004, 2);
    setTitles(hs, hs->GetXaxis()->GetTitle(), hs->GetYaxis()->GetTitle(), 0.06, 1.1, 1.5);
    
    hm->Scale(1./hm->GetSumOfWeights()); 
    hs->Scale(1./hs->GetSumOfWeights()); 
    hs->SetMaximum(1.1*(hm->GetMaximum() > hs->GetMaximum()? hm->GetMaximum(): hs->GetMaximum()));

    c0->cd();
    gPad->SetLogy(logyF[i]);
    shrinkPad(0.15, 0.2);
    hs->DrawCopy("hist");
    hm->DrawCopy("samehist");

    legg = new TLegend(x, y, x+0.15, y+0.1);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
    legge = legg->AddEntry(hs, Form("Signal"), "f"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hm, Form("Background"), "f"); legge->SetTextColor(kBlack);
    legg->Draw();

    pl->DrawLine(loF[i], 0, loF[i], hs->GetMaximum());
    pl->DrawLine(hiF[i], 0, hiF[i], hs->GetMaximum());

    c0->SaveAs(Form("%s/dist/dist-before-f%i%s.eps", outDir, i, histosF[i].Data()));

    c2->cd(i+1);
    gPad->SetLogy(logyF[i]);
    shrinkPad(0.15, 0.2);
    hs->DrawCopy("hist");
    hm->DrawCopy("samehist");

    legg = new TLegend(x, y, x+0.15, y+0.1);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
    legge = legg->AddEntry(hs, Form("Signal"), "f"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hm, Form("Background"), "f"); legge->SetTextColor(kBlack);
    legg->Draw();

    pl->DrawLine(loF[i], 0, loF[i], hs->GetMaximum());
    pl->DrawLine(hiF[i], 0, hiF[i], hs->GetMaximum());
  }
  
  c2->SaveAs(Form("%s/dist/dist-before-fact.eps", outDir));
  
  //  delete c2;
  delete hs;
  delete hm;
}




//===========================================================================================
// -- Upper limit calculation
//===========================================================================================

double anaBmm::calculateUpperLimit() {

  double nBs  = fLumiD[0] * 500.*1.e9 * 0.107 * 2.;

  double ekin = 45. / (1.e5 * (500./55.e3) * 0.107 * 2.);
  double scaleUL = ekin * nBs * 0.015 ; // this is: eff_kin * eff_total * N_Bs 
  cout << " --> scale UL (CSA07) = " << scaleUL << endl;

  ekin = 0.00028 /( (500./55.e3) * 0.107 * 2.);
  scaleUL = ekin * nBs * 0.017 ;
  cout << " --> scale UL (Spring07) = " << scaleUL << endl;

  ekin = 5000. / (24.6e6 * (500./55.e3) * 0.107 * 2.); //same as ekin = 1000. / (4.8e6 * (500./55.e3) * 0.107 * 2.); 
  scaleUL = ekin * nBs * 0.016 ;
  cout << " --> scale UL (ORCA - vCE high stat) = " << scaleUL << endl;

  ekin = 1000. / (4.8e6 * (500./55.e3) * 0.107 * 2.); 
  scaleUL = ekin * nBs * 0.016 ;
  cout << " --> scale UL (ORCA - vUL low stat) = " << scaleUL << endl;

  // -- scale = eff x Lumi x vis. XS / BR(Bs->mu mu)
  scaleUL = fEsg * fLumiD[0] * fvXsS[sgIndex] / 3.35e-9;

  cout << " --> new scale UL = " << scaleUL << endl;
  
  // -- BF < N_UL(Nobs) / scaleUL  =  N_UL(Nobs) / (epsilon * N_Bs)
  double atlasUL = expUL(7, 20, 12);
  double nUL = expUL(fNsg, fNbg, fNbgE);

  double expectedUL = nUL/scaleUL;


  cout << "Nsg: " << fNsg << " +/- " << fNsgE
       << " Nbg: " << fNbg << " +/- " << fNbgE
       << endl;

  double Scp = scp(fNsg, fNbg, fNbgE, 0.);

  for (int i = 2; i < 50.; i += 2) {
    // verbose version
//     cout << "i = " << i << "  "; 
//     scp(i*fNsg, i*fNbg, i*fNbgE, 0.,1);
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

  return expectedUL;
}


// ----------------------------------------------------------------------
void anaBmm::normalizedUpperLimit() {

  double nBs  = fLumiD[0] * 500.*1.e9 * 0.107 * 2.;
  double ekin = 1000. / (4.8e6 * (500./55.e3) * 0.107 * 2.); // == 0.107 ?

  // double nBs_kin = 55.e3 * (1000. / 4.8e6) * fLumiD[0];

  // -- BF < N_UL(Nobs) / scaleUL  =  N_UL(Nobs) / (epsilon * N_Bs)
  double atlasUL = expUL(7, 20, 12);
  double scaleUL = ekin * nBs * fEsg ; // this is: eff_kin * eff_total * N_Bs 


  // -- BF <  N_UL(n_sg + n_bg) eff_Bplus f_u / (N_Bplus eff_B0 f_s ) * BR(Bplus)
  double nUL = expUL(fNsg, fNbg, fNbgE);

  double fu    = 0.398;
  double fs    = 0.104;
  double fu_fs = fu/fs;

  double accBu   = 1.;          double accBs = 1.;
  double trgBu   = 1.;          double trgBs = 1.;
  double anaBu   = fEsg_norm;   double anaBs = fEsg;

  double eBu_eBs = (accBu * trgBu * anaBu)/ (accBs * trgBs * anaBs);

  double BR_bjk = 5.98e-5;

  double expectedUL = (nUL/fNsg_norm) * fu_fs * eBu_eBs * BR_bjk;


  cout << "Nsg: " << fNsg << " +/- " << fNsgE
       << " Nbg: " << fNbg << " +/- " << fNbgE
       << " Nnm: " << fNsg_norm << " +/- " << fNsgE_norm
       << endl;

  // double Scp = scp(fNsg, fNbg, fNbgE, 0.);

  for (int i = 2; i < 50.; i += 2) {
    // verbose version
//     cout << "i = " << i << "  "; 
//     scp(i*fNsg, i*fNbg, i*fNbgE, 0.,1);
    scp(i*fNsg, i*fNbg, i*fNbgE, 0.);
  }


  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% Upper limit (with normalization)" << endl;
  OUT  << Form("\\vdef{NormalizedUpperLimit}  {\\ensuremath{{%s} } }", (texForm31(expectedUL)).Data()) << endl;
  OUT.close();
}
  

// ----------------------------------------------------------------------
double anaBmm::expUL(double s0, double b0, double bsig) {

  double beta  = 0.9;
  double e0    = 1.; 
  double esig  = 0.25; 
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


//===========================================================================================
// -- Signal reconstruction
//===========================================================================================

void anaBmm::breco(int o, const char *hist) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TH1D *h1;
  TH1D *hs;
  TH1D *hr;

  double mean(0.), sigma(0.);

  if ( o > -1 ) {

    // -- Signal (double gauss)
    // =========================_
    if ( !strcmp(fTypeS[o].Data(), "sg") || !strcmp(fTypeS[o].Data(), "mysg") ) {
    
      h1 = (TH1D*)(fS[o]->Get(hist))->Clone();
      h1->Scale(fLumiD[0]/fLumiS[o]);
    
      c0->Clear();
      shrinkPad(0.15, 0.15);
      setFilledHist(h1, kBlack, kYellow, 1000);    
    
      f1->SetParameters(h1->GetMaximum()*0.8, h1->GetMean(), 0.5*h1->GetRMS(), 
			h1->GetMaximum()*0.2, h1->GetMean(), 3.*h1->GetRMS());

      setTitles(h1, "m_{#mu#mu} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
      h1->SetMaximum(1.4*h1->GetMaximum());
    
      h1->DrawCopy("hist");
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
    
      writeFitPar(f1, o, mean, sigma, 6);
    }
  
    // --  norm. channel (double gauss + pol1)
    // =====================================
    if ( !strcmp(fTypeS[o].Data(), "nsg") || !strcmp(fTypeS[o].Data(), "mynsg") ) {
    
      h1 = (TH1D*)(fS[o]->Get(hist))->Clone();
      h1->Scale(fLumiD[0]/fLumiS[o]);
    
      c0->Clear();
      shrinkPad(0.15, 0.15);
      setFilledHist(h1, kBlack, kYellow, 1000);   
    
      f6->SetParameters(h1->GetMaximum()*0.8, 5.3, 0.01, 
			h1->GetMaximum()*0.2, 5.3, 0.03,
			0.8*h1->GetBinContent(1), 1.);
    
    
      setTitles(h1, "m_{#mu#mu K} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
      h1->SetMaximum(1.4*h1->GetMaximum());
    
      h1->DrawCopy("hist");
      h1->Fit(f6, "0");
      f6->DrawCopy("same");
    
      mean = TMath::Sqrt((f6->GetParameter(0)*f6->GetParameter(0)*f6->GetParameter(1)*f6->GetParameter(1)
			  + f6->GetParameter(3)*f6->GetParameter(3)*f6->GetParameter(4)*f6->GetParameter(4))
			 /
			 (f6->GetParameter(0)*f6->GetParameter(0) + f6->GetParameter(3)*f6->GetParameter(3))
			 );
    
      sigma = TMath::Sqrt((f6->GetParameter(0)*f6->GetParameter(0)*f6->GetParameter(2)*f6->GetParameter(2)
			   + f6->GetParameter(3)*f6->GetParameter(3)*f6->GetParameter(5)*f6->GetParameter(5))
			  /
			  (f6->GetParameter(0)*f6->GetParameter(0) + f6->GetParameter(3)*f6->GetParameter(3))
			  );
    
      writeFitPar(f6, o, mean, sigma, 3);
    }
    
    
    tl->SetTextColor(kBlack);  
    tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
  
    tl->DrawLatex(0.16, 0.85, Form("#mu: %5.3f#pm%5.3f GeV", mean, 0.001));
    tl->DrawLatex(0.16, 0.79, Form("#sigma: %5.3f#pm%5.3f GeV", sigma, 0.001));
  
    if (h1->GetSumOfWeights()<0.001) {
      tl->DrawLatex(0.16, 0.73, Form("N: %4.2e", h1->GetSumOfWeights()));
    } else if (h1->GetSumOfWeights() < 1.) {
      tl->DrawLatex(0.16, 0.73, Form("N: %5.3f", h1->GetSumOfWeights()));
    } else if (h1->GetSumOfWeights() < 1000.) {
      tl->DrawLatex(0.16, 0.73, Form("N: %5.1f", h1->GetSumOfWeights()));
    } else {
      tl->DrawLatex(0.16, 0.73, Form("N: %4.2e", h1->GetSumOfWeights()));
    }

    tl->SetTextSize(0.06); tl->SetTextColor(kBlue);
    tl->DrawLatex(0.56, 0.65, Form("%s", fSignLeggS[o].Data()));
  
//     if ( o == 0 ) {
//       tl->DrawLatex(0.48, 0.72, "ORCA/OSCAR (priv.)");
//     } else if ( o == 1 ) {
//       tl->DrawLatex(0.50, 0.72, "official MC");
//     } else if ( o == 2 ) {
//       tl->DrawLatex(0.52, 0.72, "CSA07");
//     } else if ( o == 3 ) {
//       tl->DrawLatex(0.52, 0.72, "Spring07");
//     } else if ( o == 4 ) {
//       tl->DrawLatex(0.56, 0.72, "perfect align.");
//     } else if ( o == 5 ) {
//       tl->DrawLatex(0.52, 0.72, "short-term align.");
//     } else {  
//       tl->DrawLatex(0.52, 0.72, Form("MC%i",o));
//     }
    
  
    c0->SaveAs(Form("%s/breco/s%d-breco-%s.eps", outDir, o, hist));
    delete h1;
  }

  // -- Rare backgrounds
  // ===================
  if ( o < 0 ) {
    
    int EColor[]     = {   2,    4,    108,  6,  107,   93,   1, 
			   2,    4,    108,  6,  107,   93,   1,  
			   2,    4,    108,  6,  107,   93,   1   }; 
  
    legg = new TLegend(0.60,0.8,0.89,0.89);
    legg->SetFillStyle(1001); legg->SetBorderSize(0); legg->SetTextSize(0.035);  legg->SetFillColor(0);
  
    char hist[200]; 
    int fit(0), log(0);

    const int nhist = nMc;
    
    for (int i = 0; i < nhist - 1; i++) {  // FIX ME: last histogram is screwd up -> takes signal instead ???????
 
      if ( !strcmp( getSubGroup(fSignM[i]).Data(), "2mId") ) { 
	sprintf(hist, "c030"); 
	fit = 1; log = 0;

      } else if ( !strcmp( getSubGroup(fSignM[i]).Data(), "mIdMu+") ) { 
	sprintf(hist, "c030"); 
	fit = 0; log = 0;

      } else if ( !strcmp(fTypeM[i].Data(), "rmc") ){ 
	sprintf(hist, "c033"); 
	fit = 0; log = 1;

      } else {

	continue;
      }

      hs = (TH1D*)(fS[sgIndex]->Get(hist))->Clone();
      if (fit) hs->GetXaxis()->SetRangeUser(4.9, 5.9);
      setFilledHist(hs, kBlack, kYellow, 3004, 2, 2);
      setTitles(hs, "m_{#mu#mu} [GeV]", "events/bin", 0.06, 1.1, 1.3);
      
      hr = (TH1D*)(fM[i]->Get(hist))->Clone();
      setHist(hr,  EColor[i]);

      c0->Clear();
      shrinkPad(0.2, 0.2);
      gPad->SetLogy(1);

      hr->Scale(fMisIdM[i]*fLumiD[0]/fLumiM[i]);
   
      if ( hr->GetSumOfWeights() > 0  && hr->GetMaximum() > 0 ) {

	hs->Scale(hr->GetSumOfWeights()/hs->GetSumOfWeights());

	gPad->SetLogy(log);

	if (log) {
	  hs->GetYaxis()->SetRangeUser(0.1*hs->GetMinimum(1.e-20), 100*hs->GetMaximum());
    	} else {
	  hs->GetYaxis()->SetRangeUser(0, 1.2*hs->GetMaximum());
	}

	hs->DrawCopy("hist");
	hr->DrawCopy("same");

	if (fit) {
	  

	  f6->SetParameters(0.8*hr->GetMaximum(), 5.3, 0.03, 
			    0.8*hr->GetBinContent(1), 1.);
	  
	  hr->Fit(f6, "0");
	  f6->DrawCopy("same");


	  mean = TMath::Sqrt((f6->GetParameter(0)*f6->GetParameter(0)*f6->GetParameter(1)*f6->GetParameter(1)
			      + f6->GetParameter(3)*f6->GetParameter(3)*f6->GetParameter(4)*f6->GetParameter(4))
			     /
			     (f6->GetParameter(0)*f6->GetParameter(0) + f6->GetParameter(3)*f6->GetParameter(3))
			     );
	  
	  sigma = TMath::Sqrt((f6->GetParameter(0)*f6->GetParameter(0)*f6->GetParameter(2)*f6->GetParameter(2)
			       + f6->GetParameter(3)*f6->GetParameter(3)*f6->GetParameter(5)*f6->GetParameter(5))
			      /
			      (f6->GetParameter(0)*f6->GetParameter(0) + f6->GetParameter(3)*f6->GetParameter(3))
			      );
	  
	  
	  tl->SetTextColor(kBlack);  
	  tl->SetNDC(kTRUE); tl->SetTextSize(0.04);
	  tl->DrawLatex(0.24, 0.85, Form("#mu: %5.3f#pm%5.3f GeV", mean, 0.001));
	  tl->DrawLatex(0.24, 0.81, Form("#sigma: %5.3f#pm%5.3f GeV", sigma, 0.001));
	  tl->DrawLatex(0.24, 0.77, Form("N: %4.2e", hr->GetSumOfWeights()));
	
	} else {

	  tl->SetTextColor(kBlack);  
	  tl->SetNDC(kTRUE); tl->SetTextSize(0.04);
	  tl->DrawLatex(0.24, 0.85, Form("N: %4.2e", hr->GetSumOfWeights()));
	}

      } else {

	hs->DrawCopy("hist");
      }

      legg->Clear();
      
      legge = legg->AddEntry(hs, "Signal (scaled)", "f"); legge->SetTextColor(kBlack);
      legge = legg->AddEntry(hr, "Background", "p");
      
      legge->SetTextColor(EColor[i]);
      legg->Draw();
      
      tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
      tl->SetTextColor(kBlack);    
      tl->DrawLatex(0.22, 0.92, Form("%s",  fSignLeggM[i].Data()));
      
      c0->SaveAs(Form("%s/breco/%s-breco.eps", outDir, fSignM[i].Data()));
      delete hs;
      delete hr;
    }
  }

  gPad->SetLogy(0);
}


// ----------------------------------------------------------------------
void anaBmm::nreco(int o, const char *hist) {


  // -- Norm. channel
  // ================

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TH1D *h1;
  TH1D *h2;
  
  h1 = (TH1D*)(fS[o]->Get(hist))->Clone();
  h2 = (TH1D*)(fS[o]->Get(hist))->Clone();

  emptyBinError(h1);
  
  h1->Scale(fLumiD[0]/fLumiS[o]);
  h2->Scale(fLumiD[0]/fLumiS[o]);
    
  
  c0->Clear();
  shrinkPad(0.15, 0.15);
  setFilledHist(h1, kBlack, kYellow, 1000);
  setFilledHist(h2, kBlack, kYellow, 1000);
  
  f5->SetParameters(0.2*h1->GetMaximum(), 5.25, 0.04, 
		     0.8*h1->GetBinContent(1), 1.);
  
  setTitles(h2, "m_{#mu#mu K} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
  h1->SetMaximum(1.4*h1->GetMaximum());
  h2->SetMaximum(1.4*h2->GetMaximum());
  
  h2->DrawCopy();
  h1->Fit(f5, "0");
  f5->DrawCopy("same");
  
  f10->SetParameters(f5->GetParameter(3), f5->GetParameter(4));
  f10->SetLineStyle(7);
  f10->SetLineColor(108);
  f10->DrawCopy("same");
  
  f11->SetParameters(f5->GetParameter(0), f5->GetParameter(1), f5->GetParameter(2));
  f11->SetLineColor(2);
  f11->SetFillColor(2);
  f11->SetFillStyle(3004);
  f11->DrawCopy("same");

  double mean  = f5->GetParameter(1);
  double sigma = f5->GetParameter(2);

  double meanE  = f5->GetParError(1);
  double sigmaE = f5->GetParError(2);
     
  double totalBG = f10->Integral(h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()))/h1->GetBinWidth(2); 
  double totalSG = (f5->Integral(fMassBs - 0.1, fMassBs + 0.1)/h1->GetBinWidth(2))
                  - (f10->Integral(fMassBs - 0.1, fMassBs + 0.1)/h1->GetBinWidth(2));
  
  tl->SetTextColor(kBlack);  
  tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
  
  tl->DrawLatex(0.16, 0.85, Form("#mu: %5.3f#pm%5.3f GeV", mean, meanE));
  tl->DrawLatex(0.16, 0.79, Form("#sigma: %5.3f#pm%5.3f GeV", sigma, sigmaE));

  tl->SetTextSize(0.05);
  tl->SetTextColor(108); 
  tl->DrawLatex(0.60, 0.61, Form("N_{bg} = %4.2e", totalBG)); tl->SetTextColor(2);  
  tl->DrawLatex(0.60, 0.54, Form("N_{sg} = %4.2e", totalSG)); tl->SetTextColor(1);  
  
  tl->SetTextSize(0.06); tl->SetTextColor(kBlue);
  tl->DrawLatex(0.40, 0.72, Form("%s", fSignLeggS[o].Data()));
  
  
  c0->SaveAs(Form("%s/breco/s%d-nreco-%s.eps", outDir, o, hist));

  delete h1;
}


// ----------------------------------------------------------------------
void anaBmm::jreco(int o) {


  // -- J/Psi peak
  // =============

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TH1D *h1;
  h1 = (TH1D*)(fS[o]->Get("j101"))->Clone();

  h1->Scale(fLumiD[0]/fLumiS[o]);

  c0->Clear();
  shrinkPad(0.15, 0.15);
  setFilledHist(h1, kBlack, kYellow, 1000);

  h1->SetMaximum(1.3*h1->GetMaximum());
  setTitles(h1, "m_{#mu#mu} [GeV]", "events/bin", 0.06, 1.1, 1.3); 

  f1->SetParameters(h1->GetMaximum()*0.8, h1->GetMean(), 0.5*h1->GetRMS(), 
		    h1->GetMaximum()*0.2, h1->GetMean(), 3.*h1->GetRMS());
 
  h1->DrawCopy("hist");
  h1->Fit(f1, "0");
  f1->DrawCopy("same");
  
  double mean = TMath::Sqrt((f1->GetParameter(0)*f1->GetParameter(0)*f1->GetParameter(1)*f1->GetParameter(1)
		      + f1->GetParameter(3)*f1->GetParameter(3)*f1->GetParameter(4)*f1->GetParameter(4))
		     /
		     (f1->GetParameter(0)*f1->GetParameter(0) + f1->GetParameter(3)*f1->GetParameter(3))
		     );
  
  double sigma = TMath::Sqrt((f1->GetParameter(0)*f1->GetParameter(0)*f1->GetParameter(2)*f1->GetParameter(2)
		       + f1->GetParameter(3)*f1->GetParameter(3)*f1->GetParameter(5)*f1->GetParameter(5))
		      /
		      (f1->GetParameter(0)*f1->GetParameter(0) + f1->GetParameter(3)*f1->GetParameter(3))
		      );

  tl->SetTextColor(kBlack);  
  tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
  
  tl->DrawLatex(0.16, 0.85, Form("#mu: %5.3f#pm%5.3f GeV", mean, 0.001));
  tl->DrawLatex(0.16, 0.79, Form("#sigma: %5.3f#pm%5.3f GeV", sigma, 0.001));

  tl->SetTextSize(0.06); tl->SetTextColor(kBlue);

  if ( o == 4 ) {
    tl->DrawLatex(0.19, 0.74, "perfect align.");
    tl->DrawLatex(0.2, 0.67, "J/#psi #rightarrow #mu #mu");
  } else if ( o == 5 ) {
    tl->DrawLatex(0.19, 0.74, "short-term align.");
    tl->DrawLatex(0.2, 0.67, "J/#psi #rightarrow #mu #mu");
  } else {  
    tl->SetTextSize(0.1);
    tl->DrawLatex(0.70, 0.82, Form("MC%i",o));
  }

  c0->SaveAs(Form("%s/breco/s%d-jreco.eps", outDir, o));

  delete h1;
}


//===========================================================================================
// -- MC Validation
//===========================================================================================

void anaBmm::mcVal(const char *hname, int mode, double x, double y) {
  
  gStyle->SetOptTitle(0);

  // -- compare files1 with files2 for signal AND equivalent for background
  TString sgfiles1[]    = { TString("csg-003") , TString("csg-003") };
  TString sgfiles2[]    = { TString("csg-001") , TString("csg-003h") };

  TString bgfiles1[]    = { TString("cbg-003") , TString("cbg-003") };
  TString bgfiles2[]    = { TString("cbg-002") , TString("cbg-003h") };

  // -- sample descriptions for sg AND bg
  TString samples1[] = { TString("CMSSW (CSA07)"), TString("CSA07") }; 
  TString samples2[] = { TString("ORCA (priv.)"), TString("Spring07") };

  TString labels[]   = { TString("orca"), TString("spring07") };

  TH1D *h0, *h1;
  int n(0), m(0);
  char name_n[200], name_m[200], hist[200];
  char f_n[200], f_m[200];

  sprintf(hist,"%s", hname);  // ??? Crashes when using hname instead of hist ???

  //  int num = sizeof(sgfiles1) - 1;  // not working
  int num = 2;

  if ( (sizeof(sgfiles2) != sizeof(sgfiles1)) || (sizeof(labels) != sizeof(sgfiles1)) ||
       (sizeof(samples1) != sizeof(sgfiles1)) || (sizeof(samples2) != sizeof(sgfiles1))||
       (sizeof(bgfiles1) != sizeof(sgfiles1)) || (sizeof(bgfiles2) != sizeof(sgfiles1)) ) {

    cout << " ====> Error: Sizes of arrays are not consistent !!!" << endl;
    return; 
  }

  int logy(0); 
  if (mode > 99) {
    mode -= 100;
    logy = 1;
  }  
  
  for (int i = 0; i < num; i++ ) {

    // -- in case of missing histograms in old samples
    if (mode == 10) {
      if ( i == 0 ) { continue; }
    }


    sprintf(name_n, "%s", samples1[i].Data());
    sprintf(name_m, "%s", samples2[i].Data());

    // --------------
    // -- Signal
    // --------------
    sprintf(f_n,"%s",sgfiles1[i].Data());
    sprintf(f_m,"%s",sgfiles2[i].Data());

    n = findIndex(f_n);
    m = findIndex(f_m);  

    if ( n > -1 && m > -1 && n < nSg && m < nSg) { 
  
      fS[n]->cd();
      h0 = (TH1D*)(fS[n]->Get(hist))->Clone();
      if (!h0) {
	cout << "Histogram with name " << hist << " does not exist" << endl;
	return;
      }

      fS[m]->cd();
      h1 = (TH1D*)(fS[m]->Get(hist))->Clone();
      if (!h1) {
	cout << "Histogram with name " << hist << " does not exist" << endl;
	return;
      }

      h0->Scale(1./h0->GetSumOfWeights());
      h1->Scale(1./h1->GetSumOfWeights());
  
      c0->Clear();
      gPad->SetLogy(logy);
      shrinkPad(0.15, 0.2);
      setHist(h0, kBlack, 20, 2);
      setHist(h1, kBlue, 24, 2);
      setTitles(h0, h0->GetXaxis()->GetTitle(), h0->GetYaxis()->GetTitle(), 0.06, 1.1, 1.8);
      h0->SetMaximum(1.2*(h1->GetMaximum() > h0->GetMaximum()? h1->GetMaximum(): h0->GetMaximum()));

      h0->DrawCopy();
      h1->DrawCopy("same");

      if (x > 0) {
	legg = new TLegend(x, y, x+0.4, y+0.15);
	legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  
	legg->SetFillColor(0); legg->SetFillStyle(1000); 
	legge = legg->AddEntry(h1, Form("%s", name_m), "p"); legge->SetTextColor(kBlue);
	legge = legg->AddEntry(h0, Form("%s", name_n), "p"); legge->SetTextColor(kBlack);
	legg->Draw();
      }

      c0->Draw(); 
      c0->Update();

      c0->SaveAs(Form("%s/mcval/sgval-%s-%s.eps", outDir, labels[i].Data(), hist));

      delete h0;
      delete h1;  

    } else {
      if ( n < 0 ) cout << " ====> Error: couldn't index of find file " << f_n << endl;
      if ( m < 0 ) cout << " ====> Error: couldn't index of find file " << f_m << endl;
    }  


    // --------------
    // -- Background
    // --------------  
    sprintf(f_n,"%s",bgfiles1[i].Data());
    sprintf(f_m,"%s",bgfiles2[i].Data());

    n = findIndex(f_n);
    m = findIndex(f_m);  

    if ( n > -1 && m > -1 && n < nMc && m < nMc) {   

      fM[n]->cd();
      h0 = (TH1D*)(fM[n]->Get(hist))->Clone();
      if (0 == h0) {
	cout << "Histogram with name " << hist << " does not exist" << endl;
	return;
      }

      fM[m]->cd();
      h1 = (TH1D*)(fM[m]->Get(hist))->Clone();
      if (0 == h1) {
	cout << "Histogram with name " << hist << " does not exist" << endl;
	return;
      }

      h0->Scale(1./h0->GetSumOfWeights());
      h1->Scale(1./h1->GetSumOfWeights());
  
      c0->Clear();
      gPad->SetLogy(logy);
      shrinkPad(0.15, 0.2);
      setHist(h0, kBlack, 20, 2);
      setHist(h1, kBlue, 24, 2);
      setTitles(h0, h0->GetXaxis()->GetTitle(), h0->GetYaxis()->GetTitle(), 0.06, 1.1, 1.8);
      h0->SetMaximum(1.2*(h1->GetMaximum() > h0->GetMaximum()? h1->GetMaximum(): h0->GetMaximum()));

      h0->DrawCopy();
      h1->DrawCopy("same");

      if (x > 0) {
	legg = new TLegend(x, y, x+0.4, y+0.1);
	legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);
	legg->SetFillColor(0); legg->SetFillStyle(1000); 
	legge = legg->AddEntry(h1, Form("%s", name_m), "p"); legge->SetTextColor(kBlue);
	legge = legg->AddEntry(h0, Form("%s", name_n), "p"); legge->SetTextColor(kBlack);
	legg->Draw();
      }

      c0->Draw(); 
      c0->Update();

      c0->SaveAs(Form("%s/mcval/bgval-%s-%s.eps", outDir, labels[i].Data(), hist));

      delete h0;
      delete h1;

    } else {
      if ( n < 0 ) cout << " ====> Error: couldn't find index of file " << f_n << endl;
      if ( m < 0 ) cout << " ====> Error: couldn't find index of file " << f_m << endl;
    }
  }
}


//===========================================================================================
// -- Distributions
//===========================================================================================
void anaBmm::showDistribution(const char *hname, int mode, double x, double y) { 

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  c0->Clear();
  shrinkPad(0.15, 0.2);
  
  //  TH1 *hm = sumHistMC(hname);
  TH1 *hm = (TH1D*)fM[bgIndex]->Get(hname);
  setHist(hm, kBlack, kBlack);
  setTitles(hm, hm->GetXaxis()->GetTitle(), hm->GetYaxis()->GetTitle(), 0.06, 1.1, 1.5);

  TH1D *hs = new TH1D(*((TH1D*)fS[sgIndex]->Get(hname)));
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
    hs->Scale(fLumiD[0]/fLumiS[sgIndex]); 
    hm->Scale(fLumiD[0]/fLumiM[bgIndex]); 
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
 
  c0->SaveAs(Form("%s/dist/dist-%s.eps", outDir, hname));

  c0->SetLogy(0);

}
  

//===========================================================================================
// -- Production prosesse (GGF, FEX, GSP)
//===========================================================================================

void anaBmm::showProcesses(int signal) {

  TString label;
  if (1 == signal) {
    fS[sgIndex]->cd();
    label = TString("sproc");
  } else {
    fM[bgIndex]->cd();
    label = TString("bproc");
  }


  const char *cuts = "pt>5 && ptl1>4";

  plotProcesses(cuts, "pt", "p_{T, B} [GeV]", 10, 0., 50., 1);
  c0->SaveAs(Form("%s/proc/%s-pt.eps", outDir, label.Data()));

  plotProcesses(cuts, "iso", "I", 20, 0., 1., 0);
  c0->SaveAs(Form("%s/proc/%s-isolation.eps", outDir, label.Data()));

  plotProcesses(cuts, "rmm", "#Delta R(#mu,#mu)", 10, 0., 2., 0);
  c0->SaveAs(Form("%s/proc/%s-rmm.eps", outDir, label.Data()));

  plotProcesses(cuts, "mass", "m_{#mu#mu} [GeV]", 20, 5., 6., 0);
  c0->SaveAs(Form("%s/proc/%s-mass.eps", outDir, label.Data()));

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


//===========================================================================================
// -- Efficiency tables
//===========================================================================================

void anaBmm::effTable(TFile *f, const char *tag) {

  // -- Lumi normalisation scaling factor
  double SF(0.), comb(0.);
  char sname[200];

  if (!strcmp(tag, "mysg")) {

    SF = fLumiD[0]/fLumiS[sgIndex];
    sprintf(sname, "s%i", sgIndex);
 
  } else if (!strcmp(tag, "mymc")) {
   
    SF = fLumiD[0]/fLumiM[bgIndex];
    sprintf(sname, "m%i", bgIndex);
  
  } else if (!strcmp(tag, "mynsg")) {
   
    SF = fLumiD[0]/fLumiS[normSgIndex];
    sprintf(sname, "s%i", normSgIndex);
  
  } else if (!strcmp(tag, "mynmc")) {
   
    SF = fLumiD[0]/fLumiM[normBgIndex];
    sprintf(sname, "m%i", normBgIndex);
  
  } else {
    
    SF = 1; comb = 1;
  }

  
  // -- Analysis Efficiency
  // ----------------------

  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% -- Efficiency" << endl;

  TH1D *h = (TH1D*)f->Get("AR1");

  // -- rare/combined backgrounds 
  if (comb) {

    TH1 *sumAR1 = sumHistMC("AR1", 2, tag);
    h = (TH1D*)sumAR1;
  }
  
  double n     = h->GetBinContent(h->FindBin(0.1));
  double norm  = h->GetBinContent(h->FindBin(0.1));  // FIXME: Should be Bin(0.1) but old files have empty Bin(0.1)
  //  double norm  = h->GetBinContent(h->FindBin(1.1));  // FIXME: Should be Bin(0.1) but old files have empty Bin(0.1)
  double nkin  = h->GetBinContent(h->FindBin(1.1));
  double nexp  = SF*norm;

  cout << " #expected events [" << fLumiD[0] << "/fb]: " << SF << " x " << norm 
       << " = " << Form("%5.1f",SF*norm) << endl;


  // -- analysis efficiencies
  // ------------------------
  // -- exp. events, from (0.1)
  //  n  = h->GetBinContent(h->FindBin(1.1));
  n  = h->GetBinContent(h->FindBin(0.1));

  double nevt  = SF*norm;
  double enorm = n/norm;
  double dnorm = dEff(int(n), int(norm));
  double ekin  = n/nkin;
  double dkin  = dEff(int(n), int(nkin));

  nevt  = SF*norm;
  OUT << Form("%s", (formatTex(nevt, "nA0" , tag)).Data()) << endl;

  // -- generator kinematics
  // ------------------------
  // -- exp. events, from (1.1)
  n  = h->GetBinContent(h->FindBin(1.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nAkin"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eAkin"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eAkinE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cAkin"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cAkinE" , tag)).Data()) << endl;


  // -- L1 trigger
  // --------------
  // -- exp. events, from (11.1)
  n  = h->GetBinContent(h->FindBin(11.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nAL1"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eAL1"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eAL1E" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cAL1"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cAL1E" , tag)).Data()) << endl;


  // -- HLT trigger
  // ---------------
  // -- exp. events, from (21.1)
  n  = h->GetBinContent(h->FindBin(21.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nAHLT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eAHLT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eAHLTE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cAHLT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cAHLTE" , tag)).Data()) << endl;


  // -- Good event (rec. cand. & PV)
  // -------------------------------
  // -- exp. events, from (51.1)
  n  = h->GetBinContent(h->FindBin(51.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nABX"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eABX"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eABXE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cABX"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cABXE" , tag)).Data()) << endl;


  // -- Lepton pT
  // -------------
  // -- exp. events, from (100.1)
  n  = h->GetBinContent(h->FindBin(100.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nALPT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eALPT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eALPTE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cALPT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cALPTE" , tag)).Data()) << endl;


  // -- Delta R(mm)
  // ---------------
  // -- exp. events, from (101.1)
  n  = h->GetBinContent(h->FindBin(101.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nARmm"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eARmm"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eARmmE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cARmm"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cARmmE" , tag)).Data()) << endl;


  // -- ???? B-candidate:  mass ????
  // --------------------------------
  // -- exp. events, from (110.1)
  n  = h->GetBinContent(h->FindBin(110.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nABmass"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eABmass"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eABmassE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cABmass"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cABmassE" , tag)).Data()) << endl;


  // -- B-candidate: pT
  // -------------------
  // -- exp. events, from (120.1)
  n  = h->GetBinContent(h->FindBin(120.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nABpT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eABpT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eABpTE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cABpT"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cABpTE" , tag)).Data()) << endl;


  // -- B-candidate: eta
  // --------------------
  // -- exp. events, from (121.1)
  n  = h->GetBinContent(h->FindBin(121.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nABeta"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eABeta"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eABetaE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cABeta"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cABetaE" , tag)).Data()) << endl;


  // -- B-candidate: cos(alpha)
  // --------------------------
  // -- exp. events, from (122.1)
  n  = h->GetBinContent(h->FindBin(122.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nABcos"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eABcos"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eABcosE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cABcos"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cABcosE" , tag)).Data()) << endl;


  // -- B-candidate: Decay length significance
  // ------------------------------------------
  // -- exp. events, from (123.1)
  n  = h->GetBinContent(h->FindBin(123.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nABlxy"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eABlxy"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eABlxyE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cABlxy"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cABlxyE" , tag)).Data()) << endl;

  double fnorm = n;
  double eff = n/norm;

  // -- B-candidate: Isolation (pT-sum of track in cone)
  // ----------------------------------------------------
  // -- exp. events, from (124.1)
  n  = h->GetBinContent(h->FindBin(124.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nABiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eABiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eABisoE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cABiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cABisoE" , tag)).Data()) << endl;


  // -- B-candidate: Vertex quality (chi2)
  // --------------------------------------
  // -- exp. events, from (125.1)
  n  = h->GetBinContent(h->FindBin(125.1));

  nevt  = SF*n;
  enorm = n/norm;
  dnorm = dEff(int(n), int(norm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nABchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eABchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eABchi2E" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cABchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cABchi2E" , tag)).Data()) << endl;


  // ===================================================================================
  // -- Signal efficiency (for error estimation)
  // ===================================================================================
  if (!strcmp(tag, "mysg")) {
    fEsg0  = n/norm;
    fEsgE0 = dEff(int(n), int(norm));
  }

  // ===================================================================================

  // -- factorizing cuts: Iso & Vertex efficiency
  // --------------------------------------------

  double pnorm = h->GetBinContent(h->FindBin(220.1));

  if ( pnorm == 0. ) {

    fillInTheRest(tag);
    return;
  }
 
  // -- Isolation cut efficiciency
  // ------------------------------
  n     = h->GetBinContent(h->FindBin(224.1));

  double eff1 = n/pnorm;

  nevt  = SF*n;
  enorm = n/pnorm;
  dnorm = dEff(int(n), int(pnorm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nAfBiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eAfBiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eAfBisoE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cAfBiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cAfBisoE" , tag)).Data()) << endl;


  // -- Vertex cut efficiciency
  // ---------------------------
  n    = h->GetBinContent(h->FindBin(225.1));

  double eff2 = n/pnorm;

  nevt  = SF*n;
  enorm = n/pnorm;
  dnorm = dEff(int(n), int(pnorm));
  ekin  = n/nkin;
  dkin  = dEff(int(n), int(nkin));

  OUT << Form("%s", (formatTex(nevt,  "nAfBchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(enorm, "eAfBchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dnorm, "eAfBchi2E" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(ekin,  "cAfBchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(dkin,  "cAfBchi2E" , tag)).Data()) << endl;


  // ===================================================================================
  // -- Factorizing cuts efficiency
  // ===================================================================================

  // -- fact. Isolation
  // -------------------
  double nfact   = SF*fnorm*eff1;
  double nfactE  = SF*TMath::Sqrt(fnorm)*eff1;
  double efact   = nfact/nexp;  // = eff*eff1 for single channel
  double efactE  = (fEsgE0/fEsg0)*efact;

  if (comb) {

    nfact      = h->GetBinContent(h->FindBin(424.1));
    nfactE     = -9999; 
    efact      = nfact/nexp;
    efactE     = (fEsgE0/fEsg0)*efact;
  }
 
  OUT << Form("%s", (formatTex(nfact,  "nTotEffIso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(nfactE, "nTotEffIsoE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(efact,  "eTotEffIso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(efactE, "eTotEffIsoE" , tag)).Data()) << endl;


  // -- fact. Vertex 
  // ---------------- 
  nfact   = SF*fnorm*eff2;
  nfactE  = SF*TMath::Sqrt(fnorm)*eff2;
  efact   = nfact/nexp;  // = eff*eff2 for single channel
  efactE  = (fEsgE0/fEsg0)*efact;

  if (comb) {

    nfact      = h->GetBinContent(h->FindBin(426.1));
    nfactE     = -9999; 
    efact      = nfact/nexp;
    efactE     = (fEsgE0/fEsg0)*efact;
  }
 
  OUT << Form("%s", (formatTex(nfact,  "nTotEffChi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(nfactE, "nTotEffChi2E" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(efact,  "eTotEffChi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(efactE, "eTotEffChi2E" , tag)).Data()) << endl;


  // -- all Cuts
  // ------------
  nfact   = SF*fnorm*eff1*eff2;
  nfactE  = SF*TMath::Sqrt(fnorm)*eff1*eff2;
  efact   = nfact/nexp;  // = eff*eff2 for single channel
  efactE  = (fEsgE0/fEsg0)*efact;

  if (comb) {

    nfact      = h->GetBinContent(h->FindBin(428.1));
    nfactE     = -9999; 
    efact      = nfact/nexp;
    efactE     = (fEsgE0/fEsg0)*efact;
  }
 
  OUT << Form("%s", (formatTex(nfact,  "nAllCutsFact"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(nfactE, "nAllCutsFactE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(efact,  "eAllCutsFact"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(efactE, "eAllCutsFactE" , tag)).Data()) << endl;


  // ===================================================================================
  // -- Mass reduction in 100 MeV window
  // ===================================================================================

  double massRed0  = massReduction("c030", tag);
  double massRed1  = massReduction("c130", tag);
  double massRed2  = massReduction("c230", tag);
  double massRed5  = massReduction("c530", tag);

  double massRed   = massReduction("c330", tag);


  // --- FIX ???
  if ( !strcmp(tag, "mysg") || !strcmp(tag, "sg") ) {

//     massRed = massRed5;
  }

  // -- mass window
  // --------------

  double nmass   = massRed*nfact;
  double nmassE  = 1.6*massRed*nfactE;   // XXXX FIXME XXXX
  double emass   = nmass/nexp;
  double emassE  = (fEsgE0/fEsg0)*emass; // XXXX FIXME XXXX

  
  if (comb) {

    massRed = massRed5; 

    nmass   = massRed*nfact;
    nmassE  = -9999;                // XXXX FIXME XXXX
    emass   = efact*nmass/nfact;
    emassE  = (fEsgE0/fEsg0)*emass; // XXXX FIXME XXXX
  }
  
  OUT << Form("%s", (formatTex(nmass,  "nMassAllCutsFact"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(nmassE, "nMassAllCutsFactE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(emass,  "eMassAllCutsFact"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(emassE, "eMassAllCutsFactE" , tag)).Data()) << endl;

  OUT  << Form("\\vdef{massReduction:%s}  {\\ensuremath{{%3.2f} } }", tag, massRed) << endl;


  // ===================================================================================
  // -- Event reduction
  // ===================================================================================

  if (!strcmp(tag, "mysg")) {

    fNsg  = nmass;
    fNsgE = nmassE;
  }

  if (!strcmp(tag, "mymc")) {

    fNbg  = nmass;
    fNbgE = nmassE;
  }

  if (!strcmp(tag, "mynsg")) {

    fNsg_norm  = nmass;
    fNsgE_norm = nmassE;
  }

  if (!strcmp(tag, "mynmc")) {

    fNbg_norm  = nmass;
    fNbgE_norm = nmassE;
  }

  if (!strcmp(tag, "r0")) {

    fNrbg  = nmass;
    fNrbgE = nmassE;
  }


  // ===================================================================================
  // -- Factorizing efficiency of all cuts
  // ===================================================================================

  if (!strcmp(tag, "mysg")) {
    fEsg  = emass;
    fEsgE = emassE;
  }

  if (!strcmp(tag, "mynsg")) {
    fEsg_norm  = emass;
    fEsgE_norm = emassE;
  }


  // ===================================================================================
  // -- Statistics of the different candidate selection modes
  // ===================================================================================

  double n_cuts(0), n_fact(0);
  double e_iso(0.), e_vtx(0.);

  if (!strcmp(tag, "mysg") || !strcmp(tag, "mymc") || !strcmp(tag, "mynsg") || !strcmp(tag, "mynmc")   ) {

    for (int i= 1; i < 5; i++) {

    
      h = (TH1D*)f->Get(Form("OR1_%i",i));

      n_cuts  = h->GetBinContent(h->FindBin(126.1));
      // for bg:  n_fact  = massRed*e_iso*e_vtx*h->GetBinContent(h->FindBin(125.1)); ????

      if (h->GetBinContent(h->FindBin(220.1)) > 0) {
    
	e_iso   = h->GetBinContent(h->FindBin(224.1))/h->GetBinContent(h->FindBin(220.1));
	e_vtx   = h->GetBinContent(h->FindBin(225.1))/h->GetBinContent(h->FindBin(220.1));
	n_fact  = e_iso*e_vtx*h->GetBinContent(h->FindBin(210.1));
	// for bg:  n_fact  = massRed*e_iso*e_vtx*h->GetBinContent(h->FindBin(123.1)); ????
    
      } else {
      
	e_iso = 0., e_vtx = 0.; n_fact = 0.;

      }

      nevt  = SF*n_cuts;
      enorm = n_cuts/norm;
      dnorm = dEff(int(n_cuts), int(norm));

      OUT << Form("%s", (formatTex(n_cuts, Form("tsel%i",i),  tag)).Data()) << endl;
      OUT << Form("%s", (formatTex(nevt,   Form("nsel%i",i),  tag)).Data()) << endl;
      OUT << Form("%s", (formatTex(enorm,  Form("esel%i",i),  tag)).Data()) << endl;
      OUT << Form("%s", (formatTex(dnorm,  Form("eselE%i",i), tag)).Data()) << endl;

      nevt  = SF*n_fact;
      enorm = n_fact/norm;
      dnorm = dEff(int(n_fact), int(norm));

      OUT << Form("%s", (formatTex(n_fact, Form("tselFact%i",i),  tag)).Data()) << endl;
      OUT << Form("%s", (formatTex(nevt,   Form("nselFact%i",i),  tag)).Data()) << endl;
      OUT << Form("%s", (formatTex(enorm,  Form("eselFact%i",i),  tag)).Data()) << endl;
      OUT << Form("%s", (formatTex(dnorm,  Form("eselFactE%i",i), tag)).Data()) << endl;

    }
  }

  OUT.close();
}

// ---------------------------------------------------------------------------------------------

double anaBmm::massReduction(const char *hist, const char *sel, double window) {

  if (window < 0.) {
    window = masswi;
  }

  cout << "---------------------------------------------------------------------------" << endl;

  if ( !strcmp(sel, "mymc") ) {

    TH1D *h1 = (TH1D*)fM[bgIndex]->Get(hist)->Clone();
    TH1D *h2 = (TH1D*)fM[bgIndex]->Get(hist)->Clone();

    emptyBinError(h1);

    h2->GetXaxis()->SetRangeUser(4.8, 6.0);
    h1->GetXaxis()->SetRangeUser(4.8, 6.0);
    h1->Fit("pol1");

    TF1 *f1 = (TF1*)h1->GetFunction("pol1");

    h2->Draw();
    f1->Draw("same");
    c0->SaveAs(Form("%s/massReduction/massReduction-%s-%s.eps", outDir, sel, hist));
    
    double total = f1->Integral(h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()))/h1->GetBinWidth(2); 
    double massW = f1->Integral(fMassBs - window, fMassBs + window)/h1->GetBinWidth(2);
    double massRed = massW/total;
    
    double histT =  h1->GetSumOfWeights();
    double histW =  h1->Integral(h1->FindBin(fMassBs - window), h1->FindBin(fMassBs + window));
    double histRed = histW/histT;

    cout << endl << " --> " << hist
	 << " massRed " << massRed << " =  massW (" << massW << ") / total ("  << total << ")" << endl;

    cout << " --> " << hist
	 << " histRed " << histRed << " =  histW (" << histW << ") / histT ("  << histT << ")" << endl;

    cout << " ... in mass window +/- 100 MeV around " << fMassBs  << endl<< endl; 

    h1->GetListOfFunctions()->Clear();
  
    cout << "return mass reduction " << massRed << endl;
    return massRed;

  } else if ( !strcmp(sel, "mynmc") ) {

    TH1D *h1 = (TH1D*)fM[normBgIndex]->Get(hist)->Clone();
    TH1D *h2 = (TH1D*)fM[normBgIndex]->Get(hist)->Clone();

    emptyBinError(h1);

    h2->GetXaxis()->SetRangeUser(4.8, 6.0);
    h1->GetXaxis()->SetRangeUser(4.8, 6.0);
    h1->Fit("pol1");

    TF1 *f1 = (TF1*)h1->GetFunction("pol1");

    h2->Draw();
    f1->Draw("same");
    c0->SaveAs(Form("%s/massReduction/massReduction-%s-%s.eps", outDir, sel, hist));

    double total = f1->Integral(h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()))/h1->GetBinWidth(2); 
    double massW = f1->Integral(fMassBp - window, fMassBp + window)/h1->GetBinWidth(2);
    double massRed = massW/total;
    
    double histT =  h1->GetSumOfWeights();
    double histW =  h1->Integral(h1->FindBin(fMassBp - window), h1->FindBin(fMassBp + window));
    double histRed = histW/histT;

    cout << endl << " --> " << hist
	 << " massRed " << massRed << " =  massW (" << massW << ") / total ("  << total << ")" << endl;

    cout << " --> " << hist
	 << " histRed " << histRed << " =  histW (" << histW << ") / histT ("  << histT << ")" << endl;

    cout << " ... in mass window +/- 100 MeV around " << fMassBp  << endl<< endl; 

    h1->GetListOfFunctions()->Clear();
  
    cout << "return mass reduction " << massRed << endl;
    return massRed;

  } else if ( !strcmp(sel, "mysg") ) {

    TH1D *h1 = (TH1D*)fS[sgIndex]->Get(hist)->Clone();
    emptyBinError(h1);

    h1->Scale(fLumiD[0]/fLumiS[sgIndex]);
    
    h1->GetXaxis()->SetRangeUser(4.8, 6.);

    f1->SetParameters(h1->GetMaximum()*0.8, h1->GetMean(), 0.5*h1->GetRMS(), 
		      h1->GetMaximum()*0.2, h1->GetMean(), 3.*h1->GetRMS());
      

    shrinkPad(0.15, 0.15);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    setFilledHist(h1, kBlack, kYellow, 1000); 
    setTitles(h1, "m_{#mu#mu} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
    h1->SetMaximum(1.4*h1->GetMaximum());
    
    h1->DrawCopy("hist");
    h1->Fit(f1, "0");
    f1->DrawCopy("same");
    
    double mean = TMath::Sqrt((f1->GetParameter(0)*f1->GetParameter(0)*f1->GetParameter(1)*f1->GetParameter(1)
			  + f1->GetParameter(3)*f1->GetParameter(3)*f1->GetParameter(4)*f1->GetParameter(4))
			 /
			 (f1->GetParameter(0)*f1->GetParameter(0) + f1->GetParameter(3)*f1->GetParameter(3))
			 );

    double  sigma = TMath::Sqrt((f1->GetParameter(0)*f1->GetParameter(0)*f1->GetParameter(2)*f1->GetParameter(2)
			   + f1->GetParameter(3)*f1->GetParameter(3)*f1->GetParameter(5)*f1->GetParameter(5))
			  /
			  (f1->GetParameter(0)*f1->GetParameter(0) + f1->GetParameter(3)*f1->GetParameter(3))
			  );

    
    tl->SetTextColor(kBlack);  
    tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
  
    tl->DrawLatex(0.16, 0.85, Form("#mu: %5.3f#pm%5.3f GeV", mean, 0.001));
    tl->DrawLatex(0.16, 0.79, Form("#sigma: %5.3f#pm%5.3f GeV", sigma, 0.001));
  
    if (h1->GetSumOfWeights()<0.001) {
      tl->DrawLatex(0.16, 0.73, Form("N: %4.2e", h1->GetSumOfWeights()));
    } else if (h1->GetSumOfWeights() < 1.) {
      tl->DrawLatex(0.16, 0.73, Form("N: %5.3f", h1->GetSumOfWeights()));
    } else if (h1->GetSumOfWeights() < 1000.) {
      tl->DrawLatex(0.16, 0.73, Form("N: %5.1f", h1->GetSumOfWeights()));
    } else {
      tl->DrawLatex(0.16, 0.73, Form("N: %4.2e", h1->GetSumOfWeights()));
    }

    tl->SetTextSize(0.06); tl->SetTextColor(kBlue);
    tl->DrawLatex(0.56, 0.65, Form("%s", fSignLeggS[sgIndex].Data()));

    c0->SaveAs(Form("%s/massReduction/massReduction-%s-%s.eps", outDir, sel, hist));

    double total = f1->Integral(h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()))/h1->GetBinWidth(2); 
    double massW = f1->Integral(fMassBs - window, fMassBs + window)/h1->GetBinWidth(2);
    double massRed = massW/total;
    
    double histT =  h1->GetSumOfWeights();
    double histW =  h1->Integral(h1->FindBin(fMassBs - window), h1->FindBin(fMassBs + window));
    double histRed = histW/histT;

    
    cout << endl << " --> " << hist
	 << " massRed " << massRed << " =  massW (" << massW << ") / total ("  << total << ")" << endl;

    cout << " --> " << hist
	 << " histRed " << histRed << " =  histW (" << histW << ") / histT ("  << histT << ")" << endl;

    cout << " Peak at " << mean << ", sigma " << sigma << endl; 
    cout << " ...off by " << Form("%4.3f", fMassBs - mean) << " from m_Bs0 " << fMassBs << endl<< endl; 
    
    h1->GetListOfFunctions()->Clear();
    
    cout << "return mass reduction " << massRed << endl;

    return massRed;

  } else if ( !strcmp(sel, "mynsg") ) {

    TH1D *h1 = (TH1D*)fS[normSgIndex]->Get(hist)->Clone();

    emptyBinError(h1);

    h1->GetXaxis()->SetRangeUser(4.8, 6.0);

    f6->SetParameters(0.8*h1->GetMaximum(), 5.3, 0.01, 
		      0.8*h1->GetBinContent(1), 1.);

    h1->DrawCopy();

    h1->Fit(f6, "0");
    f6->DrawCopy("same");

    c0->SaveAs(Form("%s/massReduction/massReduction-%s-%s.eps", outDir, sel, hist));
    
    double total = f6->Integral(h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()))/h1->GetBinWidth(2); 
    double massW = f6->Integral(fMassBp - window, fMassBp + window)/h1->GetBinWidth(2);
    double massRed = massW/total;
    
    double histT =  h1->GetSumOfWeights();
    double histW =  h1->Integral(h1->FindBin(fMassBp - window), h1->FindBin(fMassBp + window));
    double histRed = histW/histT;
    
    double mean = TMath::Sqrt((f6->GetParameter(0)*f6->GetParameter(0)*f6->GetParameter(1)*f6->GetParameter(1)
			  + f6->GetParameter(3)*f6->GetParameter(3)*f6->GetParameter(4)*f6->GetParameter(4))
			 /
			 (f6->GetParameter(0)*f6->GetParameter(0) + f6->GetParameter(3)*f6->GetParameter(3))
			 );

    double  sigma = TMath::Sqrt((f6->GetParameter(0)*f6->GetParameter(0)*f6->GetParameter(2)*f6->GetParameter(2)
			   + f6->GetParameter(3)*f6->GetParameter(3)*f6->GetParameter(5)*f6->GetParameter(5))
			  /
			  (f6->GetParameter(0)*f6->GetParameter(0) + f6->GetParameter(3)*f6->GetParameter(3))
			  );

    cout << endl << " --> " << hist
	 << " massRed " << massRed << " =  massW (" << massW << ") / total ("  << total << ")" << endl;

    cout << " --> " << hist
	 << " histRed " << histRed << " =  histW (" << histW << ") / histT ("  << histT << ")" << endl;

    cout << " Peak at " << mean << ", sigma " << sigma << endl; 
    cout << " ...off by " << Form("%4.3f", fMassBp - mean) << " from m_Bs0 " << fMassBp << endl<< endl; 
    
    h1->GetListOfFunctions()->Clear();
    
    cout << "return mass reduction " << massRed << endl;
    return massRed;

  } else {    // ----------- Combined backgrounds ------------------
    
    double histT(0.), histW(0.), histRed(0.);
    double fnorm_Ch(0.), eff1_Ch(0.), eff2_Ch(0.);
  
    int nhist = nMc;
    TH1D *h1[nhist];

    for (int i = 1; i < nhist; ++i) {

      int accept = checkIf(i, sel);
      
      if ( accept ) {
      
	h1[i]  = (TH1D*)fM[i]->Get(hist)->Clone();
	
	channelEff(fM[i], fnorm_Ch, eff1_Ch, eff2_Ch);
	h1[i]->Scale((fMisIdM[i]*fLumiD[0]/fLumiM[i])*eff1_Ch*eff2_Ch);
	histT +=  h1[i]->GetSumOfWeights();
	histW +=  h1[i]->Integral(h1[i]->FindBin(fMassBs - window), h1[i]->FindBin(fMassBs + window));
      }
    }

    if ( histT ) { histRed = histW/histT; }

    
    cout << endl << " --> " << hist << " *** sel " << sel << " ***"
	 << " histRed " << histRed << " =  histW (" << histW << ") / histT ("  << histT << ")" << endl;

    cout << " ... in mass window +/- 100 MeV around " << fMassBs  << endl<< endl; 
  
    cout << "return hist. reduction " << endl;
    return histRed;
  }
}


//===========================================================================================
// -- Overlays
//===========================================================================================

void anaBmm::bgOverlay(const char *hist, const int npers) {

  int EColor[]     = {   13,    2,    4,    108,  6,  107,   93,   1, 
			 13,    2,    4,    108,  6,  107,   93,   1   }; 
  
  double min(9999.), max(0.), maxSoW(0.);
  double fnorm_Ch, eff1_Ch, eff2_Ch;  

  // -- signal
  TH1D *hs = (TH1D*)fS[sgIndex]->Get(hist)->Clone();

  if ( !strcmp(hist,"c530") ) {

    channelEff(fS[sgIndex], fnorm_Ch, eff1_Ch, eff2_Ch);
    hs->Scale((fLumiD[0]/fLumiS[sgIndex])*eff1_Ch*eff2_Ch);

  } else {
    
    hs->Scale(fLumiD[0]/fLumiS[sgIndex]);
  } 

  setFilledHist(hs, kBlack, kYellow, 1001, 2);
  setTitles(hs, "m_{#mu#mu} [GeV]", "events/bin", 0.06, 1.1, 1.3); 

  if (hs->GetMinimum(1.e-20) < min ) { min = hs->GetMinimum(1.e-20);      }
  if (hs->GetMaximum()       > max ) { max = hs->GetMaximum();            }
 
 
  // -- background
  TH1D *hb = (TH1D*)fM[bgIndex]->Get(hist)->Clone();

  if ( !strcmp(hist,"c530") ) {

    channelEff(fM[bgIndex], fnorm_Ch, eff1_Ch, eff2_Ch);
    hb->Scale((fLumiD[0]/fLumiM[bgIndex])*eff1_Ch*eff2_Ch);
  
  } else {
    
    hb->Scale(fLumiD[0]/fLumiM[bgIndex]);
  }

  hb->SetLineStyle(kDashed);
  setFilledHist(hb, 13, 13, 3005, 2);

  if (hb->GetMinimum(1.e-20) < min ) { min = hb->GetMinimum(1.e-20);      }
  if (hb->GetMaximum()       > max ) { max = hb->GetMaximum();            }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetTopMargin(0.25);
  gPad->SetRightMargin(0.2);

  gPad->SetLogy(1);

  // -- Other backgrounds
  const int nhist = nMc;
  const int sets = int(nhist/npers)+1;

  TH1D *bg[nhist];

  for (int i = 0; i < nhist; ++i) {

    bg[i] = (TH1D*)fM[i]->Get(hist)->Clone();
    
    if ( i == bgIndex || !strcmp(fTypeM[i], "nmc") ) { continue; }

    if ( !strcmp(hist,"c530") ) {

      channelEff(fM[i], fnorm_Ch, eff1_Ch, eff2_Ch);
      bg[i]->Scale((fMisIdM[i]*fLumiD[0]/fLumiM[i])*eff1_Ch*eff2_Ch);
    
    } else {

      bg[i]->Scale(fMisIdM[i]*fLumiD[0]/fLumiM[i]);
    }

    if (bg[i]->GetSumOfWeights()  > maxSoW )  { maxSoW = bg[i]->GetSumOfWeights();  }
    if (bg[i]->GetMinimum(1.e-20) < min )     { min = bg[i]->GetMinimum(1.e-20); }
    if (bg[i]->GetMaximum()       > max )     { max = bg[i]->GetMaximum();       }
  }

  // -- Draw & Legends
  legg = new TLegend(0.49,0.62,0.99,0.99);
  legg->SetFillStyle(1001); legg->SetBorderSize(0); legg->SetTextSize(0.035);  legg->SetFillColor(10); 

  int i(0), cont (0), note(0);

  for (int k = 0; (k < sets) && (i < nhist); k++ ) {
    
    if ( note ) {

      legge = legg->AddEntry(hs, Form("%s", fSignLeggS[sgIndex].Data()), "f");
      legge = legg->AddEntry(hb, Form("%s", fSignLeggM[bgIndex].Data()), "f");
    
    } else {

      legge = legg->AddEntry(hs, Form("%s (%4.1f)", fSignLeggS[sgIndex].Data(), hs->GetSumOfWeights()), "f");
      
      if ( hb->GetSumOfWeights() > 1000 ) {
	
	legge = legg->AddEntry(hb, Form("%s (%4.1e)",fSignLeggM[bgIndex].Data(), hb->GetSumOfWeights()), "f");
      }
      else if ( hb->GetSumOfWeights() < 0.1 ) {
	
	legge = legg->AddEntry(hb, Form("%s (%4.1e)",fSignLeggM[bgIndex].Data(), hb->GetSumOfWeights()), "f");
      }
      else if ( hb->GetSumOfWeights() ) {
	
	legge = legg->AddEntry(hb, Form("%s (%4.1f)",fSignLeggM[bgIndex].Data(), hb->GetSumOfWeights()), "f");
      }
    }
 
  
    // -- Draw plots
    hs->GetYaxis()->SetRangeUser(0.1*min, 100*max);
 
    hs->DrawCopy("hist"); 
    hb->DrawCopy("histsame");    
    
    cont = 1;
    for (i = npers*k+0; (cont < npers+1) && (i < nhist); ++i) {
   
      if ( i == bgIndex || !strcmp(fTypeM[i], "nmc") ) { continue; }
      
      if ( bg[i]->GetSumOfWeights() > 0 && bg[i]->GetMaximum() > 0 ) {
	
	cont++;

	setHist(bg[i], EColor[cont]);

	if ( note ) {
	  
	  legge = legg->AddEntry(bg[i], Form("%s", fSignLeggM[i].Data()), "p");

	} else {

	  if ( bg[i]->GetSumOfWeights() > 1000 ) {
	    legge = legg->AddEntry(bg[i], Form("%s (%4.1e)", fSignLeggM[i].Data(), bg[i]->GetSumOfWeights()), "p"); 
	  }
	  else if (bg[i]->GetSumOfWeights() < 0.1 ) {
	    legge = legg->AddEntry(bg[i], Form("%s (%4.1e)", fSignLeggM[i].Data(), bg[i]->GetSumOfWeights()), "p"); 
	  }
	  else {
	    legge = legg->AddEntry(bg[i], Form("%s (%4.1f)", fSignLeggM[i].Data(), bg[i]->GetSumOfWeights()), "p"); 
	  }
	}
	
	legge->SetTextColor(kBlack);  legge->SetFillStyle(1001);  legge->SetFillColor(10); 
		
	bg[i]->DrawCopy("same");

      }
    }
    
    legg->Draw();
    c0->SaveAs(Form("%s/bgOverlay-%i-%s.eps", outDir, k, hist));
    legg->Clear();
  }

  gPad->SetLogy(0);

}

// ----------------------------------------------------------------------------------

void anaBmm::nbgOverlay(const char *hist, const int npers) {

  int EColor[]     = {   13,    2,    4,    108,  6,  107,   93,   1, 
			 13,    2,    4,    108,  6,  107,   93,   1   }; 
  
  double min(9999.), max(0.), maxSoW(0.);
  double fnorm_Ch, eff1_Ch, eff2_Ch;  

  // -- signal
  TH1D *hs = (TH1D*)fS[normSgIndex]->Get(hist)->Clone();

  if ( !strcmp(hist,"c530") ) {

    channelEff(fS[normSgIndex], fnorm_Ch, eff1_Ch, eff2_Ch);
    hs->Scale((fLumiD[0]/fLumiS[normSgIndex])*eff1_Ch*eff2_Ch);

  } else {
    
    hs->Scale(fLumiD[0]/fLumiS[normSgIndex]);
  } 

  setFilledHist(hs, kBlack, kYellow, 1001, 2);
  setTitles(hs, "m_{#mu#mu} [GeV]", "events/bin", 0.06, 1.1, 1.3); 

  if (hs->GetMinimum(1.e-20) < min ) { min = hs->GetMinimum(1.e-20);      }
  if (hs->GetMaximum()       > max ) { max = hs->GetMaximum();            }
 
 
  // -- background
  TH1D *hb = (TH1D*)fM[normBgIndex]->Get(hist)->Clone();

  if ( !strcmp(hist,"c530") ) {

    channelEff(fM[normBgIndex], fnorm_Ch, eff1_Ch, eff2_Ch);
    hb->Scale((fLumiD[0]/fLumiM[normBgIndex])*eff1_Ch*eff2_Ch);
  
  } else {
    
    hb->Scale(fLumiD[0]/fLumiM[normBgIndex]);
  }

  hb->SetLineStyle(kDashed);
  setFilledHist(hb, 13, 13, 3005, 2);

  if (hb->GetMinimum(1.e-20) < min ) { min = hb->GetMinimum(1.e-20);      }
  if (hb->GetMaximum()       > max ) { max = hb->GetMaximum();            }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetTopMargin(0.25);
  gPad->SetRightMargin(0.2);

  gPad->SetLogy(1);

  // -- Other backgrounds
  const int nhist = nMc;
  const int sets = int(nhist/npers)+1;

  TH1D *bg[nhist];

  for (int i = 0; i < nhist; ++i) {

    bg[i] = (TH1D*)fM[i]->Get(hist)->Clone();
    
    if ( i == normBgIndex  || !strcmp(fTypeM[i], "mc") || !strcmp(fTypeM[i], "rmc") ) { continue; }

    if ( !strcmp(hist,"c530") ) {

      channelEff(fM[i], fnorm_Ch, eff1_Ch, eff2_Ch);
      bg[i]->Scale((fMisIdM[i]*fLumiD[0]/fLumiM[i])*eff1_Ch*eff2_Ch);
    
    } else {

      bg[i]->Scale(fMisIdM[i]*fLumiD[0]/fLumiM[i]);
    }

    if (bg[i]->GetSumOfWeights()  > maxSoW )  { maxSoW = bg[i]->GetSumOfWeights();  }
    if (bg[i]->GetMinimum(1.e-20) < min )     { min = bg[i]->GetMinimum(1.e-20); }
    if (bg[i]->GetMaximum()       > max )     { max = bg[i]->GetMaximum();       }
  }

  // -- Draw & Legends
  legg = new TLegend(0.49,0.62,0.99,0.99);
  legg->SetFillStyle(1001); legg->SetBorderSize(0); legg->SetTextSize(0.035);  legg->SetFillColor(10); 

  int i(0), cont (0), note(0);

  for (int k = 0; (k < sets) && (i < nhist); k++ ) {
    
    if ( note ) {

      legge = legg->AddEntry(hs, Form("%s", fSignLeggS[normSgIndex].Data()), "f");
      legge = legg->AddEntry(hb, Form("%s", fSignLeggM[normBgIndex].Data()), "f");
    
    } else {

      legge = legg->AddEntry(hs, Form("%s (%4.1f)", fSignLeggS[normSgIndex].Data(), hs->GetSumOfWeights()), "f");
      
      if ( hb->GetSumOfWeights() > 1000 ) {
	
	legge = legg->AddEntry(hb, Form("%s (%4.1e)",fSignLeggM[normBgIndex].Data(), hb->GetSumOfWeights()), "f");
      }
      else if ( hb->GetSumOfWeights() < 0.1 ) {
	
	legge = legg->AddEntry(hb, Form("%s (%4.1e)",fSignLeggM[normBgIndex].Data(), hb->GetSumOfWeights()), "f");
      }
      else if ( hb->GetSumOfWeights() ) {
	
	legge = legg->AddEntry(hb, Form("%s (%4.1f)",fSignLeggM[normBgIndex].Data(), hb->GetSumOfWeights()), "f");
      }
    }
 
  
    // -- Draw plots
    hs->GetYaxis()->SetRangeUser(0.1*min, 100*max);
 
    hs->DrawCopy("hist"); 
    hb->DrawCopy("histsame");    
    
    cont = 1;
    for (i = npers*k+0; (cont < npers+1) && (i < nhist); ++i) {
   
      if ( i == normBgIndex || !strcmp(fTypeM[i], "mc") || !strcmp(fTypeM[i], "rmc") ) { continue; }
      
      if ( bg[i]->GetSumOfWeights() > 0 && bg[i]->GetMaximum() > 0 ) {
	
	cont++;

	setHist(bg[i], EColor[cont]);

	if ( note ) {
	  
	  legge = legg->AddEntry(bg[i], Form("%s", fSignLeggM[i].Data()), "p");

	} else {

	  if ( bg[i]->GetSumOfWeights() > 1000 ) {
	    legge = legg->AddEntry(bg[i], Form("%s (%4.1e)", fSignLeggM[i].Data(), bg[i]->GetSumOfWeights()), "p"); 
	  }
	  else if (bg[i]->GetSumOfWeights() < 0.1 ) {
	    legge = legg->AddEntry(bg[i], Form("%s (%4.1e)", fSignLeggM[i].Data(), bg[i]->GetSumOfWeights()), "p"); 
	  }
	  else {
	    legge = legg->AddEntry(bg[i], Form("%s (%4.1f)", fSignLeggM[i].Data(), bg[i]->GetSumOfWeights()), "p"); 
	  }
	}
	
	legge->SetTextColor(kBlack);  legge->SetFillStyle(1001);  legge->SetFillColor(10); 
		
	bg[i]->DrawCopy("same");

      }
    }
    
    legg->Draw();
    c0->SaveAs(Form("%s/nbgOverlay-%i-%s.eps", outDir, k, hist));
    legg->Clear();
  }

  gPad->SetLogy(0);

}


//===========================================================================================
// -- Adding histogram - SPECIAL!
//===========================================================================================

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

  TH1D *h2 = (TH1D*)fM[bgIndex]->Get(hname); 
  if (h2) {
    h1 = new TH1D(*h2);
    h1->SetName(Form("MC:%s", hname)); 
    h1->Reset(); 
  }
    
  const int nhist = nMc;

  for (int i = 1; i < nhist - 1; ++i) {  // FIX ME: last histogram is screwd up ???????

    // -- Add (weighted) histogram
    int accept = checkIf(i, selection);

    if ( accept ) {

      if (fM[i]) {
	
	h2 = (TH1D*)fM[i]->Get(hname);  
	
	if (mode == 0) {
	  if (TMath::IsNaN(h2->GetSumOfWeights())) {
	    h2->Reset();
	    cout << "anaBmm::sumHistMC> ***** Problems with histogram " << hname 
		 << " from file " << fM[i]->GetName() << endl;
	    
	  } else {
	    
	    h1->Add(h2, fLumiD[0]/fLumiM[i]);
	  }
	}
	
	if (mode == 1) {

	}

	if (mode == 2) {
	  
	  if (TMath::IsNaN(h2->GetSumOfWeights())) {
	    h2->Reset();
	    cout << "anaBmm::sumHistMC> ***** Problems with histogram " << hname 
		 << " from file " << fM[i]->GetName() << endl;
	    
	  } else {
	    // cout << " summing " << fSignM[i] << " to " << selection << endl;
	    sumHistMC_Add(i, hname, h1, h2);
	    //h1->Add(h2,  (fMisIdM[i]*fLumiD[0]/fLumiM[i]));
	    
	  }
	}     
      }
    } 
  }

  cout << endl << endl;
  return h1; 

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
	
      } else if ( m == hist2->FindBin(426.1) ) {

	hist1->AddBinContent(m, SF_Ch*fnorm_Ch*eff2_Ch);
      
      } else if ( m == hist2->FindBin(428.1) ) {

	hist1->AddBinContent(m, SF_Ch*fnorm_Ch*eff1_Ch*eff2_Ch);
      } 
    }
  }
}

// -----------------------------------------------------
void anaBmm::channelEff(TFile *f, double &fnorm_Ch, double &eff1_Ch, double &eff2_Ch) {
  
  TH1D *tAR1 = (TH1D*)f->Get("AR1")->Clone();
  fnorm_Ch = 0.; eff1_Ch = 0.; eff2_Ch = 0.;

  fnorm_Ch = tAR1->GetBinContent(tAR1->FindBin(123.1));
  
  if ( tAR1->GetBinContent(tAR1->FindBin(220.1)) ) {
    
    eff1_Ch  = tAR1->GetBinContent(tAR1->FindBin(224.1))/
              (tAR1->GetBinContent(tAR1->FindBin(220.1)));
    
    eff2_Ch  = tAR1->GetBinContent(tAR1->FindBin(225.1))/
              (tAR1->GetBinContent(tAR1->FindBin(220.1)));
  }
}


//===========================================================================================
// -- Muon mis-identification
//===========================================================================================

void anaBmm::fakeMuons(const char *prod) {

  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.06);
  tl->SetTextColor(kBlack);
  
  int i = -1;
  char f_n[200];

  //  const int nhist = 8;
  double mis(0.), eff(0.), tot(0.);
  
  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% -- Muon identification "<< endl;
  
  // --- Mis-ID. histograms for Pions, Kaons, Protons & Muons
 
  TH1D *h;

  TH1D *Pi_pT, *fakePi_pT, *Pi_eta, *fakePi_eta;
  TH1D *K_pT, *fakeK_pT, *K_eta, *fakeK_eta;
  TH1D *P_pT, *fakeP_pT, *P_eta, *fakeP_eta;
  TH1D *Mu_pT, *effMu_pT, *Mu_eta, *effMu_eta;

  TH2D *Pi, *fakePi, *K, *fakeK, *P, *fakeP, *Mu, *effMu;

  // -- How to do this more efficient ???
  int nsg(0), nbg(0);
  TString sgfiles[20], bgfiles[20];
  // -- CSA07
  if ( !strcmp(prod,"csa07") ) {

    nsg = 2;
    sgfiles[0] = TString("csg-003");
    sgfiles[1] = TString("csg-004");

    nbg = 5;
    bgfiles[0] = TString("cbg-003");
    bgfiles[1] = TString("cbg-004");
    bgfiles[2] = TString("cbg-005");
    bgfiles[3] = TString("cbg-006");
    bgfiles[4] = TString("cbg-007");
  } 
  // -- Spring07
  else if ( !strcmp(prod,"spring07") ) {
    

    nsg = 2;
    sgfiles[0] = TString("csg-003h");
    sgfiles[1] = TString("csg-004h");

    nbg = 3;
    bgfiles[0] = TString("cbg-003h");
    bgfiles[1] = TString("cbg-004h");
    bgfiles[2] = TString("cbg-005h");
  }

  // -- Spring07
  else if ( !strcmp(prod,"orca") ) {

    nsg = 2;
    sgfiles[0] = TString("csg-001");
    sgfiles[1] = TString("csg-002");

    nbg = 1;
    bgfiles[0] = TString("cbg-002");
    
  } else {

    cout << " ====> Error: No production found, for " << prod << endl;
    return;
  }

  // -- from signal
  sprintf(f_n,"%s",sgfiles[0].Data());
  int index = findIndex(f_n);

  if ( index < 0 ) {
    cout << " ====> Error: couldn't find index of file " << f_n << endl;
    return;
  }

  cout << " Muon efficiencies / fake rates from ---> " <<  fFileS[index].Data() 
       << " (" << fSignS[index].Data() << ")" << endl;

  h = (TH1D*)fS[index]->Get("MisID");
  
  Pi_pT = (TH1D*)fS[index]->Get("m100");    
  fakePi_pT = (TH1D*)fS[index]->Get("m101");      
  Pi_eta = (TH1D*)fS[index]->Get("m110");
  fakePi_eta = (TH1D*)fS[index]->Get("m111");
  
  K_pT = (TH1D*)fS[index]->Get("m200");    
  fakeK_pT = (TH1D*)fS[index]->Get("m201");      
  K_eta = (TH1D*)fS[index]->Get("m210");
  fakeK_eta = (TH1D*)fS[index]->Get("m211");
  
  P_pT = (TH1D*)fS[index]->Get("m300");    
  fakeP_pT = (TH1D*)fS[index]->Get("m301");      
  P_eta = (TH1D*)fS[index]->Get("m310");
  fakeP_eta = (TH1D*)fS[index]->Get("m311");
  
  Mu_pT = (TH1D*)fS[index]->Get("m400");    
  effMu_pT = (TH1D*)fS[index]->Get("m401");      
  Mu_eta = (TH1D*)fS[index]->Get("m410");
  effMu_eta = (TH1D*)fS[index]->Get("m411");
  
  Pi = (TH2D*)fS[index]->Get("M100");
  fakePi = (TH2D*)fS[index]->Get("M101"); 
  K = (TH2D*)fS[index]->Get("M200");
  fakeK = (TH2D*)fS[index]->Get("M201");
  P = (TH2D*)fS[index]->Get("M300");
  fakeP = (TH2D*)fS[index]->Get("M301");
  Mu = (TH2D*)fS[index]->Get("M400");
  effMu = (TH2D*)fS[index]->Get("M401");  

  for ( int j = 1; j < nsg; j++ ) {

    sprintf(f_n,"%s",sgfiles[j].Data());
    i = findIndex(f_n);

    if ( i < 0 ) {
      cout << " ====> Error: couldn't find index of file " << f_n << endl;
      continue;
    }

    cout << " Muon efficiencies / fake rates from ---> " <<  fFileS[i].Data() 
       << " (" << fSignS[i].Data() << ")" << endl;
 
    h->Add((TH1D*)fS[i]->Get("MisID"));

    Pi_pT->Add((TH1D*)fS[i]->Get("m100"));    
    fakePi_pT->Add((TH1D*)fS[i]->Get("m101"));      
    Pi_eta->Add((TH1D*)fS[i]->Get("m110"));
    fakePi_eta->Add((TH1D*)fS[i]->Get("m111"));
    
    K_pT->Add((TH1D*)fS[i]->Get("m200"));
    fakeK_pT->Add((TH1D*)fS[i]->Get("m201"));
    K_eta->Add((TH1D*)fS[i]->Get("m210"));
    fakeK_eta->Add((TH1D*)fS[i]->Get("m211"));
    
    P_pT->Add((TH1D*)fS[i]->Get("m300"));    
    fakeP_pT->Add((TH1D*)fS[i]->Get("m301"));      
    P_eta->Add((TH1D*)fS[i]->Get("m310"));
    fakeP_eta->Add((TH1D*)fS[i]->Get("m311"));
    
    Mu_pT->Add((TH1D*)fS[i]->Get("m400"));    
    effMu_pT->Add((TH1D*)fS[i]->Get("m401"));      
    Mu_eta->Add((TH1D*)fS[i]->Get("m410"));
    effMu_eta->Add((TH1D*)fS[i]->Get("m411"));

    Pi->Add((TH2D*)fS[i]->Get("M100"));
    fakePi->Add((TH2D*)fS[i]->Get("M101"));
    K->Add((TH2D*)fS[i]->Get("M200"));
    fakeK->Add((TH2D*)fS[i]->Get("M201"));
    P->Add((TH2D*)fS[i]->Get("M300"));
    fakeP->Add((TH2D*)fS[i]->Get("M301"));
    Mu->Add((TH2D*)fS[i]->Get("M400"));
    effMu->Add((TH2D*)fS[i]->Get("M401"));

  }


  // -- from background
  for ( int j = 0; j < nbg; j++ ) {

    sprintf(f_n,"%s",bgfiles[j].Data());
    i = findIndex(f_n);

    if ( i < 0 ) {
      cout << " ====> Error: couldn't find index of file " << f_n << endl;
      continue;
    }

    cout << " Muon efficiencies / fake rates from ---> " <<  fFileM[i].Data() 
       << " (" << fSignM[i].Data() << ")" << endl;

    h->Add((TH1D*)fM[i]->Get("MisID"));

    Pi_pT->Add((TH1D*)fM[i]->Get("m100"));    
    fakePi_pT->Add((TH1D*)fM[i]->Get("m101"));      
    Pi_eta->Add((TH1D*)fM[i]->Get("m110"));
    fakePi_eta->Add((TH1D*)fM[i]->Get("m111"));
    
    K_pT->Add((TH1D*)fM[i]->Get("m200"));
    fakeK_pT->Add((TH1D*)fM[i]->Get("m201"));
    K_eta->Add((TH1D*)fM[i]->Get("m210"));
    fakeK_eta->Add((TH1D*)fM[i]->Get("m211"));
    
    P_pT->Add((TH1D*)fM[i]->Get("m300"));    
    fakeP_pT->Add((TH1D*)fM[i]->Get("m301"));      
    P_eta->Add((TH1D*)fM[i]->Get("m310"));
    fakeP_eta->Add((TH1D*)fM[i]->Get("m311"));
    
    Mu_pT->Add((TH1D*)fM[i]->Get("m400"));    
    effMu_pT->Add((TH1D*)fM[i]->Get("m401"));      
    Mu_eta->Add((TH1D*)fM[i]->Get("m410"));
    effMu_eta->Add((TH1D*)fM[i]->Get("m411"));

    Pi->Add((TH2D*)fM[i]->Get("M100"));
    fakePi->Add((TH2D*)fM[i]->Get("M101"));
    K->Add((TH2D*)fM[i]->Get("M200"));
    fakeK->Add((TH2D*)fM[i]->Get("M201"));
    P->Add((TH2D*)fM[i]->Get("M300"));
    fakeP->Add((TH2D*)fM[i]->Get("M301"));
    Mu->Add((TH2D*)fM[i]->Get("M400"));
    effMu->Add((TH2D*)fM[i]->Get("M401"));
      
  }

  // ... fake muons from pions
  tot = 0.;  mis = 0.; double frPions(0.);
  tot = h->GetBinContent(h->FindBin(0.1));
  mis = h->GetBinContent(h->FindBin(1.1));
  if ( tot != 0 ) { frPions = mis/tot; }

  // ... fake muons from kaons
  tot = 0.;  mis = 0.; double frKaons(0.);
  tot = h->GetBinContent(h->FindBin(10.1));
  mis = h->GetBinContent(h->FindBin(11.1));
  if ( tot != 0 ) { frKaons = mis/tot; }

  // ... fake muons from protons
  tot = 0.;  mis = 0.; double frProtons(0.);
  tot = h->GetBinContent(h->FindBin(20.1));
  mis = h->GetBinContent(h->FindBin(21.1));
  if ( tot != 0 ) { frProtons = mis/tot; }

  // ... not identified Muons
  tot = 0.;  eff = 0.; double effMuons(0.);
  tot = h->GetBinContent(h->FindBin(30.1));
  mis = h->GetBinContent(h->FindBin(31.1));
  if ( tot != 0 ) { effMuons = eff/tot; }


  // --- Draw combined histograms

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetLogy(0); 
  shrinkPad(0.15, 0.15);

  // -- 1D histograms
  // - pT
  fPiMid = DivideHisto(fakePi_pT, Pi_pT);
  fPiMid->GetYaxis()->SetRangeUser(0,0.021);
  fPiMid->GetYaxis()->SetNdivisions(6, kTRUE);
  setTitles(fPiMid, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  fPiMid->Draw();
  c0->SaveAs(Form("%s/muonID/%s/fakerate_pionsPt.eps", outDir, prod));

  // - eta
  TH1D* MisID = DivideHisto(fakePi_eta, Pi_eta);
  MisID->GetYaxis()->SetRangeUser(0,0.021);
  setTitles(MisID, "#eta", "#varepsilon", 0.06, 1.1, 1.1);
  MisID->Draw();
  c0->SaveAs(Form("%s/muonID/%s/fakerate_pionsEta.eps",  outDir, prod));


  // -- Kaons
  // - pT
  fKaMid = DivideHisto(fakeK_pT, K_pT);
  fKaMid->GetYaxis()->SetRangeUser(0,0.021); 
  fKaMid->GetYaxis()->SetNdivisions(6, kTRUE);
  setTitles(fKaMid, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  fKaMid->Draw();
  c0->SaveAs(Form("%s/muonID/%s/fakerate_kaonsPt.eps",  outDir, prod));

  // - eta
  MisID = DivideHisto(fakeK_eta, K_eta);
  MisID->GetYaxis()->SetRangeUser(0,0.021);
  setTitles(MisID, "#eta", "#varepsilon", 0.06, 1.1, 1.1);
  MisID->Draw();
  c0->SaveAs(Form("%s/muonID/%s/fakerate_kaonsEta.eps",  outDir, prod));


  // -- Protons
  // - pT
  fProtMid = DivideHisto(fakeP_pT, P_pT);
  fProtMid->GetYaxis()->SetRangeUser(0,0.021);
  fProtMid->GetYaxis()->SetNdivisions(6, kTRUE);
  setTitles(fProtMid, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  fProtMid->Draw();
  c0->SaveAs(Form("%s/muonID/%s/fakerate_protonsPt.eps",  outDir, prod));

  // - eta
  MisID = DivideHisto(fakeP_eta, P_eta);
  MisID->GetYaxis()->SetRangeUser(0,0.021);
  setTitles(MisID, "#eta", "#varepsilon", 0.06, 1.1, 1.1);
  MisID->Draw();
  c0->SaveAs(Form("%s/muonID/%s/fakerate_protonsEta.eps",  outDir, prod));


  // -- Muons
  // - pT
  fMuEff = DivideHisto(effMu_pT, Mu_pT);
  fMuEff->GetYaxis()->SetRangeUser(0,1.1);
  setTitles(fMuEff, "p_{T}", "#varepsilon", 0.06, 1.1, 1.1);
  //  fMuEff->GetYaxis()->SetNdivisions(6, kTRUE);
  fMuEff->Draw();
  c0->SaveAs(Form("%s/muonID/%s/efficiency_muonsPt.eps",  outDir, prod));

  // - eta
  MisID = DivideHisto(effMu_eta, Mu_eta);
  MisID->GetYaxis()->SetRangeUser(0,1.1);
  setTitles(MisID, "#eta", "#varepsilon", 0.06, 1.1, 1.1);
  MisID->Draw();
  c0->SaveAs(Form("%s/muonID/%s/efficiency_muonsEta.eps",  outDir, prod));


  // -- 2D histograms

  fPiMid2 = DivideHisto2(fakePi, Pi);
  setTitles2(fPiMid2, "p_{T}", "#eta", 0.06, 1.1, 1.1);
  fPiMid2->Draw("colz");
  c0->SaveAs(Form("%s/muonID/%s/fakerate_pionsEtaPt.eps",  outDir, prod));

  fKaMid2 = DivideHisto2(fakeK, K);
  setTitles2(fKaMid2, "p_{T}", "#eta", 0.06, 1.1, 1.1);
  fKaMid2->Draw("colz");
  c0->SaveAs(Form("%s/muonID/%s/fakerate_kaonsEtaPt.eps",  outDir, prod));

  fProtMid2 = DivideHisto2(fakeP, P);
  setTitles2(fProtMid2, "p_{T}", "#eta", 0.06, 1.1, 1.1);
  fProtMid2->Draw("colz");
  c0->SaveAs(Form("%s/muonID/%s/fakerate_protonsEtaPt.eps",  outDir, prod));

  fMuEff2 = DivideHisto2(effMu, Mu);
  setTitles2(fMuEff2, "p_{T}", "#eta", 0.06, 1.1, 1.1); 
  fMuEff2->GetZaxis()->SetRangeUser(0.,1.);
  fMuEff2->GetZaxis()->SetNdivisions(6, kTRUE);
  fMuEff2->Draw("colz");
  c0->SaveAs(Form("%s/muonID/%s/efficiency_muonsEtaPt.eps",  outDir, prod));
  
  OUT << Form("\\vdef{effIdMuonsSim}    {\\ensuremath{ {%4.3E } } }", effMuons)   << endl;
  OUT << Form("\\vdef{misIdPionsSim}    {\\ensuremath{ {%4.3E } } }", frPions)   << endl;
  OUT << Form("\\vdef{misIdKaonsSim}    {\\ensuremath{ {%4.3E } } }", frKaons)   << endl;
  OUT << Form("\\vdef{misIdProtonsSim}  {\\ensuremath{ {%4.3E } } }", frProtons) << endl;
 
  OUT << Form("\\vdef{effIdMuons}    {\\ensuremath{ {%4.3E } } }", fMu)   << endl;
  OUT << Form("\\vdef{misIdPions}    {\\ensuremath{ {%4.3E } } }", fPi)   << endl;
  OUT << Form("\\vdef{misIdKaons}    {\\ensuremath{ {%4.3E } } }", fKa)   << endl;
  OUT << Form("\\vdef{misIdProtons}  {\\ensuremath{ {%4.3E } } }", fProt) << endl;

  OUT.close();
}


//===========================================================================================
// -- Dump stuff
//===========================================================================================

void anaBmm::dumpFiles() {

  cout << endl << "Assuming data luminosity " << Form("%2.1f", fLumiD[0]) << " /fb" << endl << endl;

  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% -- X-sections, NEVT, and lumi" << endl;
  OUT << Form("%s", (formatTex(fLumiD[0], "Lumi" , "d0")).Data()) << endl;

  for (int i = 0; i < 10; ++i) {

    if (fS[i]) {

      fNexpS[i] = fLumiD[0]*fvXsS[i];

       cout << Form("Signal MC[%d]     visible cross-section: %3.2e fb, \
       Nevt: %8.0f -> Lumi: %3.2e /fb   Signature: %s Title: %s File: %s"
		    , i, fvXsS[i], fNevtS[i], fLumiS[i], fSignS[i].Data(), fSignTexS[i].Data(), fFileS[i].Data()) 
	    << endl;

      OUT << Form("\\vdef{channel:s%i} {\\ensuremath{ {%s } } }", i, fSignTexS[i].Data()) << endl;
      OUT << Form("\\vdef{vNevt:s%i} {\\ensuremath{ {%4.0f } } }", i, fNevtS[i]) << endl;
      OUT << Form("\\vdef{vLumi:s%i} {\\ensuremath{ {%s } } }", i, (texForm(fLumiS[i])).Data()) << endl;

      if ( i == sgIndex ) {
	OUT << Form("\\vdef{vXs:s%i} {\\ensuremath{ {%4.1f } } }", i, fvXsS[i]) << endl;
	OUT << Form("\\vdef{vNexp:s%i} {\\ensuremath{ {%4.1f } } }", i, fNexpS[i]) << endl;
      } else {
	OUT << Form("\\vdef{vXs:s%i} {\\ensuremath{ {%s } } }", i, (texForm(fvXsS[i])).Data()) << endl;
	//OUT << Form("%s", formatTex(fNexpS[i], "vNexp", Form("s%i", i) )) << endl;
	OUT << Form("\\vdef{vNexp:s%i} {\\ensuremath{ {%s } } }", i, (texForm(fNexpS[i])).Data()) << endl;
      }
    }
  }
 
  for (int i = 0; i < 30; ++i) {

    if (fM[i]) {

      fMisIdM[i] = getMisID(fSignM[i]);
      fNexpM[i] = fLumiD[0]*fvXsM[i]*fMisIdM[i];


      cout << Form("Background MC[%d] visible cross-section: %3.2e fb, \
      Nevt: %8.0f -> lumi: %3.2e /fb   Signature: %s Title: %s File: %s"
		   , i, fvXsM[i], fNevtM[i], fLumiM[i], fSignM[i].Data(), fSignTexM[i].Data(), fFileM[i].Data()) 
	   << endl;

      OUT << Form("\\vdef{channel:m%i} {\\ensuremath{ {%s } } }", i, fSignTexM[i].Data()) << endl;

      OUT << Form("\\vdef{channel:m%i} {\\ensuremath{ {%s } } }", i, fSignTexM[i].Data()) << endl;
      OUT << Form("\\vdef{vNevt:m%i} {\\ensuremath{ {%4.0f } } }", i, fNevtM[i]) << endl;
      OUT << Form("\\vdef{vLumi:m%i} {\\ensuremath{ {%s } } }", i, (texForm(fLumiM[i])).Data()) << endl;


      OUT << Form("\\vdef{vXs:m%i} {\\ensuremath{ {%s } } }", i, (texForm(fvXsM[i])).Data()) << endl;
      OUT << Form("%s", formatTex(fNexpM[i], "vNexp", Form("m%i", i)).Data()) << endl;
      //OUT << Form("\\vdef{vNexp:m%i} {\\ensuremath{ {%s } } }", i, (texForm(fNexpM[i])).Data()) << endl;
    }
  }
  
  OUT.close();
}


// ----------------------------------------------------------------------
void anaBmm::dumpCuts() {

  TH1D *hS = (TH1D*)fS[sgIndex]->Get("hcuts");
  TH1D *hM = (TH1D*)fM[bgIndex]->Get("hcuts");

  ofstream OUT(fNumbersFileName, ios::app);
  OUT << "% ----------------------------------------------------------------------" << endl;
  OUT << "% -- Cuts" << endl;

  char label[100];
  sprintf(label, "%s", hS->GetXaxis()->GetBinLabel(1));

  double sVal, mVal;
  int show = 0, problem = 0;

  // cout << " -- Cuts signal (mc)"  << endl;

  for (int i = 1; i < hS->GetNbinsX(); ++i) {
    if (strcmp(hS->GetXaxis()->GetBinLabel(i), "")) {

      sVal = hS->GetBinContent(i);
      mVal = hM->GetBinContent(i);

      // cout << hS->GetXaxis()->GetBinLabel(i) << " =  " << sVal << " (" <<  mVal  << ")" << endl;

      if (sVal != mVal) {
	cout << endl << "  Cut bin " << i << ": Signal  =  " << sVal << ", BG = " << mVal << " for " ;
	show = 1; problem = 1;
      } else {
	show = 0;
      }
	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "p_{T}(B_{s}) [GeV]")) {
	  OUT  << Form("\\vdef{cut:%s:ptbs}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  ptbs = sVal;
	  optimizedCut[4] = TString(Form("pt>%5.4f", sVal));
	  if ( show ) { cout << "pT(Bs) " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "p_{T}^{min}(l) [GeV]")) {
	  OUT  << Form("\\vdef{cut:%s:ptlo}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  ptmulo = sVal;
	  optimizedCut[0] = TString(Form("ptl1>%5.4f", sVal));
	  if ( show ) { cout << "pT_min " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "p_{T}^{max}(l) [GeV]")) {
	  OUT  << Form("\\vdef{cut:%s:pthi}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  if ( show ) { cout << "pT_max " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "R_{#mu#mu}^{min}")) {
	  OUT  << Form("\\vdef{cut:%s:rmmlo}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  rmmlo = sVal;
	  optimizedCut[2] = TString(Form("rmm>%5.4f", sVal));
	  if ( show ) { cout << "Rmm_min " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "R_{#mu#mu}^{max}")) {
	  OUT  << Form("\\vdef{cut:%s:rmmhi}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  rmmhi = sVal;
	  optimizedCut[3] = TString(Form("rmm<%5.4f", sVal));
	  if ( show ) { cout << "Rmm_max " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "#eta_{B}^{min}")) {
	  OUT  << Form("\\vdef{cut:%s:etalo}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  etalo = sVal;
	  optimizedCut[1] = TString(Form("eta>%5.4f", sVal));
	  optimizedCut[4] = TString(Form("etal0>%5.4f&&etal1<%5.4f", sVal, sVal));
	  if ( show ) { cout << "eta(B)_min " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "#eta_{B}^{max}")) {
	  OUT  << Form("\\vdef{cut:%s:etahi}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  etahi = sVal;
	  if ( show ) { cout << "eta(B)_max " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "TIP(l) [cm]")) {
	  OUT  << Form("\\vdef{cut:%s:tip}      {\\ensuremath{%5.3f } } ", label, sVal) << endl;
	  if ( show ) { cout << "TIP " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "#chi^2")) {
	  OUT  << Form("\\vdef{cut:%s:chi2}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  vtxhi = sVal;
	  optimizedCut[11] = TString(Form("chi2<%5.4f", sVal));
	  if ( show ) { cout << "chi2 " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "l_{3d} [cm]")) {
	  OUT  << Form("\\vdef{cut:%s:l3d}    {\\ensuremath{%5.3f } } ", label, sVal) << endl;
	  if ( show ) { cout << "L3D " << sVal << endl; }
	}


	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "cos(#alpha)")) {
	  OUT  << Form("\\vdef{cut:%s:cosalpha}    {\\ensuremath{%5.4f } } ", label, sVal) << endl;
	  coslo = sVal;
	  optimizedCut[8] = TString(Form("cosa>%5.4f", sVal));
	  optimizedCut[6] = TString(Form("cosa3>%5.4f", sVal));
	  if ( show ) { cout << "cos(alpha) " << sVal << endl; }
	}


	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "cos(#alpha3)")) {
	  OUT  << Form("\\vdef{cut:%s:cosalpha3}    {\\ensuremath{%5.4f } } ", label, sVal) << endl;
	  cos3lo = sVal;
// 	  optimizedCut[6] = TString(Form("cosa3>%5.4f", sVal));
	  if ( show ) { cout << "cos(alpha3) " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "l_{xy} [cm]")) {
	  OUT  << Form("\\vdef{cut:%s:lxy}    {\\ensuremath{%5.3f } } ", label, sVal) << endl;
	  if ( show ) { cout << "LXY " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "l_{3D}/#sigma_{3D}")) {
	  OUT  << Form("\\vdef{cut:%s:l3d/s3d}    {\\ensuremath{%3.1f } } ", label, sVal) << endl;
	  l3dlo = sVal;
	  optimizedCut[7] = TString(Form("l3d/s3d>%5.4f", sVal));
	  if ( show ) { cout << "L3D/S3D " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "l_{xy}/#sigma_{xy}")) {
	  OUT  << Form("\\vdef{cut:%s:lxy/sxy}    {\\ensuremath{%3.1f } } ", label, sVal) << endl;
	  lxylo = sVal;
	  optimizedCut[9] = TString(Form("lxy/sxy>%5.4f", sVal));
	  if ( show ) { cout << "LXY/SXY " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "I_{veto}")) {
	  OUT  << Form("\\vdef{cut:%s:isoveto}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  if ( show ) { cout << "I_veto " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "R_{I}")) {
	  OUT  << Form("\\vdef{cut:%s:isocone}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  if ( show ) { cout << "IsoCone " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "I")) {
	  OUT  << Form("\\vdef{cut:%s:isolation}    {\\ensuremath{%5.3f } } ", label, sVal) << endl;
	  isolo = sVal;
	  optimizedCut[10] = TString(Form("iso>%5.4f", sVal));
	  if ( show ) { cout << "Isolation " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "BMMSEL")) {
	  OUT  << Form("\\vdef{cut:%s:bmmsel}    {\\ensuremath{%i } } ", label, sVal) << endl;
	  if ( show ) { cout << "BMMSEL " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "SUBSEL")) {
	  OUT  << Form("\\vdef{cut:%s:subsel}    {\\ensuremath{%i } } ", label, sVal) << endl;
	  if ( show ) { cout << "SUBSEL " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "m_{min}^{#mu #mu}")) {
	  OUT  << Form("\\vdef{cut:%s:masslo}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  masslo = sVal;
	  if ( show ) { cout << "mass_min " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "m_{max}^{#mu #mu}")) {
	  OUT  << Form("\\vdef{cut:%s:masshi}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  masshi = sVal;
	  if ( show ) { cout << "mass_max " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "iso. p_{T}^{min}")) {
	  OUT  << Form("\\vdef{cut:%s:isoptmin}    {\\ensuremath{%5.1f } } ", label, sVal) << endl;
	  if ( show ) { cout << "pT_iso " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "L1 available")) {
	  OUT  << Form("\\vdef{cut:%s:setl1}    {\\ensuremath{%i } } ", label, sVal) << endl;
	  if ( show ) { cout << "L1 " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "HLT available")) {
	  OUT  << Form("\\vdef{cut:%s:sethlt}    {\\ensuremath{%i } } ", label, sVal) << endl;
	  if ( show ) { cout << "HLT " << sVal << endl; }
	}

	if (!strcmp(hS->GetXaxis()->GetBinLabel(i), "#Delta m{#mu #mu}")) {
	  OUT  << Form("\\vdef{cut:%s:masswi}    {\\ensuremath{%5.0f } } ", label, sVal) << endl;
	  masswi = 0.001*sVal;
	  if ( show ) { cout << "deltaMass " << sVal << endl; }
	}
	//      }
    }
  }

  if ( problem ) cout << "====> Error: Signal and BG MC run with different cut values !!!" << endl;
  OUT.close();
}

// ----------------------------------------------------------------------
void anaBmm::writeFitPar(TF1 *f, int o, double mean, double sigma, int npar) {

    ofstream OUT(fNumbersFileName, ios::app);

    OUT << "% ----------------------------------------------------------------------" << endl;
    OUT << "% -- BRECO mass peak numbers" << endl;
    OUT << Form("\\vdef{mBgMean:s%i} {\\ensuremath {%4.1f } }", o, 1000.*mean) << endl;
    OUT << Form("\\vdef{mBgSigma:s%i}{\\ensuremath {%4.1f } }", o, 1000.*sigma) << endl;
    
    OUT << Form("\\vdef{mBg1n:s%i} {\\ensuremath {%4.2f } }", o, f->GetParameter(0)) << endl;
    OUT << Form("\\vdef{mBg1nE:s%i}{\\ensuremath {%4.2f } }", o, f->GetParError(0)) << endl;
    
    OUT << Form("\\vdef{mBg1m:s%i}  {\\ensuremath {%4.1f } }", o, 1000.*f->GetParameter(1)) << endl;
    OUT << Form("\\vdef{mBg1mE:s%i} {\\ensuremath {%4.1f } }", o, 1000.*f->GetParError(1)) << endl;
    
    OUT << Form("\\vdef{mBg1s:s%i}  {\\ensuremath {%4.1f } }", o, 1000.*f->GetParameter(2)) << endl;
    OUT << Form("\\vdef{mBg1sE:s%i} {\\ensuremath {%4.1f } }", o, 1000.*f->GetParError(2)) << endl;
    

    if ( npar > 3 ) {

      OUT << Form("\\vdef{mBg2n:s%i} {\\ensuremath {%4.2f } }", o, f->GetParameter(3)) << endl;
      OUT << Form("\\vdef{mBg2nE:s%i}{\\ensuremath {%4.2f } }", o, f->GetParError(3)) << endl;
      
      OUT << Form("\\vdef{mBg2m:s%i}  {\\ensuremath {%4.1f } }", o, 1000.*f->GetParameter(4)) << endl;
      OUT << Form("\\vdef{mBg2mE:s%i} {\\ensuremath {%4.1f } }", o, 1000.*f->GetParError(4)) << endl;
      
      OUT << Form("\\vdef{mBg2s:s%i}  {\\ensuremath {%4.1f } }", o, 1000.*f->GetParameter(5)) << endl;
      OUT << Form("\\vdef{mBg2sE:s%i} {\\ensuremath {%4.1f } }", o, 1000.*f->GetParError(5)) << endl;
    }

    OUT.close();
}

// ----------------------------------------------------------------------
void anaBmm::fillInTheRest(const char *tag) { 

  ofstream OUT(fNumbersFileName, ios::app);

  OUT << Form("%s", (formatTex(0.,  "nMassAllCutsFact"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "nMassAllCutsFactE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0.,  "eMassAllCutsFact"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "eMassAllCutsFactE" , tag)).Data()) << endl;

  OUT << Form("%s", (formatTex(0.,  "nAllCutsFact"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "nAllCutsFactE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0.,  "eAllCutsFact"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "eAllCutsFactE" , tag)).Data()) << endl;
 
  OUT << Form("%s", (formatTex(0.,  "nTotEffChi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "nTotEffChi2E" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0.,  "eTotEffChi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "eTotEffChi2E" , tag)).Data()) << endl;
 
  OUT << Form("%s", (formatTex(0.,  "nTotEffIso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "nTotEffIsoE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0.,  "eTotEffIso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "eTotEffIsoE" , tag)).Data()) << endl;

  OUT << Form("%s", (formatTex(0.,  "nAfBchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "eAfBchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "eAfBchi2E" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0.,  "cAfBchi2"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0.,  "cAfBchi2E" , tag)).Data()) << endl;

  OUT << Form("%s", (formatTex(0.,  "nAfBiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "eAfBiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0., "eAfBisoE" , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0.,  "cAfBiso"  , tag)).Data()) << endl;
  OUT << Form("%s", (formatTex(0.,  "cAfBisoE" , tag)).Data()) << endl;

  OUT.close();
}


//===========================================================================================
// -- Optimization
//===========================================================================================

void anaBmm::overlay(const char *var, const char *cuts, double min, double max, int nbin) {


  fS[sgIndex]->cd(); 
  TH1D *hhSG = (TH1D*)gROOT->FindObject("hhSG"); 
  if (hhSG) delete hhSG;
  hhSG = new TH1D("hhSG", "", nbin, min, max);  

  TTree *s = (TTree*)gFile->Get("events");
  s->Draw(Form("%s>>hhSG",var), Form("%s", cuts), "goff"); 

  fM[bgIndex]->cd(); 
  TH1D *hhBG = (TH1D*)gROOT->FindObject("hhBG"); 
  if (hhBG) delete hhBG;
  hhBG = new TH1D("hhBG", "", nbin, min, max);  

  TTree *b = (TTree*)gFile->Get("events");
  b->Draw(Form("%s>>hhBG",var), Form("%s", cuts), "goff"); 
  
  c0->Clear();
  hhSG->Draw();
  hhBG->Draw("samee");
}


// -----------------------------------------------------------------------------------------------------
void anaBmm::loopHandOptimization(int tolerance) {


  TString extraTS[] = { TString("ptl1>"), TString("etamu"), TString("rmm>"), TString("rmm<"), 
			TString("pt>"), TString("etab"), 
			TString("cosa3>"),TString("l3d/s3d>"),TString("cosa>"),TString("lxy/sxy>"),
			TString("iso>"), TString("chi2<") } ;

  // arrays for:     pt  eta rmmlo rmmhi    ptB  eta  cosa3  l3d   cosa  lxy  iso chi2
  double mins[] = {  2., 0.5,   0.0,  0.0,     5., 0.5, 0.99, 0., 0.99,  0.,  0.8,  0.}; 
  double maxs[] = {  6., 2.4,  2.0,  2.0,    15., 2.5, 1.00,  40., 1.00,   40., 1.0, 10.};
  int bins[]    = { 20,  10,    20,    20,   20,  10,    20,  40,   20,   40,    20,  20};
  //int bins[]    = { 2,  2,    2,    2,   2,  2,    2,  2,   2,   2,    2,  2};

    
  c0->cd();
  c0->Clear();
  c0->Divide(4,3);

  char extra[200];
  for (int i = 0; i < NCUTS; i++ ) {
    
    c0->cd(i+1);
    sprintf(extra, "%s", extraTS[i].Data());
    runHandOptimization("4.8<mass&&mass<6.0&&goodL1&&goodHLT", extra, bins[i], mins[i], maxs[i], i, tolerance);

    optimizedCut[i] = TString(Form("%s %5.4f", extra, optCutValue[i]));
  }
  

  c0->SaveAs("handOptimization.pdf");

  loopOptimizedCuts2D();
  loopOptimizedCuts3D();

}

// -----------------------------------------------------------------------------------------------------
void anaBmm::runHandOptimization(const char *aoCuts, const char *extraVar, int nbin, double min, 
				 double max, int cut_nr, int tolerance) {

  double cut, maxFom(-99.), maxCut(0.);
  double lowCut(0.) , upCut(0.);
  double lowFom(0.) , upFom(0.);
  int fbin(0), iseta(0), isetaB(0);

  hopt[cut_nr]    = new TH1D(Form("fom_%s", extraVar),Form("%s", extraVar), nbin, min, max);
  
  char etaCut[200];
  TString extraVarTS = TString(Form("%s", extraVar));
  if ( extraVarTS.Contains("etamu"))  {
    iseta = 1;
  }

  if ( extraVarTS.Contains("etab"))  {
    isetaB = 1;
  }

  for (int i = 0; i < nbin; ++i) {

    cut = min + i*(max-min)/nbin; 
    
    if (iseta) { 
      sprintf(etaCut, "TMath::Abs(etal0)<%5.4f && TMath::Abs(etal1)<%5.4f", cut, cut);
      handOptimization(aoCuts, etaCut);
    } else if (isetaB) { 
      sprintf(etaCut, "TMath::Abs(eta)<%5.4f", cut);
      handOptimization(aoCuts, etaCut);
    } else {
      handOptimization(aoCuts, Form("%s %5.4f", extraVar, cut));
    }
    
    hopt[cut_nr]->SetBinContent(i+1, fFom);
    
    if (fFom > maxFom) {
      
      maxFom = fFom; 
      maxCut = cut;   
      fbin = i;
    }
  }

  if (tolerance) {
    for (int i = fbin; i > -1; --i) {

      lowCut = min + i*(max-min)/nbin; 
        
      if (iseta) { 
	sprintf(etaCut, "TMath::Abs(etal0)<%5.4f && TMath::Abs(etal1)<%5.4f", lowCut, lowCut);
	handOptimization(aoCuts, etaCut, 0);
      } else if (isetaB) { 
	sprintf(etaCut, "TMath::Abs(eta)<%5.4f", lowCut);
	handOptimization(aoCuts, etaCut, 0);
      } else {
	handOptimization(aoCuts, Form("%s %5.4f", extraVar, lowCut), 0);
      }

      if ( TMath::Abs(fFom - maxFom) > 0.01*maxFom && i > -1 ) {
	lowCut = min + (i+1)*(max-min)/nbin; 
	break;
      }

      lowFom = fFom;
    }

    for (int i = fbin; i < nbin; ++i) {

      upCut = min + i*(max-min)/nbin; 

      if (iseta) { 
	sprintf(etaCut, "TMath::Abs(etal0)<%5.4f && TMath::Abs(etal1)<%5.4f", upCut, upCut);
	handOptimization(aoCuts, etaCut, 0);
      } else if (isetaB) { 
	sprintf(etaCut, "TMath::Abs(eta)<%5.4f", upCut);
	handOptimization(aoCuts, etaCut, 0);
      } else {
	handOptimization(aoCuts, Form("%s %5.4f", extraVar, upCut), 0);
      }

      if ( TMath::Abs(fFom - maxFom) > 0.01*maxFom && i > -1 ) {
	upCut = min + (i-1)*(max-min)/nbin; 
	break;
      }

      upFom = fFom;
    }
    
    cout << "==> maximum fom: " << Form("%4.4e",maxFom) << " at cut " << extraVar << " " << maxCut 
	 << "  (" << lowCut << " - " << upCut <<" within 1%)" << endl;
  } else {
    
    cout << "==> maximum fom: " << Form("%4.4e",maxFom) << " at cut " << extraVar << " " << maxCut << endl;
  }

  hopt[cut_nr]->GetYaxis()->SetRangeUser(-0.1, 1.2*hopt[cut_nr]->GetMaximum());
  hopt[cut_nr]->Draw();
  pl->SetLineColor(2);
  pl->SetLineWidth(1);
  pl->DrawLine(maxCut, -0.1, maxCut, hopt[cut_nr]->GetMaximum());
  tl->SetTextColor(kBlack);  
  tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
  tl->DrawLatex(0.2, 0.8, Form("opt = %4.4f", maxCut)); 
  
  optCutValue[cut_nr] = maxCut;

  char filename[200];
  sprintf(filename, "handOptimization.txt"); 
  
  //   system(Form("/bin/mv %s old.%s", filename, filename));
  ofstream OUT(filename, ios::app);

  if ( tolerance ) {
    OUT << "==> maximum fom: " << Form("%4.4e",maxFom) << " at cut " << extraVar << " " << maxCut 
	<< "  (" << lowCut << " - " << upCut <<" within 1%)" << endl;
  }  else {
    OUT << "==> maximum fom: " << Form("%4.4e",maxFom) << " at cut " << extraVar << " " << maxCut << endl;
 
  }
  OUT.close();
}

// ----------------------------------------------------------------------
void anaBmm::handOptimization(const char *aoCuts, const char *extraCuts, int verbose) {
  
  // -- Fix lumi normalisation scaling factor
  double SF = fLumiD[0]/fLumiM[bgIndex];

  // -- Run on signal MC to determine efficiency 
  fS[sgIndex]->cd(); 
  TH1D *hSG = (TH1D*)gROOT->FindObject("hSG"); 
  if (!hSG) hSG = new TH1D("hSG", "", 50, 0., 10.);
  TTree *s = (TTree*)gFile->Get("events");
  s->Draw("mass>>hSG", Form("goodKinematics"), "goff"); 
  double sNorm = hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("goodL1&&goodHLT && %s", aoCuts), "goff"); 
  double s1Norm = hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("goodL1&&goodHLT && %s", extraCuts), "goff"); 
  double s2Norm  =  hSG->GetSumOfWeights(); //cout << Form("TEST: goodL1&&goodHLT && %s && %s", aoCuts, extraCuts) << endl;
  s->Draw("mass>>hSG", Form("goodL1&&goodHLT && %s && %s", aoCuts, extraCuts), "goff"); 
  double s12Norm  =  hSG->GetSumOfWeights();

  // -- Run on bg MC
  fM[bgIndex]->cd(); 
  TH1D *hBG = (TH1D*)gROOT->FindObject("hBG"); 
  if (!hBG) hBG = new TH1D("hBG", "", 50, 0., 10.);
  TTree *b = (TTree*)gFile->Get("events");
  b->Draw("mass>>hBG", Form("goodKinematics"), "goff"); 
  double bNorm   = hBG->GetSumOfWeights();
  b->Draw("mass>>hBG", Form("goodL1&&goodHLT && %s", aoCuts), "goff"); 
  double b1Norm  = hBG->GetSumOfWeights();
  b->Draw("mass>>hBG", Form("goodL1&&goodHLT && %s", extraCuts), "goff"); 
  double b2Norm  = hBG->GetSumOfWeights();
  b->Draw("mass>>hBG", Form("goodL1&&goodHLT && %s && %s", aoCuts, extraCuts), "goff"); 
  double b12Norm  =  hBG->GetSumOfWeights();
  
  double epsilon = 0.1;

  if (b12Norm < epsilon) b12Norm = epsilon;
  if ( b1Norm < epsilon) b1Norm = epsilon;
  if ( b2Norm < epsilon) b2Norm = epsilon;

//   fFom = (s12Norm/sNorm) / (1. + TMath::Sqrt(e2*b1Norm*SF)) ;

  fFom = (s12Norm/s1Norm)*(s12Norm/s1Norm) / (b12Norm/b1Norm);

  if ( verbose )  cout << "\"" << extraCuts << "\" : "<< Form("%4.4e", fFom) << " cuts: \"" << aoCuts << "\"  " << endl;

//     cout << Form("%4.4e signal: %4.0f/%4.0f = %4.4e, bg: %4.0f/%4.0f = %4.4e, (b1=%4.0f, b2= %4.0f)", 
// 			      fFom, s12Norm, s1Norm, s12Norm/s1Norm, b12Norm, b1Norm, b12Norm/b1Norm, b1Norm, b2Norm)
// 		      << " cuts: \"" << aoCuts << "\", \"" << extraCuts << "\" " << endl;
}

// -----------------------------------------------------------------------------------------------------
void anaBmm::loopOptimizedCuts2D() {

  char extra[200];
  char ao[200];

  int icut = 1;

  for (int i = 0; i < NCUTS; i++ ) {

    if ( optimizedCut[i].Contains("l3d/s3d") ) { continue; }
    if ( optimizedCut[i].Contains("cosa3") )   { continue; }
    if ( optimizedCut[i].Contains("eta") )     { continue; }

    sprintf(ao, "4.8<mass&&mass<6.0&&goodL1&&goodHLT");

    for (int j = 0; j < NCUTS; j++ ) {

      if ( i == j ) { continue; }

      if ( optimizedCut[j].Contains("l3d/s3d") ) { continue; }
      if ( optimizedCut[j].Contains("cosa3") )   { continue; }
      if ( optimizedCut[j].Contains("iso") )     { continue; }
      if ( optimizedCut[j].Contains("chi2") )    { continue; }
      if ( optimizedCut[j].Contains("eta") )     { continue; }

      sprintf(ao, "%s && %s", ao, optimizedCut[j].Data());
    }

    sprintf(extra, "%s", optimizedCut[i].Data());
    cutEfficiency(ao, extra, icut);
    icut++;
  }

  c0->cd();
  c0->Clear();
  c0->Divide(1,2);

  c0->cd(2);
  shrinkPad(0.3, 0.2);
  hoptSg->SetLineColor(kBlue);
  hoptSg->SetFillStyle(3001);
  hoptSg->GetYaxis()->SetRangeUser(0, 1.2);
  hoptSg->Draw();
  hoptBg->Draw("SAME");

  double x(0.6), y(0.75);
  legg = new TLegend(x, y, x+0.15, y+0.1);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
  legge = legg->AddEntry(hoptSg, Form("Signal"), "f"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hoptBg, Form("Background"), "f"); legge->SetTextColor(kBlack);
  legg->Draw();

  c0->cd(1);
  shrinkPad(0.3, 0.2);
  hoptSgF->SetLineColor(kBlue);
  hoptSgF->SetFillStyle(3001);
  hoptSgF->GetYaxis()->SetRangeUser(0, 1.2);
  hoptSgF->Draw();
  hoptBgF->Draw("SAME");
  legg = new TLegend(x, y, x+0.15, y+0.1);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
  legge = legg->AddEntry(hoptSg, Form("Signal 2D"), "f"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hoptBg, Form("Background 2D"), "f"); legge->SetTextColor(kBlack);
  legg->Draw();

  c0->SaveAs("optimizedCuts2D.pdf");

}


// -----------------------------------------------------------------------------------------------------
void anaBmm::loopOptimizedCuts3D() {

  char extra[200];
  char ao[200];

  int icut = 1;
  
  for (int i = 0; i < NCUTS; i++ ) {

    if ( optimizedCut[i].Contains("lxy/sxy") ) { continue; }
    if ( optimizedCut[i].Contains("cosa3") )   { continue; }
    if ( optimizedCut[i].Contains("eta") )     { continue; }

    sprintf(ao, "4.8<mass&&mass<6.0&&goodL1&&goodHLT");

    for (int j = 0; j < NCUTS; j++ ) {

      if ( i == j ) { continue; }

      if ( optimizedCut[j].Contains("lxy/sxy") ) { continue; }
      if ( optimizedCut[j].Contains("cosa3") )   { continue; }
      if ( optimizedCut[j].Contains("iso") )     { continue; }
      if ( optimizedCut[j].Contains("chi2") )    { continue; }
      if ( optimizedCut[j].Contains("eta") )     { continue; }

      sprintf(ao, "%s && %s", ao, optimizedCut[j].Data());
    }

    sprintf(extra, "%s", optimizedCut[i].Data());
    cutEfficiency(ao, extra, icut);
    icut++;

  }

  c0->cd();
  c0->Clear();
  c0->Divide(1,2);

  c0->cd(2);
  shrinkPad(0.25, 0.2);
  hoptSg->SetLineColor(kBlue);
  hoptSg->SetFillStyle(3001);
  hoptSg->GetYaxis()->SetRangeUser(0, 1.2);
  hoptSg->Draw();
  hoptBg->Draw("SAME");

  double x(0.6), y(0.75);
  legg = new TLegend(x, y, x+0.15, y+0.1);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
  legge = legg->AddEntry(hoptSg, Form("Signal 3D"), "f"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hoptBg, Form("Background 3D"), "f"); legge->SetTextColor(kBlack);
  legg->Draw();

  c0->cd(1);
  shrinkPad(0.25, 0.2);
  hoptSgF->SetLineColor(kBlue);
  hoptSgF->SetFillStyle(3001);
  hoptSgF->GetYaxis()->SetRangeUser(0, 1.2);
  hoptSgF->Draw();
  hoptBgF->Draw("SAME");
  legg = new TLegend(x, y, x+0.15, y+0.1);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); 
  legge = legg->AddEntry(hoptSg, Form("Signal 3D"), "f"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hoptBg, Form("Background 3D"), "f"); legge->SetTextColor(kBlack);
  legg->Draw();
  
  c0->SaveAs("optimizedCuts3D.pdf");
  
}


// ----------------------------------------------------------------------
void anaBmm::cutEfficiency(const char *aoCuts, const char *extraCuts, int cut_nr) {


  if (!hoptSg)   hoptSg     = new TH1D("hsgC","#varepsilon_{cut}  after cuts (w/o fact. cuts)", NCUTS, 0., NCUTS);
  if (!hoptBg)   hoptBg     = new TH1D("hbgC","#varepsilon_{cut}  after cuts (w/o fact. cuts)", NCUTS, 0., NCUTS);
  if (!hoptSgF)  hoptSgF    = new TH1D("hsgF","#varepsilon_{cut}", NCUTS, 0., NCUTS);
  if (!hoptBgF)  hoptBgF    = new TH1D("hbgF","#varepsilon_{cut}", NCUTS, 0., NCUTS);

  hoptSg->GetXaxis()->SetBinLabel(cut_nr, Form("%s", extraCuts)); hoptSg->GetXaxis()->LabelsOption("v");
  hoptBg->GetXaxis()->SetBinLabel(cut_nr, Form("%s", extraCuts)); hoptBg->GetXaxis()->LabelsOption("v");
  hoptSgF->GetXaxis()->SetBinLabel(cut_nr, Form("%s", extraCuts)); hoptSgF->GetXaxis()->LabelsOption("v");
  hoptBgF->GetXaxis()->SetBinLabel(cut_nr, Form("%s", extraCuts)); hoptBgF->GetXaxis()->LabelsOption("v");


  // -- Run on signal MC to determine efficiency 
  fS[sgIndex]->cd(); 
  TH1D *hSG = (TH1D*)gROOT->FindObject("hSG"); 
  if (!hSG) hSG = new TH1D("hSG", "", 50, 0., 10.);
  TTree *s = (TTree*)gFile->Get("events");
  s->Draw("mass>>hSG", Form("goodKinematics"), "goff"); 
  double sNorm = hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("goodL1&&goodHLT && %s", aoCuts), "goff"); 
  double s1Norm = hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("goodL1&&goodHLT && %s", extraCuts), "goff"); 
  double s2Norm  =  hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("goodL1&&goodHLT && %s && %s", aoCuts, extraCuts), "goff"); 
  double s12Norm  =  hSG->GetSumOfWeights();

  // -- Run on bg MC
  fM[bgIndex]->cd(); 
  TH1D *hBG = (TH1D*)gROOT->FindObject("hBG"); 
  if (!hBG) hBG = new TH1D("hBG", "", 50, 0., 10.);
  TTree *b = (TTree*)gFile->Get("events");
  b->Draw("mass>>hBG", Form("goodKinematics"), "goff"); 
  double bNorm   = hBG->GetSumOfWeights();
  b->Draw("mass>>hBG", Form("goodL1&&goodHLT && %s", aoCuts), "goff");
  double b1Norm  = hBG->GetSumOfWeights();
  b->Draw("mass>>hBG", Form("goodL1&&goodHLT && %s", extraCuts), "goff"); 
  double b2Norm  = hBG->GetSumOfWeights();
  b->Draw("mass>>hBG", Form("goodL1&&goodHLT && %s && %s", aoCuts, extraCuts), "goff"); 
  double b12Norm  =  hBG->GetSumOfWeights();
  
  double epsilon = 0.1;

  if ( sNorm < epsilon)  sNorm = epsilon;
  if ( s1Norm < epsilon) s1Norm = epsilon;
  if ( bNorm < epsilon)  bNorm = epsilon;
  if ( b1Norm < epsilon) b1Norm = epsilon;


  hoptSg->SetBinContent(cut_nr, s12Norm/s1Norm);
  hoptSgF->SetBinContent(cut_nr, s2Norm/sNorm);
  hoptBg->SetBinContent(cut_nr, b12Norm/b1Norm);
  hoptBgF->SetBinContent(cut_nr, b2Norm/bNorm);


  cout << "\"" << extraCuts << "\" : " << Form("signal-Eff = %4.4e (%4.4e),  bg-Eff = %4.4e (%4.4e)"
	       ,  s12Norm/s1Norm, s2Norm/sNorm, b12Norm/b1Norm, b2Norm/bNorm)  << endl;
  
  char filename[200];
  sprintf(filename, "handOptimizationEfficiencies.txt"); 
  
  ofstream OUT(filename, ios::app);

  OUT << "\"" << extraCuts << "\" : " << Form("signal-Eff = %4.4e (%4.4e),  bg-Eff = %4.4e (%4.4e)"
	       ,  s12Norm/s1Norm, s2Norm/sNorm, b12Norm/b1Norm, b2Norm/bNorm) << endl;

  OUT.close();
}

// =======================================================================================================

void anaBmm::loopULOptimization(const char *filename) {

  char pre[200]; char preSel[200];
  sprintf(preSel, "lxy/sxy>7");

  char aoCuts[200];
  sprintf(aoCuts, "goodL1 && goodHLT && TMath::Abs(etal0)<1.5 && TMath::Abs(etal1)<1.5 && rmm>-0.1 && rmm<100. && TMath::Abs(eta)<2.4000");

  int ncuts = 4; int nfact = 2;
  TString extraTS[] = { TString("ptl1>"), TString("pt>"), TString("l3d/s3d>"), TString("cosa>"), 
			TString("iso>"), TString("chi2<") };

  // arrays for:     ptMu  ptB  l3D   cosa    iso chi2
  double mins[] = {  4.,   5.,  22., 0.9985,  0.85,  2.}; 
  double maxs[] = {  7.,   9.,  25., 1.000,  1.00,  10.};
  int bins[]    = {  1,    1,    1,     4,     1,   8 };

  char extra1[200], extra2[200], extra3[200], extra4[200], extra5[200], extra6[200];

  char ao[200];
  char fact1[200];
  char fact2[200];
  
  TString cutsBestUL;
  double minUL(99.);
  double cut(0.);

  for (int i1 = 0; i1 < ncuts; i1++ ) {
    for (int i1bin = 0; i1bin < bins[i1]; ++i1bin) {  
      cut = mins[i1] + i1bin*(maxs[i1]-mins[i1])/bins[i1]; 
      sprintf(extra1, Form("%s %5.4f", extraTS[i1].Data(), cut));

      for (int i2 = i1+1; i2 < ncuts; i2++ ) {
	for (int i2bin = 0; i2bin < bins[i2]; ++i2bin) {  
	  cut = mins[i2] + i2bin*(maxs[i2]-mins[i2])/bins[i2]; 
	  sprintf(extra2, Form("%s %5.4f", extraTS[i2].Data(), cut));

	  for (int i3 = i2+1; i3 < ncuts; i3++ ) {
	    for (int i3bin = 0; i3bin < bins[i3]; ++i3bin) {  
	      cut = mins[i3] + i3bin*(maxs[i3]-mins[i3])/bins[i3]; 
	      sprintf(extra3, Form("%s %5.4f", extraTS[i3].Data(), cut));

	      for (int i4 = i3+1; i4 < ncuts; i4++ ) {
		for (int i4bin = 0; i4bin < bins[i4]; ++i4bin) {  
		  cut = mins[i4] + i4bin*(maxs[i4]-mins[i4])/bins[i4]; 
		  sprintf(extra4, Form("%s %5.4f", extraTS[i4].Data(), cut));
 
// 		  for (int i5 = i4+1; i5 < ncuts; i5++ ) {
// 		    for (int i5bin = 0; i5bin < bins[i5]; ++i5bin) {  
// 		      cut = mins[i5] + i5bin*(maxs[i5]-mins[i5])/bins[i5]; 
// 		      sprintf(extra5, Form("%s %5.4f", extraTS[i5].Data(), cut));

// 		      for (int i6 = i5+1; i6 < ncuts; i6++ ) {
// 			for (int i6bin = 0; i6bin < bins[i6]; ++i6bin) {  
// 			  cut = mins[i6] + i6bin*(maxs[i6]-mins[i6])/bins[i6]; 
// 			  sprintf(extra6, Form("%s %5.4f", extraTS[i6].Data(), cut));
			  
			  for (int j1 = ncuts; j1 < ncuts+nfact; j1++ ) {
			    for (int j1bin = 0; j1bin < bins[j1]; ++j1bin) {  
			      cut = mins[j1] + j1bin*(maxs[j1]-mins[j1])/bins[j1]; 
			      sprintf(fact1, Form("%s %5.4f", extraTS[j1].Data(), cut));
		
			      for (int j2 = j1+1; j2 < ncuts+nfact; j2++ ) {
				for (int j2bin = 0; j2bin < bins[j2]; ++j2bin) {  
				  
				  cut = mins[j2] + j2bin*(maxs[j2]-mins[j2])/bins[j2]; 
				  sprintf(fact2, Form("%s %5.4f", extraTS[j2].Data(), cut));

				  sprintf(ao, Form("%s && %s && %s && %s && %s", 
						   aoCuts, extra1, extra2, extra3, extra4));
// 						   aoCuts, extra1, extra2, extra3, extra4, extra5, extra6));

				  sprintf(pre, Form("%s && %s", preSel, extra1));
// 				  sprintf(pre, Form("%s && %s && %s && %s ", preSel, extra1, extra2, extra3));

    
				  ulOptimization(ao, pre, fact1, fact2, 0.06, filename); 
				  
				  if (fUL < minUL) {
				    
				    minUL = fUL; 
				    cutsBestUL = TString(Form("%s, fact1: %s, fact2:  %s, pre:  %s", ao, fact1, fact2, pre));
				  } 
				}
			      }
			    }
			  }
// 		        }
// 		      }
// 		    }
// 		  }
 		}
 	      }
 	    }
 	  }
	}
      }
    }
  }

  cout << "==> minimum UL: " << Form("%4.4e: ",minUL)  << cutsBestUL << " !!! " << endl;
  
  ofstream OUT(filename, ios::app);

  OUT << endl <<  "==> minimum UL: " << Form("%4.4e: ",minUL)  << cutsBestUL << " !!! " << endl;


  OUT.close();

}


// ----------------------------------------------------------------------
void anaBmm::ulOptimization(const char *aoCuts, const char *preCuts, const char *f1Cut, const char *f2Cut, 
			    double mCut, const char *filename){


  // -- Fix lumi normalisation scaling factor
  double sSF = fLumiD[0]/fLumiS[sgIndex];

  // -- Run on signal MC to determine efficiency 
  fS[sgIndex]->cd(); 
  TH1D *uSG = (TH1D*)gROOT->FindObject("uSG"); 
  if (!uSG) uSG = new TH1D("uSG", "", 120, 4.8, 6.0);
  TTree *s = (TTree*)gFile->Get("events");
  
  double sNorm = fNevtS[sgIndex];
  double massRedSG = 0.;

  s->Draw("mass>>uSG", Form("goodKinematics&&goodL1&&goodHLT"), "goff");
  if ( mCut > -0.5 || fSgReduction < 0.) {
    fSgReduction = massReduction("uSG", "mysg", mCut);
    massRedSG = fSgReduction;
  } else {
    massRedSG = fSgReduction;
  }

  s->Draw("mass>>uSG", Form("goodL1&&goodHLT && %s", aoCuts), "goff"); 
  double sCuts     = uSG->GetSumOfWeights();
  double sEff      = sCuts/sNorm;

  s->Draw("mass>>uSG", Form("goodL1&&goodHLT && %s", preCuts), "goff"); 
  double sNormF    =  uSG->GetSumOfWeights();

  s->Draw("mass>>uSG", Form("goodL1&&goodHLT && %s && %s", preCuts, f1Cut), "goff"); 
  double s1Norm    = uSG->GetSumOfWeights();
  double sEff1     = s1Norm/sNormF;

  s->Draw("mass>>uSG", Form("goodL1&&goodHLT && %s && %s", preCuts, f2Cut), "goff"); 
  double s2Norm    = uSG->GetSumOfWeights();
  double sEff2     = s2Norm/sNormF;

  fNsg  =  sCuts * sEff1 * sEff2 * massRedSG * sSF;
  fNsgE = 0.;
  //fNsgE = TMath::Sqrt(fNsg);

  fEsg  = sEff * sEff1 * sEff2 * massRedSG;

  // -- Fix lumi normalisation scaling factor
  double bSF = fLumiD[0]/fLumiM[bgIndex];

  // -- Run on signal MC to determine efficiency 
  fM[bgIndex]->cd(); 
  TH1D *uBG = (TH1D*)gROOT->FindObject("uBG"); 
  if (!uBG) uBG = new TH1D("uBG", "", 120, 4.8, 6.0);
  TTree *b = (TTree*)gFile->Get("events");
  
  double bNorm = fNevtM[bgIndex];
  double massRedBG = 0.;

  b->Draw("mass>>uBG", Form("goodKinematics&&goodL1&&goodHLT"), "goff");
  if ( mCut > -0.5 || fBgReduction < 0.) {
    fBgReduction = massReduction("uBG", "mymc", mCut);
    massRedBG = fBgReduction;
  } else {
    massRedBG = fBgReduction;
  }

  b->Draw("mass>>uBG", Form("goodL1&&goodHLT && %s", aoCuts), "goff"); 
  double bCuts     = uBG->GetSumOfWeights();
  double bEff      = bCuts/bNorm;

  b->Draw("mass>>uBG", Form("goodL1&&goodHLT && %s", preCuts), "goff"); 
  double bNormF    =  uBG->GetSumOfWeights();

  b->Draw("mass>>uBG", Form("goodL1&&goodHLT && %s && %s", preCuts, f1Cut), "goff"); 
  double b1Norm    = uBG->GetSumOfWeights();
  double bEff1     = b1Norm/bNormF;

  b->Draw("mass>>uBG", Form("goodL1&&goodHLT && %s && %s", preCuts, f2Cut), "goff"); 
  double b2Norm    = uBG->GetSumOfWeights();
  double bEff2     = b2Norm/bNormF;

  fNbg  =  bCuts * bEff1 * bEff2 * massRedBG * bSF;
  fNbgE =  fNbg * 1.6/(1.*TMath::Sqrt(bCuts));

  cout << "Calculate UL for n_S = " << fNsg << " and n_B = " << fNbg << endl;
  
  if ( (fNsg < 100 && fNbg < 1000) && (fNsg > 1e-20 && fNbg > 1e-20) ) {
    fUL = calculateUpperLimit();
  } else {
    cout << " ==> Cannot calculate UL for n_S = " << fNsg << " and n_B = " << fNbg << endl;
    fUL = 99999.;
  }
  
  cout << Form("UL = %4.4e", fUL)
       << " cuts: \"" << aoCuts << "\"  " 
       << " fact1: \"" << f1Cut << "\"  " 
       << " fact2: \"" << f2Cut << "\"  " 
       << endl;
  
  
  ofstream OUT(filename, ios::app);
  
  OUT <<  Form("UL = %4.4e", fUL)  << "   :    n_S = " << fNsg << ", n_B = " << fNbg << endl; 
  OUT << "  -> cuts: \"" << aoCuts << "\"  " 
      << " fact1: \"" << f1Cut << "\"  " 
      << " fact2: \"" << f2Cut << "\"  " 
      << " preCuts: \"" << preCuts << "\"  " 
      << endl;

  OUT.close();
}


// =======================================================================================================

void anaBmm::massResolutionEta() {

  char pre[200];
  double etalo(0.), etahi(100.);
  double veta[]  = { 0., 0.5, 1.0, 1.5, 2.5 }; int neta = 4;

  double ptcut(3.);
  double vpt[] = {3., 4., 5.}; const int npt = 3;

  TH1D *em[npt];

  for (int j = 0; j < npt; j++) { 

    em[j]= new TH1D(Form("em%i", j), "",  neta, veta); 
    setTitles(em[j], "#eta", "#Delta m / MeV", 0.06, 1.1, 1.2);
    ptcut = vpt[j];
    
    for (int i = 0; i < neta; i++ ) {
      
      etalo = veta[i];
      etahi = veta[i+1];
      
      sprintf(pre, "ptl1>%5.4f && TMath::Abs(etal0)>%5.4f && TMath::Abs(etal1)>%5.4f && TMath::Abs(etal0)<%5.4f && TMath::Abs(etal1)<%5.4f", ptcut,  etalo, etalo, etahi, etahi);
      
            
      // -- Run on signal MC to determine efficiency 
      fS[sgIndex]->cd(); 
      TH1D *mSG = (TH1D*)gROOT->FindObject("mSG"); 
      if (!mSG) mSG = new TH1D("mSG", "", 120, 4.8, 6.0);
      TTree *s = (TTree*)gFile->Get("events");
      
      s->Draw("mass>>mSG", Form("goodL1&&goodHLT && %s", pre), "goff"); 
      
      emptyBinError(mSG);
      
      mSG->Scale(fLumiD[0]/fLumiS[sgIndex]);
      
      mSG->GetXaxis()->SetRangeUser(4.8, 6.);
    
      f1->SetParameters(mSG->GetMaximum()*0.8, mSG->GetMean(), 0.5*mSG->GetRMS(), 
			mSG->GetMaximum()*0.2, mSG->GetMean(), 3.*mSG->GetRMS());
    
    
      shrinkPad(0.15, 0.15);
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      setFilledHist(mSG, kBlack, kYellow, 1000); 
      setTitles(mSG, "m_{#mu#mu} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
      mSG->GetYaxis()->SetRangeUser(0., 1.4*mSG->GetMaximum());
    
      mSG->DrawCopy("hist");
      mSG->Fit(f1, "0");
      f1->DrawCopy("same");
    
      double mean = TMath::Sqrt((f1->GetParameter(0)*f1->GetParameter(0)*f1->GetParameter(1)*f1->GetParameter(1)
				 + f1->GetParameter(3)*f1->GetParameter(3)*f1->GetParameter(4)*f1->GetParameter(4))
				/
				(f1->GetParameter(0)*f1->GetParameter(0) + f1->GetParameter(3)*f1->GetParameter(3))
				);
    
      double  sigma = TMath::Sqrt((f1->GetParameter(0)*f1->GetParameter(0)*f1->GetParameter(2)*f1->GetParameter(2)
				   + f1->GetParameter(3)*f1->GetParameter(3)*f1->GetParameter(5)*f1->GetParameter(5))
				  /
				  (f1->GetParameter(0)*f1->GetParameter(0) + f1->GetParameter(3)*f1->GetParameter(3))
				  );

      double  meanE  = getDGError(f1->GetParameter(1), f1->GetParError(1), 
				  f1->GetParameter(4), f1->GetParError(4),
				  f1->GetParameter(0), f1->GetParError(0),
				  f1->GetParameter(3), f1->GetParError(3));
    
      double  sigmaE = getDGError(f1->GetParameter(2), f1->GetParError(2), 
				  f1->GetParameter(5), f1->GetParError(5),
				  f1->GetParameter(0), f1->GetParError(0),
				  f1->GetParameter(3), f1->GetParError(3));
    
      em[j]->SetBinContent(i+1, 1000*sigma);
      em[j]->SetBinError(i+1, 1000*sigmaE);

      tl->SetTextColor(kBlack);  
      tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
    
      tl->DrawLatex(0.16, 0.85, Form("#mu: %5.3f#pm%5.3f GeV", mean, meanE));
      tl->DrawLatex(0.16, 0.79, Form("#sigma: %5.3f#pm%5.3f GeV", sigma, sigmaE));
      tl->SetTextSize(0.04);
      tl->SetTextColor(kRed);
      tl->DrawLatex(0.22, 0.60, Form("p_{T} > %1.0f", ptcut));
      tl->DrawLatex(0.18, 0.55, Form("%2.1f < #eta < %2.1f", etalo, etahi)); tl->SetTextSize(0.06); tl->SetTextColor(kBlack);
  
      if (mSG->GetSumOfWeights()<0.001) {
	tl->DrawLatex(0.16, 0.73, Form("N: %4.2e", mSG->GetSumOfWeights()));
      } else if (mSG->GetSumOfWeights() < 1.) {
	tl->DrawLatex(0.16, 0.73, Form("N: %5.3f", mSG->GetSumOfWeights()));
      } else if (mSG->GetSumOfWeights() < 1000.) {
	tl->DrawLatex(0.16, 0.73, Form("N: %5.1f", mSG->GetSumOfWeights()));
      } else {
	tl->DrawLatex(0.16, 0.73, Form("N: %4.2e", mSG->GetSumOfWeights()));
      }
    
      tl->SetTextSize(0.06); tl->SetTextColor(kBlue);
      tl->DrawLatex(0.56, 0.65, Form("%s", fSignLeggS[sgIndex].Data()));
    
      c0->SaveAs(Form("%s/massReduction/massResolution-eta%i-pt%i.eps", outDir, i, j));
    
    }
  }


  for ( int j = npt - 1; j > -1; j-- ) {
      
    em[j]->SetMarkerColor(npt-j);
    em[j]->SetLineColor(npt-j);
    em[j]->SetMarkerSize(1);
    tl->SetTextColor(npt-j);
    tl->SetTextSize(0.05);
    
    if (j==npt-1) { 
      em[j]->GetYaxis()->SetRangeUser(0., 250.);
      em[j]->DrawCopy("e");
    } else {
      em[j]->DrawCopy("samee");
    }

    tl->DrawLatex(0.2, 0.7+j*0.08, Form("p_{T}^{#mu} > %4.0f MeV", vpt[j]));
  }

  tl->SetTextColor(kBlack);
  c0->SaveAs(Form("%s/massReduction/massResolution.eps", outDir));
}



// =======================================================================================================

void anaBmm::massReductionEta() {

  char pre[200];
  double etalo(0.), etahi(100.);
  double veta[]  = { 0., 0.5, 1.0, 1.5, 2.5 }; const int neta = 4;
  double vres[]  = { 0.06, 0.08, 0.1 , 0.1  }; 

  double window(0.1);
  double vmass[] = {0.06, 0.08, 0.1 }; const int nmass = 3;


  TH1D *em[nmass];
  TH1D *emf = new TH1D("emf", "",  neta, veta);

  for (int j = 0; j < nmass; j++) { 

    em[j]= new TH1D(Form("em%i", j), "",  neta, veta); 
    setTitles(em[j], "#eta", "f", 0.06, 1.1, 1.2);
    window = vmass[j];

    for (int i = 0; i < neta; i++ ) {
      
      etalo = veta[i];
      etahi = veta[i+1];
      
      sprintf(pre, "ptl1>%5.4f && TMath::Abs(etal0)>%5.4f && TMath::Abs(etal1)>%5.4f && TMath::Abs(etal0)<%5.4f && TMath::Abs(etal1)<%5.4f", ptmulo,  etalo, etalo, etahi, etahi);
      
      
      // -- Run on signal MC to determine efficiency       
      fM[bgIndex]->cd(); 
      TH1D *mBG1 = (TH1D*)gROOT->FindObject("mBG1"); 
      if (!mBG1) mBG1 = new TH1D("mBG1", "", 120, 4.8, 6.0);
      TH1D *mBG2 = (TH1D*)gROOT->FindObject("mBG2"); 
      if (!mBG2) mBG2 = new TH1D("mBG2", "", 120, 4.8, 6.0);
      TTree *b = (TTree*)gFile->Get("events");
      
      b->Draw("mass>>mBG1", Form("goodL1&&goodHLT && %s", pre), "goff"); 
      b->Draw("mass>>mBG2", Form("goodL1&&goodHLT && %s", pre), "goff"); 
      
      emptyBinError(mBG1);
      
      mBG2->GetXaxis()->SetRangeUser(4.8, 6.0);
      mBG1->GetXaxis()->SetRangeUser(4.8, 6.0);
      mBG1->Fit("pol1");
      
      TF1 *f1 = (TF1*)mBG1->GetFunction("pol1");
      
      
      shrinkPad(0.15, 0.15);
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      setTitles(mBG2, "m_{#mu#mu} [GeV]", "events/bin", 0.06, 1.1, 1.3); 
      
      mBG2->Draw("e");
      f1->Draw("same");
      
  
      double total = f1->Integral(mBG1->GetBinLowEdge(1), mBG1->GetBinLowEdge(mBG1->GetNbinsX()))/mBG1->GetBinWidth(2); 
      double massW = f1->Integral(fMassBs - window, fMassBs + window)/mBG1->GetBinWidth(2);
      double massRed   = massW/total;
      double massRedE  = f1->GetParError(1)*massW/total;
      
      em[j]->SetBinContent(i+1, massRed);
      em[j]->SetBinError(i+1, massRedE);
      
      if ( TMath::Abs(vres[i] - vmass[j]) < 0.00001 ) {

	emf->SetBinContent(i+1, massRed);
      }


      tl->SetTextColor(kBlack);  
      tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
      
      tl->DrawLatex(0.16, 0.85, Form("f: %5.3f#pm%5.3f", massRed, massRedE));
      tl->SetTextSize(0.04);
      tl->SetTextColor(kRed);
      tl->DrawLatex(0.18, 0.60, Form("m_{#mu#mu} #pm %3.0f", 1000*window));
      tl->DrawLatex(0.18, 0.55, Form("%2.1f < #eta < %2.1f", etalo, etahi)); tl->SetTextSize(0.06); tl->SetTextColor(kBlack);
      
      tl->SetTextSize(0.06); tl->SetTextColor(kBlue);
      tl->DrawLatex(0.56, 0.75, Form("%s", fSignLeggM[bgIndex].Data()));
      
      c0->SaveAs(Form("%s/massReduction/massReduction-eta%i-mass%i.eps", outDir, i, j));
      
    }
  }

  emf->SetLineWidth(2);
  emf->GetYaxis()->SetRangeUser(0., 0.25);
  emf->DrawCopy("l");
  tl->SetTextColor(kBlack);
  tl->SetTextSize(0.05);
  tl->DrawLatex(0.2, 0.15+nmass*0.05, "#Delta m ~ 2 #sigma");

  for ( int j = nmass - 1; j > -1; j-- ) {
    
    em[j]->SetMarkerColor(j+2);
    em[j]->SetLineColor(j+2);
    em[j]->SetMarkerSize(1);
    tl->SetTextColor(j+2);

    em[j]->DrawCopy("samee");
    
    tl->DrawLatex(0.5, 0.17+j*0.05, Form("%4.0f MeV", vmass[j]*1000));
  }

  emf->DrawCopy("samel");
  tl->SetTextColor(kBlack);
  
  c0->SaveAs(Form("%s/massReduction/massReduction.eps", outDir));
}
  
  // =======================================================================================================

void anaBmm::optimizerMassEta() {

  double veta[]  = { 0., 0.5, 1.5, 2.5 }; int neta = 3;
  double vmass[] = { 0.1, 0.08, 0.06 };   int nmass = 3;
  double vmass2[] = { 50, 70, 90, 110 };

  TString cutsBestUL;
  double minUL(9999.);

  TH2D *hmeU = new TH2D("hmeU", "",  neta, veta, nmass, vmass2); 
  TH2D *hmeS = new TH2D("hmeS", "",  neta, veta, nmass, vmass2); 
  TH2D *hmeB = new TH2D("hmeB", "",  neta, veta, nmass, vmass2); 
  TH2D *hmeE = new TH2D("hmeE", "",  neta, veta, nmass, vmass2); 

  setTitles2(hmeU, "#eta","mass", 0.06, 1.1, 1.2);
  setTitles2(hmeS, "#eta","mass", 0.06, 1.1, 1.2);
  setTitles2(hmeB, "#eta","mass", 0.06, 1.1, 1.2);
  setTitles2(hmeE, "#eta","mass", 0.06, 1.1, 1.2);

  for (int i = 0; i < neta; i++ ) {
    for (int j = 0; j < nmass; j++ ) {

       cout << "quickAna(" << ptmulo  << " ," <<  veta[i]  << " ," <<  veta[i+1]  << " ," <<  rmmlo  << " ," <<  rmmhi  << " ," <<  ptbs  << " ," <<  etahi  << " ," <<  coslo  << " , -100. ," <<  lxylo  << " ," <<  l3dlo  << " ," <<  vtxhi  << " ," <<  isolo  << " ," <<  vmass[j] << ");" << endl;
      
      quickAna(ptmulo, veta[i], veta[i+1], rmmlo, rmmhi, ptbs, etahi, coslo, -100., lxylo, l3dlo, vtxhi, isolo, vmass[j]);

      hmeU->SetBinContent(i+1, j+1, fUL);
      hmeS->SetBinContent(i+1, j+1, fNsg);
      hmeB->SetBinContent(i+1, j+1, fNbg);
      hmeE->SetBinContent(i+1, j+1, fEsg);

      if (fUL < minUL) {
	
	minUL = fUL; 
	cutsBestUL = TString(Form("%4.1f < |eta| < %4.1f, mass +/- %3.0f MeV", veta[i], veta[i+1], 1000*vmass[j]));
      } 

    }
  }

  c0->cd();
  c0->Clear();
  c0->Divide(2,2);
  gStyle->SetPaintTextFormat("4.2e");
  c0->cd(1);  shrinkPad(0.2, 0.2, 0.2, 0.1); gStyle->SetOptStat(0); gPad->SetLogz(1);
  hmeU->GetZaxis()->SetRangeUser(1.e-8, 2.e-8); hmeU->SetMarkerSize(2);
  hmeU->DrawCopy("colztext45");
  c0->cd(2); 
  shrinkPad(0.2, 0.2, 0.2, 0.1); gStyle->SetOptStat(0);  
  hmeE->GetZaxis()->SetRangeUser(0., 0.01); hmeE->SetMarkerSize(2);
  hmeE->DrawCopy("colztext45");
  c0->cd(3); shrinkPad(0.2, 0.2, 0.2, 0.1); gStyle->SetOptStat(0); 
  hmeS->GetZaxis()->SetRangeUser(0., 6.); hmeS->SetMarkerSize(2);
  hmeS->DrawCopy("colztext45");
  c0->cd(4);  shrinkPad(0.2, 0.2, 0.2, 0.1); gStyle->SetOptStat(0); 
  hmeB->GetZaxis()->SetRangeUser(0., 10.); hmeB->SetMarkerSize(2);
  hmeB->DrawCopy("colztext45");
 
  tl->SetTextColor(kBlack);  
  tl->SetNDC(kTRUE); tl->SetTextSize(0.06);
  c0->cd(1);
  tl->DrawLatex(0., 0.94, Form("Best UL = %4.2e", minUL)); tl->SetTextSize(0.04);
  tl->DrawLatex(0.46, 0.95, Form("(%s)", cutsBestUL.Data())); tl->SetTextSize(0.07);
  c0->cd(2);
  tl->DrawLatex(0.25, 0.94, "#varepsilon_{sg}");
  c0->cd(3);
  tl->DrawLatex(0.25, 0.94, "N_{sg}");
  c0->cd(4);
  tl->DrawLatex(0.25, 0.94, "N_{bg}");

  c0->SaveAs("massEtaOptimizer.pdf");

}


// =======================================================================================================

void anaBmm::quickAna(double vptmu, double vetamulo, double vetamuhi, double vrmmlo, double vrmmhi, 
		      double vptb, double vetab, double vcosa, double vcosa3, double vlxysxy, double vl3ds3d, 
		      double vchi2, double viso, 
		      double vmass, double vpre) {

  char pre[200];
  //  sprintf(pre, "ptl1>%5.4f && rmm>%5.4f && rmm<%5.4f && lxy/sxy>%5.4f", vptmu, vrmmlo, vrmmhi, vpre);
  sprintf(pre, "ptl1>%5.4f && rmm>%5.4f && rmm<%5.4f && lxy/sxy>%5.4f && TMath::Abs(etal0)>%5.4f && TMath::Abs(etal1)>%5.4f && TMath::Abs(etal0)<%5.4f && TMath::Abs(etal1)<%5.4f", vptmu, vrmmlo, vrmmhi, vpre, vetamulo, vetamulo, vetamuhi, vetamuhi);

  char fact1[200];  char f1[200]; char tf1[200];
  sprintf(f1,  "chi2<%5.4f", vchi2);
  sprintf(tf1,  "$\\chi^{2} < %5.1f$", vchi2);
  
  char fact2[200];  char f2[200]; char tf2[200];
  sprintf(f2, "iso>%5.4f", viso);
  sprintf(tf2,  "$I < %4.3f$", viso);

  char cuts[500], fcuts[500];
  const int ncuts = 12; 
  TString cutsTS[] = { TString(Form("ptl1>%5.4f",vptmu)), 
		       TString(Form("TMath::Abs(etal0)>%5.4f && TMath::Abs(etal1)>%5.4f && TMath::Abs(etal0)<%5.4f && TMath::Abs(etal1)<%5.4f", vetamulo, vetamulo, vetamuhi, vetamuhi)), 
		       TString(Form("rmm>%5.4f",vrmmlo)), TString(Form("rmm<%5.4f", vrmmhi)), 
		       TString(Form("pt>%5.4f", vptb)), TString(Form("TMath::Abs(eta)<%5.4f", vetab)), 
		       TString(Form("cosa>%5.4f", vcosa)), TString(Form("cosa3>%5.4f", vcosa3)), 
		       TString(Form("lxy/sxy>%5.4f", vlxysxy)), TString(Form("l3d/s3d>%5.4f", vl3ds3d)),
		       TString(Form("chi2<%5.4f", vchi2)), TString(Form("iso>%5.4f", viso))  };

  TString tableTS[] = { TString(Form("$p_T > %5.1f$",vptmu)), 
			TString(Form("$%5.1f < |\\eta_\\mu| < %5.1f$ ", vetamulo, vetamuhi)), 
			TString(Form("$\\Delta R > %5.1f$",vrmmlo)), TString(Form("$\\Delta R < %5.1f$", vrmmhi)), 
			TString(Form("$p_T(B_{s}) > %5.1f$", vptb)), TString(Form("$\\eta_{B} < %5.1f$", vetab)), 
			TString(Form("$cos \\alpha > %5.4f$", vcosa)), TString(Form("$cos \\alpha_{3D} > %5.4f$", vcosa3)), 
			TString(Form("$l_{xy}/s_{xy} > %5.1f$", vlxysxy)), TString(Form("$l_{3D}/s_{3D}>%5.1f$", vl3ds3d)), 
			TString(Form("$\\chi^2 < %5.1f$", vchi2)), TString(Form("$I > %5.3f$", viso)) };
  
  double bEff[ncuts];
  double sNevt[ncuts];
  double bNevt[ncuts];

  double sNorm(-1.),  bNorm(-1.);
  double sNormF(-1.), bNormF(-1.);
  			
  sprintf(cuts, "goodL1&&goodHLT");
  sprintf(fcuts, "goodL1&&goodHLT");

  sNorm = fNevtS[sgIndex];
  bNorm = fNevtM[bgIndex];

  sNormF = getCutEffS(cuts, pre);
  bNormF = getCutEffB(cuts, pre);

  sprintf(fact1, "%s && %s", f1, cuts);
  double sEff1 =  getCutEffS(fact1, pre)/(1.*sNormF);
  double bEff1 =  getCutEffB(fact1, pre)/(1.*bNormF);
  cout << "fact: " << f1  << "     e_S = " << sEff1 << "     e_B = " << bEff1 << endl;
  
  sprintf(fact2, "%s && %s", f2, cuts);
  double sEff2 =  getCutEffS(fact2, pre)/(1.*sNormF);
  double bEff2 =  getCutEffB(fact2, pre)/(1.*bNormF);
  cout << "fact: " << f2  << "     e_S = " << sEff2 << "     e_B = " << bEff2 << endl;

  char filename[200];
  sprintf(filename, "cuts_table.tex"); 
  
  ofstream OUT(filename, ios::app);
  double sEff[ncuts];

  OUT << "  & " <<  " $N_S$" << "  &"  << "$\\varepsilon_S$" 
      << "  & " << " $N_B$ " << "  & " << "$\\varepsilon_B$" << " \\\\" << endl;
  OUT << "\\hline" << endl;
  

  for (int i = 0; i < ncuts; i++ ) {

    if (i<NCUTS) optimizedCut[i] = TString(Form("%s", cutsTS[i].Data()));

    sprintf(cuts, "%s && %s", cuts, cutsTS[i].Data());
    
    sNevt[i] = getCutEffS(cuts);
    bNevt[i] = getCutEffB(cuts);

    sEff[i] = sNevt[i]/(1.*sNorm);
    bEff[i] = bNevt[i]/(1.*bNorm);

    cout << "cut: " << cutsTS[i].Data() 
	 << "     e_S = " << sEff[i] << "     e_B = " << bEff[i]
	 << "     n_S = " << sNevt[i] << "     n_B = " << bNevt[i] << endl;

    OUT << tableTS[i].Data() 
	<< "  & " << sNevt[i] << "  & " << Form("%4.2e",sEff[i])
	<< "  & " << bNevt[i]  << "  & " << Form("%4.2e",bEff[i]) << " \\\\" << endl;

    if ( cutsTS[i].Contains("chi2") )  { continue; }
    if ( cutsTS[i].Contains("iso") )    { continue; }

    sprintf(fcuts, "%s", cuts);
  }

  OUT << "\\hline" << endl;

  OUT << tf1 
      << "  &-   & " << Form("%4.2e",sEff1)
      << "  &-   & " << Form("%4.2e",bEff1) << " \\\\" << endl;

  OUT << tf2 
      << "  &-   & " << Form("%4.2e",sEff2)
      << "  &-   & " << Form("%4.2e",bEff2) << " \\\\" << endl;

  OUT << "\\hline" << endl;
  
  OUT.close();

  quickUL(fcuts, pre, fact1, fact2, vmass);


}


// ----------------------------------------------------------------------
double anaBmm::getCutEffS(const char *aCuts, const char *preCuts) { 

  double nevt(-1);

  fS[sgIndex]->cd(); 
  TH1D *qSG = (TH1D*)gROOT->FindObject("qSG"); 
  if (!qSG) qSG = new TH1D("qSG", "", 120, 4.8, 6.0);
  TTree *s = (TTree*)gFile->Get("events");

  s->Draw("mass>>qSG", Form("%s && %s", aCuts, preCuts ), "goff"); 
  nevt     = qSG->GetSumOfWeights();
  
  return nevt;
  
}

// ----------------------------------------------------------------------
double anaBmm::getCutEffB(const char *aCuts, const char *preCuts) { 

  double nevt(-1);

  fM[bgIndex]->cd(); 
  TH1D *qBG = (TH1D*)gROOT->FindObject("qBG"); 
  if (!qBG) qBG = new TH1D("qBG", "", 120, 4.8, 6.0);
  TTree *b = (TTree*)gFile->Get("events");
  
  b->Draw("mass>>qBG", Form("%s && %s", aCuts, preCuts), "goff"); 
  nevt     = qBG->GetSumOfWeights();
  
  return nevt;
  
}


// ----------------------------------------------------------------------
void anaBmm::quickUL(const char *aCuts, const char *preCuts, const char *f1Cut, const char *f2Cut, double mCut) {

  // -- Fix lumi normalisation scaling factor
  double sSF = fLumiD[0]/fLumiS[sgIndex];

  // -- Run on signal MC to determine efficiency 
  fS[sgIndex]->cd(); 
  TH1D *uSG = (TH1D*)gROOT->FindObject("uSG"); 
  if (!uSG) uSG = new TH1D("uSG", "", 120, 4.8, 6.0);
  TTree *s = (TTree*)gFile->Get("events");

  double sNorm = fNevtS[sgIndex];

  s->Draw("mass>>uSG", Form("goodKinematics&&goodL1&&goodHLT"), "goff"); 

  s->Draw("mass>>uSG", Form("goodL1&&goodHLT && %s", preCuts), "goff"); 
  double sNormF    =  uSG->GetSumOfWeights();
  double massRedSG = 0.;
    if ( mCut > -0.5 || fSgReduction < 0.) {
    fSgReduction = massReduction("uSG", "mysg", mCut);
    massRedSG = fSgReduction;
  } else {
    massRedSG = fSgReduction;
  }

  s->Draw("mass>>uSG", Form("goodL1&&goodHLT && %s", aCuts), "goff"); 
  double sCuts     = uSG->GetSumOfWeights();
  double sEff      = sCuts/sNorm;

  s->Draw("mass>>uSG", Form("goodL1&&goodHLT && %s && %s", preCuts, f1Cut), "goff"); 
  double s1Norm    = uSG->GetSumOfWeights();
  double sEff1     = s1Norm/sNormF;

  s->Draw("mass>>uSG", Form("goodL1&&goodHLT && %s && %s", preCuts, f2Cut), "goff"); 
  double s2Norm    = uSG->GetSumOfWeights();
  double sEff2     = s2Norm/sNormF;

  fNsg  =  sCuts * sEff1 * sEff2 * massRedSG * sSF;
  fNsgE = 0.;
  //fNsgE = TMath::Sqrt(fNsg);

  fEsg  = sEff * sEff1 * sEff2 * massRedSG;

  // -- Fix lumi normalisation scaling factor
  double bSF = fLumiD[0]/fLumiM[bgIndex];

  // -- Run on signal MC to determine efficiency 
  fM[bgIndex]->cd(); 
  TH1D *uBG = (TH1D*)gROOT->FindObject("uBG"); 
  if (!uBG) uBG = new TH1D("uBG", "", 120, 4.8, 6.0);
  TTree *b = (TTree*)gFile->Get("events");

  double bNorm = fNevtM[bgIndex];  

  b->Draw("mass>>uBG", Form("goodKinematics&&goodL1&&goodHLT"), "goff");

  b->Draw("mass>>uBG", Form("goodL1&&goodHLT && %s", preCuts), "goff"); 
  double bNormF    =  uBG->GetSumOfWeights(); 
  double massRedBG = 0.;
  if ( mCut > -0.5 || fBgReduction < 0.) {
    fBgReduction = massReduction("uBG", "mymc", mCut);
    massRedBG = fBgReduction;
  } else {
    massRedBG = fBgReduction;
  } 

  b->Draw("mass>>uBG", Form("goodL1&&goodHLT && %s", aCuts), "goff"); 
  double bCuts     = uBG->GetSumOfWeights();
  double bEff      = bCuts/bNorm;


  b->Draw("mass>>uBG", Form("goodL1&&goodHLT && %s && %s", preCuts, f1Cut), "goff"); 
  double b1Norm    = uBG->GetSumOfWeights();
  double bEff1     = b1Norm/bNormF;

  b->Draw("mass>>uBG", Form("goodL1&&goodHLT && %s && %s", preCuts, f2Cut), "goff"); 
  double b2Norm    = uBG->GetSumOfWeights();
  double bEff2     = b2Norm/bNormF;

  double ebg  = bEff * bEff1 * bEff2 * massRedSG;
  fNbg  =  bCuts * bEff1 * bEff2 * massRedBG * bSF;
  fNbgE =  fNbg * 1.6/(1.*TMath::Sqrt(bCuts));

  cout << "Calculate UL for n_S = " << fNsg << " and n_B = " << fNbg << endl;
  
  if ( (fNsg < 100 && fNbg < 1000) && (fNsg > 1e-20 && fNbg > 1e-20) ) {
    fUL = calculateUpperLimit();
  } else {
    cout << " ==> Cannot calculate UL for n_S = " << fNsg << " and n_B = " << fNbg << endl;
    fUL = 99999.;
  }

  cout << endl << endl << "mass red.: " << mCut
	 << "     e_S = " << massRedSG << "     e_B = " << massRedBG << endl;

  cout << endl << Form(" *** UL = %4.4e *** ", fUL)
       <<  " for n_S = " << fNsg << " (e_S = " << fEsg 
       <<  "),  n_B = " << fNbg << " (e_B = " << ebg <<  ")" << endl << endl
       << " cuts: \"" << aCuts << "\"  " 
       << " fact1: \"" << f1Cut << "\"  " 
       << " fact2: \"" << f2Cut << "\"  " 
       << "\" mass +/- " << mCut << "\"  " 
       << endl << endl;

  
  char filename[200];
  sprintf(filename, "cuts_table.tex"); 
  
  ofstream OUT(filename, ios::app);

  OUT << "\\hline" << endl;
  OUT << "$m_{\\mu\\mu} \\pm " << mCut
      << " \\gev$  &-  & " << Form("%4.2e", massRedSG)
      << "  &-  & " << Form("%4.2e", massRedBG) << " \\\\" << endl;

  OUT << "\\hline" << endl;
  OUT << "\\hline" << endl;
  
  OUT << "in $ " << fLumiD[0] << " \\invfb $"
      << "  & " << Form("%4.1f",fNsg) << "  & " << Form("%4.2e",fEsg)
      << "  & " << Form("%4.1f",fNbg)  << "  & " << Form("%4.2e", ebg) << " \\\\" << endl;

  OUT << "\\hline" << endl;
  OUT << "\\hline" << endl;

  OUT << "  \\multicolumn{5}{|c|}{\\alert{{\\bf" << Form(" %4.2e",fUL) << " (90\\\% C.L.)}}}  \\\\* " << endl;


  OUT.close();

}


// =======================================================================================================

// ---------------------------------------
// -- loop for stupid loop over all  cuts
// -- ptl1, chi2, lxy/sxy, iso, cosa
void anaBmm::loopOptimization(double pt) {

  double pTlo;
  TH1D *h = (TH1D*)fS[sgIndex]->Get("hcuts");
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
  double cut1(0.), cut2(0.),   cut3(0.),   cut4(0.); 
  char scut1[200], scut2[200], scut3[200], scut4[200];
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

  OUT.close();
}


// ----------------------------------------------------------------------
void anaBmm::optimizerNumbers(const char *precuts, const char *line, double &eff, double &sg, double &bg) {

  // -- Run on signal MC to determine efficiency 
  fS[sgIndex]->cd(); 
  TH1D *hSG = (TH1D*)gROOT->FindObject("hSG"); 
  if (!hSG) hSG = new TH1D("hSG", "", 50, 5., 6.);
  TTree *s = (TTree*)gFile->Get("events");
  s->Draw("mass>>hSG", Form("%s", precuts), "goff"); 
  double norm = hSG->GetSumOfWeights();
  s->Draw("mass>>hSG", Form("%s && %s", precuts, line), "goff"); 
  sg =  hSG->GetSumOfWeights();
  eff = sg/norm;

  // -- And now MC for background expectation
  fM[bgIndex]->cd(); 
  TH1D *hBG = (TH1D*)gROOT->FindObject("hBG"); 
  if (!hBG) hBG = new TH1D("hBG", "", 50, 5., 6.);
  TTree *b = (TTree*)gFile->Get("events");
  b->Draw("mass>>hBG", Form("%s && %s", precuts, line), "goff"); 
  bg =  hBG->GetSumOfWeights();

  //   cout << " sg: " << sg << " and norm: " << norm << "bg: " << bg << endl;

}


//===========================================================================================
// -- Utilities
//===========================================================================================

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
void anaBmm::setFilledHist(TH1 *h, Int_t color, Int_t fillcolor, Int_t fillstyle, Int_t line_width, Int_t line_style) {
  // Note: 3004, 3005 are crosshatches
  // ----- 1000       is solid
  //       kYellow    comes out gray on bw printers
  h->SetLineColor(color);     h->SetLineWidth(line_width);   h->SetLineStyle(line_style);    
  h->SetFillStyle(fillstyle); h->SetFillColor(fillcolor);
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
TString anaBmm::formatTex(double n, const char *name, const char *tag) {
  
  TString out("-");

  if ( isnan(n) ) {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{NaN } } }", name, tag);
  } else if ( n > 1.e10) {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name, tag, (texForm2(n)).Data());
  } else if ( n > 1.e4) {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name, tag, (texForm(n)).Data());
  } else if ( n > 100. ) {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{%4.0f } } }", name, tag, n);
  } else if ( n > 1. ) {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{%4.1f } } }", name, tag, n);
  } else if ( n > 1.e-2) {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{%4.2f } } }", name, tag, n);
  } else if ( n > 1.e-3) {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name, tag, n);
  } else if ( n > 1.e-9 ) {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name, tag, (texForm(n)).Data());
  } else if ( n > 1.e-19 ){
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name, tag, (texForm2(n)).Data());
  } else {
    out.Form("\\vdef{%s:%s}   {\\ensuremath{{0.0 } } }", name, tag);
  }
  
  return out;
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
  } else if ( exp == 1 ) {
    factor = TString("(x 10)"); 
  } else if ( exp > -100 ) { 
    factor = TString(Form("(10^{%i})", exp));
  } else {
    factor = TString("");
  }

  return factor;

}

// ----------------------------------------------------------------------
double anaBmm::getDGError(double a, double ae, double b, double be, double c, double ce, double d, double de) {

  double error = TMath::Sqrt(2*(a*c*c)*(a*c*c)*(ae*ae)/((c*c+d*d)*(c*c+d*d)) +
			       (b*d*d)*(b*d*d)*(be*be)/((c*c+d*d)*(c*c+d*d)) +
			       (c*a*a*(d*d-c*c))*(c*a*a*(d*d-c*c))*(ce*ce)/((c*c+d*d)*(c*c+d*d)*(c*c+d*d)*(c*c+d*d)) +
			       (d*b*b*(c*c-d*d))*(d*b*b*(c*c-d*d))*(de*de)/((c*c+d*d)*(c*c+d*d)*(c*c+d*d)*(c*c+d*d)) 
			     );

  return error;
}


// ----------------------------------------------------------------------
void anaBmm::emptyBinError(TH1 *h) {

  for (int i = 0; i < h->GetNbinsX(); i++ ) {
    if ( h->GetBinContent(i) == 0 ) {
      h->SetBinError(i, 1);
    }
  }
}


//===========================================================================================
// -- Other
//===========================================================================================

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
int anaBmm::wait() {
  cout << " Continue [<RET>|q]?  "; 
  char x;
  x = getchar();
  if ((x == 'q') || (x == 'Q')) return 1;
  return 0;
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


//===========================================================================================
// -- Assignments
//===========================================================================================

int anaBmm::checkIf(int mc, const char *sel) {
  
  int accept(0);
  // -- only accepts rare backgrounds (signal)
  if ( !strcmp(sel,"r0") ) {
    if ( !strcmp(fTypeM[mc].Data(), "rmc") ) {
      accept = 1;
    }
    
  // -- only accepts backgrounds & rare backgrounds (signal)
  } else if ( !strcmp(sel,"c0") ) {
    accept = 1;
    if ( !strcmp(fTypeM[mc].Data(), "nmc") ) {
      accept = 0;
    }
    
  // -- sum-up sub-groups of rare backgrounds (2mu+, 2mId, etc)
  } else {
    if ( !strcmp(getSubGroup(fSignM[mc]).Data(), sel) ) {
      accept = 1;	    
    }
  }

  return accept;
}

// -----------------------------------------------------

int anaBmm::findIndex(const char *filename) {


  int index(-1);


  for (int i = 0; i < nSg; ++i) {
    if ( !strcmp(filename, fFileS[i].Data()) ) {
      index = i; goto rtn;
    }
  }

  for (int i = 0; i < nMc; ++i) {
    if ( !strcmp(filename, fFileM[i].Data()) ) {
      index = i; goto rtn;
    }
  }

  for (int i = 0; i < nDa; ++i) {
    if ( !strcmp(filename, fFileD[i].Data()) ) {
      index = i; goto rtn;
    }
  }

  rtn:;
  return index;
}


// -----------------------------------------------------

void anaBmm::getSignature(TString signIn, TString &signOut, TString &signOut2) {
 
 signOut  = signIn;
 signOut2 = signIn;
 
 TString sign[]       = {  TString("bsmumu")
			 , TString("bpjpsikp")
			 , TString("bbbar")

			 , TString("BB")
			 , TString("CC")
			 , TString("NonPrompt")

			 , TString("bd2pi")
			 , TString("bdpik")
			 , TString("bdpimunu")
			 , TString("bdmumupi0")

			 , TString("bs2pi")
			 , TString("bskk")
			 , TString("bskpi")
			 , TString("bskmunu")
			 , TString("bsmumug")
			 , TString("bsmumupi0")

			 , TString("bu3munu")

			 , TString("bc3munu")
			 , TString("bcjpmunu")

			 , TString("lbpk")
			 , TString("lbppi")

			 , TString("qcd")

                        };

 TString titles[]    = {   TString("\\bsmm")
			 , TString("\\bupsimmk")
			 , TString("\\bbbar")

			 , TString("\\bbmumuX")
			 , TString("\\ccmumuX")
			 , TString("\\bpsimmX")

			 , TString("\\bdpipi")
			 , TString("\\bdpik")
			 , TString("\\bdpimunu")
			 , TString("\\bdmumupz")

			 , TString("\\bspipi")
			 , TString("\\bskk")
			 , TString("\\bspik ")
			 , TString("\\bskmunu")
			 , TString("\\bsmumug")
			 , TString("\\bsmumupz")

			 , TString("\\butrmunu")

			 , TString("\\bctrmunu")
			 , TString("\\bcpsimunu")

			 , TString("\\lbpk")
			 , TString("\\lbppi")

			 , TString("M(hh)")
                         };


 TString titles2[]    = {  TString("B_{s}^{0} #rightarrow #mu^{+} #mu^{-}")
			 , TString("B^{#pm} #rightarrow J/#psi K^{#pm}")
			 , TString("b#bar{b}")

			 , TString("b#bar{b} #rightarrow #mu^{+} #mu^{-} X")
			 , TString("c#bar{c} #rightarrow #mu^{+} #mu^{-} X")
			 , TString("b #rightarrow J/#psi (#rightarrow #mu^{+} #mu^{-}) X")

			 , TString("B^{0} #rightarrow #pi^{+} #pi^{-}")
			 , TString("B^{0} #rightarrow #pi^{-} K^{+}")
			 , TString("B^{0} #rightarrow #pi^{-} #mu^{+} #nu")
			 , TString("B^{0} #rightarrow #mu^{+} #mu^{-} #pi^{0}")

			 , TString("B_{s}^{0} #rightarrow #pi^{+} #pi^{-}")
			 , TString("B_{s}^{0} #rightarrow K^{+} K^{-}")
			 , TString("B_{s}^{0} #rightarrow K^{-} #pi^{+}")
			 , TString("B_{s}^{0} #rightarrow K^{-} #mu^{+} #nu")
			 , TString("B_{s}^{0} #rightarrow #mu^{+} #mu^{-} #gamma")
			 , TString("B_{s}^{0} #rightarrow #mu^{+} #mu^{-} #pi^{0}")

			 , TString("B ^{+} #rightarrow #mu^{+} #mu^{-} #mu^{+} #nu")

			 , TString("B_{c}^{+} #rightarrow #mu^{+} #mu^{-} #mu^{+} #nu")
			 , TString("B_{c}^{+} #rightarrow J/#psi #mu^{+} #nu")

			 , TString("#Lambda_{b} #rightarrow p K^{-}")
			 , TString("#Lambda_{b} #rightarrow p #pi^{-}")

			 , TString("M(hh)")

                         };
 


 for (int i = 0; i < 22; i++) {

   if ( !strcmp(signIn.Data(), sign[i].Data()) ) {
   
     signOut  = TString(titles[i]);
     signOut2 = TString(titles2[i]);
     
   }
 }

}

// ----------------------------------------------------------------------
TString anaBmm::getSubGroup(TString signIn) {


  TString rtype;

  TString sign[]      = {  TString("bsmumu")
			 , TString("bpjpsikp")
			 , TString("bbbar")

			 , TString("BB")
			 , TString("CC")
			 , TString("NonPrompt")

			 , TString("bd2pi")
			 , TString("bdpik")
			 , TString("bdpimunu")
			 , TString("bdmumupi0")

			 , TString("bs2pi")
			 , TString("bskk")
			 , TString("bskpi")
			 , TString("bskmunu")
			 , TString("bsmumug")
			 , TString("bsmumupi0")

			 , TString("bu3munu")

			 , TString("bc3munu")
			 , TString("bcjpmunu")

			 , TString("lbpk")
			 , TString("lbppi")

			 , TString("qcd")

                        };

  TString type[]      = {   TString("2mu")
			  , TString("2mu")
			  , TString("2mu")

			  , TString("2mu")
			  , TString("2mu")
			  , TString("2mu")
		   
			  , TString("2mId")
			  , TString("2mId")
			  , TString("mIdMu+")
			  , TString("2mu+")
		   
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
 

 for (int i = 0; i < 21; i++) {

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

			 , TString("BB")
			 , TString("CC")
			 , TString("NonPrompt")

			 , TString("bd2pi")
			 , TString("bdpik")
			 , TString("bdpimunu")
			 , TString("bdmumupi0")

			 , TString("bs2pi")
			 , TString("bskk")
			 , TString("bskpi")
			 , TString("bskmunu")
			 , TString("bsmumug")
			 , TString("bsmumupi0")

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

			  , fMu*fMu
			  , fMu*fMu
			  , fMu*fMu
		   
			  , fPi*fPi
			  , fPi*fKa
			  , fPi*fMu
			  , fMu*fMu
	   
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

 
 for (int i = 0; i < 21; i++) {

   if ( !strcmp(signIn.Data(), sign[i].Data()) ) {
    
     id = eff[i];
  
   }
 }
 
 return id;

}

