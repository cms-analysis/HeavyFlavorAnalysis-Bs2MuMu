#include "initFunc.hh"

#include "TROOT.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>


ClassImp(initFunc)

using std::string;
using std::cout;
using std::endl;


// ----------------------------------------------------------------------
double f_expo(double *x, double *par) {
  return par[0]*TMath::Exp(-x[0]*par[1]);
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
double f_gauss2(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> fraction in second gaussian
  // par[4] -> mean of second gaussian
  // par[5] -> sigma of second gaussian
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
double f_gauss2c(double *x, double *par) {
  // constrained to have the same mean in the second gaussian
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> fraction in second gaussian
  // par[4] -> sigma of second gaussian
  Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.); 
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[4] > 0.) {
    arg2 = (x[0] - par[1]) / par[4];
    fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
  }
  Double_t fitval = fitval1 + fitval2;
  return fitval;
}


// ----------------------------------------------------------------------
double f_pol0(double *x, double *par) {
  return par[0]; 
}


// ----------------------------------------------------------------------
double f_pol1(double *x, double *par) {
  return par[0] + par[1]*x[0]; 
}


// ----------------------------------------------------------------------
double f_pol2(double *x, double *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0]; 
}


// ----------------------------------------------------------------------
double f_pol0_BsBlind(double *x, double *par) {
  // FIXME fixed limits!!!
  if (x[0] >= 5.2 && x[0] <= 5.45) { 
    TF1::RejectPoint();
    return 0;
  }
  return par[0]; 
}


// ----------------------------------------------------------------------
double f_pol1_BsBlind(double *x, double *par) {
  // FIXME fixed limits!!!
  if (x[0] >= 5.2 && x[0] <= 5.45) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0]; 
}


// ----------------------------------------------------------------------
// pol1 and gauss 
double f_pol1_gauss(double *x, double *par) {
  //   par[0] = normalization of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = par 0 of pol1
  //   par[4] = par 1 of pol1
  return  (f_pol1(x, &par[3]) + f_gauss(x, &par[0]));
}


// ----------------------------------------------------------------------
// pol1 and Gauss 
double f_pol1_Gauss(double *x, double *par) {
  //   par[0] = area of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = par 0 of pol1
  //   par[4] = par 1 of pol1
  return  (f_pol1(x, &par[3]) + f_Gauss(x, &par[0]));
}

// ----------------------------------------------------------------------
// pol1 and gauss2 
double f_pol1_gauss2c(double *x, double *par) {
  //   par[0] = norm of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = fraction in second gaussian
  //   par[4] = sigma of gaussian
  //   par[5] = par 0 of pol1
  //   par[6] = par 1 of pol1
  return  (f_pol1(x, &par[5]) + f_gauss2c(x, &par[0]));
}



// ----------------------------------------------------------------------
initFunc::initFunc() {
  cout << "ctor initFunc" << endl;
  fLo = 99.; 
  fHi = -99.;
}


// ----------------------------------------------------------------------
initFunc::~initFunc() {
  cout << "dtor initFunc" << endl;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol1(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0; 
  }
  TF1 *f = new TF1("f1", f_pol1, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 2);
  f->SetParNames("offset", "slope"); 			   
  //   cout << "Created f1 from " << h->GetBinLowEdge(1) << " to " << h->GetBinLowEdge(h->GetNbinsX()+1) << endl;
  
  int lbin(1), hbin(h->GetNbinsX()), EDG(4), NB(EDG+1); 
  if (fLo < fHi) {
    lbin = h->FindBin(fLo); 
    hbin = h->FindBin(fHi); 
  }

  //   cout << "lbin = " << lbin << " hbin = " << hbin << endl;
  
  double ylo = h->Integral(lbin, lbin+EDG)/NB; 
  double yhi = h->Integral(hbin-EDG, hbin)/NB;
  double dx = h->GetBinLowEdge(hbin) - h->GetBinLowEdge(lbin);

  //   cout << "ylo = " << ylo << " yhi = " << yhi << " dx = " << dx << endl;

  double p1 = (yhi-ylo)/dx;
  double p0 = ylo - p1*h->GetBinLowEdge(1);
  f->SetParameters(p0, p1); 
  //   cout << "  initialized to " << p0 << " and " << p1 << endl;
  //   cout << "  -> Integrals:  func = " 
  //        << f->Integral(h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1))/h->GetBinWidth(1)
  //        << " histogram = " 
  //        << h->GetSumOfWeights() 
  //        << endl;
  return f; 
}


// ----------------------------------------------------------------------
TF1* initFunc::pol1BsBlind(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0; 
  }
  TF1 *f = new TF1("f1", f_pol1_BsBlind, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 2);
  f->SetParNames("offset", "slope"); 			   
  //   cout << "Created f1 from " << h->GetBinLowEdge(1) << " to " << h->GetBinLowEdge(h->GetNbinsX()+1) << endl;
  
  int lbin(1), hbin(h->GetNbinsX()), EDG(4), NB(EDG+1); 
  if (fLo < fHi) {
    lbin = h->FindBin(fLo); 
    hbin = h->FindBin(fHi); 
  }

  //   cout << "lbin = " << lbin << " hbin = " << hbin << endl;
  
  double ylo = h->Integral(lbin, lbin+EDG)/NB; 
  double yhi = h->Integral(hbin-EDG, hbin)/NB;
  double dx = h->GetBinLowEdge(hbin) - h->GetBinLowEdge(lbin);

  //   cout << "ylo = " << ylo << " yhi = " << yhi << " dx = " << dx << endl;

  double p1 = (yhi-ylo)/dx;
  double p0 = ylo - p1*h->GetBinLowEdge(1);
  f->SetParameters(p0, p1); 
  //   cout << "  initialized to " << p0 << " and " << p1 << endl;
  //   cout << "  -> Integrals:  func = " 
  //        << f->Integral(h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1))/h->GetBinWidth(1)
  //        << " histogram = " 
  //        << h->GetSumOfWeights() 
  //        << endl;
  return f; 
}


// ----------------------------------------------------------------------
TF1* initFunc::pol0(TH1 *h) {

  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0; 
  }
  TF1 *f = new TF1("f1", f_pol0, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()), 1);
  f->SetParNames("constant"); 			   

  int lbin(1), hbin(h->GetNbinsX()), EDG(4), NB(EDG+1); 
  if (fLo < fHi) {
    lbin = h->FindBin(fLo); 
    hbin = h->FindBin(fHi); 
  }
  
  double ylo = h->Integral(lbin, lbin+EDG)/NB; 
  double yhi = h->Integral(hbin-EDG, hbin)/NB;
  double p0 = 0.5*(yhi+ylo);
  f->SetParameter(0, p0); 
    
  return f; 

}


// ----------------------------------------------------------------------
TF1* initFunc::pol0BsBlind(TH1 *h) {

  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0; 
  }
  TF1 *f = new TF1("f1", f_pol0_BsBlind, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()), 1);
  f->SetParNames("constant"); 			   

  int lbin(1), hbin(h->GetNbinsX()), EDG(4), NB(EDG+1); 
  if (fLo < fHi) {
    lbin = h->FindBin(fLo); 
    hbin = h->FindBin(fHi); 
  }
  
  double ylo = h->Integral(lbin, lbin+EDG)/NB; 
  double yhi = h->Integral(hbin-EDG, hbin)/NB;
  double p0 = 0.5*(yhi+ylo);
  f->SetParameter(0, p0); 
    
  return f; 

}


// ----------------------------------------------------------------------
TF1* initFunc::pol1gauss(TH1 *h) {

  TF1 *f = new TF1("f1", f_pol1_gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()), 5);
  f->SetParNames("normalization", "peak", "sigma", "constant", "slope"); 			   

  int lbin(1), hbin(h->GetNbinsX()); 
  if (fLo < fHi) {
    lbin = h->FindBin(fLo); 
    hbin = h->FindBin(fHi); 
  }
  
  double p0, p1; 
  initPol1(p0, p1, h);
  double A   = 0.5*p1*(fHi*fHi - fLo*fLo) + p0*(fHi - fLo);

  double g0 = (h->Integral(lbin, hbin) - A/h->GetBinWidth(1))*h->GetBinWidth(1);  
  double g1 = 5.3;
  double g2 = 0.04;
  
  f->SetParameters(g0, g1, g2, p0, p1); 
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7); 
  f->ReleaseParameter(1);     f->SetParLimits(1, 5.2, 5.45); 
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.010, 0.080); 

  cout << "g0 = " << g0 << " g1 = " << g1 << " g2 = " << g2 << " p0 = " << p0 << " p1 = " << p1 << endl;
  return f; 

}


// ----------------------------------------------------------------------
TF1* initFunc::pol1Gauss(TH1 *h) {

  TF1 *f = new TF1("f1", f_pol1_Gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()), 5);
  f->SetParNames("area", "peak", "sigma", "constant", "slope"); 			   

  int lbin(1), hbin(h->GetNbinsX()); 
  if (fLo < fHi) {
    lbin = h->FindBin(fLo); 
    hbin = h->FindBin(fHi); 
  }

  cout << "fLo: " << fLo << " fHi: " << fHi << " lbin: " << lbin << " hbin: " << hbin << endl;
  
  double p0, p1; 
  initPol1(p0, p1, h);
  double A   = 0.5*p1*(fHi*fHi - fLo*fLo) + p0*(fHi - fLo);

  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);  
  double g1 = 5.3;
  double g2 = 0.04;

  cout << "g0 = " << g0 << " g1 = " << g1 << " g2 = " << g2 << " p0 = " << p0 << " p1 = " << p1 << endl;
  f->SetParameters(g0, g1, g2, p0, p1); 
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7); 
  f->ReleaseParameter(1);     f->SetParLimits(1, 5.2, 5.45); 
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.010, 0.080); 

  return f; 

}


// ----------------------------------------------------------------------
TF1* initFunc::pol1gauss2c(TH1 *h, double peak, double sigma) {

  TF1 *f = (TF1*)gROOT->FindObject("f1_pol1_gauss2c"); 
  if (f) delete f; 
  f = new TF1("f1_pol1_gauss2c", f_pol1_gauss2c, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()), 7);
  f->SetParNames("norm", "peak", "sigma", "fraction", "sigma2", "constant", "slope"); 			   
  f->SetLineColor(kBlue); 
  f->SetLineWidth(2); 

  int lbin(1), hbin(h->GetNbinsX()); 
  if (fLo < fHi) {
    lbin = h->FindBin(fLo); 
    hbin = h->FindBin(fHi); 
  }

  cout << "fLo: " << fLo << " fHi: " << fHi << " lbin: " << lbin << " hbin: " << hbin 
       << " sigma = " << sigma << " peak = " << peak 
       << endl;
  
  double p0, p1; 
  initPol1(p0, p1, h);
  cout << "p0: " << p0 << " p1: " << p1 << endl;
  double A   = 0.5*p1*(fHi*fHi - fLo*fLo) + p0*(fHi - fLo);

  double sqrt2pi = 2.506628275;
  double gInt    = h->Integral(lbin, hbin) - A;  
  cout << "A: " << A << endl;
  cout << "gInt: " << gInt << endl;
  cout << "h->Integral(): " << h->Integral(lbin, hbin) << endl;

  double gaussN  = gInt/(2.*sqrt2pi*sigma)*h->GetBinWidth(1);

  cout << "initFunc> gaussN = " << gaussN << " peak = " << peak << " sigma = " << sigma << " p0 = " << p0 << " p1 = " << p1 << endl;
  f->SetParameters(gaussN, peak, sigma, 0.2, 1.8*sigma, p0, p1); 
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7); 
  f->ReleaseParameter(1);     f->SetParLimits(1, peak-0.1, peak+0.1); 
  f->ReleaseParameter(2);     f->SetParLimits(2, sigma*0.4, sigma*1.3); 
  f->ReleaseParameter(3);     f->SetParLimits(3, 0.01, 2.0); 
  f->ReleaseParameter(4);     f->SetParLimits(4, sigma*1.2, sigma*10.0); 

  return f; 


}


// ----------------------------------------------------------------------
void initFunc::initPol1(double &p0, double &p1, TH1 *h) {

  int EDG(4), NB(EDG+1); 
  int lbin(1), hbin(h->GetNbinsX()); 
  if (fLo < fHi) {
    lbin = h->FindBin(fLo); 
    hbin = h->FindBin(fHi); 
  }
  
  double dx  = fHi - fLo;
  double ylo = h->Integral(lbin, lbin+EDG)/NB; 
  double yhi = h->Integral(hbin-EDG, hbin)/NB;

  cout << "ylo: " << ylo << " yhi: " << yhi << endl;

  p1  = (yhi-ylo)/dx;
  p0  = ylo - p1*fLo;
}

