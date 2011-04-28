#include "AnalysisDistribution.hh"

#include <iostream>
#include <iomanip>

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "initFunc.hh"

#include "TMath.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

ClassImp(AnalysisDistribution)

using std::string;
using std::cout;
using std::endl;


// ----------------------------------------------------------------------
// expo and Gauss 
double f_epaG(double *x, double *par) {
  //   par[0] = area of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = par 0 of expo
  //   par[4] = par 1 of expo
  //   par[5] = par 0 of pol1
  //   par[6] = par 1 of pol1
  return  (f_p1(x, &par[5]) + f_expo(x, &par[3]) + f_Gauss(x, &par[0]));
}

// ----------------------------------------------------------------------
// expo and double Gauss 
double f_ea2G(double *x, double *par) {
  // par[0] -> area
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> fraction in second gaussian
  // par[4] -> mean
  // par[5] -> sigma
  // par[6] = par 0 of expo
  // par[7] = par 1 of expo
  return  (f_expo(x, &par[6]) + f_2G(x, &par[0]));
}

// ----------------------------------------------------------------------
AnalysisDistribution::AnalysisDistribution(const char *name, const char *title, int nbins, double lo, double hi) {
  int NBINS(120);
  fSigLo = fSigHi = fBg1Lo = fBg1Hi = fBg2Lo = fBg2Hi = 0.0;
  fMassLo = 4.8;
  fMassHi = 6.0;
  fMassPeak = -1.; 
  fMassSigma = -1.; 

  fVerbose = 0; 

  fpIF = new initFunc(); 

  string massbin[3]; 
  massbin[0] = "signal";
  massbin[1] = "sideband";
  massbin[2] = "all";
  for (int i = 0; i < 3; ++i) {
    hSi[i] = new TH1D(Form("%sSi%d", name, i), Form("%s, single cut, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hSi[i]->SetXTitle(title); 
      }
  for (int i = 0; i < 3; ++i) {
    hAo[i] = new TH1D(Form("%sAo%d", name, i), Form("%s, all other cuts, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hAo[i]->SetXTitle(title); 
  }
  for (int i = 0; i < 3; ++i) {
    hNm[i] = new TH1D(Form("%sNm%d", name, i), Form("%s, n-1 cuts, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hNm[i]->SetXTitle(title); 
  }
  for (int i = 0; i < 3; ++i) {
    hCu[i] = new TH1D(Form("%sCu%d", name, i), Form("%s, cumulative cuts, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hCu[i]->SetXTitle(title); 
  }
  for (int i = 0; i < 3; ++i) {
    hHLT[i] = new TH1D(Form("%sHLT%d", name, i), Form("%s, after HLT, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hHLT[i]->SetXTitle(title); 
  }
  for (int i = 0; i < 3; ++i) {
    hPresel[i] = new TH1D(Form("%sPresel%d", name, i), Form("%s, after Presel, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hPresel[i]->SetXTitle(title); 
  }

  hMassSi    = new TH1D(Form("%sMassSi", name), Form("%sMassSi", name), NBINS, fMassLo, fMassHi); 
  hMassAo    = new TH1D(Form("%sMassAo", name), Form("%sMassAo", name), NBINS, fMassLo, fMassHi); 
  hMassNm    = new TH1D(Form("%sMassNm", name), Form("%sMassNm", name), NBINS, fMassLo, fMassHi); 
  hMassCu    = new TH1D(Form("%sMassCu", name), Form("%sMassCu", name), NBINS, fMassLo, fMassHi); 
  hMassHLT   = new TH1D(Form("%sMassHLT", name), Form("%sMassHLT", name), NBINS, fMassLo, fMassHi); 
  hMassPresel= new TH1D(Form("%sMassPresel", name), Form("%sMassPresel", name), NBINS, fMassLo, fMassHi); 

  hMassAll   = new TH1D(Form("%sMassAll", name), Form("%sMassALL", name), NBINS, fMassLo, fMassHi); 
  hMassBGL   = new TH1D(Form("%sMassBGL", name), Form("%sMassBGL", name), NBINS, fMassLo, fMassHi); 
  hMassSG    = new TH1D(Form("%sMassSG", name), Form("%sMassSG", name), NBINS, fMassLo, fMassHi);  
  hMassBGH   = new TH1D(Form("%sMassBGH", name), Form("%sMassGBH", name), NBINS, fMassLo, fMassHi); 


  fF0 = new TF1("ADf0", f_p1, 0., 7., 2);         fF0->SetParNames("Offset", "Slope"); 			   
  fF1 = new TF1("ADf1", f_expo, 0., 7., 2);       fF1->SetParNames("Norm", "Exp"); 			   

  fP1  = new TF1("ADp1", f_p1, 0., 7., 2);        fP1->SetParNames("Offset", "Slope"); 			   
  fPG1 = new TF1("ADpg1", f_p1aG, 0., 7., 5);     fPG1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 
  fEG1 = new TF1("ADeg1", f_eaG, 0., 7., 5);      fEG1->SetParNames("Area", "Peak", "Sigma", "Norm", "Exp"); 
  fEG2 = new TF1("ADeg2", f_ea2G, 0., 7., 8);     fEG2->SetParNames("Area", "Peak", "Sigma", "Fraction2", "Peak2", "Sigma2", "Norm", "Exp"); 
  fEPG = new TF1("ADepg", f_epaG, 0., 7., 7);     fEPG->SetParNames("Area", "Peak", "Sigma", "Norm", "Exp", "Offset", "Slope");             
  
}



// ----------------------------------------------------------------------
AnalysisDistribution::AnalysisDistribution(const char *name, double SigLo, double SigHi, double Bg1Lo, double Bg1Hi, double Bg2Lo, double Bg2Hi) {

  fMassLo = 4.8;
  fMassHi = 6.0;
  fpIF = new initFunc(); 

  hMassSi    = (TH1D*)gDirectory->Get(Form("%sMassSi", name)); 
  hMassAo    = (TH1D*)gDirectory->Get(Form("%sMassAo", name)); 
  hMassNm    = (TH1D*)gDirectory->Get(Form("%sMassNm", name)); 
  hMassCu    = (TH1D*)gDirectory->Get(Form("%sMassCu", name)); 
  hMassHLT   = (TH1D*)gDirectory->Get(Form("%sMassHLT", name)); 
  hMassPresel= (TH1D*)gDirectory->Get(Form("%sMassPresel", name)); 


  hMassBGL    = (TH1D*)gDirectory->Get(Form("%sMassBGL", name)); 
  hMassBGH    = (TH1D*)gDirectory->Get(Form("%sMassBGH", name)); 
  hMassSG     = (TH1D*)gDirectory->Get(Form("%sMassSG", name)); 

  for (int i = 0; i < 3; ++i) {
    hSi[i]    = (TH1D*)gDirectory->Get(Form("%sSi%d", name, i)); 
    hAo[i]    = (TH1D*)gDirectory->Get(Form("%sAo%d", name, i)); 
    hNm[i]    = (TH1D*)gDirectory->Get(Form("%sNm%d", name, i)); 
    hCu[i]    = (TH1D*)gDirectory->Get(Form("%sCu%d", name, i)); 
    hHLT[i]   = (TH1D*)gDirectory->Get(Form("%sHLT%d", name, i)); 
    hPresel[i]= (TH1D*)gDirectory->Get(Form("%sPresel%d", name, i)); 
  }

  fSigLo = SigLo; 
  fSigHi = SigHi; 
  fBg1Lo = Bg1Lo; 
  fBg1Hi = Bg1Hi; 
  fBg2Lo = Bg2Lo;
  fBg2Hi = Bg2Hi;

  fF0 = new TF1("ADf0", f_p1, 0., 7., 2);         fF0->SetParNames("Offset", "Slope"); 			   
  fF1 = new TF1("ADf1", f_expo, 0., 7., 2);  	  fF1->SetParNames("Norm", "Exp"); 			   

  fP1  = new TF1("ADp1", f_p1, 0., 7., 2);        fP1->SetParNames("Offset", "Slope"); 			   				    
  fPG1 = new TF1("ADpg1", f_p1aG, 0., 7., 5);     fPG1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 			    
  fEG1 = new TF1("ADeg1", f_eaG, 0., 7., 5);      fEG1->SetParNames("Area", "Peak", "Sigma", "Norm", "Exp"); 				    
  fEG2 = new TF1("ADeg2", f_ea2G, 0., 7., 8);     fEG2->SetParNames("Area", "Peak", "Sigma", "Fraction2", "Peak2", "Sigma2", "Norm", "Exp");
  fEPG = new TF1("ADepg", f_epaG, 0., 7., 7);     fEPG->SetParNames("Area", "Peak", "Sigma", "Norm", "Exp", "Offset", "Slope");             
  
}



// ----------------------------------------------------------------------
AnalysisDistribution::~AnalysisDistribution() {
  if (fF0) delete fF0; 
  if (fF1) delete fF1; 
  if (fP1) delete fP1; 
  if (fPG1) delete fPG1; 
  if (fEG1) delete fEG1; 
}

// ----------------------------------------------------------------------
void AnalysisDistribution::setSigWindow(double lo, double hi) {
  fSigLo = lo;  
  fSigHi = hi;
}

// ----------------------------------------------------------------------
void AnalysisDistribution::setBg1Window(double lo, double hi) {
  fBg1Lo = lo;  
  fBg1Hi = hi;
}

// ----------------------------------------------------------------------
void AnalysisDistribution::setBg2Window(double lo, double hi) {
  fBg2Lo = lo;  
  fBg2Hi = hi;
}

// ----------------------------------------------------------------------
void AnalysisDistribution::setAnalysisCuts(AnalysisCuts *p, const char *cutname) {
  fpAnaCuts = p; 
  fCutName  = cutname;
  fCutIdx   = fpAnaCuts->getIndex(fCutName.c_str()); 
  if (fCutIdx < 0) {
    cout << "xx> AnalysisDistribution: ERROR " 
	 << fCutName << " not found in AnalysisCuts" 
	 << endl;
  }
  fHLTIdx   = fpAnaCuts->getIndex("fGoodHLT"); 
} 


// ----------------------------------------------------------------------
double AnalysisDistribution::fitMass(TH1 *h1, double &error, int mode) {

  if (0 == h1) {
    cout << "no histogram provided. STOP." << endl;
    return -1.;
  }

  TF1 *f1; 
  if (0 == mode) {
    // -- just count the number of events in the histogram
    double n = h1->GetSumOfWeights(); 
    error = TMath::Sqrt(n); 
    h1->Draw(); 
    return n; 
  } else if (1 == mode) {
    // -- blinded signal box 
    f1 = fP1; 
    setFunctionParameters(f1, h1, mode); 
    h1->Fit(f1, "", "", fMassLo, fMassHi); 
    double p0  = f1->GetParameter(0); 
    double p0E2= f1->GetParError(0)*f1->GetParError(0); 
    double p1  = f1->GetParameter(1); 
    double p1E2= f1->GetParError(1)*f1->GetParError(1); 
    double dx  = fSigHi - fSigLo;
    double d2x = fSigHi*fSigHi - fSigLo*fSigLo;
    double yield = f1->Integral(fSigLo, fSigHi)/h1->GetBinWidth(1);
    double yieldAnalytical = (p0*dx + 0.5*d2x*p1)/h1->GetBinWidth(1);
    error = TMath::Sqrt((dx*dx*p0E2 + 0.25*d2x*d2x*p1E2)/h1->GetBinWidth(1));
    cout << "yield from TF1 integral: " << yield << " vs. " << yieldAnalytical << " +/- " << error << endl;
    return yield;
  } else if (10 == mode) {
    // -- One Gaussian plus pol1
    f1 = fPG1;
    setFunctionParameters(f1, h1, mode); 
    h1->Fit(f1, "", "", fMassLo, fMassHi); 
    error = f1->GetParError(0)/h1->GetBinWidth(1); 
    return f1->GetParameter(0)/h1->GetBinWidth(1); 
  } else if (11 == mode) {
    // -- Double Gaussian with pol1
    fpIF->fLo = fMassLo; 
    fpIF->fHi = fMassHi;
    double peak = (fMassPeak>0.?fMassPeak:5.3);
    double sigma = (fMassSigma>0.?fMassSigma:0.04);
    f1 = fpIF->pol1gauss2c(h1, peak, sigma);
    TFitResultPtr r;
    r = h1->Fit(f1, "ls", "", fMassLo, fMassHi); 
    double hlimit(0.), llimit(0.); 
    f1->GetParLimits(3, llimit, hlimit);
    double relError = f1->GetParError(3)/f1->GetParameter(3); 
    double limit = TMath::Abs(f1->GetParameter(3) - llimit); 
    if (TMath::Abs(f1->GetParameter(3) - hlimit) < limit) limit = TMath::Abs(f1->GetParameter(3) - hlimit); 
    if (limit < relError || relError > 0.1) {
      cout << "%%%%% REFITTING %%%%% " << endl;
      f1->SetParameters(f1->GetParameters()); 
      f1->SetParameter(0, h1->GetMaximum()); 
      f1->FixParameter(3, 0.);
      f1->FixParameter(4, 0.);
      r = h1->Fit(f1, "ls", "", fMassLo, fMassHi); 
    }
    f1->SetParameter(5, 0.); 
    f1->SetParameter(6, 0.); 
    peak   = f1->GetParameter(1); 
    sigma  = TMath::Sqrt(f1->CentralMoment(2, fMassLo, fMassHi)); 
    double yield  = f1->Integral(peak-3.*sigma, peak+3.*sigma)/h1->GetBinWidth(1); 
    cout << "integral:       " << yield << endl;
    double ierror = f1->IntegralError(peak-3.*sigma, peak+3.*sigma, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/h1->GetBinWidth(1); 
    cout << "integral error: " << ierror << endl;
    double yieldE = f1->GetParError(0)/f1->GetParameter(0);
    error = yieldE*yield; 
    cout << "yield error: " << error << endl;
    error = ierror; 
    return yield; 
  } else if (12 == mode) {
    // -- One Gaussian plus expo
    f1 = fEG1;
    setFunctionParameters(f1, h1, mode); 
    h1->Fit(f1, "", "", fMassLo, fMassHi); 
    error = f1->GetParError(0)/h1->GetBinWidth(1); 
    return f1->GetParameter(0)/h1->GetBinWidth(1); 
  } else if (13 == mode) {
    // -- one Gaussian plus expo plus pol1
    f1 = fEPG;
    setFunctionParameters(f1, h1, mode); 
    h1->Fit(f1, "", "", fMassLo, fMassHi); 
    error = f1->GetParError(0)/h1->GetBinWidth(1); 
    return f1->GetParameter(0)/h1->GetBinWidth(1); 
  } else {
    return -2.;
  }

  return -3.; 
}


// ----------------------------------------------------------------------
// -- mode: 
//           1 pol1+Gauss
//           2 expo+Gauss
void AnalysisDistribution::setFunctionParameters(TF1 *f1, TH1 *h, int mode) {
  const int EDG(4), NB(EDG+1); 
  double ylo(0.), yhi(0.), dx(0.);
  int lo(0), hi(0); 
  double xlo(0.), xhi(0.); 
  double p0(0.), p1(0.);

  int lbin = h->FindBin(fMassLo); 
  int hbin = h->FindBin(fMassHi); 
  
  ylo = h->Integral(lbin, lbin+EDG)/NB; 
  yhi = h->Integral(hbin-EDG, hbin)/NB;
  dx = h->GetBinLowEdge(hbin) - h->GetBinLowEdge(lbin);
  
  // -- signal parameters
  double g0(0.), g1(5.3), g2(0.05); 
  xlo = g1-3.*g2;
  xhi = g1+3.*g2;
  lo = h->FindBin(xlo); 
  hi = h->FindBin(xhi);

  if (1 == mode) {
    // -- pol1 only
    p1 = (yhi-ylo)/dx;
    p0 = ylo - p1*h->GetBinLowEdge(1);
    f1->SetParameters(p0, p1); 
  } else if (10 == mode) {
    p1 = (yhi-ylo)/dx;
    p0 = ylo - p1*h->GetBinLowEdge(1);
    fF0->SetParameters(p0, p1); 
    g0 = (h->Integral(lo, hi) - fF0->Integral(xlo, xhi)/h->GetBinWidth(1))*h->GetBinWidth(1);  
    
    f1->SetParameters(g0, g1, g2, p0, p1); 
    f1->ReleaseParameter(0);     f1->SetParLimits(0, 0., 1.e7); 
    f1->ReleaseParameter(1);     f1->SetParLimits(1, 5.2, 5.45); 
    f1->ReleaseParameter(2);     f1->SetParLimits(2, 0.010, 0.080); 
  } else if (11 == mode) {
    p1 = (TMath::Log(yhi) - TMath::Log(ylo))/(-dx*NB); 
    p0 = ylo*TMath::Exp(p1*h->GetBinLowEdge(1)); 
    fF1->SetParameters(p0, p1); 
    g0 = (h->Integral(lo, hi) - fF1->Integral(xlo, xhi)/h->GetBinWidth(1))*h->GetBinWidth(1);  
    //     cout << "==> hist: " << h->Integral(lo, hi) 
    // 	 << " func: " << fF1->Integral(xlo, xhi)/h->GetBinWidth(1)
    // 	 << " g0: " << g0 
    // 	 << endl;
    f1->SetParameters(g0, g1, g2, p0, p1); 
    f1->ReleaseParameter(0);  f1->SetParLimits(0, 0., 1.e7); 
    f1->ReleaseParameter(1);  f1->SetParLimits(1, 5.2, 5.45); 
    f1->ReleaseParameter(2);  f1->SetParLimits(2, 0.010, 0.060); 
  } else if (12 == mode) {
    p1 = (TMath::Log(yhi) - TMath::Log(ylo))/(-dx*NB); 
    p0 = ylo*TMath::Exp(p1*h->GetBinLowEdge(1)); 
    fF1->SetParameters(p0, p1); 
    g0 = (h->Integral(lo, hi) - fF1->Integral(xlo, xhi)/h->GetBinWidth(1))*h->GetBinWidth(1);  

    f1->SetParameters(0.8*g0, g1, g2, 0.2, g1, 1.5*g2, p0, p1); 
    f1->ReleaseParameter(0);  f1->SetParLimits(0, 0., 1.e7); 
    f1->ReleaseParameter(1);  f1->SetParLimits(1, 5.2, 5.45); 
    f1->ReleaseParameter(2);  f1->SetParLimits(2, 0.010, 0.080); 
    f1->ReleaseParameter(3);  f1->SetParLimits(3, 0., 0.9); 
    f1->ReleaseParameter(4);  f1->SetParLimits(4, 5.2, 5.45); 
    f1->ReleaseParameter(5);  f1->SetParLimits(5, 0.020, 0.150); 
  } else if (13 == mode) {
    double offset = (yhi-ylo)/dx; 
    p1 = (TMath::Log(yhi) - TMath::Log(ylo))/(-dx*NB); 
    p0 = ylo*TMath::Exp(p1*h->GetBinLowEdge(1)); 
    fF1->SetParameters(p0, p1); 
    g0 = (h->Integral(lo, hi) - fF1->Integral(xlo, xhi)/h->GetBinWidth(1))*h->GetBinWidth(1);  

    f1->SetParameters(g0, g1, g2, p0, p1, offset, 1.); 
    f1->ReleaseParameter(0);  f1->SetParLimits(0, 0., 1.e7); 
    f1->ReleaseParameter(1);  f1->SetParLimits(1, 5.2, 5.45); 
    f1->ReleaseParameter(2);  f1->SetParLimits(2, 0.010, 0.080); 
  }
}


// ----------------------------------------------------------------------  
TH1D* AnalysisDistribution::sbsDistribution(const char *variable, const char *cut) {

  TCanvas *c0;
  int DISPLAY(1); 
  if (DISPLAY) {
    gStyle->SetOptTitle(1);
    c0 = (TCanvas*)gROOT->FindObject("c1"); 
    if (c0) {
      delete c0; 
    }
    c0 = new TCanvas("c1"); 
    c0->Clear(); 
    c0->Divide(2,3);
  }

  // -- this is not really necessary, could use the class members instead
  TH1D *hm = (TH1D*)gDirectory->Get(Form("%sMass%s", variable, cut));
  TH1D *h0 = (TH1D*)gDirectory->Get(Form("%s%s0", variable, cut));
  TH1D *h1 = (TH1D*)gDirectory->Get(Form("%s%s1", variable, cut));

  if (DISPLAY) {
    c0->cd(1);
    h0->Draw();
    c0->cd(2);
    h1->Draw();
    c0->cd(3); 
    hMassBGL->SetMinimum(0.);
    hMassBGL->SetMaximum(1.);
    hMassBGL->Draw();
    hMassBGH->Draw("same");
    hMassSG->Draw("same");
    
    c0->cd(4); 
  }

  double l0   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  double l1   = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
  double s0   = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
  double s1   = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
  double u0   = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
  double u1   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);

  // -- fit mass distribution 
  fMassPeak = 0.5*(s0+s1);
  fMassSigma= 0.2*(s1-s0);
  fMassLo   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  fMassHi   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
  fpIF->fLo = fMassLo;
  fpIF->fHi = fMassHi;

  cout << "fMass: " << fMassLo << " .. " << fMassHi << ", fMassPeak = " << fMassPeak << ", fMassSigma = " << fMassSigma << endl;

  double peak  = (fMassPeak>0.?fMassPeak:5.3);
  double sigma = (fMassSigma>0.?fMassSigma:0.04);
  TF1 *f1 = fpIF->pol1gauss2c(hm, peak, sigma);
  hm->SetMinimum(0.);

  TFitResultPtr r;
  r = hm->Fit(f1, "ls", "", fMassLo, fMassHi); 
  double hlimit(0.), llimit(0.); 
  f1->GetParLimits(3, llimit, hlimit);
  double relError = f1->GetParError(3)/f1->GetParameter(3); 
  double limit = TMath::Abs(f1->GetParameter(3) - llimit); 
  if (TMath::Abs(f1->GetParameter(3) - hlimit) < limit) limit = TMath::Abs(f1->GetParameter(3) - hlimit); 
  if (limit < relError || relError > 0.1) {
    cout << "%%%%% REFITTING %%%%% " << endl;
    f1->SetParameters(f1->GetParameters()); 
    f1->SetParameter(0, hm->GetMaximum()); 
    f1->FixParameter(3, 0.);
    f1->FixParameter(4, 0.);
    r = hm->Fit(f1, "ls", "", fMassLo, fMassHi); 
  }
  hm->DrawCopy();

  //   hm->Fit(f1, "l", "", fMassLo, fMassHi); 
  //   hm->DrawCopy(); 

  // -- compute integrals
  TF1 *fpol1 = (TF1*)gROOT->FindObject("fpol1"); 
  if (fpol1) delete fpol1; 
  fpol1 = new TF1("fpol1", f_p1, fMassLo, fMassHi, 2);
  fpol1->SetParameters(f1->GetParameter(5), f1->GetParameter(6)); 
  if (DISPLAY) {
    c0->cd(5); 
    fpol1->Draw();
  }

  double bgl = fpol1->Integral(l0, l1); 
  double sg  = fpol1->Integral(s0, s1); 
  double bgh = fpol1->Integral(u0, u1); 

  cout << "bgl (" << l0 << " .. " << l1 << ") = " << bgl << endl;
  cout << "sg  (" << s0 << " .. " << s1 << ") = " << sg << endl;
  cout << "bgh (" << u0 << " .. " << u1 << ") = " << bgh << endl;
  
  TH1D *h = (TH1D*)h0->Clone(Form("sbs_%s%s", variable, cut)); 
  h->Sumw2(); 
  h->Add(h0, h1, 1., -sg/(bgl+bgh)); 

  if (DISPLAY) {
    c0->cd(6); 
    h->Draw();
  }

  return h; 
}
 



// ----------------------------------------------------------------------
void AnalysisDistribution::fill(double value, double mass) {
  int mBin(-1); 

  // -- these histograms are just to KNOW afterwards what the mass windows were exactly
  hMassAll->Fill(mass); 
  if ((fSigLo < mass) && (mass < fSigHi)) {
    hMassSG->Fill(mass); 
    mBin = 0; 
  }
  if ((fBg1Lo < mass) && (mass < fBg1Hi)) {
    hMassBGL->Fill(mass); 
    mBin = 1; 
  }
  if ((fBg2Lo < mass) && (mass < fBg2Hi)) {
    hMassBGH->Fill(mass); 
    mBin = 1; 
  }

  if (fVerbose > 0) {
    cout << "value: " << value 
	 << " mass: " << mass 
	 << " fCutIdx: " << fCutIdx 
	 << " nm: " << fpAnaCuts->nMinus1CutsTrue(fCutIdx)
	 << " si: " << fpAnaCuts->singleCutTrue(fCutIdx)
	 << endl;
  }
  
  if (fpAnaCuts->singleCutTrue(fCutIdx)) {
    if (mBin > -1) hSi[mBin]->Fill(value); 
    hSi[2]->Fill(value); 
    hMassSi->Fill(mass);
  }
  if (fpAnaCuts->nMinus1CutsTrue(fCutIdx)) {
    if (mBin > -1) hNm[mBin]->Fill(value); 
    hNm[2]->Fill(value); 
    hMassNm->Fill(mass);
  }
  if (fpAnaCuts->cumulativeCutTrue(fCutIdx)) {
    if (mBin > -1) hCu[mBin]->Fill(value); 
    hCu[2]->Fill(value); 
    hMassCu->Fill(mass);
  }
  if (fpAnaCuts->allOtherCutsTrue(fCutIdx)) {
    if (mBin > -1) hAo[mBin]->Fill(value); 
    hAo[2]->Fill(value); 
    hMassAo->Fill(mass);
  }
  if (fpAnaCuts->singleCutTrue(fHLTIdx)) {
    if (mBin > -1) hHLT[mBin]->Fill(value); 
    hHLT[2]->Fill(value); 
    hMassHLT->Fill(mass);
  }

  if (true == *fpPreselCutTrue) {
    if (mBin > -1) hPresel[mBin]->Fill(value); 
    hPresel[2]->Fill(value); 
    hMassPresel->Fill(mass);
  } 

}

