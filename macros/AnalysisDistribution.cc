#include "AnalysisDistribution.hh"

#include <iostream>
#include <iomanip>

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"

#include "TMath.h"
#include "TDirectory.h"

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
AnalysisDistribution::AnalysisDistribution(const char *name, const char *title, int nbins, double lo, double hi) {
  fSigLo = fSigHi = fBg1Lo = fBg1Hi = fBg2Lo = fBg2Hi = 0.0;
  fMassLo = 4.8;
  fMassHi = 6.0;

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

  int NBINS(30); 
  hMassSi = new TH1D(Form("%sMassSi", name), Form("%sMassSi", name), NBINS, fMassLo, fMassHi); 
  hMassAo = new TH1D(Form("%sMassAo", name), Form("%sMassAo", name), NBINS, fMassLo, fMassHi); 
  hMassNm = new TH1D(Form("%sMassNm", name), Form("%sMassNm", name), NBINS, fMassLo, fMassHi); 
  hMassCu = new TH1D(Form("%sMassCu", name), Form("%sMassCu", name), NBINS, fMassLo, fMassHi); 
  hMassHLT= new TH1D(Form("%sMassHLT", name), Form("%sMassHLT", name), NBINS, fMassLo, fMassHi); 

  fF0 = new TF1("ADf0", f_p1, 0., 6., 2);  
  fF1 = new TF1("ADf1", f_expo, 0., 6., 2);  

  fP1  = new TF1("ADp1", f_p1, 0., 6., 2);        fP1->SetParNames("Offset", "Slope"); 			   
  fPG1 = new TF1("ADpg1", f_p1aG, 0., 6., 5);     fPG1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 
  fEG1 = new TF1("ADeg1", f_eaG, 0., 6., 5);      fEG1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 
  fEG2 = new TF1("ADeg2", f_ea2G, 0., 6., 8);     fEG2->SetParNames("Area", "Peak", "Sigma", "Fraction2", "Peak2", "Sigma2", "Offset", "Slope"); 

  fEPG = new TF1("ADepg", f_epaG, 0., 6., 7);     fEPG->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope", "Offset", "Slope"); 
  
}



// ----------------------------------------------------------------------
AnalysisDistribution::AnalysisDistribution(const char *name, double SigLo, double SigHi, double Bg1Lo, double Bg1Hi, double Bg2Lo, double Bg2Hi) {

  fMassLo = 4.8;
  fMassHi = 6.0;

  hMassSi = (TH1D*)gDirectory->Get(Form("%sMassSi", name)); 
  hMassAo = (TH1D*)gDirectory->Get(Form("%sMassAo", name)); 
  hMassNm = (TH1D*)gDirectory->Get(Form("%sMassNm", name)); 
  hMassCu = (TH1D*)gDirectory->Get(Form("%sMassCu", name)); 
  hMassHLT= (TH1D*)gDirectory->Get(Form("%sMassHLT", name)); 

  for (int i = 0; i < 3; ++i) {
    hSi[i] = (TH1D*)gDirectory->Get(Form("%sSi%d", name, i)); 
    hAo[i] = (TH1D*)gDirectory->Get(Form("%sAo%d", name, i)); 
    hNm[i] = (TH1D*)gDirectory->Get(Form("%sNm%d", name, i)); 
    hCu[i] = (TH1D*)gDirectory->Get(Form("%sCu%d", name, i)); 
    hHLT[i]= (TH1D*)gDirectory->Get(Form("%sHLT%d", name, i)); 
  }

  fSigLo = SigLo; 
  fSigHi = SigHi; 
  fBg1Lo = Bg1Lo; 
  fBg1Hi = Bg1Hi; 
  fBg2Lo = Bg2Lo;
  fBg2Hi = Bg2Hi;

  fF0 = new TF1("ADf0", f_p1, 0., 6., 2);  
  fF1 = new TF1("ADf1", f_expo, 0., 6., 2);  

  fP1  = new TF1("ADp1", f_p1, 0., 6., 2);        fP1->SetParNames("Offset", "Slope"); 			   
  fPG1 = new TF1("ADpg1", f_p1aG, 0., 6., 5);     fPG1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 
  fEG1 = new TF1("ADeg1", f_eaG, 0., 6., 5);      fEG1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 
  fEG2 = new TF1("ADeg2", f_ea2G, 0., 6., 8);     fEG2->SetParNames("Area", "Peak", "Sigma", "Fraction2", "Peak2", "Sigma2", "Offset", "Slope"); 
  fEPG = new TF1("ADepg", f_epaG, 0., 6., 7);     fEPG->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope", "Offset", "Slope"); 
  
}



// ----------------------------------------------------------------------
AnalysisDistribution::~AnalysisDistribution() {
  delete fF0; 
  delete fF1; 
  delete fP1; 
  delete fPG1; 
  delete fEG1; 
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
    // -- One Gaussian plus expo
    f1 = fEG1;
    setFunctionParameters(f1, h1, mode); 
    h1->Fit(f1, "", "", fMassLo, fMassHi); 
    error = f1->GetParError(0)/h1->GetBinWidth(1); 
    return f1->GetParameter(0)/h1->GetBinWidth(1); 
  } else if (12 == mode) {
    // -- Two Gaussian plus expo
    f1 = fEG2;
    setFunctionParameters(f1, h1, mode); 
    h1->Fit(f1, "", "", fMassLo, fMassHi); 
    // FIXME: Wrong yield/error for DOUBLE gaussian!
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
  double p0(0.), p1(0.);

  int lbin = h->FindBin(fMassLo); 
  int hbin = h->FindBin(fMassHi); 
  
  ylo = h->Integral(lbin, lbin+EDG)/NB; 
  yhi = h->Integral(hbin-EDG, hbin)/NB;
  dx = h->GetBinLowEdge(hbin) - h->GetBinLowEdge(lbin);
  
  // -- signal parameters
  double g0(0.), g1(5.3), g2(0.05); 
  lo = h->FindBin(g1)-3.*g2/h->GetBinWidth(1); 
  hi = h->FindBin(g1)+3.*g2/h->GetBinWidth(1);

  if (1 == mode) {
    // -- pol1 only
    p1 = (yhi-ylo)/dx;
    p0 = ylo - p1*h->GetBinLowEdge(1);
    f1->SetParameters(p0, p1); 
  } else if (10 == mode) {
    p1 = (yhi-ylo)/dx;
    p0 = ylo - p1*h->GetBinLowEdge(1);
    fF0->SetParameters(p0, p1); 
    g0 = (h->Integral(lo, hi) - fF0->Integral(lo, hi))*h->GetBinWidth(1);  
    
    f1->SetParameters(g0, g1, g2, p0, p1); 
    f1->ReleaseParameter(0);     f1->SetParLimits(0, 0., 1.e7); 
    f1->ReleaseParameter(1);     f1->SetParLimits(1, 5.2, 5.45); 
    f1->ReleaseParameter(2);     f1->SetParLimits(2, 0.010, 0.080); 
  } else if (11 == mode) {
    p1 = (TMath::Log(yhi) - TMath::Log(ylo))/(-dx*NB); 
    p0 = ylo*TMath::Exp(p1*h->GetBinLowEdge(1)); 
    fF1->SetParameters(p0, p1); 
    g0 = (h->Integral(lo, hi) - fF1->Integral(lo, hi))*h->GetBinWidth(1);  

    f1->SetParameters(g0, g1, g2, p0, p1); 
    f1->ReleaseParameter(0);  f1->SetParLimits(0, 0., 1.e7); 
    f1->ReleaseParameter(1);  f1->SetParLimits(1, 5.2, 5.45); 
    f1->ReleaseParameter(2);  f1->SetParLimits(2, 0.010, 0.080); 
  } else if (12 == mode) {
    p1 = (TMath::Log(yhi) - TMath::Log(ylo))/(-dx*NB); 
    p0 = ylo*TMath::Exp(p1*h->GetBinLowEdge(1)); 
    fF1->SetParameters(p0, p1); 
    g0 = (h->Integral(lo, hi) - fF1->Integral(lo, hi))*h->GetBinWidth(1);  

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
    g0 = (h->Integral(lo, hi) - fF1->Integral(lo, hi))*h->GetBinWidth(1);  

    f1->SetParameters(g0, g1, g2, p0, p1, offset, 1.); 
    f1->ReleaseParameter(0);  f1->SetParLimits(0, 0., 1.e7); 
    f1->ReleaseParameter(1);  f1->SetParLimits(1, 5.2, 5.45); 
    f1->ReleaseParameter(2);  f1->SetParLimits(2, 0.010, 0.080); 
  }
}

  
 



// ----------------------------------------------------------------------
void AnalysisDistribution::fill(double value, double mass) {
  int mBin(-1); 
  if ((fSigLo < mass) && (mass < fSigHi)) mBin = 0; 
  if ((fBg1Lo < mass) && (mass < fBg1Hi)) mBin = 1; 
  if ((fBg2Lo < mass) && (mass < fBg2Hi)) mBin = 1; 
  
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

}

