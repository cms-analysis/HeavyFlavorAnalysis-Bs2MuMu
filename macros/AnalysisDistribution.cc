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
AnalysisDistribution::AnalysisDistribution(const char *name, const char *title, int nbins, double lo, double hi) {
  fSigLo = fSigHi = fBg1Lo = fBg1Hi = fBg2Lo = fBg2Hi = 0.0;

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

  hMassSi = new TH1D(Form("%sMassSi", name), Form("%sMassSi", name), 60, 4.8, 6.0); 
  hMassAo = new TH1D(Form("%sMassAo", name), Form("%sMassAo", name), 60, 4.8, 6.0); 
  hMassNm = new TH1D(Form("%sMassNm", name), Form("%sMassNm", name), 60, 4.8, 6.0); 
  hMassCu = new TH1D(Form("%sMassCu", name), Form("%sMassCu", name), 60, 4.8, 6.0); 
  hMassHLT= new TH1D(Form("%sMassHLT", name), Form("%sMassHLT", name), 60, 4.8, 6.0); 
  
}



// ----------------------------------------------------------------------
AnalysisDistribution::AnalysisDistribution(const char *name, double SigLo, double SigHi, double Bg1Lo, double Bg1Hi, double Bg2Lo, double Bg2Hi) {

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
  
}



// ----------------------------------------------------------------------
AnalysisDistribution::~AnalysisDistribution() {

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

  TF1 *f1; 
  if (0 == mode) {
    // -- just count the number of events in the histogram
    double n = h1->GetSumOfWeights(); 
    error = TMath::Sqrt(n); 
    return n; 
  } else if (1 == mode) {
    // -- One Gaussian plus pol1
    f1 = new TF1("f1", f_p1aG, 0., 6., 5);
    f1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 
    setFunctionParameters(f1, h1, 0); 
  } else if(2 == mode) {
    // -- One Gaussian plus expo
    f1 = new TF1("f1", f_eaG, 0., 6., 5);
    f1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 
    setFunctionParameters(f1, h1, 1); 
  } else {
    return -2.;
  }

  h1->Fit("f1", "", ""); 
  
  error = f1->GetParError(0)/h1->GetBinWidth(1); 
  return f1->GetParameter(0)/h1->GetBinWidth(1); 

}


// ----------------------------------------------------------------------
void AnalysisDistribution::setFunctionParameters(TF1 *f1, TH1 *h, int mode) {
  const int EDG(4), NB(EDG+1); 
  double lo(0.), hi(0.), dx(0.);
  double p0(0.), p1(0.);

  // -- (background) pol1 parameters from end points
  lo = h->Integral(1, 1+EDG); 
  hi = h->Integral(h->GetNbinsX()-EDG, h->GetNbinsX());
  dx = h->GetBinLowEdge(h->GetNbinsX()+1) - h->GetBinLowEdge(1);
  p1 = (hi-lo)/NB/dx;
  p0 = lo/NB - p1*h->GetBinLowEdge(1);

  // -- signal parameters
  double g0(0.), g1(5.3), g2(0.05); 
  lo = h->FindBin(g1)-3.*g2/h->GetBinWidth(1); 
  hi = h->FindBin(g1)+3.*g2/h->GetBinWidth(1);
  TF1 f0("ADf0", f_p1, 0., 6., 2);  f0.SetParameters(p0, p1); 
  g0 = (h->Integral(lo, hi) - f0.Integral(lo, hi))*h->GetBinWidth(1);  

  if (0 == mode) {
    f1->SetParameters(g0, g1, g2, p0, p1); 
    f1->ReleaseParameter(1); 
    f1->ReleaseParameter(2); 
    f1->SetParLimits(1, 5.2, 5.45); 
    f1->SetParLimits(2, 0.015, 0.050); 
  }

  // -- (background) expo parameters from end points
  lo = h->Integral(1, 1+EDG); 
  hi = h->Integral(h->GetNbinsX()-EDG, h->GetNbinsX());
  dx = h->GetBinLowEdge(h->GetNbinsX()+1) - h->GetBinLowEdge(1);
  p1 = (hi-lo)/NB/dx;
  p0 = lo/NB - p1*h->GetBinLowEdge(1);
  


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

