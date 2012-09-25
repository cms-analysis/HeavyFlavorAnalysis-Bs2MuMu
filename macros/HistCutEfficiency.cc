#include "HistCutEfficiency.hh"

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"

#include <iostream>
#include <iomanip>

using namespace std; 

ClassImp(HistCutEfficiency)

// ----------------------------------------------------------------------
HistCutEfficiency::HistCutEfficiency(TH1 *h1) {
  fH = h1; 

  fVerbose = 0; 
  fIncludeOverflow = 1; 
}

// ----------------------------------------------------------------------
HistCutEfficiency::HistCutEfficiency(TH1 *h1, double cut, int includeOverflow) {
  fH = h1; 

  fVerbose = 0; 
  fIncludeOverflow = includeOverflow; 
  
  eff(cut); 
}


// ----------------------------------------------------------------------
HistCutEfficiency::~HistCutEfficiency() {
}


// ----------------------------------------------------------------------
void HistCutEfficiency::eff(TH1* h1, double cut, int includeOverflow) {
  fH = h1; 
  fIncludeOverflow = includeOverflow; 
  eff(cut); 
}

// ----------------------------------------------------------------------
void HistCutEfficiency::eff(TH1* h1, double cut) {
  fH = h1; 
  eff(cut); 
}


// ----------------------------------------------------------------------
void HistCutEfficiency::eff(double cut) {
//	cout << "eff: " << " cut: " << cut << endl;
  double epsilon = 1.e-4 * fH->GetBinWidth(1);
  int firstbin = 1; 
  if (fIncludeOverflow) firstbin = 0; 

  int lastbin  = fH->GetNbinsX(); 
  if (fIncludeOverflow) ++lastbin; 

  nTot   = fH->Integral(firstbin, lastbin);
  int cutbin = fH->FindBin(cut+epsilon); 
  lCut   = fH->Integral(cutbin, lastbin);
  uCut   = fH->Integral(firstbin, cutbin-1);
  if (nTot == 0) {
    hiEff = 0;
    hiErr = 0;
    loEff = 0;
    loErr = 0;
  } else {
    hiEff = lCut/nTot; 
    hiErr = dEff(static_cast<int>(lCut), static_cast<int>(nTot));
    
    loEff = uCut/nTot; 
    loErr = dEff(static_cast<int>(uCut), static_cast<int>(nTot));
  }
  if (fVerbose) {
    cout << "cut     : " << cut << endl;
    cout << "epsilon : " << epsilon << endl;
    cout << "firstbin: " << firstbin << endl;
    cout << "lastbin : " << lastbin << endl;
    cout << "cutbin  : " << cutbin << endl;
    cout << "ntot    : " << nTot << endl;
    cout << "ucut    : " << uCut << endl;
    cout << "loEff   : " << loEff << endl;
    cout << "lcut    : " << lCut << endl;
    cout << "hiEff   : " << hiEff << endl;
  }
}


// ----------------------------------------------------------------------
void HistCutEfficiency::eff(TH1* h1, double cut1, double cut2, int includeOverflow) {
  fH = h1; 
  fIncludeOverflow = includeOverflow; 
  eff(cut1, cut2); 
}

// ----------------------------------------------------------------------
void HistCutEfficiency::eff(TH1* h1, double cut1, double cut2) {
  cout << "eff0: " << " cut1: " << cut1 << " " << cut2 << endl;
  fH = h1; 
  eff(cut1, cut2); 
}


// ----------------------------------------------------------------------
void HistCutEfficiency::eff(double cut1, double cut2) {
  cout << "eff1: " << " cut1: " << cut1 << " " << cut2 << endl;
  double epsilon = 1.e-4 * fH->GetBinWidth(1);
  int firstbin = 1; 
  if (fIncludeOverflow) firstbin = 0; 

  int lastbin  = fH->GetNbinsX(); 
  if (fIncludeOverflow) ++lastbin; 
  
  double ntot   = fH->Integral(firstbin, lastbin);
  
  int cut1bin = fH->FindBin(cut1+epsilon); 
  int cut2bin = fH->FindBin(cut2-epsilon); 
  
  double cut   = fH->Integral(cut1bin, cut2bin);
  
  inEff = cut/ntot; 
  inErr = dEff(static_cast<int>(cut), static_cast<int>(ntot));

  outEff = (ntot-cut)/ntot; 
  outErr = dEff(static_cast<int>(ntot-cut), static_cast<int>(ntot));
  
  if (fVerbose) {
    cout << "cut1    : " << cut1 << endl;
    cout << "cut2    : " << cut2 << endl;
    cout << "epsilon : " << epsilon << endl;
    cout << "firstbin: " << firstbin << endl;
    cout << "lastbin : " << lastbin << endl;
    cout << "cut1bin : " << cut1bin << endl;
    cout << "cut2bin : " << cut2bin << endl;
    cout << "ntot    : " << ntot << endl;
    cout << "cut     : " << cut << endl;
    cout << "inEff   : " << inEff << endl;    
    cout << "outEff  : " << outEff << endl;
  }
}


