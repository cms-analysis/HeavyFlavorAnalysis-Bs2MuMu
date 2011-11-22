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
  int firstbin = 1; 
  if (fIncludeOverflow) firstbin = 0; 

  int lastbin  = fH->GetNbinsX(); 
  if (fIncludeOverflow) ++lastbin; 

  double ntot   = fH->Integral(firstbin, lastbin);

  int cutbin = fH->FindBin(cut); 

  double lcut   = fH->Integral(cutbin, lastbin);
  
  hiEff = lcut/ntot; 
  hiErr = dEff(static_cast<int>(lcut), static_cast<int>(ntot));
  
  double ucut   = fH->Integral(firstbin, cutbin-1);

  if (fVerbose) {
    cout << "firstbin: " << firstbin << endl;
    cout << "lastbin : " << lastbin << endl;
    cout << "cutbin  : " << cutbin << endl;
    cout << "ntot    : " << ntot << endl;
    cout << "ucut    : " << ucut << endl;
    cout << "lcut    : " << lcut << endl;
  }

  loEff = ucut/ntot; 
  loErr = dEff(static_cast<int>(ucut), static_cast<int>(ntot));
}


