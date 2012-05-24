#ifndef HISTCUTEFFICIENCY
#define HISTCUTEFFICIENCY

#include "TH1.h"

#include <iostream>
#include <fstream>




class HistCutEfficiency: public TObject {

public:

  HistCutEfficiency(TH1 *h); 
  HistCutEfficiency(TH1 *h1, double cut, int includeOverflow = 1);
  ~HistCutEfficiency(); 

  void eff(double cut);
  void eff(TH1* h1, double cut);
  void eff(TH1* h1, double cut, int includeOverFlow);

  void eff(double cut1, double cut2);
  void eff(TH1* h1, double cut1, double cut2);
  void eff(TH1* h1, double cut, double cut2, int includeOverFlow);

  TH1 *fH; 
  int fIncludeOverflow; 

  double nTot, lCut, uCut; // to hold total, number ABOVE the cut, number BELOW the cut
  double loEff, loErr; // efficiency of integral BELOW the cut
  double hiEff, hiErr; // efficiency of integral ABOVE the cut
  double inEff, inErr; // efficiency of integral inside the cuts
  double outEff, outErr; // efficiency of integral outside the cuts

  int fVerbose; 

  ClassDef(HistCutEfficiency,1) //Testing HistCutEfficiency

}; 

#endif
