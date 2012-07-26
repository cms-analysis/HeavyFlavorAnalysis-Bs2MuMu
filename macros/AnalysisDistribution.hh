#ifndef ANALYSISDISTRIBUTION
#define ANALYSISDISTRIBUTION

#include "TString.h"
#include "TObject.h"
#include "TH1.h"
#include "TF1.h"

#include <iostream>
#include <fstream>
#include <string>

#include "AnalysisCuts.hh"
#include "initFunc.hh"

class AnalysisDistribution: public TObject {

public:

  AnalysisDistribution(const char *name, const char *title, int nbins, double lo, double hi); 
  AnalysisDistribution(const char *name, double SigLo=5.2, double SigHi=5.5, 
		       double Bg1Lo=5.0, double Bg1Hi=5.2, double Bg2Lo=5.5, double Bg2Hi=5.7); 
  ~AnalysisDistribution(); 

  void setAnalysisCuts(AnalysisCuts *p, const char *cutname);
  void setSigWindow(double lo, double hi);
  void setBg1Window(double lo, double hi);
  void setBg2Window(double lo, double hi);

  void fill(double value, double m);

  double fitMass(TH1 *h, double &error, int mode = 0); 
  // -- default version: pol1 + Gauss
  TH1D*  sbsDistribution(const char *variable, const char *cut); 
  // -- expo + Gauss
  TH1D* sbsDistributionExpoGauss(const char *variable, const char *cut);
  // -- expo + error function + Gauss
  TH1D* sbsDistributionExpoErrGauss(const char *variable, const char *cut, double preco=5.1);
  // -- pol1 + error function + Gauss
  TH1D* sbsDistributionPol1ErrGauss(const char *variable, const char *cut, double preco=5.1);
  
  void   setPreselCut(bool *p) {fpPreselCutTrue = p;} 

  std::string fCutName; 
  int fCutIdx, fHLTIdx; 

  double fMassLo, fMassHi, fMassPeak, fMassSigma;

  double fSigLo, fSigHi; 
  double fBg1Lo, fBg1Hi; 
  double fBg2Lo, fBg2Hi; 
 
  AnalysisCuts *fpAnaCuts; 

  bool *fpPreselCutTrue; 
  
  TH1D *hSi[3], *hAo[3], *hCu[3], *hNm[3], *hHLT[3], *hPresel[3]; 

  TH1D *hMassSi, *hMassAo, *hMassCu, *hMassNm, *hMassHLT, *hMassPresel; 

  TH1D *hMassAll, *hMassBGL, *hMassSG, *hMassBGH; 

//   TF1 *fF0, *fF1; 
//   TF1 *fP1, *fPG1, *fEG1, *fEG2, *fEPG; 

  int fVerbose;
  std::string fDirectory;
  std::string fControlPlotsFileName;
  
  initFunc *fpIF; 

  ClassDef(AnalysisDistribution,1) //Testing AnalysisDistribution
}; 

#endif
