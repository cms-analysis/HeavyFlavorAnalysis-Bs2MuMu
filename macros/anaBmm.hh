#ifndef ANABMM
#define ANABMM

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TRolke.h"
#include "TVirtualPad.h"  // access to gPad

#include <iostream>
#include <fstream>
#include <vector>

#include "initFunc.hh"

struct numbers {
  std::string name;
  double effGenFilter, effGenFilterE; 
  double fitYield, fitYieldE;
  double genFileYield, genYield, recoYield, muidYield, trigYield, candYield, ana0Yield, anaNmcYield, anaYield;
  double acc, accE;
  double effMuidMC, effMuidMCE, effTrigMC, effTrigMCE;
  double effMuidPid, effMuidPidE, effTrigPid, effTrigPidE;
  double effCand, effCandE; 
  double effAna, effAnaE; 
  double effTot, effTotE; 
};



class anaBmm: public TObject {

public:

  anaBmm(const char *files, const char *dir = ".", int mode = 0);
  ~anaBmm();

  // -- initialization and setup
  // ---------------------------
  void init(const char *files, const char *dir, int mode);
  void loadFiles(const char *files);
  TFile* loadFile(std::string file, std::string type);

  // -- main methods
  // --------------
  void makeAll(int channels = 3);
  void dumpCutNames();
  void allEffTables();
  void computeNormUL();
  void computeCsBF();
  void effTable(std::string mode);
  void breco(TH1D *h); 
  void acceptanceAndPreselection(numbers &a);
  void effTree(int mode);
  TH1* loopTree(int mode);

  void plotVar(const char *plotstring, const char *cuts, const char *options = "");
  void testUL(const char *cuts);
  void optimizeCut(const char *cut, double lo, double hi, const char *otherCuts); 
  void optimizeUL(int nruns = 10); 
  void normYield(TH1 *h, int mode, double lo = 5.15, double hi=5.5);
  void csYield(TH1 *h, int mode, double lo = 5.25, double hi=5.6);
  void bgBlind(TH1 *h, int mode = 2, double lo = 4.5, double hi = 6.5); 
  double barlow(int nobs, double bg = 0., double bgE = 0., double sE = 0.);
  void rolkeM3();
  void rolkeM3(int x, double bm, double em, double sde, double sdb);
  void printNumbers(numbers &a);


  // -- Utilities and helper methods
  // -------------------------------
  void setErrors(TH1D *h);
  std::string formatTex(double n, std::string name, std::string tag);
  void replaceAll(std::string &s, std::string a, std::string b);
  int  wait();
  void makeCanvas(int i = 3);

  // -- Files for Signal and Normalization modes in data and MC
#define MAXFILES 20
  TFile *fpData[MAXFILES], *fpMc[MAXFILES];
  double fDataLumi[MAXFILES], fMcLumi[MAXFILES];
  int fNData, fNMc;
  std::string fDataString[MAXFILES], fMcString[MAXFILES];
  int fSgData, fSgMc, fNoData, fNoMc;  // the indices for the default files
  int fCsData, fCsMc; // control sample Bs -> J/psi phi

  int fShow; 
  TString fFile; 

  // -- output histograms
  TFile *fHistFile; 

  // -- functions
  TF1 *f0, *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9; 

  int fMode; 

  // -- Display utilities
  int fFont; 
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5;
  TLatex *tl; 
  TBox *box;
  TArrow *pa; 
  TLine *pl; 
  TLegend *legg;
  TLegendEntry *legge;

  char line[200];
  int fFixFit; 

  std::string fDirectory;
  std::string fNumbersFileName;
  std::string fSample; 
  std::string fSuffix; 


  initFunc  *fpFunc;
  TRolke    *fpRolke;

  // -- cuts
  double M2PT, ISO1, CHI2, FLS3D, ALPHA;
  
  // -- analysis numbers
  double fSigLo, fSigHi;
  double fBgLo, fBgHi;
  double fPeak, fWidth;
  double fNormLo, fNormHi;
  double fCsLo, fCsHi; 

  double fNul, fUL;
  int    fNobs, fNobsExp;
  double fBgExp, fBgExpE; 
  double fBgHist, fBgHistE; 
  double fAcc, fAccE, fAccNum;
  double fEff, fEffE;
  double fEffAna, fEffAnaE;
  double fEffPresel, fEffPreselE;
  double fEffMuID, fEffMuIDE;
  double fEffTrig, fEffTrigE;
  double fEffTot, fEffTotE;

  double fSigGenFilter, fSigGenFilterE; 
  double fNormGenFilter, fNormGenFilterE; 
  double fCsGenFilter, fCsGenFilterE; 

  // -- normalization numbers
  double fNormSig, fNormSigE;
  double fNormAcc, fNormAccE, fNormAccNum;
  double fNormEff, fNormEffE;
  double fNormEffAna, fNormEffAnaE;
  double fNormEffPresel, fNormEffPreselE;
  double fNormEffMuID, fNormEffMuIDE;
  double fNormEffTrig, fNormEffTrigE;
  double fNormEffTot, fNormEffTotE;

  // -- control sample numbers
  double fCsSig, fCsSigE;
  double fCsAcc, fCsAccE, fCsAccNum;
  double fCsEff, fCsEffE;
  double fCsEffAna, fCsEffAnaE;
  double fCsEffPresel, fCsEffPreselE;
  double fCsEffMuID, fCsEffMuIDE;
  double fCsEffTrig, fCsEffTrigE;
  double fCsEffTot, fCsEffTotE;


  numbers fNumbersSig, fNumbersNorm, fNumbersCS; 

  double fBF, fu, fs;
  double fMassLo, fMassHi;

  ClassDef(anaBmm,1) //Testing anaBmm
};


#endif

