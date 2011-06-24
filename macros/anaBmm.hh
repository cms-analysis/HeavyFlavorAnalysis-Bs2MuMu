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
#include <map>

#include "initFunc.hh"

struct numbers {
  std::string name;
  int index;
  double effGenFilter, effGenFilterE; 
  double fitYield, fitYieldE;
  double genFileYield, genYield, genChanYield, recoYield, chanYield, muidYield, trigYield, candYield;
  double genFileYieldE, genYieldE, genChanYieldE, recoYieldE, chanYieldE, muidYieldE, trigYieldE, candYieldE;
  double ana0Yield,  anaYield,  anaMuonYield,  anaTriggerYield,  anaWmcYield;
  double ana0YieldE, anaYieldE, anaMuonYieldE, anaTriggerYieldE, anaWmcYieldE;
  double cFrac, cFracE, acc, accE, accChan, accChanE;
  double effChan, effChanE; 
  double accMuidMC, accMuidMCE, accTrigMC, accTrigMCE;
  double effMuidMC, effMuidMCE, effTrigMC, effTrigMCE;
  double effMuidPid, effMuidPidE, effTrigPid, effTrigPidE;
  double effCand, effCandE; 
  double effAna, effAnaE; 
  double effTot, effTotE, aEffProdMC, aEffProdMCE, effProdMC, effProdMCE, effProdPid, effProdPidE; 
  double effTotChan, effTotChanE; 
  double prodGenYield, combGenYield, chanGenYield; // eps*A corrected
  // -- signal stuff
  double expSignal;
  double bgObs, bgBsExp, bgBsExpE, bgBdExp, bgBdExpE; 
  double bsObs, bdObs; 
  double bsRare, bsRareE, bdRare, bdRareE; 
  double pss, pssE, pdd, pddE;
  double psd, psdE, pds, pdsE;
  double mBdLo, mBdHi, mBsLo, mBsHi;
};


struct cuts {
  int index; 
  double mBdLo, mBdHi, mBsLo, mBsHi;
  double etaMin, etaMax, pt; 
  double m1pt, m2pt, m1eta, m2eta;
  double iso1, chi2dof, alpha, fls3d, docatrk; 
};



class anaBmm: public TObject {

public:

  anaBmm(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~anaBmm();

  // -- initialization and setup
  // ---------------------------
  void init(const char *files, const char *cuts, const char *dir, int mode);
  void loadFiles(const char *files);
  TFile* loadFile(std::string file, std::string type);

  // -- main methods
  // --------------
  void makeAll(int channels = 3);

  void accEffFromEffTree(std::string fname, numbers &a, cuts &b, int proc = -1);
  void filterEff(int mode);

  void computeNormUL();
  void computeCsBF();

  void allEffTables();
  void effTable(std::string mode, const char *region = "A");
  double computeDelta(const char *s1, const char *s2, int relative = 1);

  void sbsDistributionOverlay(std::string file1, std::string file2, const char *selection="Ao", const char *region = ""); 
  void sbsDistributionOverlaySameFile(std::string file1, std::string r1 = "B_iso5", 
				      std::string r2 = "E_iso5", std::string sel = "Ao", std::string L1="barrel", std::string L2="endcap"); 

  void breco(TH1D *h); 
  void effTree(int mode);
  TH1* loopTree(int mode, int proc = -1);

  void plotVar(const char *plotstring, const char *cuts, const char *options = "");
  void testSimpleUL(const char *cuts);
  void testUL(int ichan = 0);
  void optimizeULs(int nruns, int seed = -1);
  void bestUL(const char *fname, int chan = 0, int nselections = 10); 
  void testEff(const char *additionalCuts, 
	       const char *basicCuts="gmuid&&hlt&&gtqual&&m1pt>3&&m2pt>2.&&acos(cosa)<0.04&&iso1>0.7&&fls3d>13&&chi2<1.5&&pt>4");
  void optimizeCut(const char *cut, double lo, double hi, const char *otherCuts); 
  void playUL(int nruns = 10); 

  void triggerSignal(const char *cuts = "iso1>0.6&&fls3d>4"); 
  void triggerNorm(const char *cuts = "iso1>0.6&&fls3d>4&&chi2/dof<2"); 

  void normYield(TH1 *h, int mode, double lo = 5.15, double hi=5.5);
  void csYield(TH1 *h, int mode, double lo = 5.25, double hi=5.6);
  void bgBlind(TH1 *h, int mode = 2, double lo = 4.5, double hi = 6.5); 

  void rareBg(); 
  void isoMean(const char *var = "iso1", const char *cuts="gmuid&&hlt&&acos(cosa)<0.05&&chi2/dof<2&&fls3d>10", int maxPvN = 15);
  void isoProcess(const char *var, const char *cuts, int loop = 1);
  void plotWithCut(const char *var, const char *cuts, double cut, const char *title, double hmin, double hmax); 
  void puEff(const char *var, double cut=0.75, const char *ylabel="#epsilon(I>0.75)", const char* file="NoData", const char *selection="Ao"); 

  // -- Utilities and helper methods
  // -------------------------------
  void readCuts(const char *filename); 
  void printCuts(ostream &OUT); 
  void initNumbers(numbers *a); 
  void printNumbers(numbers &a, ostream &OUT);
  void printUlcalcNumbers();
  void printCsBFNumbers();
  void texUlcalcNumbers(const char *filename, const char *output); 
  int  detChan(double m1eta, double m2eta);

  void dumpSamples();
  void dumpCutNames();
  void setErrors(TH1D *h);
  std::string formatTex(double n, std::string name, int digits);
  std::string formatTex(double n, std::string name, std::string tag);
  void replaceAll(std::string &s, std::string a, std::string b);
  int  wait();
  void makeCanvas(int i = 3);
  void newLegend(double x1, double y1, double x2, double y2);

  double barlow(int nobs, double bg = 0., double bgE = 0., double sE = 0.);
  void rolkeM3();
  void rolkeM3(int x, double bm, double em, double sde, double sdb);

  // -- Files for Signal and Normalization modes in data and MC
#define MAXFILES 30
  TFile *fpData[MAXFILES], *fpMc[MAXFILES];
  double fDataLumi[MAXFILES], fMcLumi[MAXFILES];
  int fNData, fNMc;
  std::string fDataString[MAXFILES], fMcString[MAXFILES];
  int fSgData, fSgMc, fBdMc, fNoData, fNoMc;  // the indices for the default files
  int fCsData, fCsMc; // control sample Bs -> J/psi phi

  int fShow; 
  TString fFile; 

  // -- output histograms
  TFile *fHistFile; 
  ofstream fOUT, fTEX; 

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
  std::string fUlcalcFileName;
  std::string fSample; 
  std::string fSuffix; 

  std::map<std::string, TFile*> fF; 
  std::map<std::string, double> fLumi, fFilterEff; 
  std::map<std::string, std::string> fName; 

  std::vector<cuts*> fCuts; 
  std::vector<TH1D*> fhMassAcc;
  std::vector<TH1D*> fhMassChan;
  std::vector<TH1D*> fhMassAbsNoCuts;
  std::vector<TH1D*> fhMassNoCuts;
  std::vector<TH1D*> fhMassWithAnaCuts;
  std::vector<TH1D*> fhMassWithMuonCuts; 
  std::vector<TH1D*> fhMassWithTriggerCuts; 
  std::vector<TH1D*> fhMassWithAllCuts;
  std::vector<TH1D*> fhMassWithMassCuts;
  std::vector<TH1D*> fhMassWithAnaCutsManyBins;
  std::vector<TH1D*> fhMassWithMuonCutsManyBins; 
  std::vector<TH1D*> fhMassWithTriggerCutsManyBins; 
  std::vector<TH1D*> fhMassWithAllCutsManyBins; 
  std::vector<TH1D*> fhMassWithMassCutsManyBins; 
  std::vector<TH1D*> fhMuId;
  std::vector<TH1D*> fhMuTr;
  std::vector<TH1D*> fh0PidTrigger, fh1PidTrigger, fh0PidMuID, fh1PidMuID;
  std::vector<TH1D*> fh0MCTrigger, fh1MCTrigger, fh0MCMuID, fh1MCMuID;


  unsigned int fNchan;
  int fChan, fComb, fOver; 

  initFunc  *fpFunc;
  TRolke    *fpRolke;

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
  double fBgHistExp, fBgHistExpE; 
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


  //  numbers fNumbersSig, fNumbersNorm, fNumbersCS; 
  std::vector<numbers*> fNumbersBs, fNumbersBd, fNumbersNorm, fNumbersCS; 

  double fBF, fu, fs;
  double fMassLo, fMassHi;

  bool fDoPrint;

  ClassDef(anaBmm,1) //Testing anaBmm
};


#endif

