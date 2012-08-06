#ifndef PLOTCLASS
#define PLOTCLASS

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
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include "../macros/AnalysisCuts.hh"

// -- TMVA related
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
struct readerData {
  float pt, eta, m1eta, m2eta, m1pt, m2pt;
  float fls3d, alpha, maxdoca, pvip, pvips, iso, docatrk, chi2dof, closetrk; 
  float m;
};


#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../macros/initFunc.hh"

struct RedTreeData {
  int run, evt, ls, tm, pr, procid, pvn, rr;
  bool json, hlt, cb;
  double bdt, bdt2, pvw8;
  //  npv

  bool gmuid, gmupt, gmueta, gtqual, gtpt, gteta;
  double w8mu, w8tr;

  double pvlip, pvlips, pvlip2, pvlips2, pvip, pvips;

  int q, type;
  double pt, eta, phi, tau, m, cm, cosa, alpha, iso;
  int isotrk, closetrk; 
  double chi2, dof, pchi2dof, fls3d, fl3d, flxy, fl3dE, flsxy, docatrk, docatrkbdt, maxdoca, lip, lipE, tip, tipE;
  
  double osiso, osreliso, osmpt, osmptrel, osmdr; 

  int m1q, m2q;  
  double m1pt, m1eta, m1phi, m1ip, m1chi2;
  double m2pt, m2eta, m2phi, m2ip, m2chi2;
  double kpt, keta, kphi; 
  double k1pt, k1eta, k1phi, k2pt, k2eta, k2phi; 

  bool m1id, m1gt, m2id, m2gt, kgt; 
  int m1pix, m1bpix, m1bpixl1, m1pv; 
  int m2pix, m2bpix, m2bpixl1, m2pv; 

  double mudist, mudeltar; 
  double g1pt, g2pt, g3pt, g4pt, g1eta, g2eta, g3eta, g4eta, gtau;
  double t1pt, t1eta, t2pt, t2eta, t3pt, t3eta, t4pt, t4eta; 

  double mpsi, mkk;
  double psipt, psieta, psiphi, phipt, phieta, phiphi, dr;
  double md0, dm, ptd0; 

  double hm1pt, hm1eta, hm1phi, hm2pt, hm2eta, hm2phi; 
};


struct cuts {
  int index; 
  std::string xmlFile;
  double mBdLo, mBdHi, mBsLo, mBsHi;
  double etaMin, etaMax, pt; 
  double m1pt, m2pt, m1eta, m2eta;
  double iso, chi2dof, alpha, fls3d, docatrk; 
  // -- new variables
  double closetrk, pvlip, pvlips; 
  double bdt; 
  double maxdoca, pvlip2, pvlips2;
  double pvip, pvips; 
};

struct numbers {
  std::string name;
  int index;
  double effPtReco, effPtRecoE;
  double fitYield, fitYieldE;
  double fitYieldC, fitYieldCE;
  //  double genFileYield, genYield, genChanYield, recoYield, chanYield, muidYield, trigYield, candYield;
  double genAccFileYield, genYield, genFileYield, genAccYield;
  double effGenFilter, effGenFilterE; // only for the non-genAccFileYield!
  double recoYield, muidYield, trigYield, candYield;
  //  double genFileYieldE, genYieldE, genChanYieldE, recoYieldE, chanYieldE, muidYieldE, trigYieldE, candYieldE;
  double genFileYieldE, genYieldE, recoYieldE, muidYieldE, trigYieldE, candYieldE;
  double ana0Yield,  anaYield,  anaMuonYield,  anaTriggerYield,  anaWmcYield;
  double ana0YieldE, anaYieldE, anaMuonYieldE, anaTriggerYieldE, anaWmcYieldE;
  //  double cFrac, cFracE, acc, accE, accChan, accChanE;
  double acc, accE;
  double effChan, effChanE; 
  double effMuidMC, effMuidMCE, effTrigMC, effTrigMCE;
  double effMuidTNP, effMuidTNPE, effTrigTNP, effTrigTNPE;
  double effMuidTNPMC, effMuidTNPMCE, effTrigTNPMC, effTrigTNPMCE;
  double effCand, effCandE, effCandTE; 
  double effAna, effAnaE; 
  double effTot, effTotE, aEffProdMC, aEffProdMCE, effProdMC, effProdMCE, effProdTNP, effProdTNPE; 
  double effTotChan, effTotChanE; 
  double prodGenYield, combGenYield, chanGenYield; // eps*A corrected
  // -- signal stuff
  double expSignal;
  double bgObs, bgBsExp, bgBsExpE, bgBsExpTE, bgBdExp, bgBdExpE, bgBdExpTE;  // observed sideband bg and expected backgrounds in Bs/Bs signal windows
  double bsObs, bdObs; // observed events in Bs/B0 signal windows
  double bsExpObs, bdExpObs, bsExpObsE, bdExpObsE, bsExpObsTE, bdExpObsTE; // expected observation in Bs/B0 signal windows
  double tauBs, tauBsE, tauBd, tauBdE; 
  double offHiRare, offHiRareE, offLoRare, offLoRareE, bsRare, bsRareE, bdRare, bdRareE; 
  double offHi, offLo;
  double bsNoScaled, bsNoScaledE, bdNoScaled, bdNoScaledE;
  double pss, pssE, pssTE, pdd, pddE, pddTE;
  double psd, psdE, psdTE, pds, pdsE, pdsTE;
  double mBdLo, mBdHi, mBsLo, mBsHi;
  // -- total errors: quad sum of stat and syst errors
  double accTE, effAnaTE, effTotTE; 
  double tauTE; 
  double fitYieldTE;
  double effMuidTE, effMuidTNPTE, effMuidMCTE, effMuidTNPMCTE, effTrigTE, effTrigMCTE, effTrigTNPTE, effTrigTNPMCTE; 
};

// ----------------------------------------------------------------------
class plotClass: public TObject {

public:

  plotClass(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotClass();

  virtual void cd(const char *file) {fF[file]->cd();}
  virtual numbers* getNumbersNo(int i) {return fNumbersNo[i];}

  // -- initialization and setup
  // ---------------------------
  virtual void init(const char *files, const char *cuts, const char *dir, int mode);
  virtual void loadFiles(const char *files);
  virtual TFile* loadFile(std::string file);

  // -- main methods
  // --------------
  virtual void makeAll(int verbose);

  void loopTree(int mode, int proc = -1);
  void accEffFromEffTree(std::string fname, std::string dname, numbers &a, cuts &b, int proc);
  void filterEfficiency(std::string fname, std::string name);
  void normYield(TH1 *h, int mode, double lo = 5.15, double hi=5.5, double preco=-1.);
  void csYield(TH1 *h, int mode, double lo = 5.25, double hi=5.6, double preco=-1.);
  void dpYield(TH1 *h, int mode, double lo=5.1, double hi=5.6, int bdtc=-10);
  void bgBlind(TH1 *h, int mode = 2, double lo = 4.5, double hi = 6.5); 

  void singleEventPrintout(std::string suffix, std::string st, int ievt);
  void printNumbers(numbers &a, ostream &OUT);
  void initNumbers(numbers *a); 
  int  detChan(double m1eta, double m2eta);


  virtual void dumpSamples();
  virtual void dumpCutNames(const char *h);
  virtual void readCuts(const char *filename); 
  virtual void printCuts(ostream &OUT); 
  virtual double recalcMass(double m1, double m2);

  virtual void setErrors(TH1D *h);
  virtual void stamp(double x1, std::string text1, double x2, std::string text2);
  virtual std::string scientificTex(double n, double nE, std::string name, double base = 1.e-2, int digits = 2);
  virtual std::string formatTex(double n, std::string name, int digits, int sgn = 0);
  virtual std::string formatTex(double n, std::string name, std::string tag);
  virtual void drawArrow(double height, int mode = 0, double y = 0.1);
  virtual void drawBox(int mode, double hi = 0.5, int ylo = 0.01);
  virtual void replaceAll(std::string &s, std::string a, std::string b);
  virtual void makeCanvas(int i = 3);
  virtual void newLegend(double x1, double y1, double x2, double y2, std::string title = "");


  std::string fFiles; 
  // -- Files for Signal and Normalization modes in data and MC
  std::map<std::string, TFile*> fF; 
  std::map<std::string, double> fLumi, fFilterEff, fBF, fBFE, fProdR; 
  std::map<std::string, std::string> fName; 

  std::map<std::string, double> fNgen;

  int fShow; 
  std::string fFile; 
  std::string fStampString, fStampCms;

  // -- output histograms and numbers
  TFile *fHistFile; 
  ofstream fOUT, fTEX; 

  bool fDoPrint; 
  int fVerbose; 

  // -- cuts 
  std::string fCutsFileName;
  std::vector<cuts*> fCuts; 

  // -- functions
  TF1 *f0, *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9; 

  int fMode; 
  double fPreco; 

  // -- Display utilities
  int fFont; 
  double fSize; 
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5;
  TLatex *tl; 
  TBox *box;
  TArrow *pa; 
  TLine *pl; 
  TLegend *legg;
  TLegendEntry *legge;

  char line[200];

  std::string fDirectory;
  std::string fNumbersFileName;
  std::string fSample; 
  std::string fSuffix; 

  std::string fRareName; 

  unsigned int fNchan;
  int fChan, fComb, fOver; 

  initFunc  *fpFunc;

  double fMassLo, fMassHi
    , fBgLo, fBgHi
    , fSgLo, fSgHi
    , fNoLo, fNoHi
    , fCsLo, fCsHi
    ;


  std::vector<TH1D*> fhMassAcc;
  std::vector<TH1D*> fhMassChan;
  std::vector<TH1D*> fhMassAbsNoCuts;
  std::vector<TH1D*> fhMassNoCuts;
  std::vector<TH1D*> fhMassNoCutsManyBins;
  std::vector<TH1D*> fhMassWithAnaCuts;
  std::vector<TH1D*> fhMassWithMuonCuts; 
  std::vector<TH1D*> fhMassWithTriggerCuts; 
  std::vector<TH1D*> fhMassWithAllCuts;
  std::vector<TH1D*> fhMassWithAllCutsBlind;
  std::vector<TH1D*> fhNorm, fhNormC, fhDstarPi;
  std::vector<TH1D*> fhMassWithMassCuts;
  std::vector<TH1D*> fhMassWithAnaCutsManyBins;
  std::vector<TH1D*> fhMassWithMuonCutsManyBins; 
  std::vector<TH1D*> fhMassWithTriggerCutsManyBins; 
  std::vector<TH1D*> fhMassWithAllCutsManyBins; 
  std::vector<TH1D*> fhMassWithMassCutsManyBins; 
  std::vector<TH1D*> fhMuId, fhMuIdMC;
  std::vector<TH1D*> fhMuTr, fhMuTrMC;
  std::vector<TH1D*> fh0TNPTrigger, fh1TNPTrigger, fh0TNPMuID, fh1TNPMuID;
  std::vector<TH1D*> fh0TNPMCTrigger, fh1TNPMCTrigger, fh0TNPMCMuID, fh1TNPMCMuID;
  std::vector<TH1D*> fh0MCTrigger, fh1MCTrigger, fh0MCMuID, fh1MCMuID;

  TH1D *hMassPiPi, *hMassKPi, *hMassPiK; 


  std::vector<numbers*> fNumbersBs, fNumbersBd, fNumbersNo, fNumbersCs, fNumbersBla; 

  bool fDoPrintSingleEvent;
  int  fPrintSingleEvt, fPrintSingleRun;
  bool fDoUseBDT;
  bool fDoApplyCowboyVeto, fDoApplyCowboyVetoAlsoInSignal; 
  bool fInvertedIso;
  bool fNormProcessed; 

  double fBgExp, fBgExpE; 
  double fBgHist, fBgHistE, fBgHistLo, fBgHistHi; 
  double fBgHistExp, fBgHistExpE;

  double fNoSig, fNoSigE; 
  double fCsSig, fCsSigE; 
  double fDpSig, fDpSigE; 

  double fu, fs, fsfu, fsfuE;


  // -- stuff to run over the tree from any derived class (using loopFunction() to hook in specifics)
  virtual TTree* getTree(std::string mode);
  virtual void setupTree(TTree *t, std::string mode); 
  virtual void loopOverTree(TTree *t, std::string mode, int function, int nevts = -1); 
  virtual void candAnalysis(int mode);
  virtual void loopFunction(int function) {std::cout << "replace me" << std::endl;} 
  virtual void calcBDT();
  struct RedTreeData fb; 
  std::vector<TMVA::Reader*> fReaderEven; 
  std::vector<TMVA::Reader*> fReaderOdd; 
  bool fIsMC, fIsSignal;
  double fBDT; 
  readerData frd; 
  double MASSMIN, MASSMAX, 
    SIGBOXMIN, SIGBOXMAX, 
    BGLBOXMIN, BGLBOXMAX, 
    BGHBOXMIN, BGHBOXMAX
    ; 

  bool fPreselection,
    fWideMass, 
    fGoodHLT, fGoodMuonsID, fGoodMuonsPt, fGoodMuonsEta, fGoodTracks, fGoodTracksPt, fGoodTracksEta, fGoodQ,
    fGoodPvAveW8, fGoodIp, fGoodIpS, fGoodMaxDoca, fGoodPt, fGoodEta, fGoodAlpha, fGoodFLS, fGoodChi2, fGoodIso, fGoodCloseTrack, 
    fGoodDocaTrk, fGoodLastCut; 

  AnalysisCuts fAnaCuts; 

  ClassDef(plotClass,1) //Testing plotClass
};

TMVA::Reader* setupReader(std::string xmlFile, readerData& rd);


#endif

