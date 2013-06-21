#ifndef PLOTCLASS
#define PLOTCLASS

#include <stdio.h>
#include <stdlib.h>
#include <cstdarg>

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

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../macros/AnalysisCuts.hh"

#include "../macros2/preselection.hh"

// -- TMVA related
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

struct readerData {
  float pt, eta, m1eta, m2eta, m1pt, m2pt;
  float fls3d, alpha, maxdoca, pvip, pvips, iso, docatrk, chi2dof, closetrk; 
  float m;
  float closetrks1, closetrks2, closetrks3;
  float m1iso, m2iso; 
  float pvdchi2, othervtx;
  float pvlip2, pvlips2;
};


#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../macros/initFunc.hh"
#include "RedTreeData.hh"

struct cuts {
  int index; 
  std::string xmlFile;
  double mBdLo, mBdHi, mBsLo, mBsHi;
  double etaMin, etaMax, pt; 
  double m1pt, m2pt, m1eta, m2eta;
  double iso, chi2dof, alpha, fls3d, docatrk; 
  // -- new variables
  double closetrk, pvlip, pvlips; 
  double bdtpt, bdt, bdtMax; 
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
  double recoYield, muidYield, trigYield, chanYield, candYield;
  //  double genFileYieldE, genYieldE, genChanYieldE, recoYieldE, chanYieldE, muidYieldE, trigYieldE, candYieldE;
  double genFileYieldE, genYieldE, recoYieldE, chanYieldE, muidYieldE, trigYieldE, candYieldE;
  double absNoCutsYield, ana0Yield,  anaYield,  anaMuonYield,  anaTriggerYield,  anaWmcYield;
  double ana0YieldE, anaYieldE, anaMuonYieldE, anaTriggerYieldE, anaWmcYieldE;
  //  double cFrac, cFracE, acc, accE, accChan, accChanE;
  double acc, accE;
  double effChan, effChanE; 
  double effMuidMC, effMuidMCE, effTrigMC, effTrigMCE;
  double effMuidTNP, effMuidTNPE, effTrigTNP, effTrigTNPE;
  double effMuidTNPMC, effMuidTNPMCE, effTrigTNPMC, effTrigTNPMCE;
  double effRhoMuidTNPMC, effRhoMuidTNPMCE, effRhoTrigTNPMC, effRhoTrigTNPMCE;
  double effRhoMuidTNP, effRhoMuidTNPE, effRhoTrigTNP, effRhoTrigTNPE;
  double rhoMuidTNPMC, rhoMuidTNPMCE, rhoTrigTNPMC, rhoTrigTNPMCE;
  double effCand, effCandE, effCandTE; 
  double effAna, effAnaE; 
  double effTot, effTotE, aEffProdMC, aEffProdMCE, effProdMC, effProdMCE, effProdTNP, effProdTNPE, effProdTNPMC, effProdTNPMCE; 
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
  double pls, plsE, plsTE, pld, pldE, pldTE;
  double phs, phsE, phsTE, phd, phdE, phdTE;
  double psd, psdE, psdTE, pds, pdsE, pdsTE;
  double mBdLo, mBdHi, mBsLo, mBsHi;
  // -- total errors: quad sum of stat and syst errors
  double accTE, effAnaTE, effTotTE; 
  double tauTE; 
  double fitYieldTE;
  double effMuidTE, effMuidTNPTE, effMuidMCTE, effMuidTNPMCTE, effTrigTE, effTrigMCTE, effTrigTNPTE, effTrigTNPMCTE, 
    effRhoMuidTNPTE, effRhoMuidTNPMCTE, effRhoTrigTNPTE, effRhoTrigTNPMCTE; 
  // -- new numbers
  double fBgPeakLo,   fBgPeakHi,   fBgPeakBs,   fBgPeakBd;    // (B+ normalized peaking background)    
  double fBgPeakLoE1, fBgPeakHiE1, fBgPeakBsE1, fBgPeakBdE1;  
  double fBgPeakLoE2, fBgPeakHiE2, fBgPeakBsE2, fBgPeakBdE2;  
  double fBgRslLo,    fBgRslHi,    fBgRslBs,    fBgRslBd;     // (B+ normalized rare sl background)    
  double fBgRslLoE1,  fBgRslHiE1,  fBgRslBsE1,  fBgRslBdE1;   
  double fBgRslLoE2,  fBgRslHiE2,  fBgRslBsE2,  fBgRslBdE2;   
  double fBgRareLo,   fBgRareHi,   fBgRareBs,   fBgRareBd;    // (B+ normalzied rare sl plus peaking background)   
  double fBgRareLoE1, fBgRareHiE1, fBgRareBsE1, fBgRareBdE1;  
  double fBgRareLoE2, fBgRareHiE2, fBgRareBsE2, fBgRareBdE2;  
  double fBgRslsLo,   fBgRslsHi,   fBgRslsBs,   fBgRslsBd;    // (scaled rare sl background)    
  double fBgRslsLoE1, fBgRslsHiE1, fBgRslsBsE1, fBgRslsBdE1;  
  double fBgRslsLoE2, fBgRslsHiE2, fBgRslsBsE2, fBgRslsBdE2;  
  double fBgCombLo,   fBgCombHi,   fBgCombBs,   fBgCombBd;    // combinatorial background
  double fBgCombLoE1, fBgCombHiE1, fBgCombBsE1, fBgCombBdE1;                            
  double fBgCombLoE2, fBgCombHiE2, fBgCombBsE2, fBgCombBdE2;                            
  double fBgNonpLo,   fBgNonpHi,   fBgNonpBs,   fBgNonpBd;    // (scaled non-peaking background)
  double fBgNonpLoE1, fBgNonpHiE1, fBgNonpBsE1, fBgNonpBdE1;                            
  double fBgNonpLoE2, fBgNonpHiE2, fBgNonpBsE2, fBgNonpBdE2;                            
  double fBgTotLo,    fBgTotHi,    fBgTotBs,    fBgTotBd;    // (sum of the above)      
  double fBgTotLoE1,  fBgTotHiE1,  fBgTotBsE1,  fBgTotBdE1;                             
  double fBgTotLoE2,  fBgTotHiE2,  fBgTotBsE2,  fBgTotBdE2;                             
  double fSgLo,       fSgHi,       fSgBs,       fSgBd;       // (Bs signal yield)       
  double fSgLoE1,     fSgHiE1,     fSgBsE1,     fSgBdE1;                                
  double fSgLoE2,     fSgHiE2,     fSgBsE2,     fSgBdE2;                                
  double fBdLo,       fBdHi,       fBdBs,       fBdBd;       // (Bd signal yield)       
  double fBdLoE1,     fBdHiE1,     fBdBsE1,     fBdBdE1;                                
  double fBdLoE2,     fBdHiE2,     fBdBsE2,     fBdBdE2;                                
  double fSgTot,      fBdTot,      fNoTot,      fCsTot;      // (total yields)
  double fSgTotE1,    fBdTotE1,    fNoTotE1,    fCsTotE1; 
  double fSgTotE2,    fBdTotE2,    fNoTotE2,    fCsTotE2; 
  double fFitSg,      fFitBd,      fFitNo,      fFitNoC,     fFitCs,    fFitCsC;   // (fitted yields) 
  double fFitSgE1,    fFitBdE1,    fFitNoE1,    fFitNoCE1,   fFitCsE1,  fFitCsCE1;                    
  double fFitSgE2,    fFitBdE2,    fFitNoE2,    fFitNoCE2,   fFitCsE2,  fFitCsCE2;                    
  double fSgAndBgLo,  fSgAndBgHi,  fSgAndBgBs,  fSgAndBgBd;  // (sum of signal and backgrounds)       
  double fSgAndBgLoE1,fSgAndBgHiE1,fSgAndBgBsE1,fSgAndBgBdE1;                                         
  double fSgAndBgLoE2,fSgAndBgHiE2,fSgAndBgBsE2,fSgAndBgBdE2;                                         
  double fObsLo,      fObsHi,      fObsBs,      fObsBd;      // (observed counts)

};

// ----------------------------------------------------------------------
class plotClass: public TObject {

public:

  plotClass(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default", int mode = 0);
  virtual ~plotClass();

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

  void accEffFromEffTree(std::string fname, std::string dname, numbers &a, cuts &b, int proc);
  void accEffFromEffTreeBac(std::string fname, std::string dname, numbers &a, cuts &b, int proc);
  void filterEfficiency(std::string fname, std::string name);
  void normYield(TH1 *h, int mode, double lo = 5.15, double hi=5.5, double preco=-1.);
  void csYield(TH1 *h, int mode, double lo = 5.25, double hi=5.6, double preco=-1.);
  void dpYield(TH1 *h, int mode, double lo=5.1, double hi=5.6, int bdtc=-10);
  void bgBlind(TH1 *h, int mode = 2, double lo = 4.5, double hi = 6.5); 

  // A modified class which includes the landau peak 14/9/12
  void normYield2(TH1 *h, int mode, double lo = 5.15, double hi=5.5, double preco=-1.);
  // A modified class which includes 2 related gaussians  28/9/12
  void csYield2(TH1 *h, int mode, double lo = 5.25, double hi=5.6, double fraction = -1., double preco=-1.);

  void singleEventPrintout(std::string suffix, std::string st, int ievt);
  void printNumbers(numbers &a, ostream &OUT);
  void initNumbers(numbers *a, bool initAll = true); 
  int  detChan(double m1eta, double m2eta);
  void checkAgainstDuplicates(std::string mode); 
  void reduceTree(TTree *t);
  double quadraticSum(int n, ...); 

  virtual void dumpSamples();
  virtual void dumpCutNames(const char *h);
  virtual void readCuts(const char *filename); 
  virtual void printCuts(ostream &OUT); 
  virtual double recalcMass(double m1, double m2);

  virtual void saveHist(TH1* h, std::string name);
  virtual void setErrors(TH1D *h);
  virtual void stamp(double x1, std::string text1, double x2, std::string text2);
  virtual std::string scientificTex(double n, double nE, const char *name, double base = 1.e-2, int digits = 2);
  virtual std::string formatTex(double n, std::string name, int digits, int sgn = 0);
  virtual std::string formatTex(double n, std::string name, std::string tag);
  virtual std::string formatTex(std::string s, std::string name);
  virtual void drawArrow(double height, int mode = 0, double y = 0.1);
  virtual void drawBox(int mode, double hi = 0.5, double ylo = 0.01);
  virtual void replaceAll(std::string &s, std::string a, std::string b);
  virtual void makeCanvas(int i = 3);
  virtual void newLegend(double x1, double y1, double x2, double y2, std::string title = "");
  virtual double getValueByLabel(TH1D *h, std::string label); 
  virtual void rmSubString(std::string &sinput, const std::string &remove);
  virtual void rmPath(string &sInput);
  virtual double getMaximum(TH1 *h1, TH1 *h2);

  std::string fFiles; 
  // -- Files for Signal and Normalization modes in data and MC
  std::map<std::string, TFile*> fF; 
  std::map<std::string, double> fLumi, fFilterEff, fBF, fBFE, fProdR; 
  std::map<std::string, std::string> fName; 

  std::map<std::string, double> fNgen;

  int fShow, fYear; 
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

  double fAccPt, fAccEtaGen, fAccEtaRec;


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


  std::vector<TH1D*> fhGenAndAccNumbers;
  std::vector<TH1D*> fhMassAbsNoCuts;
  std::vector<TH1D*> fhMassNoCuts;
  std::vector<TH1D*> fhMassNoCutsManyBins;
  std::vector<TH1D*> fhMassWithAnaCuts;
  std::vector<TH1D*> fhMassWithAnaCutsManyBins;
  std::vector<TH1D*> fhMassWithMuonCuts; 
  std::vector<TH1D*> fhMassWithMuonCutsManyBins; 
  std::vector<TH1D*> fhMassWithTriggerCuts; 
  std::vector<TH1D*> fhMassWithTriggerCutsManyBins; 
  std::vector<TH1D*> fhMassWithAllCuts;
  std::vector<TH1D*> fhMassWithAllCutsManyBins; 
  std::vector<TH1D*> fhMassWithAllCutsBlind;

  std::vector<TH1D*> fhW8MassWithAllCuts;
  std::vector<TH1D*> fhW8MassWithAllCutsManyBins; 
  std::vector<TH1D*> fhW8MassWithAllCutsBlind;

  std::vector<TH1D*> fhMassWithMassCuts;
  std::vector<TH1D*> fhMassWithMassCutsManyBins; 
  std::vector<TH1D*> fhNorm, fhNormC, fhDstarPi;
  std::vector<TH1D*> fhMuId, fhMuIdMC;
  std::vector<TH1D*> fhMuTr, fhMuTrMC;
  std::vector<TH2D*> fhBdtMass;

  std::vector<TH2D*> fhAccAll, fhAccPass; 
  std::vector<TH1D*> fhAccPtAll, fhAccPtPass, fhAccEtaAll, fhAccEtaPass; 

  std::vector<numbers*> fNumbersBs, fNumbersBd, fNumbersNo, fNumbersCs, fNumbersBla; 

  bool fDoPrintSingleEvent;
  int  fPrintSingleEvt, fPrintSingleRun;
  bool fDoUseBDT, fDoApplyMuonPtCuts;
  bool fDoApplyCowboyVeto, fDoApplyCowboyVetoAlsoInSignal; 
  bool fInvertedIso;
  bool fNormProcessed; 

  double fBsBgExp, fBsBgExpE, fBdBgExp, fBdBgExpE, fLoBgExp, fHiBgExp; // total expected bg
  double fBsSlBgExp, fBsSlBgExpE, fBdSlBgExp, fBdSlBgExpE, fLoSlBgExp, fHiSlBgExp; // (scaled) rare semileptonic bg
  double fBsCoBgExp, fBsCoBgExpE, fBdCoBgExp, fBdCoBgExpE, fLoCoBgExp, fHiCoBgExp; // combinatorial bg
  double fBgHist, fBgHistE, fBgHistLo, fBgHistHi; 
  double fBgHistExp, fBgHistExpE;

  double fNoSig, fNoSigE; 
  double fCsSig, fCsSigE; 
  double fDpSig, fDpSigE; 
  
  double fCsKstFrac;
  double fNoErrTurnon;
  
  double fu, fs, fsfu, fsfuE;


  // -- stuff to run over the tree from any derived class (using loopFunction() to hook in specifics)
  virtual TTree* getTree(std::string mode);
  virtual void setupTree(TTree *t, std::string mode); 
  virtual void loopOverTree(TTree *t, std::string mode, int function, int nevts = -1, int nstart = 0); 
  virtual void candAnalysis(int mode);
  virtual void loopFunction(int function, int mode = 0) {std::cout << "replace me" << std::endl;} 
  
  virtual TMVA::Reader* setupReader(std::string xmlFile, readerData &rd);
  virtual void calcBDT();
  
  struct RedTreeData fb; 
  //  std::vector<TMVA::Reader*> fReaderEvents0; 
  //  std::vector<TMVA::Reader*> fReaderEvents1; 
  //  std::vector<TMVA::Reader*> fReaderEvents2; 
  TMVA::Reader* fReaderEvents0[2]; 
  TMVA::Reader* fReaderEvents1[2]; 
  TMVA::Reader* fReaderEvents2[2]; 
  bool fIsMC, fIsSignal;
  double fBDT; 
  readerData frd; 
  double MASSMIN, MASSMAX, SIGBOXMIN, SIGBOXMAX, BGLBOXMIN, BGLBOXMAX, BGHBOXMIN, BGHBOXMAX; 
  
  bool fGoodAcceptance, fPreselection, fWideMass, fGoodHLT, fGoodMuonsID, 
       fGoodBdtPt, fGoodMuonsPt, fGoodMuonsEta, fGoodTracks, fGoodTracksPt, fGoodTracksEta; 
  bool fGoodQ, fGoodPvAveW8, fGoodLip, fGoodLipS, fGoodIp, fGoodIpS, fGoodMaxDoca,
       fGoodPt, fGoodEta, fGoodAlpha, fGoodFLS, fGoodChi2, fGoodIso;
  bool fGoodCloseTrack, fGoodDocaTrk, fGoodJpsiCuts, fGoodLastCut; 

  double fW8, fW8MisId, fW8MmuID, fW8Mtrig, fW8DmuID, fW8Dtrig;

  int fRunMin, fRunMax; // if you want to look at a specific run range
  
  AnalysisCuts fAnaCuts; 
  
  std::string fSetup; 
  bool fSaveSmallTree, fSaveLargerTree; 
  bool fIsCowboy;

  PidTable *fptT1, *fptT2; 
  PidTable *fptM; 

  PidTable *fptT1MC, *fptT2MC;
  PidTable *fptMMC; 

  // -- split into seagull and cowboys
  PidTable *fptSgT1, *fptSgT2;
  PidTable *fptSgM; 

  PidTable *fptSgT1MC, *fptSgT2MC;
  PidTable *fptSgMMC; 

  PidTable *fptCbT1, *fptCbT2;
  PidTable *fptCbM; 

  PidTable *fptCbT1MC, *fptCbT2MC;
  PidTable *fptCbMMC; 
  
  PidTable *fptFakePosKaons, *fptFakePosPions, *fptFakePosProtons;
  PidTable *fptFakeNegKaons, *fptFakeNegPions, *fptFakeNegProtons;

  double fEpsilon; 

  ClassDef(plotClass,1) //Testing plotClass
};

//TMVA::Reader* setupReader(std::string xmlFile, readerData& rd);


#endif

