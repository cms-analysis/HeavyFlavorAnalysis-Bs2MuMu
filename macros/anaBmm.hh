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
#include "TH2.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TVirtualPad.h"  // access to gPad

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <vector>

class anaBmm: public TObject {

public:

  anaBmm(const char *files = "nada");

  // -- initialization and setup
  // ---------------------------
  void init(const char *files);
  void loadFiles();
  void loadFiles(const char *filename);
  void loadDa(const char *filename, double vxs, const char *sign, const char *type);
  void loadMc(const char *filename, double vxs, const char *sign, const char *type);
  void loadSg(const char *filename, double vxs, const char *sign, const char *type);
  void dumpFiles();

  // -- main methods
  // --------------
  void makeAllPlots();
  void jreco(int o = 4); 
  void breco(int o = 2, const char *hist = "c030"); 
  void nreco(int o = 5, const char *hist = "c030");
  void bgOverlay(const char *s = "c030", const int npers = 6);
  void nbgOverlay(const char *s = "c030", const int npers = 6);
  void dumpCuts();
  void effTables();
  void effTable(TFile *f, const char *tag); 
  void showDistributions(int offset = 0, int wait = 0); 
  void showDistribution(const char *hname, int mode = 0, double x = -9999., double y = -9999.);
  void calculateUpperLimit(); 
  void normalizedUpperLimit();
  double massReduction(const char *s = "c030", const char *tag = "m0");

  void sgRecos();
  void combiRecos();
  void bgOverlays();
  void fakeMuons();

  void signalPlots();

  void plotAll(const char *s = "c030");
 

  // -- temporary 'working' methods
  void loopOptimization(double pt = 4.); 
  void optimizerNumbers(const char *precuts, const char *line, double &eff, double &sg, double &bg);
  void runOptimization(const char *aoCuts, const char *extraVar, int nbin, double min, double max);
  void handOptimization(const char *aoCuts, const char *extraCuts); 


  // -- Special studies
  // ------------------
  void showProcesses(int signal = 1);
  void plotProcesses(const char *cut, const char *var, const char *axis, int bin, double lo, double hi, int legend = 1);

  void mcValidation(int wait = 0); 
  void mcVal(const char *hname); 

  // -- Utilities and helper methods
  // -------------------------------
  double  expUL(double sig, double bg);
  TH1* sumHistMC(const char *hist, int draw = 0, const char *sel = "c0");
  void sumHistMC_Add(int ch, const char *hist, TH1 *hist1, TH1 *hist2);
  void channelEff(TFile *f, double &fnorm_Ch, double &eff1_Ch, double &eff2_Ch);
  TH1D* DivideHisto(TH1D *hist1, TH1D *hist2);
  TH2D* DivideHisto2(TH2D *hist1, TH2D *hist2);
  void writeFitPar(TF1 *f, int o, double mean, double sigma, int npar);
 
  int  wait();

  void makeCanvas(int i = 3);
  void shrinkPad(double b = 0.1, double l = 0.1, double r = 0.1, double t = 0.1);
  void setTitles(TH1 *h, const char *sx, const char *sy, 
		 float size = 0.05, float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 132);
  void setTitles2(TH2 *h, const char *sx, const char *sy, 
		 float size = 0.05, float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 132);
  void setHist(TH1 *h, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);
  void emptyBinError(TH1 *h);
  void setFilledHist(TH1 *h, int lcol = kBlack, int fcol = kYellow, int fstyle = 1000, int width = 1, int style = 1);

  TString texForm(double e);
  TString texForm2(double e);
  TString texForm31(double e);
  TString formatTex(double n, const char *name, const char *tag);
  TString scaleFactor(int exp);
  double dBinomial(double n, double N);
  double dEff(int in, int iN);
  double nul(int nobs); 
  double barlow(int nobs, double bg = 0., double bgE = 0., double sE = 0.);
  void   getSignature(TString in, TString &out, TString &out2);
  TString getSubGroup(TString in);
  int checkIf(int mc, const char *sel);
  double getMisID(TString in);

  // -- Files for Signal, Data, Monte Carlo, and Control Samples
  int nSg, nMc, nDa;
  TFile *fS[10], *fD[10], *fM[30], *fCS[10];

  int fShow; 
  TString fFile; 
  TString fNumbersFileName;


private:

  // -- for normalization
  double fvXsS[10],  fvXsD[10],  fvXsM[30];  // this is entered by hand from NSEL/NGEN
  double fLumiS[10], fLumiD[10], fLumiM[30]; // this is computed in loadFiles()
  double fNevtS[10], fNevtD[10], fNevtM[30]; // this is filled in loadFiles()
  double fNexpS[10], fNexpD[10], fNexpM[30]; // this is filled in loadFiles()
  double fMisIdM[30];                         // this is filled in muonMisId()
  TString fSignS[10], fSignD[10], fSignM[30]; // this is filled in loadFiles()
  TString fTypeS[10], fTypeD[10], fTypeM[30]; // this is filled in loadFiles()
  TString fSignTexS[10], fSignTexD[10], fSignTexM[30]; // this is filled in loadFiles()
  TString fSignLeggS[10], fSignLeggD[10], fSignLeggM[30]; // this is filled in loadFiles()

  // -- misc
  double fFom;

  double fMassBs, fMassBp;

  // -- Upper Limit determination
  double fEsg0, fEsgE0;       // efficiency and error (w/o fact. in 5-6 GeV window)
  double fEsg, fEsgE;         // efficiency and error
  double fNsg, fNsgE;         // expected number of signal events 
  double fNbg, fNbgE;         // expected number of background events 
  double fNrbg, fNrbgE;       // expected number of rare background events

  double fEsg_norm, fEsgE_norm;         // efficiency and error (norm. channel)
  double fNsg_norm, fNsgE_norm;         // expected number of signal events (norm. channel)
  double fNbg_norm, fNbgE_norm;         // expected number of background events (norm. channel)

  int sgIndex, bgIndex, normSgIndex, normBgIndex; 
  
  double fMu
    , fPi
    , fKa
    , fProt
    ;
  
  double fBp_max
    ,fBs_max;

  TH1D *fMuEff
    , *fPiMid
    , *fKaMid
    , *fProtMid
    ;

  TH2D *fMuEff2
    , *fPiMid2
    , *fKaMid2
    , *fProtMid2
    ;
  
  // -- functions
  TF1 *f0, *f1, *f2, *f3, *f4, *f5, *f6; 
  TF1 *f10, *f11;

  // -- totals counters for background table
  double fBgTotal, fBgTotalE;
  double fSgTotal, fSgTotalE;

  // -- Display utilities
  int fFont; 
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5;
  TLatex *tl; 
  TBox *box;
  TArrow *pa; 
  TLine *pl; 
  TLegend *legg;
  TLegendEntry *legge;

  char inDir[200];
  char outDir[200];
  char line[200];

  char fSuffix[100]; 

  ClassDef(anaBmm,1) //Testing anaBmm
};


#endif

