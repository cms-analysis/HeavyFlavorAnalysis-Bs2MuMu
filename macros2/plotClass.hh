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
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TVirtualPad.h"  // access to gPad

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../macros/initFunc.hh"

struct cuts {
  int index; 
  double mBdLo, mBdHi, mBsLo, mBsHi;
  double etaMin, etaMax, pt; 
  double m1pt, m2pt, m1eta, m2eta;
  double iso, chi2dof, alpha, fls3d, docatrk; 
};


class plotClass: public TObject {

public:

  plotClass(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotClass();

  virtual void cd(const char *file) {fF[file]->cd();}

  // -- initialization and setup
  // ---------------------------
  virtual void init(const char *files, const char *cuts, const char *dir, int mode);
  virtual void loadFiles(const char *files);
  virtual TFile* loadFile(std::string file);

  // -- main methods
  // --------------
  virtual void makeAll(int verbose = 0);

  virtual void dumpSamples();
  virtual void dumpCutNames(const char *h);
  virtual void readCuts(const char *filename); 
  virtual void printCuts(ostream &OUT); 

  virtual void setErrors(TH1D *h);
  virtual void stamp(double x1, std::string text1, double x2, std::string text2);
  virtual std::string scientificTex(double n, double nE, std::string name, double base = 1.e-2, int digits = 2);
  virtual std::string formatTex(double n, std::string name, int digits);
  virtual std::string formatTex(double n, std::string name, std::string tag);
  virtual void drawArrow(double height, int mode = 0, int color = kBlue);
  virtual void drawBox(int mode, double hi = 0.5, int ylo = 0.01);
  virtual void replaceAll(std::string &s, std::string a, std::string b);
  virtual void makeCanvas(int i = 3);
  virtual void newLegend(double x1, double y1, double x2, double y2);

  // -- Files for Signal and Normalization modes in data and MC
  std::map<std::string, TFile*> fF; 
  std::map<std::string, double> fLumi, fFilterEff; 
  std::map<std::string, std::string> fName; 

  int fShow; 
  std::string fFile; 

  // -- output histograms and numbers
  TFile *fHistFile; 
  ofstream fOUT, fTEX; 

  bool fDoPrint; 
  int fVerbose; 

  // -- cuts 
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

  unsigned int fNchan;
  int fChan, fComb, fOver; 

  initFunc  *fpFunc;

  double fMassLo, fMassHi
    , fBgLo, fBgHi
    , fSgLo, fSgHi
    , fNoLo, fNoHi
    , fCsLo, fCsHi
    ;

  ClassDef(plotClass,1) //Testing plotClass
};


#endif

