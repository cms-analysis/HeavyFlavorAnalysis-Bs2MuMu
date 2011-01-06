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
#include "TVirtualPad.h"  // access to gPad

#include <iostream>
#include <fstream>
#include <vector>

class anaBmm: public TObject {

public:

  anaBmm(const char *files, const char *dir = ".", int mode = 0);

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
  void effTable(std::string mode);


  // -- Utilities and helper methods
  // -------------------------------
  std::string formatTex(double n, std::string name, std::string tag);
  void replaceAll(std::string &s, std::string a, std::string b);
  int  wait();
  void makeCanvas(int i = 3);

  // -- Files for Signal and Normalization modes in data and MC
#define MAXFILES 20
  TFile *fpData[MAXFILES], *fpMc[MAXFILES];
  int fNData, fNMc;
  std::string fDataString[MAXFILES], fMcString[MAXFILES];
  int fSgData, fSgMc, fNoData, fNoMc;  // the indices for the default files

  int fShow; 
  TString fFile; 


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

  ClassDef(anaBmm,1) //Testing anaBmm
};


#endif
