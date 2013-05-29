#ifndef TMVA1
#define TMVA1

#include <iostream>
#include <fstream>
#include <vector>
#include <map>


#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"

#include "../macros2/RedTreeData.hh"
#include "../macros2/preselection.hh"

struct bdtSetup {
  int NTrees, nEventsMin, MaxDepth, NNodesMax;
  int nCuts;
  float AdaBoostBeta; 
};


struct readerData {
  float pt, eta, m1eta, m2eta, m1pt, m2pt;
  float fls3d, alpha, maxdoca, pvip, pvips, iso, docatrk, chi2dof, closetrk; 
  float closetrks1, closetrks2, closetrks3; 
  float m1iso, m2iso, pvdchi2, othervtx;
  float pvlips2, pvlip2; 
  float m;
};


struct files {
  std::string dname; 
  std::string sname; 
  std::string rname; 
};


class tmva1: public TObject {

  public:

  tmva1(int year = 2012, std::string vars = "pt:eta:alpha:fls3d:maxdoca:pvip:pvips:iso:m1iso:m2iso");
  ~tmva1();

  TCanvas* getC0();
  void train(std::string oname = "TMVA-0", std::string filename = "/scratch/ursl/tmva-trees.root", int nsg = -1, int nbg = -1);
  void apply(const char *fname = "TMVA-0");
  void analyze(const char *fname = "TMVA-0");
  void makeAll(int offset, std::string filename = "", int clean = 1);
  void make(int offset, std::string filename, int evt, int clean);

  void reAnalyze(int imin, int imax);
  void createInputFile(std::string fname, int randomSeed = -1);
  void cleanup(std::string fname);
  TH1D* getRanking(std::string fname, std::string prefix, std::string after);
  
  void setupTree(TTree *t, RedTreeData &b);
  void calcBDT();
  
  void mvas(std::string fname="TMVA");
  void redrawStats(double x, double y, const char *newname= "newstat", int color=kRed); 

  void setBDTParameters(std::string pars) {fBDTParameters = pars;}
  void setNTrees(int i) {fBdtSetup.NTrees = i;}
  void setnEventsMin(int i) {fBdtSetup.nEventsMin = i;};
  void setMaxDepth(int i) {fBdtSetup.MaxDepth = i;};
  void setnCuts(int i) {fBdtSetup.nCuts = i;};
  void setAdaBoostBeta(float f) {fBdtSetup.AdaBoostBeta = f;}; 
  void setNNodesMax(int f) {fBdtSetup.NNodesMax = f;}; 
  void newLegend(double x1, double y1, double x2, double y2, std::string title = "");
  void writeOut(TFile*, TH1*);

  void setApply0() {fApplyOn0 = true;  fApplyOn1 = false; fApplyOn2 = false;};  // apply on even events, train on odd
  void setApply1() {fApplyOn0 = false; fApplyOn1 = true;  fApplyOn2 = false;};  // apply on odd events, train on even
  void setApply2() {fApplyOn0 = false; fApplyOn1 = false; fApplyOn2 = true; };  // apply on odd events, train on even
  void setTrainAntiMuon(bool yes) {fTrainAntiMuon = yes;};
  void setChannel(int channel) {fChannel = channel;};

  TTree* createTree(struct readerData &rd);
  double bgBlind(TH1 *h, int mode, double lo, double hi); 

  void toyRun(std::string modifier, std::string vars, std::string bdtpars, int seed = 0, int nsg = 15000, int nbg = 10000);
  void createToyData(std::string sgfilename, std::string bgfilename, std::string ofilename, int seed, int nsg, int nbg);

  void analyzeTexFiles(std::string ofilename, int start, int end); 
  void readTexFile(std::string filename, std::vector<std::string> &lines);
  float parseTexLine(std::string line);

  files fInputFiles;
  bdtSetup fBdtSetup;

  bool fApplyOn0, fApplyOn1, fApplyOn2; 
  bool fTrainAntiMuon;
  int fChannel, fYear;
  double fRsigma, fLumiScale; 
  std::string fVariables, fBDTParameters; 

  RedTreeData ftd;  
  readerData frd; 
  double fBDT, fBDT0, fBDT1, fBDT2; 
  std::vector<TMVA::Reader*> fReader; 

  TLegend *legg;
  TLegendEntry *legge;
  TLatex *tl; 
  
  ClassDef(tmva1,1) //Testing tmva1
};

TMVA::Reader* setupReader(std::string xmlFile, readerData &rd);
void setTitles(TH1 *h, const char *sx, const char *sy, 
	       float size = 0.05, float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 42);
void setHist(TH1 *h, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);
void shrinkPad(double b = 0.1, double l = 0.1, double r = 0.1, double t = 0.1);

#endif

