#ifndef TMVA2
#define TMVA2

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
  float m;
};


struct files {
  std::string dname; 
  std::string sname; 
  std::string rname; 
};


class tmva2: public TObject {

  public:

  tmva2();
  ~tmva2();

  TCanvas* getC0();
  void train(std::string oname = "TMVA-0", std::string filename = "/scratch/ursl/tmva-trees.root");
  void apply(const char *fname = "TMVA-0");
  void analyze(const char *fname = "TMVA-0");
  void makeAll(int offset, std::string filename = "/scratch/ursl/tmva-trees.root", int clean = 1);
  void make(int offset, std::string filename, int evt, int clean);

  void reAnalyze(int imin, int imax);
  void createInputFile(std::string fname);
  void cleanup(std::string fname);
  TH1D* getRanking(std::string fname, std::string prefix, std::string after);
  
  void setupTree(TTree *t, RedTreeData &b);
  void calcBDT();
  
  void mvas(std::string fname="TMVA");
  void redrawStats(double x, double y, const char *newname= "newstat", int color=kRed); 

  void setNTrees(int i) {fBdtSetup.NTrees = i;}
  void setnEventsMin(int i) {fBdtSetup.nEventsMin = i;};
  void setMaxDepth(int i) {fBdtSetup.MaxDepth = i;};
  void setnCuts(int i) {fBdtSetup.nCuts = i;};
  void setAdaBoostBeta(float f) {fBdtSetup.AdaBoostBeta = f;}; 
  void setNNodesMax(float f) {fBdtSetup.NNodesMax = f;}; 
  void newLegend(double x1, double y1, double x2, double y2, std::string title = "");
  void writeOut(TFile*, TH1*);

  void setApply0() {fApplyOn0 = true;  fApplyOn1 = false; fApplyOn2 = false;};  // apply on even events, train on odd
  void setApply1() {fApplyOn0 = false; fApplyOn1 = true;  fApplyOn2 = false;};  // apply on odd events, train on even
  void setApply2() {fApplyOn0 = false; fApplyOn1 = false; fApplyOn2 = true; };  // apply on odd events, train on even
  void setTrainAntiMuon(bool yes) {fTrainAntiMuon = yes;};
  void setChannel(int channel) {fChannel = channel;};

  void toyRuns(std::string ifilename, int nruns = 10); 
  void createToyData(std::string ifilename, std::string ofilename, int seed = 100, int nsg=20000, int nbg=25000);
  void trainOnToyData(std::string iname = "/scratch/ursl/tmva-toy-100.root", std::string oname = "toy-100"); 
  TTree* createTree(struct readerData &rd);


  files fInputFiles;
  bdtSetup fBdtSetup;

  bool fApplyOn0, fApplyOn1, fApplyOn2; 
  bool fTrainAntiMuon;
  int fChannel;
  double fRsigma; 

  RedTreeData ftd;  
  readerData frd; 
  double fBDT, fBDT0, fBDT1, fBDT2; 
  std::vector<TMVA::Reader*> fReader; 

  TLegend *legg;
  TLegendEntry *legge;
  TLatex *tl; 
  
  ClassDef(tmva2,1) //Testing tmva2
};

TMVA::Reader* setupReader(std::string xmlFile, readerData &rd);
void setTitles(TH1 *h, const char *sx, const char *sy, 
	       float size = 0.05, float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 42);
void setHist(TH1 *h, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);

#endif

