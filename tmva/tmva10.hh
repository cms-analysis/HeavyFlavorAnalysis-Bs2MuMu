#ifndef TMVA10
#define TMVA10

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

struct bdtSetup {
  int NTrees, nEventsMin, MaxDepth, NNodesMax;
  int nCuts;
  float AdaBoostBeta; 
};


struct readerData {
  // -- tight muon variables
  float glbNChi2;
  int   validMuonHits, nMatchedStations, validPixelHits, trkLayerWithHits;
  // -- Luca's additional variables
  float trkValidFract; 
  float pt, eta; 
  float segComp, chi2LocMom, chi2LocPos, glbTrackProb;
  float NTrkVHits, NTrkEHitsOut;
  // -- Mario's and other additional variables
  float kink;
  float dpt, dptrel, deta, dphi, dr; 
};


class tmva10: public TObject {

  public:

  tmva10(std::string vars = "all");
  ~tmva10();

  TCanvas* getC0();

  void createInputFile(std::string ofile, 
		       std::string sgfile = "/shome/ursl/scratch/muonidtree-bstomumu.root", 
		       std::string bgfile = "/shome/ursl/scratch/muonidtree-hadrons.root");

  void train(std::string oname = "TMVA-muid-0", std::string filename = "/scratch/ursl/tmva-trees.root");

  TH1D* getRanking(std::string fname, std::string prefix, std::string after);
  
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

  bdtSetup fBdtSetup;

  std::string fVariables, fBDTParameters; 

  TLegend *legg;
  TLegendEntry *legge;
  TLatex *tl; 
  
  ClassDef(tmva10,1) //Testing tmva10
};

#endif

