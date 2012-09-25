#ifndef PLOTREDUCEDOVERLAYS
#define PLOTREDUCEDOVERLAYS

#include "plotClass.hh"
#include "plotBDT.hh"

#include <TDirectory.h>
#include <TH1.h>
#include <TTree.h>

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "../macros/AnalysisDistribution.hh"


class plotReducedOverlays: public plotClass {

public:

  plotReducedOverlays(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default");
  ~plotReducedOverlays();

  virtual void loopFunction(int mode); 
  virtual void loopFunction1(); 


  void makeAll(std::string selection = "Presel"); 
  void makeSampleOverlay(std::string sample1, std::string sample2, std::string channel, std::string selection);
  void makeSample(std::string sample1, std::string selection);
  void makeOverlay(std::string sample1, std::string sample2, std::string channel, std::string selection);
  void allSystematics();

  void systematics(std::string sample1, std::string selection, int chan = 0);
  void bookDistributions(std::string sample);
  void fillDistributions();
  void sbsDistributions(std::string sample, std::string selection = "Presel", std::string what = "");
  void overlay(std::string sample1, std::string sample2, std::string selection = "Presel", std::string what = ""); 
  
  AnalysisDistribution* bookDistribution(std::string hn, std::string ht, std::string hc, int nbins, double lo, double hi); 

  //  void loop(TTree *t, string mode, int offset);

  TDirectory *fHistDir; 

#define NAD 2
  AnalysisDistribution   
  *fpMuonsEta[NAD], *fpMuon1Pt[NAD], *fpMuon2Pt[NAD]
    , *fpPt[NAD], *fpP[NAD], *fpEta[NAD] 
    , *fpAlpha[NAD]
    , *fpIso[NAD], *fpCloseTrk[NAD], *fpDocaTrk[NAD]   
    , *fpChi2Dof[NAD], *fpPChi2Dof[NAD] 
    , *fpFLS3d[NAD], *fpFL3d[NAD], *fpFL3dE[NAD] 
    , *fpMaxDoca[NAD], *fpIp[NAD], *fpIpS[NAD]
    , *fpBDT[NAD]   
    , *fpPvZ[NAD], *fpPvN[NAD], *fpPvAveW8[NAD]
    ;
  int fOffset;
  int fMode;

  vector<std::string> fDoList; 

  ClassDef(plotReducedOverlays,1) //Testing plotReducedOverlays

};



#endif

