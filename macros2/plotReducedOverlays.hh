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

  plotReducedOverlays(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default", int mode = 0);
  ~plotReducedOverlays();

  virtual void loopFunction(int function, int mode); 
  virtual void loopFunction1(int mode); 


  void makeAll(std::string selection = "Presel"); 
  void makeSampleOverlay(std::string sample1, std::string sample2, std::string channel, std::string selection);
  void makeSample(std::string sample1, std::string selection, std::string = "B");

  void makeOverlay(std::string sample1, std::string sample2, std::string channel, std::string selection);
  void allSystematics();
  void systematics(std::string sample1, std::string selection, int chan = 0);

  // -- 
  void bookDistributions(std::string sample);
  void fillDistributions();
  void sbsDistributions(std::string sample, std::string selection = "Presel", std::string what = "");
  void overlay(std::string sample1, std::string sample2, std::string selection = "Presel", std::string what = ""); 
  
  AnalysisDistribution* bookDistribution(std::string hn, std::string ht, std::string hc, int nbins, double lo, double hi); 
  AnalysisDistribution* bookSpecialDistribution(string hn, string ht, string hc, int nbins, double lo, double hi, bool *presel);

  //  void loop(TTree *t, string mode, int offset);

  TDirectory *fHistDir; 

#define NAD 2
  AnalysisDistribution   
  *fpMuonsEta[NAD], *fpMuon1Pt[NAD], *fpMuon2Pt[NAD]
    , *fpPt[NAD], *fpP[NAD], *fpPz[NAD], *fpEta[NAD] 
    , *fpAlpha[NAD]
    , *fpIso[NAD], *fpCloseTrk[NAD], *fpDocaTrk[NAD]   
    , *fpChi2Dof[NAD], *fpPChi2Dof[NAD] 
    , *fpFLS3d[NAD], *fpFL3d[NAD], *fpFL3dE[NAD] 
    , *fpMaxDoca[NAD], *fpIp[NAD], *fpIpS[NAD]
    , *fpBDT[NAD]   
    , *fpBDTSel0[NAD], *fpBDTSel1[NAD], *fpBDTSel2[NAD]
    , *fpPvZ[NAD], *fpPvN[NAD], *fpPvAveW8[NAD]
    ;
  int fOffset;
  int fMode;

  vector<std::string> fDoList; 


  bool fSel0, fSel1, fSel2; 

  ClassDef(plotReducedOverlays,1) //Testing plotReducedOverlays

};



#endif

