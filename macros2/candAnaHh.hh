#ifndef CANDANAHH_H
#define CANDANAHH_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include "candAna.hh"

class candAnaHh : public candAna {
  
public:
  candAnaHh(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAnaHh();

  void        candAnalysis();
  void        hhAnalysis();
  void        efficiencyCalculation();
  
  void        processType(); 
  void        genMatch(); 
  void        recoMatch(); 
  void        candMatch(); 
  
  void        bookHist();
  void        readCuts(string filename, int dump);
  
  void        evtAnalysis(TAna01Event *evt);
  bool        anaMC(TAna01Event *evt);

  //void        moreBasicCuts();

  int         truthMatch(TAnaCand *pC, int verbose = 0); 
  void        dumpHFTruthCand(TAnaCand *pC); 
  void        dumpHFHhCand(TAnaCand *pC); 

private:
  TTree * tree;
  int fcands, ftm[10];
  float fmds[10], fmdz[10];
  float ffls3d[10],fchi2[10],falpha[10],fm[10],fdr[10];
  float fpt[10],fdoca[10],fweight[10],fptpi1[10],fptpi2[10];
  int fclose[10];
  float fiso[10], fperp1[10], fperp2[10];
  float fm1[10],fm2[10],fm3[10],fm4[10];

  // Additional variables 
  int MUON_VETO;
  double HH_MLO, HH_MHI;


};

#endif
