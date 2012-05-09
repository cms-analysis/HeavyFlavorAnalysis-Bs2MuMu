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

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaJet.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"

#include "candAna.hh"
#include "bmm2Reader.hh"


class candAnaHh : public candAna {
  
public:
  candAnaHh(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAnaHh();

  void        evtAnalysis(TAna01Event *evt);
  bool        anaMC(TAna01Event *evt);
  void        candAnalysis();
  void        moreBasicCuts();

  int         truthMatch(TAnaCand *pC, int verbose = 0); 
  void        dumpHFTruthCand(TAnaCand *pC); 
  void        dumpHFHhCand(TAnaCand *pC); 

  void        readCuts(string filename, int dump);
  
  void        bookHist();

private:
  TTree * tree;
  int fcands, ftm[10];
  float fmds[10], fmdz[10];
  float ffls3d[10],fchi2[10],falpha[10],fm[10],fdr[10];
  float fpt[10],fdoca[10],fweight[10],fptpi1[10],fptpi2[10];
};

#endif