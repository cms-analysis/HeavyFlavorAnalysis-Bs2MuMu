#ifndef CANDANADSTAR_H
#define CANDANADSTAR_H

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

class candAnaDstar : public candAna {
  
public:
  candAnaDstar(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAnaDstar();

  //void        evtAnalysis(TAna01Event *evt);
  bool        anaMC(TAna01Event *evt);
  void        candAnalysis();
  void        moreBasicCuts();

  int         truthMatch(TAnaCand *pC, int verbose = 0); 
  void        dumpHFTruthCand(TAnaCand *pC); 
  void        dumpHFDstarCand(TAnaCand *pC); 

  void        readCuts(string filename, int dump);
  
  void        bookHist();

 

private:
  TTree * tree;
  int ftm, fnclose;
  bool fmuid1, fmuid2, fmumat1, fmumat2;
  float fmds, fmdz;
  float ffls3d,fchi2,falpha,falpha2,fqpis,fdr;
  float fpt,fptdz,fptpis,fptpi,fptk;
  float fpvd, fiso;
  float feta, fetapi, fetak;
  float fchipi, fchik;
  float mudr1, mudr2;
};

#endif
