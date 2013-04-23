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
  bool        anaMC();
  void        candAnalysis();
  void        moreBasicCuts();

  int         truthMatch(TAnaCand *pC, int verbose = 0); 
  void        dumpHFTruthCand(TAnaCand *pC); 
  void        dumpHFDstarCand(TAnaCand *pC); 
  void        readCuts(string filename, int dump);
  void        bookHist();
  void        genMatch();
  void        recoMatch();
  void        candMatch();
  void        efficiencyCalculation();
  int         doTest(TAnaCand *pC, int mode =-1); 
  int         getJpsi(int &idx1, int &idx2);

  double      doTriggerMatchingTest(TAnaTrack *pt, int trig = 0); // match a single track to HLT

private:
  TTree * tree;
  int ftm, fnclose;
  bool fmuid1, fmuid2, fmumat1, fmumat2;
  float fmds, fmdz;
  float ffls3d,fchi2,falpha,fqpis,fdr;
  float fpt,fptdz,fptpis,fptpi,fptk;
  float fpvd, fiso;
  float feta, fetapi, fetak;
  float fchipi, fchik;
  float mudr1, mudr2;  
  float match1dr, match2dr;
  int fnumHltMuons, fnumHltPureMuons; // for number of hlt matched muons
  bool fmuidmva1, fmuidmva2; // MVA muon id
  double fmva1, fmva2; // MVA
  // test
  float match1dr1, match2dr1,match1dr2, match2dr2,match1dr3, match2dr3,match1dr4, match2dr4;
  vector<int> hltMatchedMuons;
  bool fMatchedPion,fMatchedKaon,fRejectPion,fRejectKaon;
  int foundJpsis;
};

#endif
