#ifndef CANDANABU2JPSIK_H
#define CANDANABU2JPSIK_H

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

class candAnaBu2JpsiK : public candAna {
  
public:
  candAnaBu2JpsiK(bmm2Reader *pReader, std::string name, std::string cutsFile);
  virtual ~candAnaBu2JpsiK();

  void        candAnalysis();
  void        efficiencyCalculation();
  void        moreBasicCuts();
  void        moreReducedTree(TTree *);

  void        genMatch(); 
  void        recoMatch(); 
  void        candMatch(); 
  
  void        readCuts(string filename, int dump);
  void        bookHist();
  void        fillCandidateHistograms(int offset);

  
  int          JPSITYPE; 
  double       JPSIMASSLO, JPSIMASSHI;

  double       fKaonPt, fKaonEta, fKaonPhi;
  double       fKPtGen, fKEtaGen;
  int          fKGenID; 
  double       fKaonPtNrf, fKaonEtaNrf;
  int          fKaonTkQuality;
  double       fJpsiMass, fJpsiPt, fJpsiEta, fJpsiPhi;

  bool         fGoodJpsiMass;
  bool         fKa1Missid, fKa1MuMatch;
  bool         fKa1Missid2, fKa1MuMatch2; 
  float        fKa1MuMatchR, fKa1MuMatchR2, fKa1MuMatchR3, fKa1MuMatchR4, fKa1MuMatchR5, fKa1MuMatchR6;



  // -- TM
  int          fGenK1Tmi; 
  int          fRecK1Tmi; 

  // -- effTree
  float fETk1pt, fETk1eta, fETg3pt, fETg3eta;
  int   fETk1q; 
  bool  fETk1gt;

};

#endif
