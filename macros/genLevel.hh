#ifndef GENLEVEL_H
#define GENLEVEL_H

#include <iostream>
#include <vector>
#include <utility>

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

#include "treeReader01.hh"

#define PIPMASS 0.13957
#define ELMASS  0.000511
#define MUMASS  0.10567


class genLevel : public treeReader01 {

public:
  genLevel(TChain *tree, TString evtClassName);
  ~genLevel();

  void         bookHist();
  void         startAnalysis();
  void         eventProcessing();
  void         endAnalysis();
  void         fillHist();
  void         readCuts(TString filename, int dump = 1);
  void         initVariables();
  int          muonType(TGenCand *pCand);
  void         bbbarCrossSection();
  void         printBdecays(); 

  int          NTOTAL; 
  double       XSECTION;

};

#endif
