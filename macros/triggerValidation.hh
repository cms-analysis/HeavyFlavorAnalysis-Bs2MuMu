#ifndef TRIGGERVALIDATION_H
#define TRIGGERVALIDATION_H

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

using namespace std;

#define DR      57.29577951
#define PIPMASS 0.13957
#define ELMASS  0.000511
#define MUMASS  0.10567


class triggerValidation : public treeReader01 {

public:
  triggerValidation(TChain *tree, TString evtClassName);
  ~triggerValidation();

  void         bookHist();
  void         startAnalysis();
  void         eventProcessing();
  void         fillHist();
  void         readCuts(TString filename, int dump = 1);
  void         initVariables();
  void         closeHistFile();

  void         l1tValidation();
  void         lttValidation();
  void         hltValidation();

  int          fHLTWasRun[256], fHLTResult[256], fHLTWasRunResult[256];
  TString      fHLTNames[256];

  int          fL1TWasRun[128], fL1TResult[128], fL1TWasRunResult[128];
  TString      fL1TNames[128];

  int          fLTTWasRun[128], fLTTResult[128], fLTTWasRunResult[128];
  TString      fLTTNames[128];

};

#endif
