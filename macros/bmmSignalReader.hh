#ifndef BMMSIGNALREADER_H
#define BMMSIGNALREADER_H

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

#include "bmmReader.hh"

class bmmSignalReader : public bmmReader {

public:
  bmmSignalReader(TChain *tree, TString evtClassName);
  ~bmmSignalReader();

  void         bookHist();
  void         startAnalysis();
  void         eventProcessing();
  void         fillHist();
  void         readCuts(TString filename, int dump = 1);
  void         initVariables();
  void         candidateSelection(int mode = 0); 

};

#endif
