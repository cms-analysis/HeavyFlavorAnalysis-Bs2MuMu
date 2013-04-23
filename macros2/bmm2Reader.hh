#ifndef BMM2READER_H
#define BMM2READER_H

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

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/JSON.hh"

#include "../macros/treeReader01.hh"

#define DR      57.29577951

class PidTable; 
class MuScleFitCorrector; 
class candAna;

class bmm2Reader : public treeReader01 {

public:
  bmm2Reader(TChain *tree, TString evtClassName);
  ~bmm2Reader();

  virtual void   startAnalysis();
  virtual void   eventProcessing();
  virtual void   readCuts(TString filename, int dump = 1);
  virtual void   bookHist();

  virtual void   processType();
  virtual void   setYear(int year) {fYear = year;}

  std::vector<candAna*> lCandAnalysis;

//   // -- PidTables
//   PidTable *fpMuonID;
//   PidTable *fpMuonTr, *fpMuonTr1, *fpMuonTr2;
  
//   PidTable *ptSgMUID, *ptCbMUID; 
//   PidTable *ptSgMUT1, *ptCbMUT1; 
//   PidTable *ptSgMUT2, *ptCbMUT2; 

  MuScleFitCorrector *msc; 

  int fProcessType; 
  int fYear; 
};

#endif
