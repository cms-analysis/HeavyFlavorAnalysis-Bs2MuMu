#ifndef BMMREADER_H
#define BMMREADER_H

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

#define DR      57.29577951

class bmmReader : public treeReader01 {

public:
  bmmReader(TChain *tree, TString evtClassName);
  ~bmmReader();

  void         bookHist();
  void         startAnalysis();
  void         eventProcessing();
  void         fillHist();
  void         readCuts(TString filename, int dump = 1);
  void         initVariables();

  virtual void MCKinematics();  
  virtual void L1TSelection();  
  virtual void HLTSelection();  

  virtual void pvStudy();  

  virtual void trackSelection();  
  virtual void muonSelection();  
  virtual void candidateSelection(int mode = 0);  // 0 = closest in r-phi

  // -- Cut values
  double 
      BSPTLO
    , BSETALO
    , BSETAHI   
    , MUPTLO
    , MUPTHI
    , MUETALO
    , MUETAHI   
    ;
  int TYPE;

  // -- Variables
  TAnaCand    *fpCand; 

  double      fCandPt, fCandMass;

  bool        fGoodMCKinematics, fGoodL1, fGoodHLT, fGoodEvent;
  bool        fGoodMuonsID, fGoodMuonsPT;
  bool        fGoodTracks, fGoodTracksPT;
  bool        fGoodCandPT;

};

#endif
