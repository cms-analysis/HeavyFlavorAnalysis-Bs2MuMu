#ifndef MUCHARMREADER_H
#define MUCHARMREADER_H

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


class muCharmReader : public treeReader01 {

public:
  muCharmReader(TChain *tree, TString evtClassName);
  ~muCharmReader();

  void         bookHist();
  void         startAnalysis();
  void         eventProcessing();
  void         fillHist();
  void         readCuts(TString filename, int dump = 1);
  void         initVariables();

  virtual void MCKinematics();  
  virtual void L1TSelection();  
  virtual void HLTSelection();  
  virtual void trackSelection();  
  virtual void muonSelection();  
  virtual void candidateSelection(int mode = 0);  // 0 = closest in r-phi
  virtual void fillTMCand(TAnaCand *pCand, int type);

  // -- Cut values
  double 
      CHARMPTLO
    , CHARMETALO
    , CHARMETAHI   
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
