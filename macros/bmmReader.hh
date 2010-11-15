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

  virtual void bookHist();
  virtual void startAnalysis();
  virtual void eventProcessing();
  virtual void fillHist();
  virtual void readCuts(TString filename, int dump = 1);
  virtual void initVariables();
  virtual void insertCand(TAnaCand* pCand);

  virtual void pvStudy();  

  virtual void MCKinematics();  
  virtual void L1TSelection();  
  virtual void HLTSelection();  
  virtual void trackSelection();  
  virtual void muonSelection();  
  virtual void candidateSelection(int mode = 0);  // 0 = closest in r-phi
  virtual void fillCandidateVariables(); 


  // -- Cut values
  double 
      CANDPTLO
    , CANDETALO
    , CANDETAHI   
    , TRACKPTLO
    , TRACKPTHI
    , TRACKTIP
    , TRACKLIP
    , MUPTLO
    , MUPTHI
    , MUETALO
    , MUETAHI   
    ;
  int TYPE, MUID, TRACKQUALITY;

  // -- Variables
  bool                    fGoodMCKinematics;
  bool                    fGoodL1, fGoodHLT;

  // -- vectors with the candidates of the specified type
  std::vector<TAnaCand *> fCands;  
  TAnaCand               *fpCand;       // the 'best' candidate

  std::vector<bool>       fGoodMuonsID, fGoodMuonsPT;
  std::vector<bool>       fGoodTracks, fGoodTracksPT;
  std::vector<bool>       fGoodCandPT;

  bool                    fGoodEvent;

  // -- variables for reduced tree, they are from fpCand
  double                  fCandPt, fCandM; 

  double       SIGBOXMIN, SIGBOXMAX; 
  double       BGLBOXMIN, BGLBOXMAX; 
  double       BGHBOXMIN, BGHBOXMAX; 
};

#endif
