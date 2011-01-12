#ifndef BMMBS2JPSIPHIREADER_H
#define BMMBS2JPSIPHIREADER_H

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

class bmmBs2JpsiPhiReader : public bmmReader {

public:
  bmmBs2JpsiPhiReader(TChain *tree, TString evtClassName);
  ~bmmBs2JpsiPhiReader();

  void    moreBasicCuts(); 
  void    bookHist();
  void    startAnalysis();
  void    eventProcessing();
  void    readCuts(TString filename, int dump = 1);
  void    initVariables();
  void    MCKinematics();
  void    candidateSelection(int mode = 0); 
  void    fillCandidateVariables();
  void    fillCandidateHistograms(); 
  int     tmCand(TAnaCand *pC);
  void    insertCand(TAnaCand* pCand);
  
  double       fKa1Pt, fKa1Eta, fKa1Phi;
  double       fKa2Pt, fKa2Eta, fKa2Phi;

  // -- Additional variables and cuts for Bs -> J/psi phi
  double            MKKLO, MKKHI, DELTAR;
  double            fDeltaR, fMKK;
  bool              fGoodDeltaR, fGoodMKK;

  AnalysisDistribution   *fpMKK, *fpDeltaR; 

};

#endif
