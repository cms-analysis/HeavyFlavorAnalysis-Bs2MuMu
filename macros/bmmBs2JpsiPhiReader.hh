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
  void    efficiencyCalculation();
  void    candidateSelection(int mode = 0); 
  void    fillCandidateVariables();
  void    fillCandidateHistograms(); 
  int     tmCand(TAnaCand *pC);
  int     tmCand2(TAnaCand *pC);
  void    insertCand(TAnaCand* pCand);

  void    genMatch();
  void    recoMatch();
  void    candMatch();
  
  double       fKa1Pt, fKa1Eta, fKa1Phi;
  double       fKa2Pt, fKa2Eta, fKa2Phi;

  double       fKa1PtNrf, fKa1EtaNrf;
  double       fKa2PtNrf, fKa2EtaNrf;

  double       fKa1PtGen, fKa1EtaGen, fKa2PtGen, fKa2EtaGen;
  int          fKa1TkQuality, fKa2TkQuality;

  // -- TM 
  int                     fGenK1Tmi, fGenK2Tmi; 
  int                     fRecK1Tmi, fRecK2Tmi; 
  
  // -- effTree
  float fETk1pt, fETk1eta, fETg3pt, fETg3eta;
  float fETk2pt, fETk2eta, fETg4pt, fETg4eta;
  int   fETk1q,  fETk2q; 
  bool  fETk1gt, fETk2gt;

  // -- Additional variables and cuts for Bs -> J/psi phi
  double            MKKLO, MKKHI, DELTAR;
  double            fDeltaR, fMKK;
  bool              fGoodDeltaR, fGoodMKK;

  AnalysisDistribution   *fpMKK, *fpDeltaR; 

};

#endif
