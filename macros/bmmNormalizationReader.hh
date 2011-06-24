#ifndef BMMNORMALIZATIONREADER_H
#define BMMNORMALIZATIONREADER_H

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

class bmmNormalizationReader : public bmmReader {

public:
  bmmNormalizationReader(TChain *tree, TString evtClassName);
  ~bmmNormalizationReader();

  void    moreBasicCuts();
  void    bookHist();
  void    startAnalysis();
  void    eventProcessing();
  void    fillCandidateHistograms(int offset = 0);
  void    readCuts(TString filename, int dump = 1);
  void    initVariables();
  void    MCKinematics();
  void    efficiencyCalculation();
  void    candidateSelection(int mode = 0); 
  void    fillCandidateVariables(); 
  int     tmCand(TAnaCand *pC);
  int     tmCand2(TAnaCand *pC);
  int     partialReco(TAnaCand *pCand);
  int     fromB(TGenCand *pCand);

  void    genMatch();
  void    recoMatch();
  void    candMatch();

  double       fKaonPt, fKaonEta, fKaonPhi;
  double       fKPtGen, fKEtaGen;
  double       fKaonPtNrf, fKaonEtaNrf;
  int          fKaonTkQuality;
  double       fJpsiPt, fJpsiEta, fJpsiPhi;

  // -- TM
  int          fGenK1Tmi; 
  int          fRecK1Tmi; 

  // -- effTree
  float fETk1pt, fETk1eta, fETg3pt, fETg3eta;
  int   fETk1q; 
  bool  fETk1gt;

  AnalysisDistribution *fpKaonPt[NAD], *fpKaonEta[NAD], *fpMpsi[NAD], *fpPsiPt[NAD], *fpPsiEta[NAD];

};

#endif
