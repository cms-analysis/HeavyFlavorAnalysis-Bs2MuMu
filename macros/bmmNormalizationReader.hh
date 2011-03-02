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
  void    fillCandidateHistograms();
  void    readCuts(TString filename, int dump = 1);
  void    initVariables();
  void    MCKinematics();
  void    efficiencyCalculation();
  void    candidateSelection(int mode = 0); 
  void    fillCandidateVariables(); 
  int     tmCand(TAnaCand *pC);
  int     tmCand2(TAnaCand *pC);

  void    genMatch();
  void    recoMatch();
  void    candMatch();

  double       fKaonPt, fKaonEta, fKaonPhi;
  double       fKPtGen, fKEtaGen;
  double       fKaonPtNrf, fKaonEtaNrf;
  int          fKaonTkQuality;

  // -- TM
  int          fGenK1Tmi; 
  int          fRecK1Tmi; 

};

#endif
