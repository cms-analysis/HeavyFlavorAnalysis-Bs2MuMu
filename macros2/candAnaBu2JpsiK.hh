#ifndef CANDANABU2JPSIK_H
#define CANDANABU2JPSIK_H

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

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"

#include "candAna.hh"

class candAnaBu2JpsiK : public candAna {
  
public:
  candAnaBu2JpsiK(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAnaBu2JpsiK();

  void        candAnalysis();
  void        efficiencyCalculation();
  void        moreBasicCuts();

  void        genMatch(); 
  void        recoMatch(); 
  void        candMatch(); 
  
  void        readCuts(string filename, int dump);
  void        bookHist();
  void        fillCandidateHistograms(int offset);

  
  int          JPSITYPE; 
  double       JPSIMASSLO, JPSIMASSHI;

  double       fKaonPt, fKaonEta, fKaonPhi;
  double       fKPtGen, fKEtaGen;
  double       fKaonPtNrf, fKaonEtaNrf;
  int          fKaonTkQuality;
  double       fJpsiMass, fJpsiPt, fJpsiEta, fJpsiPhi;

  bool         fGoodJpsiMass;
  bool         fKa1Missid, fKa1MuMatch;

  // -- TM
  int          fGenK1Tmi; 
  int          fRecK1Tmi; 

  // -- effTree
  float fETk1pt, fETk1eta, fETg3pt, fETg3eta;
  int   fETk1q; 
  bool  fETk1gt;


  // -- AnalysisDistributionsa
  AnalysisDistribution *fpKaonPt[NAD], *fpKaonEta[NAD], *fpMpsi[NAD], *fpPsiPt[NAD], *fpPsiEta[NAD];

};

#endif
