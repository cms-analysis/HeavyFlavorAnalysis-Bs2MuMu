#ifndef CANDANABS2JPSIPHI_H
#define CANDANABS2JPSIPHI_H

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


#include "candAna.hh"

class candAnaBs2JpsiPhi : public candAna {
  
public:
  candAnaBs2JpsiPhi(bmm2Reader *pReader, std::string name, std::string cutsFile);
  virtual ~candAnaBs2JpsiPhi();

  void        candAnalysis();
  void        efficiencyCalculation();
  void        moreBasicCuts();
  void        moreReducedTree(TTree *);

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


  double       fKa1Pt, fKa1Eta, fKa1Phi;
  double       fKa2Pt, fKa2Eta, fKa2Phi;
  double       fPhiPt, fPhiEta, fPhiPhi; 
  bool         fKa1Missid, fKa2Missid, fKa1MuMatch, fKa2MuMatch;
  // tests d.k.   
  bool         fKa1Missid2, fKa1MuMatch2; 
  float        fKa1MuMatchR, fKa1MuMatchR2, fKa1MuMatchR3, fKa1MuMatchR4, fKa1MuMatchR5;
  bool         fKa2Missid2, fKa2MuMatch2; 
  float        fKa2MuMatchR, fKa2MuMatchR2, fKa2MuMatchR3, fKa2MuMatchR4, fKa2MuMatchR5;


  double       fKa1PtNrf, fKa1EtaNrf;
  double       fKa2PtNrf, fKa2EtaNrf;

  double       fKa1PtGen, fKa1EtaGen, fKa2PtGen, fKa2EtaGen;
  int          fKa1GenID, fKa2GenID; 
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
  int               PHITYPE; 
  double            MKKLO, MKKHI, DELTAR;
  double            fDeltaR, fMKK;
  bool              fGoodDeltaR, fGoodMKK;


};

#endif
