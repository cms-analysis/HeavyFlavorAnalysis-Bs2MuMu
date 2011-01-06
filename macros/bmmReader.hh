#ifndef BMMREADER_H
#define BMMREADER_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>

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

#include "treeReader01.hh"
#include "AnalysisCuts.hh"
#include "AnalysisDistribution.hh"

#define DR      57.29577951

class bmmReader : public treeReader01 {

public:
  bmmReader(TChain *tree, TString evtClassName);
  ~bmmReader();

  virtual void   startAnalysis();
  virtual void   basicCuts();
  virtual void   moreBasicCuts();
  virtual void   candidateCuts();
  virtual void   moreCandidateCuts();

  virtual void   eventProcessing();
  virtual void   MCKinematics();  
  virtual void   L1TSelection();  
  virtual void   HLTSelection();  
  virtual void   trackSelection();  
  virtual void   muonSelection();  
  virtual void   bookHist();
  virtual void   fillHist();
  virtual void   readCuts(TString filename, int dump = 1);
  virtual void   readFile(std::string filename, std::vector<std::string> &lines);
  virtual void   initVariables();

  virtual void   studyL1T(); 

  virtual AnalysisDistribution* bookDistribution(const char *hn, const char *ht, const char *hc, int nbins, double lo, double hi);
  virtual void                  candidateSelection(int mode = 0);  // 0 = closest in r-phi
  virtual void                  fillCandidateVariables(); 
  virtual void                  fillCandidateHistograms();
  virtual void                  insertCand(TAnaCand*);
  virtual int                   tmCand(TAnaCand*);
  virtual double                isoClassic(TAnaCand*); 
  virtual double                isoClassicOnePv(TAnaCand*); 
  virtual int                   checkCut(const char *, TH1D *); 


  // -- Cut values
  double 
      CANDPTLO
    , CANDETALO
    , CANDETAHI   
    , CANDCOSALPHA
    , CANDFLS3D
    , CANDFLSXY
    , CANDVTXCHI2
    , CANDISOLATION
    , CANDDOCATRK
    , TRACKPTLO
    , TRACKPTHI
    , TRACKETALO
    , TRACKETAHI
    , TRACKTIP
    , TRACKLIP
    , MUPTLO
    , MUPTHI
    , MUETALO
    , MUETAHI   
    , MUIP
    ;
  int TYPE, SELMODE, MUIDMASK, MUIDRESULT, TRACKQUALITY;
  std::vector<std::string> HLTPath, L1TPath; 

  bool fL1TMu0, fL1TMu3;
  bool fHLTMu0, fHLTMu3;
  

  // -- vectors with the candidates of the specified type
  std::vector<TAnaCand *> fCands;  
  TAnaCand               *fpCand;       // the 'best' candidate
  int                     fCandIdx; 

  std::vector<bool>       fvGoodTracks, fvGoodTracksPt, fvGoodTracksEta;
  std::vector<bool>       fvGoodMuonsID, fvGoodMuonsPt, fvGoodMuonsEta;
  std::vector<bool>       fvGoodCand, fvGoodCandPt;

  bool                    fGoodEvent;

  // -- variables for reduced tree, they are from fpCand
  int                     fCandTM; 
  double                  fMu1Pt, fMu1Eta, fMu1Phi, fMu2Pt, fMu2Eta, fMu2Phi;
  double                  fMu1W8Mu, fMu1W8Tr, fMu2W8Mu, fMu2W8Tr; 
  double                  fPvX, fPvY, fPvZ; 
  double                  fCandPt, fCandEta, fCandPhi, fCandM, fCandW8Tr, fCandW8Mu; 
  double                  fCandCosA, fCandIso, fCandIso1, fCandChi2, fCandDof, fCandProb, fCandFLS3d, fCandFLSxy; 
  double                  fCandDocaTrk, fMu1IP, fMu2IP; 

  double       SIGBOXMIN, SIGBOXMAX; 
  double       BGLBOXMIN, BGLBOXMAX; 
  double       BGHBOXMIN, BGHBOXMAX; 

  // -- Cut Variables
  bool                    fWideMass; 
  bool                    fGoodMCKinematics;
  bool                    fGoodL1T, fGoodHLT;
  bool                    fGoodTracks, fGoodTracksPt, fGoodTracksEta, fGoodMuonsID, fGoodMuonsPt, fGoodMuonsEta; 
  bool                    fGoodPt, fGoodEta, fGoodCosA, fGoodIso, fGoodChi2, fGoodFLS; 
  bool                    fGoodDocaTrk, fGoodIP; 

  AnalysisCuts fAnaCuts; 

  // -- Analysis distributions
  AnalysisDistribution   *fpAllEvents, *fpHLT, *fpPvZ, *fpTracksPt, 
    *fpMuonsID, *fpMuonsPt, *fpMuonsEta, 
    *fpPt, *fpEta, 
    *fpCosA, *fpCosA0, 
    *fpIso, *fpIso1, 
    *fpDoca, *fpIP,
    *fpChi2, *fpChi2Dof, *fpProb, 
    *fpFLS3d, *fpFLSxy, 
    *fpDocaTrk, *fpIP1, *fpIP2;   

  // -- PidTables
  PidTable *fpMuonID;
  PidTable *fpMuonTr;

};

#endif
