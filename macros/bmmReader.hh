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
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/JSON.hh"

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
  virtual void   efficiencyCalculation();
  virtual void   L1TSelection();  
  virtual void   HLTSelection();  
  virtual void   trackSelection();  
  virtual void   muonSelection();  
  virtual void   bookHist();
  virtual void   fillHist();
  virtual void   readCuts(TString filename, int dump = 1);
  virtual void   readFile(std::string filename, std::vector<std::string> &lines);
  virtual void   initVariables();
  virtual bool   muonID(TAnaTrack *pT);

  virtual void   genMatch();
  virtual void   recoMatch();
  virtual void   candMatch();
  virtual bool   evtFoundInCN(int evt);


  virtual void   studyL1T(); 

  virtual AnalysisDistribution* bookDistribution(const char *hn, const char *ht, const char *hc, int nbins, double lo, double hi);
  virtual void                  candidateSelection(int mode = 0);  // 0 = closest in r-phi
  virtual void                  fillCandidateVariables(); 
  virtual void                  fillCandidateHistograms();
  virtual void                  insertCand(TAnaCand*);
  virtual int                   tmCand(TAnaCand*);
  virtual int                   tmCand2(TAnaCand*);
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
    , JPSIMASSLO
    , JPSIMASSHI
    ;
  int TYPE, SELMODE, MUIDMASK, MUIDRESULT, TRACKQUALITY, JPSITYPE;
  std::vector<std::string> HLTPath, L1TPath; 
  std::string JSONFILE;

  bool fL1TMu0, fL1TMu3;
  bool fHLTMu0, fHLTMu3;

  bool fGoodEffCand;

  // -- vectors with the candidates of the specified type
  std::vector<TAnaCand *> fCands;  
  TAnaCand               *fpCand;       // the 'best' candidate
  int                     fCandIdx; 

  std::vector<bool>       fvGoodTracks, fvGoodTracksPt, fvGoodTracksEta;
  std::vector<bool>       fvGoodMuonsID, fvGoodMuonsPt, fvGoodMuonsEta;
  std::vector<bool>       fvGoodCand, fvGoodCandPt;

  bool                    fGoodEvent;

  // -- TM
  int                     fGenM1Tmi, fGenM2Tmi, fNGenPhotons; 
  int                     fRecM1Tmi, fRecM2Tmi; 
  int                     fCandTmi; 
 
  // -- variables for reduced tree, they are from fpCand
  int                     fJSON;
  int                     fCandTM, fCandType; 
  int                     fMu1TkQuality, fMu2TkQuality, fMu1Q, fMu2Q, fCandQ;
  bool                    fMu1Id, fMu2Id;
  double                  fHltMu1Pt, fHltMu1Eta, fHltMu1Phi, fHltMu2Pt, fHltMu2Eta, fHltMu2Phi;
  double                  fMu1Pt, fMu1Eta, fMu1Phi, fMu2Pt, fMu2Eta, fMu2Phi;
  double                  fMu1PtGen, fMu2PtGen, fMu1EtaGen, fMu2EtaGen;
  double                  fMu1PtNrf, fMu2PtNrf, fMu1EtaNrf, fMu2EtaNrf; // "now refitted"
  int                     fMu1Pix, fMu1BPix, fMu1BPixL1, fMu2Pix, fMu2BPix, fMu2BPixL1;
  double                  fMu1W8Mu, fMu1W8Tr, fMu2W8Mu, fMu2W8Tr; 
  double                  fPvX, fPvY, fPvZ; 
  double                  fJpsiMass;
  double                  fCandPt, fCandEta, fCandPhi, fCandM, fCandW8Tr, fCandW8Mu; 
  double                  fCandCosA, fCandIso, fCandIso1, fCandChi2, fCandDof, fCandProb, fCandFLS3d, fCandFLSxy; 
  double                  fCandDocaTrk, fMu1IP, fMu2IP, fCandPvTip, fCandPvTipE, fCandPvLip, fCandPvLipE; 

  double       MASSMIN,   MASSMAX; 
  double       SIGBOXMIN, SIGBOXMAX; 
  double       BGLBOXMIN, BGLBOXMAX; 
  double       BGHBOXMIN, BGHBOXMAX; 

  // -- Cut Variables
  bool                    fWideMass; 
  bool                    fGoodMCKinematics;
  bool                    fGoodL1T, fGoodHLT;
  bool                    fGoodTracks, fGoodTracksPt, fGoodTracksEta, fGoodMuonsID, fGoodMuonsPt, fGoodMuonsEta; 
  bool                    fGoodJpsiMass;
  bool                    fGoodQ, fGoodPt, fGoodEta, fGoodCosA, fGoodIso, fGoodChi2, fGoodFLS; 
  bool                    fGoodDocaTrk, fGoodIP; 

  bool                    fPreselection; 

  AnalysisCuts fAnaCuts; 

  // -- Analysis distributions
  AnalysisDistribution   *fpAllEvents, *fpHLT, *fpPvZ,  
    *fpTracksQual, *fpTracksPt,  *fpTracksEta, 
    *fpMuonsID, *fpMuonsPt, *fpMuonsEta, *fpMuon1Pt, *fpMuon2Pt, *fpMuon1Eta, *fpMuon2Eta,
    *fpMpsi,
    *fpQ, *fpPt, *fpEta, 
    *fpCosA, *fpCosA0, *fpAlpha,
    *fpIso, *fpIso1, 
    *fpDoca, *fpIP,
    *fpChi2, *fpChi2Dof, *fpProb, 
    *fpFLS3d, *fpFLSxy, 
    *fpDocaTrk, *fpIP1, *fpIP2;   

  // -- another reduced tree
  TTree       *fEffTree;
  bool fETm1gt, fETm2gt, fETm1id, fETm2id;
  int fETm1q, fETm2q; 
  float fETcandMass;
  float fETm1pt, fETm1eta, fETg1pt, fETg1eta;
  float fETm2pt, fETm2eta, fETg2pt, fETg2eta;

  // -- PidTables
  PidTable *fpMuonID;
  PidTable *fpMuonTr;

  JSON *fpJSON; 

  std::vector<int> fEventVector;

};

#endif
