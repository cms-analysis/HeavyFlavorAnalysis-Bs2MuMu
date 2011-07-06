#ifndef BMMREADER_H
#define BMMREADER_H

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

  virtual int    partialReco(TAnaCand *pCand);
  virtual void   processType();
  virtual void   genMatch();
  virtual void   recoMatch();
  virtual void   candMatch();
  virtual bool   evtFoundInCN(int evt);


  virtual void   studyL1T(); 

  virtual AnalysisDistribution* bookDistribution(const char *hn, const char *ht, const char *hc, int nbins, double lo, double hi);
  virtual void                  candidateSelection(int mode = 0);  // 0 = closest in r-phi
  virtual void                  fillCandidateVariables(); 
  virtual void                  fillCandidateHistograms(int offset);
  virtual void                  insertCand(TAnaCand*);
  virtual int                   tmCand(TAnaCand*);
  virtual int                   tmCand2(TAnaCand*);
  virtual double                isoClassic(TAnaCand*); 
  virtual double                isoClassicOnePv(TAnaCand*, double r = 1.0, double ptmin = 0.9); 
  virtual double                isoClassicWithDOCA(TAnaCand*, float dca, double r = 1.0, double ptmin = 0.9); 
  virtual double                isoWithDOCA(TAnaCand*, float dca, double r = 1.0, double ptmin = 0.9); 
  virtual int                   checkCut(const char *, TH1D *); 
  virtual std::string           splitTrigRange(string hlt, int &r1, int &r2);


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
  int TYPE, SELMODE, MUIDMASK, MUIDRESULT, TRACKQUALITY, JPSITYPE, TRUTHCAND, IGNORETRIGGER;
  std::vector<std::string> HLTPath, L1TPath; 
  std::map<std::string, int> HLTRangeMin;
  std::map<std::string, int> HLTRangeMax;

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
  int                     fGenBTmi; 
  int                     fGenM1Tmi, fGenM2Tmi, fNGenPhotons; 
  int                     fRecM1Tmi, fRecM2Tmi; 
  int                     fCandTmi; 
  int                     fGenBpartial; 
  int                     fProcessType;
 
  // -- variables for reduced tree, they are from fpCand
  bool                    fJSON;
  int                     fCandTM, fCandType; 
  int                     fMu1TkQuality, fMu2TkQuality, fMu1Q, fMu2Q, fMu1Chi2, fMu2Chi2, fCandQ;
  bool                    fMu1Id, fMu2Id;
  double                  fHltMu1Pt, fHltMu1Eta, fHltMu1Phi, fHltMu2Pt, fHltMu2Eta, fHltMu2Phi;
  double                  fMu1Pt, fMu1Eta, fMu1Phi, fMu2Pt, fMu2Eta, fMu2Phi;
  double                  fMu1PtGen, fMu2PtGen, fMu1EtaGen, fMu2EtaGen;
  double                  fMu1PtNrf, fMu2PtNrf, fMu1EtaNrf, fMu2EtaNrf; // "now refitted"
  int                     fMu1Pix, fMu1BPix, fMu1BPixL1, fMu2Pix, fMu2BPix, fMu2BPixL1;
  double                  fMu1W8Mu, fMu1W8Tr, fMu2W8Mu, fMu2W8Tr; 
  double                  fPvX, fPvY, fPvZ, fPvNtrk; 
  int                     fPvN;
  double                  fJpsiMass;
  double                  fCandPt, fCandEta, fCandPhi, fCandM, fCandM2, fCandW8Tr, fCandW8Mu; 
  double                  fCandCosA, fCandA, fCandChi2, fCandDof, fCandProb, fCandFL3d, fCandFL3dE, fCandFLS3d, fCandFLSxy; 
  double                  fCandIso, fCandIso1, fCandIso2, fCandIso3, fCandIso4, fCandIso5;
  double                  fIsoR05Pt03, fIsoR05Pt05, fIsoR05Pt07, fIsoR05Pt09, fIsoR05Pt11;
  double                  fIsoR07Pt03, fIsoR07Pt05, fIsoR07Pt07, fIsoR07Pt09, fIsoR07Pt11;
  double                  fIsoR10Pt03, fIsoR10Pt05, fIsoR10Pt07, fIsoR10Pt09, fIsoR10Pt11;
  int                     fCandItrk, fCandI0trk, fCandI4trk; 
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

  bool fBarrel; 
  AnalysisCuts fAnaCuts; 

  // -- Analysis distributions
#define NAD 18

  std::map<string, int> fRegion;

  AnalysisDistribution   *fpAllEvents[NAD], *fpHLT[NAD], *fpPvZ[NAD], *fpPvN[NAD], *fpPvNtrk[NAD],  
    *fpTracksQual[NAD], *fpTracksPt[NAD],  *fpTracksEta[NAD], 
    *fpMuonsID[NAD], *fpMuonsPt[NAD], *fpMuonsEta[NAD], *fpMuon1Pt[NAD], *fpMuon2Pt[NAD], *fpMuon1Eta[NAD], *fpMuon2Eta[NAD],
    *fpQ[NAD], *fpPt[NAD], *fpEta[NAD], 
    *fpCosA[NAD], *fpCosA0[NAD], *fpAlpha[NAD],
    *fpIso[NAD], *fpIso1[NAD], *fpIso2[NAD], *fpIso3[NAD], *fpIso4[NAD], *fpIso5[NAD], *fpI0trk[NAD], *fpI4trk[NAD], 
    *fpIsoPv1[NAD], *fpIsoPv2[NAD], *fpIsoPv3[NAD], *fpIsoPv4[NAD], *fpIsoPv5[NAD], *fpIsoPv6[NAD],  
    *fpIso1Pv1[NAD], *fpIso1Pv2[NAD], *fpIso1Pv3[NAD], *fpIso1Pv4[NAD], *fpIso1Pv5[NAD], *fpIso1Pv6[NAD],  
    *fpIso4Pv1[NAD], *fpIso4Pv2[NAD], *fpIso4Pv3[NAD], *fpIso4Pv4[NAD], *fpIso4Pv5[NAD], *fpIso4Pv6[NAD],  
    *fpIso5Pv1[NAD], *fpIso5Pv2[NAD], *fpIso5Pv3[NAD], *fpIso5Pv4[NAD], *fpIso5Pv5[NAD], *fpIso5Pv6[NAD],  
    *fpFLS3dPv1[NAD], *fpFLS3dPv2[NAD], *fpFLS3dPv3[NAD], *fpFLS3dPv4[NAD], *fpFLS3dPv5[NAD], *fpFLS3dPv6[NAD],  
    *fpFLSxyPv1[NAD], *fpFLSxyPv2[NAD], *fpFLSxyPv3[NAD], *fpFLSxyPv4[NAD], *fpFLSxyPv5[NAD], *fpFLSxyPv6[NAD],  
    *fpAlphaPv1[NAD], *fpAlphaPv2[NAD], *fpAlphaPv3[NAD], *fpAlphaPv4[NAD], *fpAlphaPv5[NAD], *fpAlphaPv6[NAD],  

    *fpIsoR05Pt03[NAD], *fpIsoR05Pt05[NAD], *fpIsoR05Pt07[NAD], *fpIsoR05Pt09[NAD], *fpIsoR05Pt11[NAD],
    *fpIsoR07Pt03[NAD], *fpIsoR07Pt05[NAD], *fpIsoR07Pt07[NAD], *fpIsoR07Pt09[NAD], *fpIsoR07Pt11[NAD],
    *fpIsoR10Pt03[NAD], *fpIsoR10Pt05[NAD], *fpIsoR10Pt07[NAD], *fpIsoR10Pt09[NAD], *fpIsoR10Pt11[NAD],

    *fpDoca[NAD], *fpIP[NAD],
    *fpChi2[NAD], *fpChi2Dof[NAD], *fpProb[NAD], 
    *fpFLS3d[NAD], *fpFLSxy[NAD], 
    *fpFL3d[NAD], *fpFL3dE[NAD], 
    *fpDocaTrk[NAD], *fpIP1[NAD], *fpIP2[NAD];   

  // -- another reduced tree
  TTree       *fEffTree;
  bool fETm1gt, fETm2gt, fETm1id, fETm2id;
  int fETm1q, fETm2q; 
  float fETgpt, fETgeta; 
  float fETcandMass;
  float fETm1pt, fETm1eta, fETg1pt, fETg1eta;
  float fETm2pt, fETm2eta, fETg2pt, fETg2eta;

  // -- PidTables
  PidTable *fpMuonID;
  PidTable *fpMuonTr, *fpMuonTr1, *fpMuonTr2;

  std::vector<int> fEventVector;

};

#endif
