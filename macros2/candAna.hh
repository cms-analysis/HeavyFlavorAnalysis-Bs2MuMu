#ifndef CANDANA_H
#define CANDANA_H

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
#include <TDirectory.h>
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

#include "../macros/AnalysisCuts.hh"
#include "../macros/AnalysisDistribution.hh"

#include "bmm2Reader.hh"


class ReadBDT; 
class ReadBDT2; 

struct isoNumbers {
  double iso; 
  int    pvTracks; 
  int    clTracks; 
  int    Tracks; 
};

class candAna {
  
public:
  candAna(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAna();

  virtual void        evtAnalysis(TAna01Event *evt);
  virtual void        candAnalysis();
  virtual void        efficiencyCalculation();
  
  virtual int         nearestPV(int pvIdx, double maxDist = 99.);
  virtual double      constrainedMass();
  virtual void        runRange();
  virtual void        genMatch(); 
  virtual void        recoMatch(); 
  virtual void        candMatch(); 
  virtual void        triggerSelection();
  virtual void        fillCandidateHistograms(int offset);
    
  virtual void        bookHist();
  virtual AnalysisDistribution* bookDistribution(const char *hn, const char *ht, const char *hc, int nbins, double lo, double hi);

  virtual void        basicCuts();
  virtual void        moreBasicCuts();
  virtual void        candidateCuts();
  virtual void        moreCandidateCuts();

  virtual void        readCuts(std::string fileName, int dump = 1);
  virtual void        readFile(std::string fileName, std::vector<std::string> &lines);

  virtual bool        goodTrack(TAnaTrack *pt);
  virtual bool        goodMuon(TAnaTrack *pt, int mask = 0);
  virtual bool        tightMuon(TAnaTrack *pt);

  virtual std::string splitTrigRange(std::string tl, int &r1, int &r2);

  virtual double      isoClassicWithDOCA(TAnaCand*, double dca, double r = 0.7, double ptmin = 0.9); 
  virtual int         nCloseTracks(TAnaCand*, double dca, double pt = 0.5); 

  virtual TAnaCand*   osCand(TAnaCand *pC);
  virtual double      osIsolation(TAnaCand *pC, double r = 1.0, double ptmin = 0.9); 
  virtual int         osMuon(TAnaCand *pC, double r = 1.0); 

  virtual void        isolationStudy(double doca);
  virtual void        bookIsoPlots();  
  virtual void        fillIsoPlots();

  std::string fName; 
  std::string fCutFile; 
  TDirectory *fHistDir; 
  bmm2Reader *fpReader; 
  TTree *fTree; 
  TAna01Event *fpEvt;
  TAnaCand *fpCand, *fpOsCand;
  int fCandIdx; 

  int fVerbose;
  int fIsMC;

  int fRun, fEvt, fLS;
  int fEvent; 
  int fRunRange;

  double       MASSMIN,   MASSMAX; 
  double       SIGBOXMIN, SIGBOXMAX; 
  double       BGLBOXMIN, BGLBOXMAX; 
  double       BGHBOXMIN, BGHBOXMAX; 
  
  double 
  CANDPTLO, CANDETALO, CANDETAHI
    , CANDCOSALPHA, CANDALPHA
    , CANDFLS3D, CANDFLSXY, CANDVTXCHI2
    , CANDISOLATION, CANDDOCATRK, CANDCLOSETRK
    , PVAVEW8, CANDLIP, CANDLIPS, CANDLIP2, CANDLIPS2
    , CANDDOCA, CANDIP, CANDIPS
    , TRACKPTLO, TRACKPTHI, TRACKETALO, TRACKETAHI
    , TRACKTIP, TRACKLIP
    , MUPTLO, MUPTHI
    , MUETALO, MUETAHI, MUIP
    ;
  
  int BLIND, TYPE, SELMODE, MUIDMASK, MUIDRESULT, TRACKQUALITY, TRUTHCAND, IGNORETRIGGER;

  std::map<std::string, pair<int, int> > HLTRANGE;


  bool fBarrel, fWideMass; 
  AnalysisCuts fAnaCuts; 

  // -- TM
  int     fGenBTmi; 
  int     fGenM1Tmi, fGenM2Tmi, fNGenPhotons; 
  int     fRecM1Tmi, fRecM2Tmi; 
  int     fCandTmi; 
  int     fGenBpartial; 
  int     fProcessType;
  double  fGenLifeTime; 

  // -- variables for reduced tree, they are from fpCand
  bool    fJSON, fCowboy;
  int     fCandTM, fCandType; 
  double  fCandBDT, fCandBDT2; 
  int     fMu1TkQuality, fMu2TkQuality, fMu1Q, fMu2Q, fMu1Chi2, fMu2Chi2, fCandQ, fMu1PV, fMu2PV;
  bool    fMu1Id, fMu2Id;  
  double  fMuDist, fMuDeltaR;
  double  fHltMu1Pt, fHltMu1Eta, fHltMu1Phi, fHltMu2Pt, fHltMu2Eta, fHltMu2Phi;
  double  fMu1Pt, fMu1Eta, fMu1Phi, fMu2Pt, fMu2Eta, fMu2Phi;
  double  fMu1PtGen, fMu2PtGen, fMu1EtaGen, fMu2EtaGen;
  double  fMu1PtNrf, fMu2PtNrf, fMu1EtaNrf, fMu2EtaNrf; // "now refitted"
  int     fMu1Pix, fMu1BPix, fMu1BPixL1, fMu2Pix, fMu2BPix, fMu2BPixL1;
  double  fMu1W8Mu, fMu1W8Tr, fMu2W8Mu, fMu2W8Tr; 
  double  fPvX, fPvY, fPvZ, fPvNtrk, fPvNdof, fPvAveW8; 
  int     fPvN;
  double  fCandPt, fCandP, fCandTau, fCandEta, fCandPhi, fCandM, fCandM2, fCandW8Tr, fCandW8Mu; 
  double  fCandCosA, fCandA;
  double  fCandChi2, fCandDof, fCandProb, fCandFL3d, fCandFL3dE, fCandFLS3d, fCandFLxy, fCandFLSxy, fCandDoca; 
  double  f2MChi2,   f2MDof,   f2MProb,   f2MFL3d,   f2MFL3dE,   f2MFLS3d,   f2MFLSxy; 
  double  fCandIso;
  int     fCandIsoTrk, fCandCloseTrk, fCandPvTrk, fCandI0trk, fCandI1trk, fCandI2trk; 
  double  fCandDocaTrk, fCandDocaTrkBdt, fMu1IP, fMu1IPE, fMu2IP, fMu2IPE; 
  double  fCandPvTip, fCandPvTipE, fCandPvTipS, fCandPvLip, fCandPvLipE, fCandPvLipS, fCandPvIp, fCandPvIpE, fCandPvIpS;
  double  fCandPvLip2, fCandPvLipS2, fCandPvLip12, fCandPvLipE12, fCandPvLipS12; 

  double  fOsMuonPt, fOsMuonPtRel, fOsIso, fOsRelIso, fOsMuonDeltaR;

  // -- another reduced tree
  TTree       *fEffTree;
  bool fETm1gt, fETm2gt, fETm1id, fETm2id;
  int fETm1q, fETm2q; 
  float fETgpt, fETgeta; 
  float fETcandMass;
  float fETm1pt, fETm1eta, fETg1pt, fETg1eta;
  float fETm2pt, fETm2eta, fETg2pt, fETg2eta;

  bool    fGoodEffCand; 

  TAnaTrack *fpMuon1, *fpMuon2; 

  // -- isolation study
  isoNumbers fIsoR03Pt03, fIsoR03Pt05, fIsoR03Pt07, fIsoR03Pt09, fIsoR03Pt11;
  isoNumbers fIsoR05Pt03, fIsoR05Pt05, fIsoR05Pt07, fIsoR05Pt09, fIsoR05Pt11;
  isoNumbers fIsoR07Pt03, fIsoR07Pt05, fIsoR07Pt07, fIsoR07Pt09, fIsoR07Pt11;
  isoNumbers fIsoR09Pt03, fIsoR09Pt05, fIsoR09Pt07, fIsoR09Pt09, fIsoR09Pt11;
  isoNumbers fIsoR10Pt03, fIsoR10Pt05, fIsoR10Pt07, fIsoR10Pt09, fIsoR10Pt11;
  isoNumbers fIsoR11Pt03, fIsoR11Pt05, fIsoR11Pt07, fIsoR11Pt09, fIsoR11Pt11;

  string  fHLTPath;
  bool    fGoodHLT, fGoodMuonsID, fGoodMuonsPt, fGoodMuonsEta, fGoodTracks, fGoodTracksPt, fGoodTracksEta;
  bool    fGoodPvAveW8, fGoodPvLip, fGoodPvLipS, fGoodPvLip2, fGoodPvLipS2, fGoodMaxDoca, fGoodIp, fGoodIpS; 
  bool    fGoodQ, fGoodPt, fGoodEta, fGoodCosA, fGoodAlpha, fGoodIso, fGoodCloseTrack, fGoodChi2, fGoodFLS; 
  bool    fGoodDocaTrk, fGoodLastCut; 

  bool    fPreselection; 

  // -- Analysis distributions
  std::map<std::string, int> fRegion;
#define NAD 13
  AnalysisDistribution   *fpHLT[NAD], *fpPvZ[NAD], *fpPvN[NAD], *fpPvNtrk[NAD], *fpPvAveW8[NAD]  
    , *fpTracksQual[NAD], *fpTracksPt[NAD],  *fpTracksEta[NAD] 
    , *fpMuonsID[NAD], *fpMuonsPt[NAD], *fpMuonsEta[NAD], *fpMuon1Pt[NAD], *fpMuon2Pt[NAD], *fpMuon1Eta[NAD], *fpMuon2Eta[NAD]
    , *fpPt[NAD], *fpP[NAD], *fpEta[NAD] 
    , *fpCosA[NAD], *fpAlpha[NAD]
    , *fpIso[NAD], *fpIsoTrk[NAD], *fpCloseTrk[NAD]
    , *fpChi2[NAD], *fpChi2Dof[NAD], *fpProb[NAD] 
    , *fpFLS3d[NAD], *fpFLSxy[NAD] 
    , *fpFL3d[NAD], *fpFL3dE[NAD] 
    , *fpDocaTrk[NAD]   
    , *fpBDT[NAD]   
    , *fpLip[NAD], *fpLipE[NAD], *fpLipS[NAD] 
    , *fpTip[NAD], *fpTipE[NAD], *fpTipS[NAD] 
    , *fpLip12[NAD], *fpLipE12[NAD], *fpLipS12[NAD] 
    , *fpLip2[NAD], *fpLipS2[NAD]
    , *fpMaxDoca[NAD], *fpIp[NAD], *fpIpS[NAD]
    , *fp2MChi2[NAD],  *fp2MChi2Dof[NAD], *fp2MProb[NAD] 
    , *fp2MFLS3d[NAD], *fp2MFLSxy[NAD] 
    , *fp2MFL3d[NAD],  *fp2MFL3dE[NAD] 
    , *fpOsIso[NAD],  *fpOsRelIso[NAD] 
    , *fpOsMuonPt[NAD],  *fpOsMuonDeltaR[NAD], *fpOsMuonPtRel[NAD]

    , *fpOsIsoGGF[NAD], *fpOsIsoGSP[NAD], *fpOsIsoFEX[NAD]  
    , *fpOsRelIsoGGF[NAD], *fpOsRelIsoGSP[NAD], *fpOsRelIsoFEX[NAD] 
    , *fpOsMuonPtGGF[NAD], *fpOsMuonPtGSP[NAD], *fpOsMuonPtFEX[NAD]
    , *fpOsMuonPtRelGGF[NAD], *fpOsMuonPtRelGSP[NAD], *fpOsMuonPtRelFEX[NAD]
    , *fpOsMuonDeltaRGGF[NAD], *fpOsMuonDeltaRGSP[NAD], *fpOsMuonDeltaRFEX[NAD] 
    , *fpIsoGGF[NAD], *fpIsoGSP[NAD], *fpIsoFEX[NAD]
    , *fpCloseTrkGGF[NAD], *fpCloseTrkGSP[NAD], *fpCloseTrkFEX[NAD]
    , *fpDocaTrkGGF[NAD], *fpDocaTrkGSP[NAD], *fpDocaTrkFEX[NAD]   
    , *fpPtGGF[NAD], *fpPtGSP[NAD], *fpPtFEX[NAD]   
    ;
  
  // -- Analysis distributions in bins of n(PV)
#define NADPV 15
  AnalysisDistribution   *fpNpvPvN[NADPV][NAD];
  AnalysisDistribution   *fpNpvAveW8[NADPV][NAD];
  AnalysisDistribution   *fpNpvChi2Dof[NADPV][NAD];
  AnalysisDistribution   *fpNpvProb[NADPV][NAD];
  AnalysisDistribution   *fpNpvFLS3d[NADPV][NAD];
  AnalysisDistribution   *fpNpvFLSxy[NADPV][NAD];
  AnalysisDistribution   *fpNpvAlpha[NADPV][NAD];
  AnalysisDistribution   *fpNpvDocaTrk[NADPV][NAD];
  AnalysisDistribution   *fpNpvIso[NADPV][NAD];
  AnalysisDistribution   *fpNpvIsoTrk[NADPV][NAD];
  AnalysisDistribution   *fpNpvCloseTrk[NADPV][NAD];
  AnalysisDistribution   *fpNpvLip[NADPV][NAD];
  AnalysisDistribution   *fpNpvLipS[NADPV][NAD];
  AnalysisDistribution   *fpNpvLip2[NADPV][NAD];
  AnalysisDistribution   *fpNpvLipS2[NADPV][NAD];
  AnalysisDistribution   *fpNpvMaxDoca[NADPV][NAD];
  AnalysisDistribution   *fpNpvIp[NADPV][NAD];
  AnalysisDistribution   *fpNpvIpS[NADPV][NAD];
  AnalysisDistribution   *fpNpvIso0[NADPV][NAD];
  AnalysisDistribution   *fpNpvIso1[NADPV][NAD];
  AnalysisDistribution   *fpNpvIso2[NADPV][NAD];
  AnalysisDistribution   *fpNpvIso3[NADPV][NAD];
  AnalysisDistribution   *fpNpvIso4[NADPV][NAD];
  AnalysisDistribution   *fpNpvIso5[NADPV][NAD];

  AnalysisDistribution   *fpEtaFLS3d[NADPV][NAD];

  // -- Isolation study
#define NISO 10
  AnalysisDistribution   
    *fpIsoR03Pt03[NISO], *fpIsoR03Pt05[NISO], *fpIsoR03Pt07[NISO], *fpIsoR03Pt09[NISO], *fpIsoR03Pt11[NISO],      
    *fpIsoR05Pt03[NISO], *fpIsoR05Pt05[NISO], *fpIsoR05Pt07[NISO], *fpIsoR05Pt09[NISO], *fpIsoR05Pt11[NISO],      
    *fpIsoR07Pt03[NISO], *fpIsoR07Pt05[NISO], *fpIsoR07Pt07[NISO], *fpIsoR07Pt09[NISO], *fpIsoR07Pt11[NISO],      
    *fpIsoR09Pt03[NISO], *fpIsoR09Pt05[NISO], *fpIsoR09Pt07[NISO], *fpIsoR09Pt09[NISO], *fpIsoR09Pt11[NISO],    
    *fpIsoR10Pt03[NISO], *fpIsoR10Pt05[NISO], *fpIsoR10Pt07[NISO], *fpIsoR10Pt09[NISO], *fpIsoR10Pt11[NISO],  
    *fpIsoR11Pt03[NISO], *fpIsoR11Pt05[NISO], *fpIsoR11Pt07[NISO], *fpIsoR11Pt09[NISO], *fpIsoR11Pt11[NISO];

  AnalysisDistribution   
    *fpTk0R03Pt03[NISO], *fpTk0R03Pt05[NISO], *fpTk0R03Pt07[NISO], *fpTk0R03Pt09[NISO], *fpTk0R03Pt11[NISO],      
    *fpTk0R05Pt03[NISO], *fpTk0R05Pt05[NISO], *fpTk0R05Pt07[NISO], *fpTk0R05Pt09[NISO], *fpTk0R05Pt11[NISO],      
    *fpTk0R07Pt03[NISO], *fpTk0R07Pt05[NISO], *fpTk0R07Pt07[NISO], *fpTk0R07Pt09[NISO], *fpTk0R07Pt11[NISO],      
    *fpTk0R09Pt03[NISO], *fpTk0R09Pt05[NISO], *fpTk0R09Pt07[NISO], *fpTk0R09Pt09[NISO], *fpTk0R09Pt11[NISO],    
    *fpTk0R10Pt03[NISO], *fpTk0R10Pt05[NISO], *fpTk0R10Pt07[NISO], *fpTk0R10Pt09[NISO], *fpTk0R10Pt11[NISO],  
    *fpTk0R11Pt03[NISO], *fpTk0R11Pt05[NISO], *fpTk0R11Pt07[NISO], *fpTk0R11Pt09[NISO], *fpTk0R11Pt11[NISO];

  AnalysisDistribution   
    *fpTk1R03Pt03[NISO], *fpTk1R03Pt05[NISO], *fpTk1R03Pt07[NISO], *fpTk1R03Pt09[NISO], *fpTk1R03Pt11[NISO],      
    *fpTk1R05Pt03[NISO], *fpTk1R05Pt05[NISO], *fpTk1R05Pt07[NISO], *fpTk1R05Pt09[NISO], *fpTk1R05Pt11[NISO],      
    *fpTk1R07Pt03[NISO], *fpTk1R07Pt05[NISO], *fpTk1R07Pt07[NISO], *fpTk1R07Pt09[NISO], *fpTk1R07Pt11[NISO],      
    *fpTk1R09Pt03[NISO], *fpTk1R09Pt05[NISO], *fpTk1R09Pt07[NISO], *fpTk1R09Pt09[NISO], *fpTk1R09Pt11[NISO],    
    *fpTk1R10Pt03[NISO], *fpTk1R10Pt05[NISO], *fpTk1R10Pt07[NISO], *fpTk1R10Pt09[NISO], *fpTk1R10Pt11[NISO],  
    *fpTk1R11Pt03[NISO], *fpTk1R11Pt05[NISO], *fpTk1R11Pt07[NISO], *fpTk1R11Pt09[NISO], *fpTk1R11Pt11[NISO];

  AnalysisDistribution   
    *fpTk2R03Pt03[NISO], *fpTk2R03Pt05[NISO], *fpTk2R03Pt07[NISO], *fpTk2R03Pt09[NISO], *fpTk2R03Pt11[NISO],      
    *fpTk2R05Pt03[NISO], *fpTk2R05Pt05[NISO], *fpTk2R05Pt07[NISO], *fpTk2R05Pt09[NISO], *fpTk2R05Pt11[NISO],      
    *fpTk2R07Pt03[NISO], *fpTk2R07Pt05[NISO], *fpTk2R07Pt07[NISO], *fpTk2R07Pt09[NISO], *fpTk2R07Pt11[NISO],      
    *fpTk2R09Pt03[NISO], *fpTk2R09Pt05[NISO], *fpTk2R09Pt07[NISO], *fpTk2R09Pt09[NISO], *fpTk2R09Pt11[NISO],    
    *fpTk2R10Pt03[NISO], *fpTk2R10Pt05[NISO], *fpTk2R10Pt07[NISO], *fpTk2R10Pt09[NISO], *fpTk2R10Pt11[NISO],  
    *fpTk2R11Pt03[NISO], *fpTk2R11Pt05[NISO], *fpTk2R11Pt07[NISO], *fpTk2R11Pt09[NISO], *fpTk2R11Pt11[NISO];


  ReadBDT *fBdtReader; 
  ReadBDT2 *fBdt2Reader; 
};

#endif
