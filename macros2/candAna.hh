#ifndef CANDANA_H
#define CANDANA_H

#include <string>
#include <vector>
#include <map>

#include <TH1.h>
#include <TH2.h>

#include "bmm2Reader.hh"
#include "RedTreeData.hh"
#include "preselection.hh"
#include "../macros/AnalysisCuts.hh"


struct isoNumbers {
  double iso; 
  int    pvTracks; 
  int    clTracks; 
  int    Tracks; 
};

// -- TMVA related
#include "TMVA/Reader.h"
struct readerData {
  float pt, eta, m1eta, m2eta, m1pt, m2pt;
  float fls3d, alpha, maxdoca, pvip, pvips, iso, docatrk, chi2dof, closetrk; 
  float m;
};

class TTree; 
class TDirectory; 

// ----------------------------------------------------------------------
class candAna {
  
public:
  candAna(bmm2Reader *pReader, std::string name, std::string cutsFile);
  virtual ~candAna();

  virtual void        evtAnalysis(TAna01Event *evt);
  virtual void        candAnalysis();
  virtual void        efficiencyCalculation();
  virtual void        setupReducedTree(TTree *);
  
  virtual int         nearestPV(int pvIdx, double maxDist = 99.);
  virtual void        getSigTracks(std::vector<int> &v, TAnaCand *pC);
  virtual double      constrainedMass();
  virtual void        runRange();
  virtual void        genMatch(); 
  virtual void        recoMatch(); 
  virtual void        candMatch(); 
  virtual void        triggerSelection();
  virtual void        fillCandidateHistograms(int offset);
  virtual void        fillRedTreeData();
  virtual void        replaceAll(std::string &s, std::string a, std::string b);

  virtual TMVA::Reader* setupReader(std::string xmlFile, readerData &rd);
  virtual void        calcBDT();
  virtual int         detChan(double m1eta, double m2eta);
    
  virtual void        bookHist();

  virtual void        basicCuts();
  virtual void        moreBasicCuts();
  virtual void        candidateCuts();
  virtual void        moreCandidateCuts();

  virtual void        readCuts(std::string fileName, int dump = 1);
  virtual void        readFile(std::string fileName, std::vector<std::string> &lines);

  virtual bool        highPurity(TAnaTrack *pt);
  virtual bool        tightMuon(TAnaTrack *pt);
  virtual bool        tightMuon(TSimpleTrack *pt);

  virtual std::string splitTrigRange(std::string tl, int &r1, int &r2);

  virtual double      isoClassicWithDOCA(TAnaCand*, double dca, double r = 0.7, double ptmin = 0.9); 
  virtual int         nCloseTracks(TAnaCand*, double dca, double pt = 0.5); 

  virtual TAnaCand*   osCand(TAnaCand *pC);
  virtual double      osIsolation(TAnaCand *pC, double r = 1.0, double ptmin = 0.9); 
  virtual int         osMuon(TAnaCand *pC, double r = 1.0); 
  virtual bool        doTriggerMatching(); // match the 2 muons from the dimuon to HLT
  virtual bool        doTriggerMatching(TAnaTrack *pt, bool anyTrig = false); // match a single track to HLT
  virtual void        boostGames();
  virtual double      matchToMuon(TAnaTrack *pt, bool skipSame = false); // match a single track to ALL muons
  virtual void        play(); 
  // To return the full deltaR not just a bool
  virtual double      doTriggerMatchingR(TAnaTrack *pt, bool anyTrig = false); // match a single track to HLT

  std::string fName; 
  std::string fCutFile; 
  TDirectory *fHistDir; 
  bmm2Reader *fpReader; 
  TTree *fTree, *fAmsTree; 
  TAna01Event *fpEvt;
  TAnaCand *fpCand, *fpOsCand;
  int fCandIdx; 

  int fVerbose;
  int fIsMC;

  Long64_t fRun, fEvt;
  int fLS;
  int fEvent; 
  int fRunRange;
  int fYear;

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
  double  fGenLifeTime, fGenMass; 

  // -- BDT
  std::vector<TMVA::Reader*> fReaderEvents0; 
  std::vector<TMVA::Reader*> fReaderEvents1; 
  std::vector<TMVA::Reader*> fReaderEvents2; 
  double  fBDT; 
  readerData frd; 


  // -- variables for reduced tree, they are from fpCand
  bool    fJSON, fCowboy;
  int     fCandTM, fCandType; 
  int     fMu1TkQuality, fMu2TkQuality, fMu1Q, fMu2Q, fCandQ, fMu1PV, fMu2PV;
  bool    fMu1Id, fMu2Id, fHLTmatch;  
  double  fMuDist, fMuDeltaR, fMu1Chi2, fMu2Chi2;
  double  fHltMu1Pt, fHltMu1Eta, fHltMu1Phi, fHltMu2Pt, fHltMu2Eta, fHltMu2Phi;
  double  fMu1Pt, fMu1Eta, fMu1Phi, fMu2Pt, fMu2Eta, fMu2Phi;
  double  fMu1PtGen, fMu2PtGen, fMu1EtaGen, fMu2EtaGen, fMu1PhiGen, fMu2PhiGen;
  double  fMu1PtNrf, fMu2PtNrf, fMu1EtaNrf, fMu2EtaNrf; // "now refitted"
  int     fMu1TrkLayer, fMu1Pix, fMu1BPix, fMu1BPixL1, fMu2TrkLayer, fMu2Pix, fMu2BPix, fMu2BPixL1;
  double  fMu1W8Mu, fMu1W8Tr, fMu2W8Mu, fMu2W8Tr; 
  double  fPvX, fPvY, fPvZ, fPvNtrk, fPvNdof, fPvAveW8; 
  int     fPvN;
  double  fCandPt, fCandP, fCandTau, fCandEta, fCandPhi, fCandM, fCandME, fCandM2, fCandW8Tr, fCandW8Mu; 
  double  fCandCosA, fCandA;
  double  fCandChi2, fCandDof, fCandProb, fCandFL3d, fCandFL3dE, fCandFLS3d, fCandFLxy, fCandFLSxy, fCandDoca; 
  double  f2MChi2,   f2MDof,   f2MProb,   f2MFL3d,   f2MFL3dE,   f2MFLS3d,   f2MFLSxy; 
  double  fCandIso;
  int     fCandIsoTrk, fCandCloseTrk, fCandPvTrk, fCandI0trk, fCandI1trk, fCandI2trk; 
  double  fCandDocaTrk, fCandDocaTrkBdt, fMu1IP, fMu1IPE, fMu2IP, fMu2IPE; 
  double  fCandPvTip, fCandPvTipE, fCandPvTipS, fCandPvLip, fCandPvLipE, fCandPvLipS, fCandPvIp, fCandPvIpE, fCandPvIpS;
  double  fCandPvIp3D, fCandPvIpE3D, fCandPvIpS3D; 
  double  fCandPvLip2, fCandPvLipS2, fCandPvLip12, fCandPvLipE12, fCandPvLipS12; 

  double  fOsMuonPt, fOsMuonPtRel, fOsIso, fOsRelIso, fOsMuonDeltaR;

  int     fChan; 

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

  bool    fBadEvent;
  int     fhltType; // to hold the HLT information d.k.
 
  struct RedTreeData fRTD;

};

#endif //  CANDANA_H
