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

class candAna {
  
public:
  candAna(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAna();

  virtual void        evtAnalysis(TAna01Event *evt);
  virtual void        candAnalysis();

  virtual void        processType(); 
  virtual void        genMatch(); 
  virtual void        recoMatch(); 
  virtual void        candMatch(); 
  virtual void        triggerSelection();

  virtual void        bookHist();

  virtual void        readCuts(std::string fileName, int dump = 1);
  virtual void        readFile(std::string fileName, std::vector<std::string> &lines);

  virtual bool        muonID(TAnaTrack *pT);
  virtual std::string splitTrigRange(std::string tl, int &r1, int &r2);

  virtual double      isoClassic(TAnaCand*); 
  virtual double      isoClassicOnePv(TAnaCand*, double r = 1.0, double ptmin = 0.9); 
  virtual double      isoClassicWithDOCA(TAnaCand*, float dca, double r = 1.0, double ptmin = 0.9); 
  virtual double      isoWithDOCA(TAnaCand*, float dca, double r = 1.0, double ptmin = 0.9); 

  std::string fName; 
  std::string fCutFile; 
  TDirectory *fHistDir; 
  bmm2Reader *fpReader; 
  TTree *fTree; 
  TAna01Event *fpEvt;
  TAnaCand *fpCand;

  int fVerbose;
  int fIsMC;

  int fRun, fEvt, fLS;


  double       MASSMIN,   MASSMAX; 
  double       SIGBOXMIN, SIGBOXMAX; 
  double       BGLBOXMIN, BGLBOXMAX; 
  double       BGHBOXMIN, BGHBOXMAX; 
  
  double 
  CANDPTLO, CANDETALO, CANDETAHI
    , CANDCOSALPHA, CANDALPHA
    , CANDFLS3D, CANDFLSXY, CANDVTXCHI2
    , CANDISOLATION, CANDDOCATRK
    , TRACKPTLO, TRACKPTHI, TRACKETALO, TRACKETAHI
    , TRACKTIP, TRACKLIP
    , MUPTLO, MUPTHI
    , MUETALO, MUETAHI, MUIP
    ;
  
  int BLIND, TYPE, SELMODE, MUIDMASK, MUIDRESULT, TRACKQUALITY, TRUTHCAND, IGNORETRIGGER;

  std::map<std::string, pair<int, int> > HLTRANGE;
  //  std::map<std::string, int> HLTRangeMax;


  bool fBarrel, fWideMass; 
  AnalysisCuts fAnaCuts; 

  // -- TM
  int     fGenBTmi; 
  int     fGenM1Tmi, fGenM2Tmi, fNGenPhotons; 
  int     fRecM1Tmi, fRecM2Tmi; 
  int     fCandTmi; 
  int     fGenBpartial; 
  int     fProcessType;
 
  // -- variables for reduced tree, they are from fpCand
  bool    fJSON;
  int     fCandTM, fCandType; 
  int     fMu1TkQuality, fMu2TkQuality, fMu1Q, fMu2Q, fMu1Chi2, fMu2Chi2, fCandQ;
  bool    fMu1Id, fMu2Id;
  double  fHltMu1Pt, fHltMu1Eta, fHltMu1Phi, fHltMu2Pt, fHltMu2Eta, fHltMu2Phi;
  double  fMu1Pt, fMu1Eta, fMu1Phi, fMu2Pt, fMu2Eta, fMu2Phi;
  double  fMu1PtGen, fMu2PtGen, fMu1EtaGen, fMu2EtaGen;
  double  fMu1PtNrf, fMu2PtNrf, fMu1EtaNrf, fMu2EtaNrf; // "now refitted"
  int     fMu1Pix, fMu1BPix, fMu1BPixL1, fMu2Pix, fMu2BPix, fMu2BPixL1;
  double  fMu1W8Mu, fMu1W8Tr, fMu2W8Mu, fMu2W8Tr; 
  double  fPvX, fPvY, fPvZ, fPvNtrk; 
  int     fPvN;
  double  fCandPt, fCandEta, fCandPhi, fCandM, fCandM2, fCandW8Tr, fCandW8Mu; 
  double  fCandCosA, fCandA, fCandChi2, fCandDof, fCandProb, fCandFL3d, fCandFL3dE, fCandFLS3d, fCandFLSxy; 
  double  fCandIso, fCandIso1, fCandIso2, fCandIso3, fCandIso4, fCandIso5;
  double  fIsoR05Pt03, fIsoR05Pt05, fIsoR05Pt07, fIsoR05Pt09, fIsoR05Pt11;
  double  fIsoR07Pt03, fIsoR07Pt05, fIsoR07Pt07, fIsoR07Pt09, fIsoR07Pt11;
  double  fIsoR10Pt03, fIsoR10Pt05, fIsoR10Pt07, fIsoR10Pt09, fIsoR10Pt11;
  int     fCandItrk, fCandI0trk, fCandI4trk; 
  double  fCandDocaTrk, fMu1IP, fMu2IP, fCandPvTip, fCandPvTipE, fCandPvLip, fCandPvLipE; 

  bool    fGoodHLT, fGoodMuonsID, fGoodMuonsPt, fGoodMuonsEta, fGoodTracks, fGoodTracksPt, fGoodTracksEta;
  bool    fGoodQ, fGoodPt, fGoodEta, fGoodCosA, fGoodIso, fGoodChi2, fGoodFLS; 
  bool    fGoodDocaTrk, fGoodIP; 

  bool    fPreselection; 
};

#endif
