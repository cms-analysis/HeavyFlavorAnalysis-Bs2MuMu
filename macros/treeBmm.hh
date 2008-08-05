#ifndef TREEBMM_H
#define TREEBMM_H

#define NSUBSEL 4

#include <iostream>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"


using namespace std;

#define DR      57.29577951
#define PIPMASS 0.13957
#define ELMASS  0.000511
#define MUMASS  0.10567


class treeBmm {
public:
  treeBmm(TChain *tree, TString evtClassName);
  ~treeBmm();
  void        init(TString evtClassName);

  void        readCuts(TString filename, int dump = 1, double pt = 0., double eta = 10.);
  void        decayChannel(TString signature, int dump = 1);
  void        openHistFile(TString filename);
  void        closeHistFile();
  void        chainFileName(const char *name);
  void        bookHist();
  void        setTitles(TH1 *h, const char *sx, const char *sy,
			float size = 0.06, float xoff = 1.2, float yoff = 1.4, float lsize = 0.06, int font = 132);
  void        setTitles2(TH2 *h, const char *sx, const char *sy,
			float size = 0.06, float xoff = 1.2, float yoff = 1.4, float lsize = 0.06, int font = 132);

  void        startAnalysis();
  int         loop(int nevents = 1, int start = -1);
  int         getSigCand(int cand = 531, int sel = 2, int crit = 1);
  int         getNormCand(int cand = 521, int subcand = 433, int sel = 2, int crit = 1);
  void        eventProcessing();

  void        isMC(int i) {fIsMC = i;}
  void        debugLevel(int i) {fDebug = i;}
  void        candSelEff(int i) {fOffCand = i;}
  void        triggerDecisions(int i, int j) {SETL1 = i; SETHLT = j;}
  void        dumpEvent();

  void        sameSign(int i) {fSameSign = i;}
  void        isSignal(int i)  {fIsSignal = i;}

  void        processDiscrimination();
  void        muonEfficiency();
  void        vertexResiduals();
  void        ptResiduals();

  void        initVariables();
  void        kinematicSelection(int offset = 0);
  void        L1Selection(int offset = 100);
  void        HLTSelection(int offset = 200);
  void        candidateSelection(int offset = 300);
  void        trackSelection(int offset = 400);
  void        trackProperties();
  void        candidateProperties();
  void        kaonCandidates();

  void        fillHist();
  void        fillTriggerEff();
  void        fillAnalysisEff(); 
  void        offlineEff(int bs_index, int sel); 
  void        book(int offset);
  void        book2(int offset);
  void        histogram(int offset); 
  void        signalPlots(); 
  void        fillEventStats(); 

  int         genPartType(int i);
  int         getNrOfDaughters(int i);

  void        boxes();


private:

//   void        dumpSmallTree();
//   TFile      *fFSarah;
//   TTree      *fTSarah;
//   int         fGm0, fGm1;

  TChain     *fpChain;        // pointer to the analyzed TTree or TChain
  TFile      *fpHistFile;     // for output histograms and reduced trees
  TString     fCutFile;       // contains file with the cut definitions
  TString     fChannel;       // Channel signature (used for rare background)
  TString     fChainFileName; // the name of the chain file
  int         fNentries;      // number of events in chain
  int         fEvent;         // current event number

  TAna00Event *fpEvt;
  TAnaTrack   *fpL1, *fpL2, *fpK;
  TAnaCand    *fpB, *fpJpsi;

  int         fBCI;            // Index of B candidate 
  int         fSTI[3];         // Index of signal tracks
  int         fSel, fSubSel, fOffCand;
  double      fMassB;

  std::vector<int>    fKTI;        // Index for all kaon candidates 
  std::vector<int>    fKCI;        // Index for all kaon candidte tracks 

  int fDebug   // bit encoded: {....} = track (8), B cand (4), Event (2).
    , fIsMC, fIsSignal, fSameSign;

  // -- Histogram pointers 
  TTree *fTree;


  TH1D *fNgen, *fNrec, *fNglb, *fErec, *fEglb;
  TH1D *fNgenJ, *fNdecJ, *fNgenB, *fNdecB;
  TH2D *fNR0, *fNR1, *fNR2, *fNR3, *fNR4, *fNR5, *fNR6;
  TH1D *fER1, *fAR1, *fMisID;

  TH1D *fOR1[NSUBSEL];
  int  fSubCand[NSUBSEL];

  // -- Cut values
  double PTBS
    , PTLO
    , PTHI
    , ETALO
    , ETAHI   
    , RMMLO
    , RMMHI   
    , TIPHI
    , VTXCHI
    , COSALPHA
    , LXYLO
    , LXYSXYLO
    , L3DLO
    , L3DS3DLO
    , PRESEL
    , MASSLO
    , MASSHI
    , MASSWI
    , ISOLATION
    , ISOCONE
    , ISOPTLO
    ; 

  int TYPE
    , ISOVETO
    , BMMSEL
    , SUBSEL
    , SETL1
    , SETHLT
    ;

  // - generator cuts
  double fMinPt
    , fMaxEta
    ;

  int fgB, fgJ, fgBmm, fgJmm, fgrB, fgrJ, fnB, fnJ;
  int fgMu, fnMu, frMu, fgK, frK;
  int fnSigMu, frSigMu, frSigK;
  int fnTmMu, fnTmK;
  int fSig1, fSig2, fSig3, fSigB, fSigJ;

  // -- Vertex resolution
  double fDvtxPerp
    , fDvtxPar
    , fDpVtxX
    , fDpVtxY
    , fDpVtxZ
    , fDsVtxX
    , fDsVtxY
    , fDsVtxZ
    , fFltR
    , fFltG
    , fFltRes
    , fFltR3D
    , fFltG3D
    , fFltRes3D
    ;

  // -- Event selection counters
  int fGoodKinematics
    , fGoodL1
    , fGoodHLT
    , fGoodPV
    , fGoodCand
    , fGoodEvent
    , fGoodAna
    , fGoodAnaF
    , fGoodIsoF
    , fGoodVtxF
    , fGoodTT
    , fProcessType
    ;

  int fGoodPtMu
    , fGoodRmm
    , fGoodMass
    , fGoodWindow
    , fGoodPtB
    , fGoodEtaB
    , fGoodCosa
    , fGoodLength
    , fGood3DLength
    , fGoodIso
    , fGoodVtx
    , fGoodPresel
    ;

  int fSkipEvent;

  // -- Boxes
  int fSignalBox, fLoSide, fHiSide;

  // -- Event variables
  int fNorm;

  double fEvtWeight;

  int fRunNr
    , fEvtNr
    , fD0
    , fD1
    , fD2
    , fMcL0
    , fMcL1
    , fMcK
    , fMomL0
    , fMomL1
    , fMomK
    , fGMoL0
    , fGMoL1
    , fGMoK
    , fTruthL0
    , fTruthL1
    , fTruthK
    , fTruthJ
    , fTruthB
    ;

  double fTkL0
    , fTkL1
    , fTkK
    , fLvL0
    , fLvL1
    , fLvK
    , fMuL0
    , fMuL1
    , fMuK
    ;


  // -- Secondary vertex variables
  double fL3d
    , fS3d
    , fLxy
    , fSxy
    , fChi2
    , fNdof
    , fProb
    ;

  // -- Variables for muon tracks / kaon track
  double fQL0, fQL1, fQK
    , fPtL0, fPtL1, fPtK
    , fDptL0, fDptL1, fDptK
    , fEtaL0, fEtaL1, fEtaK
    , fTipL0, fTipL1, fTipK
    ;



  // -- Variables for the genereated and reconstr. B candidate
  double fgQL0
    , fgQL1
    , fgMuL0
    , fgMuL1
    , fgPtL0
    , fgPtL1
    , fgEtaL0
    , fgEtaL1
    , fgChi2
    , fgS3d
    , fgRMM
    ;

  // -- Variables for the chosen B candidate
  double fMass
    , fPt, fGPt
    , fP
    , fEta
    , fTau
    , fTxy
    , fDphi
    , fDeta
    , fRMM
    , fIso
    , fIso05, fIso06, fIso08, fIso10, fIso11, fIso12, fIso14
    ;

  int fIsoVeto
    , fIsoVeto05, fIsoVeto06, fIsoVeto08, fIsoVeto10, fIsoVeto11, fIsoVeto12, fIsoVeto14;

  // -- Variables for the chosen Jpsi candidate
  double fMassJ
    , fPtJ, fGPtJ
    , fPJ
    , fEtaJ
    , fTauJ
    , fTxyJ
    , fRMJ1, fRMJ2, fRKJ
    ;

  // -- B and SV
  double fCosAngle
    , fCosAngle3
  ;
};

// ----------------------------------------------------------------------
inline void mk4Vector(TLorentzVector &p4, const Double_t p, const Double_t t, const Double_t f, const Double_t
 m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

#endif
