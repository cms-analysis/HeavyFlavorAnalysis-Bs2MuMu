#ifndef TREEBMM_H
#define TREEBMM_H

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

#include "..//rootio/TAna00Event.hh"
#include "../rootio/TGenCand.hh"
#include "../rootio/TAnaCand.hh"
#include "../rootio/TAnaTrack.hh"
#include "../rootio/TAnaVertex.hh"


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
  void        initVariables();   // returns number of reconstructed Bs candidates
  void        vertexResolution();
  void        eventProcessing();

  void        isMC(int i) {fIsMC = i;}
  void        debugLevel(int i) {fDebug = i;}
  void        dumpEvent();

  void        sameSign(int i) {fSameSign = i;}
  void        isSignal(int i)  {fIsSignal = i;}

  void        kinematicSelection(int offset = 0);
  void        processDiscrimination();
  void        muonMisIdRate();
  void        L1Selection(int offset = 100);
  void        cutJunk(int offset = 200);
  void        HLTSelection(int offset = 300);
  void        trackSelection(int offset = 400);
  void        pidSelection(int offset = 500);
  void        candidateSelection(int offset = 600);
  void        normSample(int offset = 700);
  void        kaonPlots();

  void        fillHist();
  void        book(int offset);
  void        book2(int offset);
  void        histogram(int offset); 
  void        fillAnalysisEff(); 

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

  TAna00Event*fpEvt; 
  TAnaTrack  *fpL1, *fpL2, *fpK;
  TAnaCand   *fpB, *fpJpsi;       
  int        fSTI[3];         // Index for signal track

  std::vector<int>    fKTI;        // Index for all kaon candidates 
  std::vector<int>    fKCI;        // Index for all kaon candidte tracks 

  int fDebug   // bit encoded: {....} = track (8), B cand (4), Event (2).
    , fIsMC, fIsSignal, fSameSign;

  // -- Histogram pointers 
  TTree *fTree;

  TH1D  *fER1, *fAR1, *fMisID;

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
    , ISOVETO
    , ISOLATION
    , ISOCONE
    ; 

  // - generator cuts
  double fMinPt
    , fMaxEta
    ;

  int TYPE;

  // -- Vertex resolution
  double fDvtxPerp
    , fDvtxPar
    , fDpVtxX
    , fDpVtxY
    , fDpVtxZ
    , fFltR
    , fFltG
    , fFltRes
    ;

  // -- Event selection counters
  int fAllEvent
    , fGoodCand
    , fGoodPV
    , fGoodTT
    , fGoodPid
    , fGoodHLTTrigger
    , fGoodEvent
    , fGoodL1Trigger
    , fGoodKinematics
    , fProcessType
    ;

  int fSkipEvent;

  // -- Boxes
  int fSignalBox, fLoSide, fHiSide;

  // -- Event variables
  int fRun
    , fD0
    , fD1
    , fD2
    , fGenMuIDL0
    , fGenMuIDL1
    , fGenMuIDK
    ;

  double fMuIDL0
    , fMuIDL1
    , fMuIDK;


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



  // -- Variables for THE B candidate
  double fMass
    , fPt, fGPt
    , fP
    , fEta
    , fTau
    , fTxy
    , fRMM
    , fIso
    , fIso05, fIso06, fIso08, fIso10, fIso12, fIso14
    ;

  int fIsoVeto
    , fIsoVeto05, fIsoVeto06, fIsoVeto08, fIsoVeto10, fIsoVeto12, fIsoVeto14;

  // -- Variables for THE Jpsi candidate
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
