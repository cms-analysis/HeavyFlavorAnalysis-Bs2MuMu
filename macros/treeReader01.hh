#ifndef TREEREADER01_H
#define TREEREADER01_H

#include <iostream>

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
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#define DR      57.29577951
#define PIPMASS 0.13957
#define ELMASS  0.000511
#define MUMASS  0.10567


class treeReader01 {
public:
  treeReader01(TChain *tree, TString evtClassName);
  virtual      ~treeReader01();
  virtual void init(TString evtClassName);

  virtual void openHistFile(TString filename);
  virtual void closeHistFile();
  virtual void bookHist();
  virtual void readCuts(TString filename, int dump = 1);

  virtual void startAnalysis();
  virtual void endAnalysis();
  virtual int  loop(int nevents = 1, int start = -1);
  virtual void eventProcessing();
  virtual void initVariables(); 
  virtual void fillHist();
  virtual bool goodRun();
  virtual void setVerbosity(int f) {fVerbose = f;}
  virtual void setMC(int f) {std::cout << Form("setMC(%d)", f) << std::endl; fIsMC = f;}
  virtual void runBlind() {std::cout << "runBlind()" << std::endl; BLIND = 1;}

protected:

  TChain      *fpChain;        // pointer to the analyzed TTree or TChain
  TFile       *fpHistFile;     // for output histograms and reduced trees
  TString      fChainFileName; // the name of the chain file
  TString      fCutFile;       // contains file with the cut definitions

  TAna01Event *fpEvt; 

  // -- Pre-filled variables
  int          fNentries;      // number of events in chain; filled in treeReader01::treeReader01()
  int          fEvent;         // current sequential event number in chain; filled in treeReader01::loop()
  int          fEvt;           // current event number; filled in treeReader01::loop()
  int          fRun;           // current run number; filled in treeReader01::loop()
  int          fLS;            // current lumi section; filled in treeReader01::loop()


  // -- Histogram pointers 
  TTree       *fTree;


  // -- Cut values
  double 
      PTLO
    , PTHI
    , ETALO
    , ETAHI   
    ;
  int TYPE;
  int BLIND; 

  int fVerbose;
  int fIsMC;
};

// ----------------------------------------------------------------------
inline void mk4Vector(TLorentzVector &p4, const Double_t p, const Double_t t, const Double_t f, const Double_t m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

#endif
