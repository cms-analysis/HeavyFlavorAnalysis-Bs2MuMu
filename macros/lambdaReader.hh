#ifndef GUARD_LAMBDAREADER_H 
#define GUARD_LAMBDAREADER_H

#include <iostream>
#include <vector>
#include <utility>

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

#include "treeReader01.hh"

#define DR      57.29577951
#define PIPMASS 0.13957
#define ELMASS  0.000511
#define MUMASS  0.10567

class decayMap
{
    public:
	typedef unsigned int pos_t;
	typedef std::map<std::string,pos_t> map_t;
	typedef std::map<pos_t,std::string> revmap_t;
	decayMap();
	pos_t getPos(std::string key);
	void printMap();
	
    private:
	map_t myMap;
	revmap_t myRevMap;
	void createRevMap();
};

class lambdaReader : public treeReader01 {

public:
  lambdaReader(TChain *tree, TString evtClassName);
  ~lambdaReader();

  void         bookHist();
  void         startAnalysis();
  void         endAnalysis();
  void         eventProcessing();
  void         fillHist();
  void         readCuts(TString filename, int dump = 1);
  void         initVariables();

  virtual void MCKinematics();  
  virtual void L1TSelection();  
  virtual void HLTSelection();  
  virtual void trackSelection();  
  virtual void muonSelection();  
  virtual void candidateSelection(int mode = 0);  // 0 = closest in r-phi
  virtual void fillTMCand(TAnaCand *pCand, int type);

  virtual void doBplus(TAnaCand *pCand);
  virtual void doDzero(TAnaCand *pCand);
  //virtual void doJpsi(TAnaCand *pCand);
  //virtual void doUpsilon(TAnaCand *pCand);


  // -- Cut values
  double 
      CHARMPTLO
    , CHARMETALO
    , CHARMETAHI   
    , KAPTLO
    , PIPTLO
    , VTXCHI2
    , MUPTLO
    , MUPTHI
    , MUETALO
    , MUETAHI   
    ;
  int TYPE, MCTYPE;
  // Cut valuas FRANK
  int cutLbCandidate;
  int cutMuId1, cutMuId2;
  double cutMjpMin, cutMjpMax;
  double cutAlphal0Max, cutMl0Max, cutD3dl0Min, cutPtl0Min;

  // -- Variables
  TAnaCand    *fpCand1, *fpCand2; 

  // -- Candidate variables
  double fmlb, fml0, fmjp; // m
  double fptlb, fptl0, fptjp; // pt
  double fetalb, fetal0, fetajp; // eta
  double fphilb, fphil0, fphijp; // phi

  double fpt1m, fpt2m, fptpr, fptpi; // kinematic variables of granddaughters
  double feta1m, feta2m, fetapr, fetapi;
  double fphi1m, fphi2m, fphipr, fphipi;
  int    fq1m, fq2m, fqpr, fqpi;
  int    fid1m, fid2m; // muon id

  // -- sig track variables
  double fSgptpr, fSgptpi; 
  double fSgetapr, fSgetapi;
  double fSgphipr, fSgphipi;
  
  double fKshypo; // mass of proton as Ks
  
  double falphalb, falphal0; // alpha
  double fmaxdocalb, fmaxdocal0, fmaxdocajp; // maxdoca
  
  double fd3lb, fd3l0, fd3jp;    // 3d distance
  double fd3Elb, fd3El0, fd3Ejp;
  
  double fdxylb, fdxyl0, fdxyjp; // 2d distance
  double fdxyElb, fdxyEl0, fdxyEjp;

  double fchi2lb, fchi2l0, fchi2jp; // quality of vtx fit
  double fndoflb, fndofl0, fndofjp;

  double fdRprpi, fdRmumu, fdRl0jp; // deltaR of pairs for convenience
  
  double fanglbl0; // angle between lb and l0 (from cands, not alpha)

  int fnDaughters, fnGrandDaughters; // generator info
  int fmapRef1Gen; //, fmapRef2Gen;
  int fIsSig;

  // for genCand stuff
  TTree* fGenTree;
  double fgmlb, fgmlbsw; // mass of lb from ppi and swapped mass hypotheses
  double fgml0, fgml0sw;
  double fgptpr, fgptpi, fgptl0;
  double fgppr, fgppi, fgpl0;
  double fgdRprpi;
  int    fgnDaughters, fgnGrandDaughters; // no of daughters from gen truth
  int    fgmapRef1Gen; // pointers to entries in decay maps
  //int    fgmapRef1Gen, fgmapRef2Gen; // pointers to entries in decay maps
  TGenCand *gcPrev;
  decayMap myDecayMap1Gen;
  decayMap myDecayMap2Gen;

};

#endif
