#ifndef GUARD_LAMBDAREADER_H 
#define GUARD_LAMBDAREADER_H

#include <iostream>
#include <vector>
#include <utility>
#include <sstream>

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


/*! Class for handling maps of decay chains

  It consists of two maps to handle requests in both ways.
  */
class decayMap
{
    public:
	typedef unsigned int pos_t;
	typedef std::map<std::string,pos_t> map_t;
	typedef std::map<pos_t,std::string> revmap_t;

	decayMap();
	pos_t getPos(std::string key); // returns the key, adds entry if not existing
	void printMap(); // prints the map to cout
	unsigned int readFile(std::string filename); // reads a file into a map
	
    private:
	map_t myMap; // main map
	revmap_t myRevMap; // reverse map, created upon need
	void createRevMap();
};

class lambdaReader : public treeReader01 {

public:
  lambdaReader(TChain *tree, TString evtClassName);
  ~lambdaReader();

  void         bookHist();
  void         bookReducedTree();
  void         startAnalysis();
  void         endAnalysis();
  void         eventProcessing();
  void         fillHist();
  void         readCuts(TString filename, int dump = 1);
  void         initVariables();
  void	       doGenLevelStuff();

  template <typename T> void setCut(T &var, std::string value)
  {
      std::istringstream iss(value);
      iss >> var;
  };

  template <typename T> void setCut(T &var, std::string value, TH1D *hcuts, int ibin, std::string title)
  {
      std::istringstream iss(value);
      iss >> var;
      hcuts->SetBinContent(ibin, var);
      hcuts->GetXaxis()->SetBinLabel(ibin, title.c_str());
  };

  // -- Cut values
  int CUTLbCandidate;
  int CUTMuId1, CUTMuId2;
  double CUTMjpMin, CUTMjpMax;
  double CUTAlphal0Max, CUTMl0Max, CUTD3dl0Min, CUTPtl0Min;
  bool CUTReadDecayMaps, CUTPrintDecayMaps;
  bool CUTgenTreeCandsOnly;
  std::string CUTDecayMap1, CUTDecayMap2;

  // -- Variables
  TAnaCand    *fpCand1, *fpCand2; 

  // -- Candidate variables
  double fmlb, fml0, fmjp; // m
  double fptlb, fptl0, fptjp; // pt
  double fplb, fpl0, fpjp; // p
  double fetalb, fetal0, fetajp; // eta
  double fphilb, fphil0, fphijp; // phi

  double frpt1m, frpt2m, frptpr, frptpi; // kinematic variables of granddaughters
  double freta1m, freta2m, fretapr, fretapi;
  double frphi1m, frphi2m, frphipr, frphipi;
  int    frq1m, frq2m, frqpr, frqpi;
  int    frid1m, frid2m; // muon id

  // -- sig track variables, S is capitalized intentionally for better distinction to the reco-variants
  double fSpt1m,  fSpt2m,  fSptpr, fSptpi;
  double fSeta1m, fSeta2m, fSetapr, fSetapi;
  double fSphi1m, fSphi2m, fSphipr, fSphipi;
  
  double fKshypo; // mass of proton as Ks
  
  double falphalb, falphal0; // alpha
  double fmaxdocalb, fmaxdocal0, fmaxdocajp; // maxdoca
  
  double fd3lb, fd3l0, fd3jp;    // 3d distance
  double fd3Elb, fd3El0, fd3Ejp;
  double fctlb, fctl0; // ctau
  double fbtlbx, fbtlby, fbtlbz; // beta vector
  double fbtl0x, fbtl0y, fbtl0z;
  
  double fdxylb, fdxyl0, fdxyjp; // 2d distance
  double fdxyElb, fdxyEl0, fdxyEjp;

  double fchi2lb, fchi2l0, fchi2jp; // quality of vtx fit
  double fndoflb, fndofl0, fndofjp;
  double fproblb, fprobl0, fprobjp;

  double fchi21m, fchi22m, fchi2pr, fchi2pi; // quality of track fit
  double fprob1m, fprob2m, fprobpr, fprobpi;
  int    fndof1m, fndof2m, fndofpr, fndofpi;
  int    fqual1m, fqual2m, fqualpr, fqualpi;

  double fdRprpi, fdRmumu, fdRl0jp; // deltaR of pairs for convenience
  
  double fanglbl0; // angle between lb and l0 (from cands, not alpha)

  int fnDaughters, fnGrandDaughters; // generator info
  int fmapRef1Gen, fmapRef2Gen;
  int fIsSig; // true if the cand block for this evt contains a signal decay
  int fIsMCmatch; // true if truth matched

  // for genCand stuff
  TTree* fGenTree;
  double fgmlb, fgmlbsw; // mass of lb from ppi and swapped mass hypotheses
  double fgml0, fgml0sw;
  double fgptpr, fgptpi, fgptmu1, fgptmu2, fgptl0;
  double fgetamu1, fgetamu2;
  double fgppr, fgppi, fgpl0;
  double fgdRprpi, fgdRmumu, fgdRl0lb;
  double fganprpi, fganmumu, fganl0lb, fganl0jp;
  double fganl0mumin, fganl0muPt;
  int    fgnDaughters, fgnGrandDaughters; // no of daughters from gen truth
  int    fgmapRef1Gen, fgmapRef2Gen; // pointers to entries in decay maps
  int    fghasCand;
  TGenCand *gcPrev;
  decayMap myDecayMap1Gen;
  decayMap myDecayMap2Gen;

};

#endif
