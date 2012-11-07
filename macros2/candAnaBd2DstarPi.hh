#ifndef CANDANABD2DSTARPI_H
#define CANDANABD2DSTARPI_H

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

#include "candAna.hh"
#include "bmm2Reader.hh"


class candAnaBd2DstarPi : public candAna {
  
public:
  candAnaBd2DstarPi(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAnaBd2DstarPi();

  void        candAnalysis();
  void        efficiencyCalculation();
  void        moreBasicCuts();

  void        genMatch(); 
  void        recoMatch(); 
  void        candMatch(); 
  
  void        readCuts(string filename, int dump);
  void        bookHist();
  void        fillCandidateHistograms(int offset);


  // -- Additional variables and cuts for Bd -> Dstar pi
  int               D0TYPE, DSTARTYPE; 
  double            MD0LO, MD0HI, DELTAR, DELTAM;
  double            fDeltaR, fMD0, fMDs, fDeltaM;
  bool              fGoodDeltaR, fGoodMD0, fGoodMDs, fGoodDeltaM;

  double fPTD0, fPTPiS, fPTPi2, fPTPi, fPTKa;


};

#endif
