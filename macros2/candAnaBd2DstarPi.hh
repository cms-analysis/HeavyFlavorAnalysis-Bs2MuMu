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

#include "candAna.hh"

class candAnaBd2DstarPi : public candAna {
  
public:
  candAnaBd2DstarPi(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAnaBd2DstarPi();

  void        candAnalysis();
  void        efficiencyCalculation();
  void        moreBasicCuts();
  void        moreReducedTree(TTree *);

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

  double fPTDs, fPTD0, fPTPiS, fPTPiB, fPTPi, fPTKa;


};

#endif
