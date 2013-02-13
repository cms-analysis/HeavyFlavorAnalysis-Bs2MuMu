#ifndef CANDANAMUMU_H
#define CANDANAMUMU_H

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

class candAnaMuMu : public candAna {
  
public:
  candAnaMuMu(bmm2Reader *pReader, std::string name, std::string cutsFile);
  virtual ~candAnaMuMu();

  void        candAnalysis();
  void        efficiencyCalculation();
  
  void        processType(); 
  void        genMatch(); 
  void        recoMatch(); 
  void        candMatch(); 
  
  void        bookHist();
  void        readCuts(string filename, int dump);
};

#endif
