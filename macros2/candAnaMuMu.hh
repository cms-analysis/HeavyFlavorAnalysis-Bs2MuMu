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

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaJet.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"

#include "candAna.hh"
#include "bmm2Reader.hh"


class candAnaMuMu : public candAna {
  
public:
  candAnaMuMu(bmm2Reader *pReader, std::string name, std::string cutsFile);
  ~candAnaMuMu();

  void        candAnalysis();
  
  void        processType(); 
  void        genMatch(); 
  void        recoMatch(); 
  void        candMatch(); 
  
  void        bookHist();
  void        readCuts(string filename, int dump);

  double fMMSpecial;
};

#endif
