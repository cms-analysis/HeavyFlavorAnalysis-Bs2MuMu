#ifndef _HFDUMPGENERATOR_h_
#define _HFDUMPGENERATOR_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


class TFile;
class TTree;
class TAna01Event;

// ----------------------------------------------------------------------
class HFDumpGenerator : public edm::EDAnalyzer {
 public:
  explicit HFDumpGenerator(const edm::ParameterSet&);
  ~HFDumpGenerator();
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int         fVerbose;
  std::string fGenCandidatesLabel, fGenEventLabel;

};

#endif
