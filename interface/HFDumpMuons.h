#ifndef _HFDUMPMUONS_h_
#define _HFDUMPMUONS_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ----------------------------------------------------------------------
class HFDumpMuons : public edm::EDAnalyzer {
 public:
  explicit HFDumpMuons(const edm::ParameterSet&);
  ~HFDumpMuons();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose; 
  edm::InputTag fMuonsLabel;

};

#endif
