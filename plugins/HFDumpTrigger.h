#ifndef GUARD_HFDUMPTRIGGER_H
#define GUARD_HFDUMPTRIGGER_H

// system include files
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
class HFDumpTrigger : public edm::EDAnalyzer {
 public:
  explicit HFDumpTrigger(const edm::ParameterSet&);
  ~HFDumpTrigger();
  
 private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run &, const edm::EventSetup &);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose;
  int           fNevt;

  std::string   fHLTProcessName; 
  edm::InputTag fL1GTReadoutRecordLabel; 
  edm::InputTag fL1GTmapLabel;
  edm::InputTag fL1MuonsLabel;
  edm::InputTag fTriggerEventLabel;
  edm::InputTag fHLTResultsLabel;
};

#endif
