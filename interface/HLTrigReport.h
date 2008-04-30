#ifndef HLTrigReport_h
#define HLTrigReport_h
 
/** \class HLTrigReport
 *
 *  
 *  This class is an EDAnalyzer implementing TrigReport (statistics
 *  printed to log file) for HL triggers
 *
 *  $Date: 2008/04/11 15:06:59 $
 *  $Revision: 1.2 $
 *
 *  \author Martin Grunewald
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/TriggerNames.h"
 

#include <iostream>
#include<vector>
#include<string>

//
// class declaration
//
 
class HLTrigReport : public edm::EDAnalyzer {
 
 public:
  explicit HLTrigReport(const edm::ParameterSet&);
  ~HLTrigReport();
  virtual void beginJob();
  virtual void endJob();
 
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
 
 private:
 
  edm::InputTag hlTriggerResults_;  // Input tag for TriggerResults
  edm::TriggerNames triggerNames_;  // TriggerNames class
 
  unsigned int  nEvents_;           // number of events processed
 
  unsigned int  nWasRun_;           // # where at least one HLT was run
  unsigned int  nAccept_;           // # of accepted events
  unsigned int  nErrors_;           // # where at least one HLT had error
 
  std::vector<unsigned int> hlWasRun_; // # where HLT[i] was run
  std::vector<unsigned int> hlAccept_; // # of events accepted by HLT[i]
  std::vector<unsigned int> hlErrors_; // # of events with error in HLT[i]
 
  std::vector<std::string>  hlNames_;  // name of each HLT algorithm
  bool init_;                          // vectors initialised or not
 
};
 
#endif //HLTrigReport_h
