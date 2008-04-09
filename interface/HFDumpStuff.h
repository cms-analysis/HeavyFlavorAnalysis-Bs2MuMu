#ifndef _HFDUMPSTUFF_h_
#define _HFDUMPSTUFF_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TFile;
class TTree;
class TAna00Event;

class TrackAssociatorBase;

// ----------------------------------------------------------------------
class HFDumpStuff : public edm::EDAnalyzer {
 public:
  explicit HFDumpStuff(const edm::ParameterSet&);
  ~HFDumpStuff();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  std::string          fGenEventScaleLabel, fPrimaryVertexLabel,
                       fCandidates1Label, fCandidates2Label, fCandidates3Label, 
                       fMETLabel, fGenMETLabel;

};

#endif
