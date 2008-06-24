#ifndef _BSDUMPSTUFF_h_
#define _BSDUMPSTUFF_h_

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
class BSDumpStuff : public edm::EDAnalyzer {
 public:
  explicit BSDumpStuff(const edm::ParameterSet&);
  ~BSDumpStuff();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int                  fNevt;
  int                  fVerbose;

  std::string          fGenEventScaleLabel;
  std::string          fPrimaryVertexLabel;

};

#endif