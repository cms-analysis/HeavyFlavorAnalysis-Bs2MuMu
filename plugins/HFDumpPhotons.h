#ifndef _HFDUMPPHOTONS_h_
#define _HFDUMPPHOTONS_h_

#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


class TFile;
class TTree;
class TAna01Event;

// ----------------------------------------------------------------------
class HFDumpPhotons : public edm::EDAnalyzer {
 public:
  explicit HFDumpPhotons(const edm::ParameterSet&);
  ~HFDumpPhotons();

 private:
  virtual void              beginJob();
  virtual void              analyze(const edm::Event&, const edm::EventSetup&);
  virtual void              endJob();
  
  edm::InputTag             fPFLabel;  
  edm::InputTag             fPhotonsLabel;
  
  int                       fVerbose, fDoTruthMatching; 
  bool                      fRunOnAOD;

};

#endif
