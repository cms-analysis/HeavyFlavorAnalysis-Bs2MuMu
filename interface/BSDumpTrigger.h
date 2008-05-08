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
class TAna00Event;


// ----------------------------------------------------------------------
class BSDumpTrigger : public edm::EDAnalyzer {
 public:
  explicit BSDumpTrigger(const edm::ParameterSet&);
  ~BSDumpTrigger();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fNevt;

  int           fVerbose;
  std::string   fparticleMap;
  std::string   fL1MuLabel;
  std::string   fL1TriggerName;
  edm::InputTag fHLTriggerLabel;
  std::string   fHLTriggerName;

  std::string   fHLTFilterObject0;
  std::string   fHLTFilterObject1;
  std::string   fHLTFilterObject2;
  std::string   fHLTFilterObject3;
  std::string   fHLTFilterObject4;

};

