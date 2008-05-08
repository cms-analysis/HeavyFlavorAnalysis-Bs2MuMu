#ifndef _BSDUMPSIGNAL_h_
#define _BSDUMPSIGNAL_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ----------------------------------------------------------------------
class BSDumpSignal : public edm::EDAnalyzer {
 public:
  explicit BSDumpSignal(const edm::ParameterSet&);
  ~BSDumpSignal();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose; 
  std::string   fTracksLabel;
  edm::InputTag fMuonsLabel;
  std::string   fPrimaryVertexLabel;
  double        fMuonPt; 
  double        fKaonPt; 
  double        fDeltaR;

};

#endif
