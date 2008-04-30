#ifndef _BMMDUMPSIGNAL_h_
#define _BMMDUMPSIGNAL_h_

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
class BmmDumpSignal : public edm::EDAnalyzer {
 public:
  explicit BmmDumpSignal(const edm::ParameterSet&);
  ~BmmDumpSignal();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose; 
  std::string   fTracksLabel, fPrimaryVertexLabel;
  edm::InputTag fMuonsLabel;
  double        fMuonPt, fKaonPt, fDeltaR;

};

#endif
