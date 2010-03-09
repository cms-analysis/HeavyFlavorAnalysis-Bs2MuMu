#ifndef _HFDUMPTRACKS_h_
#define _HFDUMPTRACKS_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

class TFile;
class TTree;
class TAna01Event;

class TrackAssociatorBase;

// ----------------------------------------------------------------------
class HFDumpTracks : public edm::EDAnalyzer {
 public:
  explicit HFDumpTracks(const edm::ParameterSet&);
  ~HFDumpTracks();
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag        fTracksLabel, fPrimaryVertexLabel,
                       fGenEventLabel, fSimTracksLabel,
                       fAssociatorLabel, fTrackingParticlesLabel;
  edm::InputTag        fMuonsLabel, fCaloMuonsLabel;

  int                  fVerbose, fDoTruthMatching; 

  reco::Vertex         fPV;

  const TrackAssociatorBase *fAssociator;
};

#endif
