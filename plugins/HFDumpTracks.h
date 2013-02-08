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
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TSimpleTrack.hh"

#define HFMAXTRACK 10000

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
  virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag        fTracksLabel, fPrimaryVertexLabel, fBeamSpotLabel,
                       fGenEventLabel, fSimTracksLabel,
                       fAssociatorLabel, fTrackingParticlesLabel;
  edm::InputTag        fMuonsLabel;

  int                  fVerbose, fDoTruthMatching;
  bool                 fDumpSimpleTracks, fDumpRecTracks;

  reco::Vertex         fPV;

  const TrackAssociatorBase *fAssociator;
  PropagateToMuon fPropMuon;

};

#endif
