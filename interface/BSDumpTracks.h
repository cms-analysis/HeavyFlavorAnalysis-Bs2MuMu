#ifndef _BSDUMPTRACKS_h_
#define _BSDUMPTRACKS_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include <iostream>
#include <string>
#include <map>
#include <set>

#include <TROOT.h>
#include <TSystem.h>

class TFile;
class TTree;
class TAna00Event;

class TrackAssociatorBase;

// ----------------------------------------------------------------------
class BSDumpTracks : public edm::EDAnalyzer {
 public:
  explicit BSDumpTracks(const edm::ParameterSet&);
  ~BSDumpTracks();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int idRecTrack(double pt, double eta, double phi, double ept = 0.2, double eeta = 0.01, double ephi = 0.01);

  const reco::TrackCollection      *theTkCollection;

  int                  fNevt;

  int                  fVerbose;
  std::string          fGenEventLabel; 
  std::string          fSimTracksLabel;
  std::string          fTrackingParticlesLabel;
  std::string          fTracksLabel;
  std::string          fAssociatorLabel;
  std::string          fL1MuLabel;


  edm::InputTag        fMuonsLabel1;
  edm::InputTag        fMuonsLabel2;
  int                  fDoTruthMatching;

  TrackAssociatorBase *fAssociator;
};

#endif
