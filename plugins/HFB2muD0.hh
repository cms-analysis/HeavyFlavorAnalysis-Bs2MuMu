#ifndef _HFB2MUD0_H_
#define _HFB2MUD0_H_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "DataFormats/VertexReco/interface/Vertex.h"

// ----------------------------------------------------------------------

class HFB2muD0 : public edm::EDAnalyzer {
 public:
  explicit HFB2muD0(const edm::ParameterSet&);
  ~HFB2muD0();

  TAnaVertex DoVertexFit(std::vector<reco::Track> &Tracks);
  bool KineCut(reco::TrackCollection::const_iterator iTrack);
  bool massCut(TLorentzVector Track, double mass, double sigma = 0.4);
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose; 
  edm::InputTag fTracksLabel, fPrimaryVertexLabel;
  edm::InputTag fMuonsLabel;

  double        fMuonPt, fProtonPt, fPionPt, fDeltaR;
  double        fTrackPt, fDeltaMD0, fDeltaMDs, fSlowPt;

  reco::Vertex  fPV;

  edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
