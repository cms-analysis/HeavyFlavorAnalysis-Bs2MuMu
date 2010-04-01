#ifndef _HFB2D0DS_H_
#define _HFB2D0DS_H_

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

class HFB2D0Ds : public edm::EDAnalyzer {
public:
  explicit HFB2D0Ds(const edm::ParameterSet&);
  ~HFB2D0Ds();
  
  TAnaVertex DoVertexFit(std::vector<reco::Track> &Tracks);
  bool KineCut(reco::TrackCollection::const_iterator iTrack);
  bool massCut(TLorentzVector Track, double mass, double sigma = 0.4);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  int           fVerbose; 
  edm::InputTag fTracksLabel, fPrimaryVertexLabel;
  
  double        fPionPt, fKaonPt;
  double        fDeltaMD0, fDeltaMDs, fSlowPt;
  double        fDeltaR, fVtxErr, fVtxProb, fVtxAngle, fVtxRDist, fVtxADist;
  
  reco::Vertex  fPV;

  edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
