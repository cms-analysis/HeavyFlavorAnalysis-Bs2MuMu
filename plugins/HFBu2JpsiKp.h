#ifndef _HFBU2JPSIKP_h_
#define _HFBU2JPSIKP_h_

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
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"


// ----------------------------------------------------------------------

class HFBu2JpsiKp : public edm::EDAnalyzer {
 public:
  explicit HFBu2JpsiKp(const edm::ParameterSet&);
  ~HFBu2JpsiKp();

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose; 
  edm::InputTag   fTracksLabel, fPrimaryVertexLabel;
  edm::InputTag fMuonsLabel;

  double        fMuonPt; 
  int           fPsiMuons; 
  double        fPsiWindow, fBuWindow;
  double        fTrackPt, fDeltaR;
  int           fVertexing, fType;

  reco::Vertex  fPV;

  edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
