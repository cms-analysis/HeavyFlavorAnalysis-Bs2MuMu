#ifndef _HFBS2JPSIPHI_h_
#define _HFBS2JPSIPHI_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"


#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"


// ----------------------------------------------------------------------

class HFBs2JpsiPhi : public edm::EDAnalyzer {
 public:
  explicit HFBs2JpsiPhi(const edm::ParameterSet&);
  ~HFBs2JpsiPhi();

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose; 
  edm::InputTag fTracksLabel, fPrimaryVertexLabel;
  edm::InputTag fMuonsLabel;

  double        fMuonPt; 
  int           fPsiMuons;
  double        fPsiWindow, fPhiWindow, fBsWindow; 
  double        fTrackPt, fDeltaR;
  int           fVertexing, fType; 

  reco::Vertex  fPV;

  edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
