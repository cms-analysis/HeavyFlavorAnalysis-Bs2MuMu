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
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


  int           doVertexFit(std::vector<reco::Track> &Tracks, RefCountedKinematicTree &bsTree);
  void          doJpsiVertexFit(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2, TAnaCand *pCand);

  int           fVerbose; 
  edm::InputTag fTracksLabel, fPrimaryVertexLabel;
  edm::InputTag fMuonsLabel;

  double        fMuonPt, fTrackPt, fDeltaR;
  int           fType; 

  reco::Vertex  fPV;

  edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
