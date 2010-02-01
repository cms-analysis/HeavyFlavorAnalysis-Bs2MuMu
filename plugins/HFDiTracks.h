#ifndef _HFDITRACKS_h_
#define _HFDITRACKS_h_

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

class HFDiTracks : public edm::EDAnalyzer {
 public:
  explicit HFDiTracks(const edm::ParameterSet&);
  ~HFDiTracks();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void         doVertexFit(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2, TAnaCand *pCand);
  void         doVertexFit(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2);
  void         fillCandAndSignal(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2, TAnaCand *pCand);
  void         fillCandAndSignal(std::vector<reco::Track> &Tracks, int iMuon1, int iMuon2);

  int           fVerbose;
  int           fVertexing; 
  std::string   fTracksLabel, fPrimaryVertexLabel;
  edm::InputTag fMuonsLabel;

  double        fMuonPt, fTrackPt, fTrackMass, fMassLow, fMassHigh;
  int           fType; 

  reco::Vertex  fPV;

  edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
