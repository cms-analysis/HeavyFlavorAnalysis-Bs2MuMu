#ifndef GUARD_HFLAMDAS_H
#define GUARD_HFLAMDAS_H

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

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

// ----------------------------------------------------------------------

class HFLambdas : public edm::EDAnalyzer {
 public:
  explicit HFLambdas(const edm::ParameterSet&);
  ~HFLambdas();

    typedef unsigned int count_t;
    typedef unsigned int index_t;
    typedef double value_t;
    typedef std::pair<index_t,index_t> duplet_t;
    typedef std::vector<std::pair<index_t, TLorentzVector> > trackList_t;
    //typedef std::pair<GlobalPoint,float> poca_t; // point of closest approach pair: point and distance

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  bool inMassWindow(const value_t& v, const value_t& m, const value_t& dm);

  int fillCand(const TLorentzVector& tlvCand, const int& trk1, const int& trk2, const int& type);
  int fillCand(const TLorentzVector& tlvCand, const TwoTrackMinimumDistance& ttmd, const int& trk1, const int& trk2, const int& type);

  //poca_t calculatePoca(const edm::Handle<edm::View<reco::Track> >& tracks, int track1, int track2);
  TwoTrackMinimumDistance calculatePoca(const edm::Handle<edm::View<reco::Track> >& tracks, int track1, int track2);

  int           fVerbose; 
  int           fMaxTracks, fMinTracks; 
  edm::InputTag fTracksLabel, fPrimaryVertexLabel;
  edm::InputTag fMuonsLabel;

  int           fUseMuon;
  double        fPhiWindow, fJPsiWindow, fL0Window; 
  double        fMuonPt, fProtonPt, fPionPt, fTrackPt, fDeltaR;
  double        fMaxDoca, fMaxVtxChi2, fMinVtxSigXY, fMinVtxSig3d, fMinCosAngle, fMinPtCand; 

  int           fType; 

  reco::Vertex  fPV;

  edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
