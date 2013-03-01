#ifndef _HFDUMPMUONS_h_
#define _HFDUMPMUONS_h_

#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaMuon.hh"


class TFile;
class TTree;
class TAna01Event;

// ----------------------------------------------------------------------
class HFDumpMuons : public edm::EDAnalyzer {
 public:
  explicit HFDumpMuons(const edm::ParameterSet&);
  ~HFDumpMuons();
  
 private:
  virtual void              beginJob();
  virtual void              beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void              analyze(const edm::Event&, const edm::EventSetup&);
  virtual void              endJob();
  void                      fillMuon(const reco::Muon& tr, int type);
  void                      fillCaloMuon(const reco::CaloMuon& tr, int type);
  void                      extrapolateTracks(); 
  bool                      doExtrapolate(double pt, double eta);
  void                      findVertex(TAnaMuon *anaMu, std::set<unsigned> *trkIcs, double *prob);
  std::vector<unsigned int> muonStatHits(const reco::Track& tr);
  edm::InputTag             fTracksLabel;
  edm::InputTag             fMuonsLabel;
  edm::InputTag             fBeamSpotLabel;
  edm::InputTag             fPrimaryVertexLabel;
  edm::InputTag             fCaloMuonsLabel;
  
  edm::Handle<edm::View<reco::Track> > *fhTracks;
  edm::ESHandle<TransientTrackBuilder> fTTB;
  const reco::MuonCollection *fMuonCollection;
  
  double                    fMaxTrackDistToStore;
  double                    fDocaVertex; // try vertexing only with tracks closer than this
  unsigned                  fKeepBest; // number of candidates to keep for iterative vertex search
  unsigned                  fMaxCandTracks; // max number of tracks for the muon candidate vertex

  int                       fVerbose, fDoTruthMatching; 
  bool                      fRunOnAOD;

  const reco::BeamSpot         *fBeamSpot;
  const reco::VertexCollection *fVertexCollection;
  PropagateToMuon           fpropM1, fpropM2;
  std::vector<xpTrack>      fXpTracks;
};

#endif
