/*
 *  HFSequentialVertexFit.h
 *
 *  Created by Christoph on 29.4.10.
 *
 */

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h" 

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDecayTree.h"

class HFSequentialVertexFit
{
 public:
  
  HFSequentialVertexFit(edm::Handle<edm::View<reco::Track> > hTracks, const TransientTrackBuilder *TTB, reco::Vertex &PV, int verbose = 0);
  virtual ~HFSequentialVertexFit();
  
  void doFit(HFDecayTree *tree);
  
  int fVerbose;
  reco::Vertex fPV;
  const TransientTrackBuilder* fpTTB;
  edm::Handle<edm::View<reco::Track> > fhTracks;
  
 private:
  void fitTree(HFDecayTree *tree);
  void saveTree(HFDecayTree *tree);

  // wrapper for the template routine bellow
  inline TAnaCand *addCandidate(HFDecayTree *tree, VertexState &wrtVertexState) {return addCand<VertexState>(tree,wrtVertexState);}
  inline TAnaCand *addCandidate(HFDecayTree *tree, reco::Vertex &wrtVertex) {return addCand<reco::Vertex>(tree,wrtVertex);}

  template<class T>
  TAnaCand *addCand(HFDecayTree *tree, T &toVertex);

  double getParticleMass(int particleID);
};
