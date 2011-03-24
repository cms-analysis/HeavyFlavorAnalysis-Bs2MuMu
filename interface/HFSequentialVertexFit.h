/*
 *  HFSequentialVertexFit.h
 *
 *  Created by Christoph on 29.4.10.
 *
 */

#ifndef GUARD_HFSEQUENTIALVERTEXFIT_H
#define GUARD_HFSEQUENTIALVERTEXFIT_H

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h" 

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFDecayTree.h"

class HFSequentialVertexFit
{
 public:
  
  HFSequentialVertexFit(edm::Handle<edm::View<reco::Track> > hTracks, const TransientTrackBuilder *TTB, edm::Handle<reco::VertexCollection> pvCollection, const MagneticField *field, int verbose = 0);
  virtual ~HFSequentialVertexFit();
  
  void doFit(HFDecayTree *tree);
  
 private:
  bool fitTree(HFDecayTree *tree);
  void saveTree(HFDecayTree *tree);
  
  double getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
  double getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles);

  // wrapper for the template routine bellow
  TAnaCand *addCandidate(HFDecayTree *tree, VertexState *wrtVertexState = NULL);
  float getParticleMass(int particleID, float *mass_sigma);
 private: // instance variables
	int fVerbose;
	const TransientTrackBuilder* fpTTB;
	edm::Handle<edm::View<reco::Track> > fhTracks;
	edm::Handle<reco::VertexCollection> fPVCollection;
	const MagneticField* magneticField;
};

#endif

