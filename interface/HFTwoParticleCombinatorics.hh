#ifndef GUARD_HFTWOPARTICLECOMBINATORICS_H
#define GUARD_HFTWOPARTICLECOMBINATORICS_H

#include <vector>
#include <utility>

#include <TLorentzVector.h>
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

typedef struct{
  int   id1;
  int   id2;
  TLorentzVector mother;
} HFTwoParticleState;

// ----------------------------------------------------------------------
class HFTwoParticleCombinatorics {
  
public:
  
  HFTwoParticleCombinatorics(int verbose = 0);
  HFTwoParticleCombinatorics(int verbose, const TransientTrackBuilder *fTTB );
  ~HFTwoParticleCombinatorics();
  void combine(std::vector<std::pair<int, int> > &combList, 
	       std::vector<std::pair<int, TLorentzVector> > &tlist1, 
	       std::vector<std::pair<int, TLorentzVector> > &tlist2, 
	       double loMass = 0.4, double hiMass = 20., int rmDuplicate = 0);

  void combine(std::vector<HFTwoParticleState > &combList, 
	       std::vector<std::pair<int, reco::Track> > &tlist1, double mass1,
	       std::vector<std::pair<int, reco::Track> > &tlist2, double mass2,
	       double loMass, double hiMass, double maxDoca, int rmDuplicate);


  const TransientTrackBuilder *ttBuilder;
  int fVerbose; 
};

#endif

