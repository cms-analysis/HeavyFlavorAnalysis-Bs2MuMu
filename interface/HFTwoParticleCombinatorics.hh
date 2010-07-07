#ifndef GUARD_HFTWOPARTICLECOMBINATORICS_H
#define GUARD_HFTWOPARTICLECOMBINATORICS_H

#include <vector>
#include <utility>

#include <TLorentzVector.h>

// ----------------------------------------------------------------------
class HFTwoParticleCombinatorics {
  
public:
  
  HFTwoParticleCombinatorics(int verbose = 0);
  ~HFTwoParticleCombinatorics();
  void combine(std::vector<std::pair<int, int> > &combList, 
	       std::vector<std::pair<int, TLorentzVector> > &tlist1, 
	       std::vector<std::pair<int, TLorentzVector> > &tlist2, 
	       double loMass = 0.4, double hiMass = 20., int rmDuplicate = 0);

  int fVerbose; 
};

#endif

