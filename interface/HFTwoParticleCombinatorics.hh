#include <vector>
#include <utility>

#include <TLorentzVector.h>

using namespace std;

// ----------------------------------------------------------------------
class HFTwoParticleCombinatorics {
  
public:
  
  HFTwoParticleCombinatorics(int verbose = 0);
  ~HFTwoParticleCombinatorics();
  void combine(vector<pair<int, int> > &combList, 
	       vector<pair<int, TLorentzVector> > &tlist1, 
	       vector<pair<int, TLorentzVector> > &tlist2, 
	       double loMass = 0.4, double hiMass = 20., int rmDuplicate = 0);

  int fVerbose; 
};
