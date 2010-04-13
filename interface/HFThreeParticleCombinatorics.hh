#include <vector>
#include <utility>

#include <TLorentzVector.h>

using namespace std;

class triplet{

public:
  triplet(int a, int b, int c) : ka(a) {       // first one is always kaon
     if(b<c){                                  // pions are ordered
           pi1=b; pi2=c;
     } else {
	   pi1=c; pi2=b;}
  };
  
   bool operator==(triplet t) {return (ka==t.ka && pi1==t.pi1 && pi2==t.pi2); };
   unsigned int ka; 
   unsigned int pi1; 
   unsigned int pi2;


};

// ----------------------------------------------------------------------
class HFThreeParticleCombinatorics {
  
public:
  
  HFThreeParticleCombinatorics(int verbose = 0);
  ~HFThreeParticleCombinatorics();
  void combine(vector<triplet > &combList, 
	       vector<pair<int, TLorentzVector> > &tlist1, 
	       vector<pair<int, TLorentzVector> > &tlist2, 
	       double loMass = 0.4, double hiMass = 20., int rmDuplicate = 0);

  int fVerbose; 
};
