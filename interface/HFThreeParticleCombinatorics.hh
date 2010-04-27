#include <vector>
#include <utility>

#include <TLorentzVector.h>

using namespace std;

class triplet{

public:
  triplet(int a, int b, int c) : kaon(a) {       // first one is always kaon
     if(b<c){                                    // pions are ordered
           pion1=b; pion2=c;
     } else {
	   pion1=c; pion2=b;}
  };
  
  bool operator==(triplet t) {return (kaon==t.ka() && pion1==t.pi1() && pion2==t.pi2()); };
  unsigned int pi1(){ return pion1;};
  unsigned int pi2(){ return pion2;};
  unsigned int ka(){ return kaon;};
  unsigned int pi3(){ return kaon;};

private:
   unsigned int kaon; 
   unsigned int pion1; 
   unsigned int pion2;
};

// ----------------------------------------------------------------------
class HFThreeParticleCombinatorics {
  
public:
  
  HFThreeParticleCombinatorics(int verbose = 0);
  ~HFThreeParticleCombinatorics();
  void combine(vector<triplet > &combList, 
	       vector<pair<int, TLorentzVector> > &tlist, 
	       double loMass = 0.2, double hiMass = 0.8, int rmDuplicate = 0);
  void combine(vector<triplet > &combList, 
	       vector<pair<int, TLorentzVector> > &tlist1, 
	       vector<pair<int, TLorentzVector> > &tlist2, 
	       double loMass = 0.4, double hiMass = 20., int rmDuplicate = 0);

  int fVerbose; 
};
