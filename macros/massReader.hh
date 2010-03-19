#include "treeReader01.hh"

using namespace std;

class massReader : public treeReader01 {
	
public:
  massReader(TChain *tree, TString evtClassName);
  ~massReader();
  
  void bookHist();
  void eventProcessing();
  
private:
  TTree *reduced_tree;
  int fCandidate;
  TVector3 fMomentum;
  TVector3 *fMomentumPtr;
  double fMass;
};
