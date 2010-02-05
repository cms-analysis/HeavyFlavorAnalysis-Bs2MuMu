#ifndef _HFTRUTHCANDIDATE_h_
#define _HFTRUTHCANDIDATE_h_

#include <memory>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace edm;
using namespace std;

class HFTruthCandidate : public EDAnalyzer {

public:

  explicit HFTruthCandidate(const ParameterSet&);
  ~HFTruthCandidate();
  
  virtual void beginJob() ;
  virtual void analyze(const Event&, const EventSetup&);
  virtual void endJob() ;

private:
  
  InputTag      fTracksLabel;
  int           fMotherID, fType, fGenType;
  vector<int>   fDaughtersID;

  int           fStableDaughters; 
  multiset<int> fDaughtersSet;
  multiset<int> fDaughtersGammaSet;
  multiset<int> fDaughtersGamma2Set;




};

#endif
