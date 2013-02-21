#ifndef _HFTRUTHCANDIDATE_h_
#define _HFTRUTHCANDIDATE_h_

#include <memory>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HFTruthCandidate : public edm::EDAnalyzer {

public:

  explicit HFTruthCandidate(const edm::ParameterSet&);
  ~HFTruthCandidate();
  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

private:
  
  edm::InputTag      fTracksLabel, fPrimaryVertexLabel, fBeamSpotLabel, fMuonsLabel;

  bool               fPartialDecayMatching;
  int                fMotherID, fType, fGenType;
  std::vector<int>   fDaughtersID;

  int                fStableDaughters; 
  std::multiset<int> fDaughtersSet;
  std::multiset<int> fDaughtersGammaSet;
  std::multiset<int> fDaughtersGamma2Set;


  double             fMaxDoca;
  int                fVerbose; 

  edm::ESHandle<TransientTrackBuilder> fTTB;

};

#endif
