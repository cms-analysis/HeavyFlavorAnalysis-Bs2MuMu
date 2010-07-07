#ifndef GUARD_HFKINEMATICVERTEXFIT_H
#define GUARD_HFKINEMATICVERTEXFIT_H

#include <vector>

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

// ----------------------------------------------------------------------
class HFKinematicVertexFit {
  
public:
  
  HFKinematicVertexFit(const TransientTrackBuilder *TTB, reco::Vertex &PV, int type = 0, int verbose = 0);
  ~HFKinematicVertexFit();

  void doJpsiFit(std::vector<reco::Track> &trackList, std::vector<int> &trackIndices, std::vector<double> &trackMasses, int type = 0);
  
  int                         fType, fVerbose; 
  reco::Vertex                fPV;   
  const TransientTrackBuilder *fpTTB;   
};

#endif

