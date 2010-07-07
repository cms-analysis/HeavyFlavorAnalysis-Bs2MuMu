#ifndef GUARD_HFKALMANVERTEXFIT_H
#define GUARD_HFKALMANVERTEXFIT_H

#include <vector>

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

// ----------------------------------------------------------------------
class HFKalmanVertexFit {
  
public:
  
  HFKalmanVertexFit(const TransientTrackBuilder *TTB, reco::Vertex &PV, int type = 0, int verbose = 0);
  ~HFKalmanVertexFit();

  // -- do a Kalman Vertex Fit
  // Note: ntracksToVertex:    how many tracks are used for vertexing
  //       trackMasses[i] < 0: a second invariant mass is computed without the tracks whose mass < 0
  void doFit(std::vector<reco::Track> &trackList, std::vector<int> &trackIndices, std::vector<double> &trackMasses, 
	     int type = 0, int ntracksToVertex = -1);

  void setNoCuts();
  
  // -- Just combine the 4-std::vectors
  void doNotFit(std::vector<reco::Track> &trackList, std::vector<int> &trackIndices, std::vector<double> &trackMasses, int type = 0);
  
  int                         fType, fVerbose; 
  reco::Vertex                fPV;   
  const TransientTrackBuilder *fpTTB;   

  double                      fMaxDoca; 
  double                      fVtxChi2;
  double                      fVtxSigXY;
  double                      fVtxSig3d;
  double                      fCosAngle;
  double                      fPtCand;
};

#endif

