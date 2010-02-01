#include <vector>

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

using namespace std;
using namespace reco;
using namespace edm;

// ----------------------------------------------------------------------
class HFKinematicVertexFit {
  
public:
  
  HFKinematicVertexFit(const TransientTrackBuilder *TTB, Vertex &PV, int type = 0, int verbose = 0);
  ~HFKinematicVertexFit();

  void doJpsiFit(vector<Track> &trackList, vector<int> &trackIndices, vector<double> &trackMasses, int type = 0);
  
  int                   fType, fVerbose; 
  Vertex                fPV;   
  const TransientTrackBuilder *fpTTB;   
};
