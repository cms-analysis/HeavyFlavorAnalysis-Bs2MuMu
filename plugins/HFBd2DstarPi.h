#ifndef _HFBD2DSTARPI_h_
#define _HFBD2DSTARPI_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"


// ----------------------------------------------------------------------

class HFBd2DstarPi : public edm::EDAnalyzer {
public:
	explicit HFBd2DstarPi(const edm::ParameterSet&);
	~HFBd2DstarPi();
	
private:
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	int fVerbose;
	edm::InputTag   fTracksLabel, fPrimaryVertexLabel;
	
	double		  fTrackPt, fSlowPionPt;
	double        fD0Window, fDeltaM;
	double        fDeltaR, fMaxDoca, fMaxD0, fMaxDz;
	int           fType;
	
	edm::ESHandle<TransientTrackBuilder> fTTB;
};

#endif
