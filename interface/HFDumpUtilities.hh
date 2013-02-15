#ifndef HFDUMPUTILITIES_H
#define HFDUMPUTILITIES_H

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TSimpleTrack.hh"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

void fillSimpleTrack(TSimpleTrack *pTrack, const reco::Track &trackView, 
		     int tidx, int mid, int gidx, const reco::VertexCollection *vc);

void fillAnaTrack(TAnaTrack *pTrack, const reco::Track &trackView, int tidx, int gidx,
		  const reco::VertexCollection *vc, const reco::MuonCollection *mc, const reco::BeamSpot *bs);

int getPv(int tidx, const reco::VertexCollection *vc);

int muonID(const reco::Muon &rm);

//void cleanupTruthMatching();
void cleanupTruthMatching(edm::Handle<edm::View<reco::Track> > &hTracks, edm::ESHandle<MagneticField> &magfield);

#endif
