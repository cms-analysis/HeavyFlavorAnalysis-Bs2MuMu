#ifndef HFDUMPUTILITIES_H
#define HFDUMPUTILITIES_H

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


void fillSimpleTrack(TSimpleTrack *pTrack, const reco::Track &trackView, 
		     int tidx, int mid, int gidx, const reco::VertexCollection *vc);

void fillAnaTrack(TAnaTrack *pTrack, const reco::Track &trackView, int tidx, int gidx,
		  const reco::VertexCollection *vc, const reco::MuonCollection *mc, const reco::BeamSpot *bs);

int getPv(int tidx, const reco::VertexCollection *vc);

int muonID(const reco::Muon &rm);

#endif
