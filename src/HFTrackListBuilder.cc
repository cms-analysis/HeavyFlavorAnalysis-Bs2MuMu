/*
 *  HFTrackListBuilder.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 10.11.12.
 *
 */

#include <algorithm>

#include "HeavyFlavorAnalysis/Bs2MuMu/interface/HFTrackListBuilder.hh"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

using namespace std;

HFTrackListBuilder::HFTrackListBuilder(edm::Handle<edm::View<reco::Track> > &hTracks, edm::Handle<reco::MuonCollection> &hMuons, int verbose) :
	fhTracks(hTracks),
	fhMuons(hMuons),
	fVerbose(verbose),
	fCallerName("HFTrackListBuilder"),
	fMaxD0(999.),
	fMaxDz(999.),
	fMinPt(-1.)
{} // HFTrackListBuilder()

std::vector<int> HFTrackListBuilder::getMuonList()
{
	reco::MuonCollection::const_iterator muonIt;
	vector<int> trackList;
	vector<int>::iterator trackIt;
	
	trackList.reserve(20); // 20 muons should be enough
	
	for (muonIt = fhMuons->begin(); muonIt != fhMuons->end(); ++muonIt) {
		int ixMu = muonIt->track().index();
		if (ixMu >= 0) trackList.push_back(ixMu);
	}
	
	if (fVerbose > 0) {
		cout << "==>" << fCallerName << "> nMuons = " << fhMuons->size() << endl;
		cout << "==>" << fCallerName << "> nMuonIndices = " << trackList.size() << endl;
	}
	
	trackIt = std::remove_if(trackList.begin(), trackList.end(), *this);
	trackList.erase(trackIt,trackList.end());
	
	return trackList;
} // getMuonList()

std::vector<int> HFTrackListBuilder::getTrackList()
{
	vector<int> trackList; // allocate capacity
	int ix;
	
	trackList.reserve(300);
	for (ix = 0; (unsigned)ix < fhTracks->size(); ix++) {
		if ( !(*this)(ix))
			trackList.push_back(ix);
	}
	
	return trackList;
} // getTrackList()

bool HFTrackListBuilder::operator()(int ix)
{
	reco::TrackBaseRef rTrackView(fhTracks,ix);
	reco::Track tTrack(*rTrackView);
	bool result = tTrack.d0() > fMaxD0 || tTrack.dz() > fMaxDz || tTrack.pt() < fMinPt;
	
	return result;
} // operator(int ix);
