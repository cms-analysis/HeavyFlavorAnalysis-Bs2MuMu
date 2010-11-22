#include "kpIsolationReader.hh"

#include <TH2D.h>

kpIsolationReader::kpIsolationReader(TChain *tree, TString evtClassName) : kpReader(tree,evtClassName)
{
	std::cout << "kpIsolationReader instantiated..." << std::endl;
}

kpIsolationReader::~kpIsolationReader() {}

void kpIsolationReader::bookHist()
{
	kpReader::bookHist(); // superclass
	
	reduced_tree->Branch("tracks_ix",fTracksIx,Form("tracks_ix[%d]/I",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ip",fTracksIP,Form("tracks_ip[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_pt",fTracksPT,Form("tracks_pt[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ptrel",fTracksPTRel,Form("tracks_ptrel[%d]/F",NBR_TRACKS_STORE));
} // bookHist()

int kpIsolationReader::loadCandidateVariables(TAnaCand *pCand)
{
	double ip,pt,ptrel;
	TAnaTrack *pTrack;
	TVector3 uVector;
	int result = kpReader::loadCandidateVariables(pCand); // load all the variables...
	if(!result) goto bail;
	
	// clean entries
	for (unsigned j = 0; j < NBR_TRACKS_STORE; j++) fTracksIx[j] = -1;
	memset(fTracksIP,0,sizeof(fTracksIP));
	memset(fTracksPT,0,sizeof(fTracksPT));
	memset(fTracksPTRel,0,sizeof(fTracksPTRel));
	
	// doesn't match the cuts
	if(!applyCut()) {
		result = 0;
		goto bail;
	}
	
	// fill the histogramms...
	for (unsigned j = 0; j < pCand->fNstTracks.size(); j++) {
		
		pTrack = fpEvt->getRecTrack(pCand->fNstTracks[j].first);
		uVector = pCand->fPlab.Unit();
		pt = pTrack->fPlab.Perp();
		ip = pCand->fNstTracks[j].second.first;
		ptrel = (pTrack->fPlab - (pTrack->fPlab * uVector) * uVector).Mag();

		if (j < NBR_TRACKS_STORE) {
			fTracksIx[j] = pCand->fNstTracks[j].first;
			fTracksIP[j] = ip;
			fTracksPT[j] = pt;
			fTracksPTRel[j] = ptrel;
		}
	}
	
bail:
	return result;
} // loadCandidateVariables()
