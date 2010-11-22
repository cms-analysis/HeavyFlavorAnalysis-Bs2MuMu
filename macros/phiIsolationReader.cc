#include "phiIsolationReader.hh"

phiIsolationReader::phiIsolationReader(TChain *tree, TString evtClassName) : phiReader(tree,evtClassName)
{} // phiIsolationReader()

phiIsolationReader::~phiIsolationReader() {}

void phiIsolationReader::bookHist()
{
	phiReader::bookHist();
		
	reduced_tree->Branch("tracks_ip",fTracksIP,Form("tracks_ip[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_pt",fTracksPT,Form("tracks_pt[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ptrel",fTracksPTRel,Form("tracks_ptrel[%d]/F",NBR_TRACKS_STORE));
} // bookHist()

int phiIsolationReader::loadCandidateVariables(TAnaCand *pCand)
{
	double ip,pt,ptrel;
	TAnaTrack *pTrack;
	TVector3 uVector;
	int result = phiReader::loadCandidateVariables(pCand); // load all the variables...
	if(!result) goto bail;
	
	// clean entries
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
			fTracksIP[j] = ip;
			fTracksPT[j] = pt;
			fTracksPTRel[j] = ptrel;
		}
	}
	
bail:
	return result;
} // loadCandidateVariables()
