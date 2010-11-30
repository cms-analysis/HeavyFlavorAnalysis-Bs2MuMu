#include "kpIsolationReader.hh"

#include <TH2D.h>

kpIsolationReader::kpIsolationReader(TChain *tree, TString evtClassName) : kpReader(tree,evtClassName),fCutNearestTrackPt(0.0)
{
	std::cout << "kpIsolationReader instantiated..." << std::endl;
}

kpIsolationReader::~kpIsolationReader() {}

void kpIsolationReader::bookHist()
{
	kpReader::bookHist(); // superclass
	
	reduced_tree->Branch("tracks_ix",fTracksIx,Form("tracks_ix[%d]/I",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ip",fTracksIP,Form("tracks_ip[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ipe",fTracksIPE,Form("tracks_ipe[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_pt",fTracksPT,Form("tracks_pt[%d]/F",NBR_TRACKS_STORE));
	reduced_tree->Branch("tracks_ptrel",fTracksPTRel,Form("tracks_ptrel[%d]/F",NBR_TRACKS_STORE));
} // bookHist()

int kpIsolationReader::loadCandidateVariables(TAnaCand *pCand)
{
	float ip,ipE,pt,ptrel;
	TAnaTrack *pTrack;
	TVector3 uVector;
	unsigned j,k;
	int result = kpReader::loadCandidateVariables(pCand); // load all the variables...
	if(!result) goto bail;
	
	// clean entries
	for (unsigned j = 0; j < NBR_TRACKS_STORE; j++) fTracksIx[j] = -1;
	memset(fTracksIP,0,sizeof(fTracksIP));
	memset(fTracksIPE,0,sizeof(fTracksIPE));
	memset(fTracksPT,0,sizeof(fTracksPT));
	memset(fTracksPTRel,0,sizeof(fTracksPTRel));
	
	// doesn't match the cuts
	if(!applyCut()) {
		result = 0;
		goto bail;
	}
	
	// fill the histogramms...
	for (j = 0, k = 0; j < pCand->fNstTracks.size(); j++) {
		
		pTrack = fpEvt->getRecTrack(pCand->fNstTracks[j].first);
		uVector = pCand->fPlab.Unit();
		pt = pTrack->fPlab.Perp();
		ip = pCand->fNstTracks[j].second.first;
		ipE = pCand->fNstTracks[j].second.second;
		ptrel = (pTrack->fPlab - (pTrack->fPlab * uVector) * uVector).Mag();
		
		if (pt > fCutNearestTrackPt && k < NBR_TRACKS_STORE) {
			fTracksIx[k] = pCand->fNstTracks[j].first;
			fTracksIP[k] = ip;
			fTracksIPE[k] = ipE;
			fTracksPT[k] = pt;
			fTracksPTRel[k] = ptrel;
			k++;
		}
	}
	
bail:
	return result;
} // loadCandidateVariables()

bool kpIsolationReader::parseCut(char *cutName, float cutLow, float cutHigh, int dump)
{
	using std::cout;
	using std::endl;
	bool parsed;
	
	parsed = (strcmp(cutName,"NEAR_PT") == 0);
	if (parsed) {
		fCutNearestTrackPt = cutLow;
		if (dump) cout << "NEAR_PT: " << fCutNearestTrackPt << endl;
		goto bail;
	}
	
	// call the super routine
	parsed = kpReader::parseCut(cutName, cutLow, cutHigh, dump);
	
bail:
	return parsed;
} // parseCut()
