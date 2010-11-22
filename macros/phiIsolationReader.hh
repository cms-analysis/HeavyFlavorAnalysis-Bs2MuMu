#ifndef PHIISOLATIONREADER_H
#define PHIISOLATIONREADER_H

#include "phiReader.hh"

#define NBR_TRACKS_STORE 20

class phiIsolationReader : public phiReader {
	
	public:
		phiIsolationReader(TChain *tree, TString evtClassName);
		virtual ~phiIsolationReader();
		
		virtual void bookHist();
	
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
	
	protected:
		float fTracksIP[NBR_TRACKS_STORE];
		float fTracksPT[NBR_TRACKS_STORE];
		float fTracksPTRel[NBR_TRACKS_STORE];
};

#endif
