#ifndef KPISOLATIONREADER_H
#define KPISOLATIONREADER_H

#include "kpReader.hh"

#define NBR_TRACKS_STORE 20

class kpIsolationReader : public kpReader {
	
	public:
		kpIsolationReader(TChain *tree, TString evtClassName);
		virtual ~kpIsolationReader();
		
		virtual void bookHist();
	
	protected:
		virtual int loadCandidateVariables(TAnaCand *pCand);
	
	protected:
		int fTracksIx[NBR_TRACKS_STORE];
		float fTracksIP[NBR_TRACKS_STORE];
		float fTracksPT[NBR_TRACKS_STORE];
		float fTracksPTRel[NBR_TRACKS_STORE];
};

#endif
