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
		virtual bool parseCut(char *cutName, float cutLow, float cutHigh, int dump = 1);
	
	protected:
		int fTracksIx[NBR_TRACKS_STORE];
		float fTracksIP[NBR_TRACKS_STORE];
		float fTracksIPE[NBR_TRACKS_STORE];
		float fTracksPT[NBR_TRACKS_STORE];
		float fTracksPTRel[NBR_TRACKS_STORE];
	
	private:
		// Additional cut variables
		float fCutNearestTrackPt;
};

#endif
