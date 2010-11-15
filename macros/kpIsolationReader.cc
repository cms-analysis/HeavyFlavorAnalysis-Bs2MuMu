#include "kpIsolationReader.hh"

#include <TH2D.h>

kpIsolationReader::kpIsolationReader(TChain *tree, TString evtClassName) : kpReader(tree,evtClassName),pt_ix(NULL),ip_ix(NULL),pt_ip(NULL),ptrel_ip(NULL)
{
	std::cout << "kpIsolationReader instantiated..." << std::endl;
}

kpIsolationReader::~kpIsolationReader() {}

void kpIsolationReader::bookHist()
{
	kpReader::bookHist(); // superclass
	
	pt_ix = new TH2D("pt_ix","pt vs ix",100,0.0,100.0,100,0.0,8.0); // first 100 indices and pt up to 8.0 GeV
	ip_ix = new TH2D("ip_ix","ip vs ix",100,0.0,100.0,100,0.0,2.0); // first 100 indices and ip up to 2 cm
	pt_ip = new TH2D("pt_ip","pt vs ip",100,0.0,2.0,100,0.0,8.0); // ip up to 2 cm and pt up to 8.0 GeV
	ptrel_ip = new TH2D("ptrel_ip","ptrel vs ip",100,0.0,2.0,100,0.0,10.0); // ip up to 2 cm and ptrel up to 10 GeV
} // bookHist()

int kpIsolationReader::loadCandidateVariables(TAnaCand *pCand)
{
	double ip,pt;
	TAnaTrack *pTrack;
	TVector3 uVector;
	int result = kpReader::loadCandidateVariables(pCand); // load all the variables...
	if(!result) goto bail;
	
	// doesn't match the cuts
	if(!applyCut()) {
		result = 0;
		goto bail;
	}
	
	// fill the histogramms...
	for (unsigned j = 0; j < pCand->fNstTracks.size(); j++) {
		
		pTrack = fpEvt->getRecTrack(pCand->fNstTracks[j].first);
		pt = pTrack->fPlab.Perp();
		ip = pCand->fNstTracks[j].second.first;
		
		pt_ix->Fill(j,pt);
		ip_ix->Fill(j,ip);
		pt_ip->Fill(ip,pt);
		
		uVector = pCand->fPlab.Unit();
		ptrel_ip->Fill(ip, (pTrack->fPlab - (pTrack->fPlab * uVector) * uVector).Mag() );
	}
	
bail:
	return result;
} // loadCandidateVariables()
