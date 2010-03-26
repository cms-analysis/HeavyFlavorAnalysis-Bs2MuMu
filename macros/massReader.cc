#include "massReader.hh"
#include <cstdlib>

#define require_true(COND,LABEL) if( !(COND) ) goto LABEL

// test a <= x < b
static inline int in_interval(int x, int a, int b)
{
	return (a <= x) && (x < b);
} // in_interval()

massReader::massReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName),reduced_tree(NULL)
{
  fMomentumPtr = &fMomentum;
} // massReader()

massReader::~massReader()
{ } // ~massReader()

void massReader::eventProcessing()
{
	int j,nc;
	TAnaCand *pCand;
	
	// Fill a reduced tree with the following types
	//	- momentum
	//	- mass
	//	- type
	nc = fpEvt->nCands();
	for (j=0; j<nc; j++) {
		pCand = fpEvt->getCand(j);
		
		// Save in the tree
		fCandidate = pCand->fType;
		fMomentum = pCand->fPlab;
		fMass = pCand->fMass;
		fTruth = checkTruth(pCand, abs(pCand->fType) % 1000);
		reduced_tree->Fill();
	}
} // massReader::eventProcessing()

void massReader::bookHist()
{
	// create the tree
	reduced_tree = new TTree("T","Candidate Mass / Eta / Charge");
	
	// and add the branches
	reduced_tree->Branch("candidate",&fCandidate,"candidate/I");
	reduced_tree->Branch("p","TVector3",&fMomentumPtr);
	reduced_tree->Branch("mass",&fMass,"mass/D");
	reduced_tree->Branch("truth",&fTruth,"truth/I");
} // massReader::bookHist()

int massReader::checkTruth(TAnaCand *cand, int truth_type)
{
	TAnaTrack *sgTrack;
	TAnaTrack *recTrack;
	TGenCand *truthParticle;
	TGenCand *trackParticle;
	int j,succes = 0;
	int nSigs,nRecs,nGens;
	
	
	nSigs = fpEvt->nSigTracks();
	nRecs = fpEvt->nRecTracks();
	nGens = fpEvt->nGenCands();
	
	// get the first track and get the originating particle of type 'truth_type' => truthParticle
	require_true(in_interval(cand->fSig1, 0, nSigs),bail);
	sgTrack = fpEvt->getSigTrack(cand->fSig1);
	
	require_true(in_interval(sgTrack->fIndex, 0, nRecs),bail);
	recTrack = fpEvt->getRecTrack(sgTrack->fIndex);
	
	require_true(in_interval(recTrack->fGenIndex, 0, nGens),bail);
	truthParticle = fpEvt->getGenCand(recTrack->fGenIndex);
	
	while (abs(truthParticle->fID) != truth_type && in_interval(truthParticle->fMom1, 0, nGens))
		truthParticle = fpEvt->getGenCand(truthParticle->fMom1);

	// check if our original is a valid.
	require_true(abs(truthParticle->fID) == truth_type,bail);
	
	// reconstruct the other tracks
	for (j = cand->fSig1+1; j <=cand->fSig2; j++) {
		
		require_true(in_interval(j, 0, nSigs),bail);
		sgTrack = fpEvt->getSigTrack(j);
		
		require_true(in_interval(sgTrack->fIndex, 0, nRecs),bail);
		recTrack = fpEvt->getRecTrack(sgTrack->fIndex);
		
		require_true(in_interval(recTrack->fGenIndex, 0, nGens),bail);
		trackParticle = fpEvt->getGenCand(recTrack->fGenIndex);
		
		while (abs(trackParticle->fID)!=truth_type && in_interval(trackParticle->fMom1, 0, nGens))
			trackParticle = fpEvt->getGenCand(trackParticle->fMom1);
		
		// check the particle
		require_true(trackParticle->fNumber==truthParticle->fNumber,bail);
	}
	
	// still here? then every track originated from the right particle
	succes = 1;
	
bail:
	return succes;
} // massReader::checkTruth()
