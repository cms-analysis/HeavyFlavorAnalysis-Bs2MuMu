#include "massReader.hh"
#include <cstdlib>

#define require_true(COND,LABEL) if( !(COND) ) goto LABEL

using namespace std;

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
	
	// Fill a reduced tree
	nc = fpEvt->nCands();
	for (j = 0; j<nc; j++) {
		
		if(loadCandidateVariables(fpEvt->getCand(j)))
			reduced_tree->Fill();
	}
} // massReader::eventProcessing()

int massReader::loadCandidateVariables(TAnaCand *pCand)
{
	TAnaCand *momCand;
	
	// Save in the tree
	fCandidate = pCand->fType;
	fMomentum = pCand->fPlab;

	fMass = pCand->fMass;
	fTruth = checkTruth(pCand, abs(pCand->fType) % 1000);
	fTwoMuon = checkMuons(pCand);
	fDxy = pCand->fVtx.fDxy;
	fDxyE = pCand->fVtx.fDxyE;
	fChi2 = pCand->fVtx.fChi2;
	fNdof = pCand->fVtx.fNdof;
	
	if (pCand->fMom >= 0) {
		momCand = fpEvt->getCand(pCand->fMom);
		fAlpha = pCand->fPlab.Angle(pCand->fVtx.fPoint - momCand->fVtx.fPoint);
	} else
		fAlpha = pCand->fPlab.Angle(pCand->fVtx.fPoint - fpEvt->bestPV()->fPoint);
	
	return 1;
} // loadCandidateVariables()

void massReader::bookHist()
{
	// create the tree
	reduced_tree = new TTree("T","Candidate Mass / Eta / Charge");
	
	// and add the branches
	reduced_tree->Branch("candidate",&fCandidate,"candidate/I");
	reduced_tree->Branch("p","TVector3",&fMomentumPtr);
	reduced_tree->Branch("mass",&fMass,"mass/D");
	reduced_tree->Branch("truth",&fTruth,"truth/I");
	reduced_tree->Branch("two_muon",&fTwoMuon,"two_muon/I");
	reduced_tree->Branch("dxy",&fDxy,"dxy/D");
	reduced_tree->Branch("dxye",&fDxyE,"dxye/D");
	reduced_tree->Branch("alpha",&fAlpha,"alpha/D");
	reduced_tree->Branch("chi2",&fChi2,"chi2/D");
	reduced_tree->Branch("Ndof",&fNdof,"Ndof/D");
} // massReader::bookHist()

int massReader::checkTruth(TAnaCand *cand, int truth_type)
{
	TAnaTrack *sgTrack;
	TAnaTrack *recTrack;
	TGenCand *truthParticle;
	TGenCand *trackParticle;
	int j,success = 0;
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
		require_true(trackParticle->fNumber == truthParticle->fNumber,bail);
	}
	
	// still here? then every track originated from the right particle
	success = 1;
	
bail:
	return success;
} // massReader::checkTruth()

int massReader::checkMuons(TAnaCand *cand)
{
	TAnaTrack *track;
	unsigned nbrMu = 0, totalMu = 0;
	int j;
	
	// return true if all muons come from the muon list
	
	// check for invalid candidates
	if (cand->fSig1 < 0)
		return 0;
	
	for (j = cand->fSig1; j <= cand->fSig2; j++) {
		
		track = fpEvt->getSigTrack(j); // signal track
		// check if this is supposed to be a muon!
		if (abs(track->fMCID) == 13) {
			totalMu++; // muon
			track = fpEvt->getRecTrack(track->fIndex); // rec track
			if (track->fMuID > 0 && (track->fMuID & 6))
				nbrMu++;
		}
	}
	
	return (nbrMu == totalMu);
} // checkMuons()
