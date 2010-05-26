#include "ksReader.hh"

using namespace std;

ksReader::ksReader(TChain *tree, TString evtClassName) :
	massReader(tree,evtClassName),unrecoverableDecays(0),recoverableCounter(0)
{
	fPlabMu1Ptr = &fPlabMu1;
	fPlabMu2Ptr = &fPlabMu2;
	fPlabPi1Ptr = &fPlabPi1;
	fPlabPi2Ptr = &fPlabPi2;
	
	// build the true decay
	trueDecay.insert(13);
	trueDecay.insert(13);
	trueDecay.insert(211);
	trueDecay.insert(211);
	trueDecay.insert(310);
	trueDecay.insert(443);
	trueDecay.insert(511);
	
	cout << "ksReader instantiated..." << endl;
} // ksReader()

ksReader::~ksReader() {
	cout << "Possible recoverable decays: " << recoverableCounter << endl;
	cout << "Unrecovered Decays: " << unrecoverableDecays << endl;
}

int ksReader::loadCandidateVariables(TAnaCand *pCand)
{
	TAnaCand *ksCand = NULL;
	TAnaTrack *track;
	int real_type;
	int algo_type;
	int result;
	int j;
	int firstMu = 1;
	int firstPi = 1;
	int type;
	
	// debugging
	static int counter = 0;
	
	// default initialization
	fMassJPsi = -1.0;
	fMassKs = -1.0;
	fAlphaKs = -1.0;
	fDxyKs = -1.0;
	fDxyeKs = -1.0;
	fChi2Ks = -1.0;

	fPlabMu1 = TVector3();
	fPlabMu2 = TVector3();
	fPlabPi1 = TVector3();
	fPlabPi2 = TVector3();
	
	real_type = pCand->fType % 1000;
	algo_type = pCand->fType / 1000;
	
	result = massReader::loadCandidateVariables(pCand);
		
	if(real_type != 511) goto bail;
	
	// tracks
	for (j = pCand->fSig1; j <= pCand->fSig2; j++) {
		
		track = fpEvt->getSigTrack(j);
		type = abs(track->fMCID);
		track = fpEvt->getRecTrack(track->fIndex);
		if (type == 13) {
			// muon
			if(firstMu)	fPlabMu1 = track->fPlab;
			else		fPlabMu2 = track->fPlab;
			firstMu = 0;
			
		} else if (type == 211) {
			// pion
			if(firstPi)	fPlabPi1 = track->fPlab;
			else		fPlabPi2 = track->fPlab;
			firstPi = 0;
		}
	}
	
	
	// Ks subvariables
	for (j = pCand->fDau1; j <= pCand->fDau2 && j>= 0; j++) {
		
		TAnaCand *tmpCand = fpEvt->getCand(j);
		if(tmpCand->fType % 1000 == 310) { // ks
			ksCand = tmpCand;
			break;
		}
	}
	if (ksCand) {
		fMassKs = ksCand->fMass;
		fAlphaKs = ksCand->fPlab.Angle(ksCand->fVtx.fPoint - pCand->fVtx.fPoint);
		fDxyKs = ksCand->fVtx.fDxy;
		fDxyeKs = ksCand->fVtx.fDxyE;
		fChi2Ks = ksCand->fVtx.fChi2;
	}
	
	// jpsi subvariables
	ksCand = NULL;
	for (j = pCand->fDau1; j <= pCand->fDau2; j++) {
		
		TAnaCand *tmpCand = fpEvt->getCand(j);
		if (tmpCand->fType % 1000 == 443) { // J/Psi
			ksCand = tmpCand;
			break;
		}
	}
	if (ksCand) {
		fMassJPsi = ksCand->fMass;
	}
	
bail:
	counter++;
	return result;
} // ksReader()

void ksReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("mass_jpsi",&fMassJPsi,"mass_jpsi/D");
	reduced_tree->Branch("alpha_ks",&fAlphaKs,"alpha_ks/D");
	reduced_tree->Branch("mass_ks",&fMassKs,"mass_ks/D");
	reduced_tree->Branch("dxy_ks",&fDxyKs,"dxy_ks/D");
	reduced_tree->Branch("dxye_ks",&fDxyeKs,"dxye_ks/D");
	reduced_tree->Branch("chi2_ks",&fChi2Ks,"chi2_ks/D");
	reduced_tree->Branch("plab_mu1","TVector3",&fPlabMu1Ptr);
	reduced_tree->Branch("plab_mu2","TVector3",&fPlabMu2Ptr);
	reduced_tree->Branch("plab_pi1","TVector3",&fPlabPi1Ptr);
	reduced_tree->Branch("plab_pi2","TVector3",&fPlabPi2Ptr);
} // bookHist()

// a more sophisticated truth matching for the B0 particles
int ksReader::checkTruth(TAnaCand *cand, int truth_type)
{
	int result = massReader::checkTruth(cand, truth_type);
	multiset<int> particles;
	TAnaTrack *track;
	TGenCand *gen;
	
	truth_type = abs(truth_type);
	
	// check more deeply the B0 candidates
	if (result && truth_type == 511) {
		
		// check they come from the right decay channel!
		track = fpEvt->getSigTrack(cand->fSig1);
		track = fpEvt->getRecTrack(track->fIndex);
		
		// get the index of the B0 GenCandidate
		gen = fpEvt->getGenCand(track->fGenIndex);
		while (abs(gen->fID) != truth_type) // this works as result == 1!
			gen = fpEvt->getGenCand(gen->fMom1);
		
		buildDecay(gen->fNumber,&particles);
		particles.erase(22); // remove Bremsstahlung
		
		result = (particles == trueDecay);
	}
	
	return result;
} // checkTruth()

// return the number of tracks and identified muons
unsigned ksReader::buildDecay(int genIx, multiset<int> *particles, unsigned *nbrMuons)
{
	TGenCand *pGen;
	TAnaTrack *pTrack;
	unsigned result = 0;
	int j,nc;
	
	pGen = fpEvt->getGenCand(genIx);
	particles->insert(abs(pGen->fID));
	
	for (j = pGen->fDau1; j <= pGen->fDau2 && j>=0; j++)
		result += buildDecay(j,particles,nbrMuons);
	
	// count the tracks!!!
	nc = fpEvt->nRecTracks();
	for (j = 0; j < nc; j++) {
		pTrack = fpEvt->getRecTrack(j);
		if (genIx == pTrack->fGenIndex) {
			result++;
			
			if (abs(pGen->fID) == 13 && pTrack->fMuID > 0 && (pTrack->fMuID & 6) && nbrMuons)
				(*nbrMuons)++;
			
			// we want at most one track added per trace
			break;
		}
	}
	
	return result;
} // ksReader::buildDecay()

void ksReader::eventProcessing()
{
	unsigned j,nc;
	int k;
	unsigned nbrMuons,nbrTracks;
	TGenCand *pGen;
	TAnaCand *pCand;
	multiset<int> particles;
	set<int> recoverableGens;
	set<int>::const_iterator it;
	
	massReader::eventProcessing(); // call the super routine to save reduced tree!
	
	// create a list of all true candidates!!
	nc = fpEvt->nGenCands();
	for (j = 0; j < nc; j++) {
		pGen = fpEvt->getGenCand(j);
		
		if (abs(pGen->fID) == 511) {
			
			// this is a B0 decay
			particles.clear();
			nbrMuons = 0;
			
			nbrTracks = buildDecay(j,&particles,&nbrMuons);
			
			// remove Bremsstrahlung
			particles.erase(22);
			if (particles == trueDecay && nbrTracks == 4) {
				
				if (nbrMuons > 0) {
					recoverableCounter++;
					recoverableGens.insert(j);
				}
			}
		}
	}
	
	// remove all the generators we actually have recovered!
	// we only check the 600XXX candidates!
	nc = fpEvt->nCands();
	for (j = 0; j < nc; j++) {
		pCand = fpEvt->getCand(j);
		
		if (pCand->fType != 600511)
			continue;
		
		k = getGenIndex(pCand,511);
		recoverableGens.erase(k);
	}
	
	// dump the stuff
	for (it = recoverableGens.begin(); it != recoverableGens.end(); ++it)
		unrecoverableDecays++;
} // eventProcessing()

int ksReader::getGenIndex(TAnaCand *pCand, int candID)
{
	TAnaTrack *pTrack;
	TGenCand *momCand,*pGen;;
	int j;

	// set the momCand out of the first signal track...
	pTrack = fpEvt->getSigTrack(pCand->fSig1);
	pTrack = fpEvt->getRecTrack(pTrack->fIndex);
	if (pTrack->fGenIndex < 0) return -1;
	momCand = fpEvt->getGenCand(pTrack->fGenIndex);

	while (abs(momCand->fID) != candID && momCand->fMom1 >= 0)
		momCand = fpEvt->getGenCand(momCand->fMom1);
		
	if (abs(momCand->fID) != candID) return -1;

	for (j = pCand->fSig1+1; j <= pCand->fSig2; j++) {
		
		pTrack = fpEvt->getSigTrack(j);
		pTrack = fpEvt->getRecTrack(pTrack->fIndex);
		
		if (pTrack->fGenIndex < 0) return -1;
		pGen = fpEvt->getGenCand(pTrack->fGenIndex);
		
		while(abs(pGen->fID) != candID && pGen->fMom1 >= 0)
			pGen = fpEvt->getGenCand(pGen->fMom1);
			
		if((abs(pGen->fID)!=candID) || (pGen->fNumber != momCand->fNumber)) return -1;
	}

	return momCand->fNumber;
} // getGenIndex()
