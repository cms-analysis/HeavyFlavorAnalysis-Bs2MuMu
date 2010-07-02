#include "kpReader.hh"

#include <utility>

using namespace std;

kpReader::kpReader(TChain *tree, TString evtClassName) : massReader(tree,evtClassName),total_counter(0),reco_counter(0)
{
	// set the pointers to save in the tree
	fPlabMu1Ptr = &fPlabMu1;
	fPlabMu2Ptr = &fPlabMu2;
	fPlabKpPtr = &fPlabKp;
	
	fTreeName = "kpReader reduced tree";
	
	// the true decay for truth matching...
	trueDecay.insert(13);	// mu
	trueDecay.insert(13);	// mu
	trueDecay.insert(321);	// Kp
	trueDecay.insert(443);	// J/Psi
	trueDecay.insert(521);	// B+
	
	cout << "kpReader instantiated..." << endl;
} // kpReader()

kpReader::~kpReader()
{
	// just dump some stuff.
	cout << "~kpReader:" << endl;
	cout << "\treco counter: " << reco_counter << endl;
	cout << "\ttotal counter: " << total_counter << endl;
} // ~kpReader()

void kpReader::eventProcessing()
{
	unsigned j,nc;
	TGenCand *pGen;
	multiset<int> current_decay;
	
	decay_indices.clear();
	
	// count the number of real decay using the generator block
	nc = fpEvt->nGenCands();
	for (j = 0; j < nc; j++) {
		pGen = fpEvt->getGenCand(j);
		if (abs(pGen->fID) == 521) {
			current_decay.clear();
			buildDecay(pGen,&current_decay);
			current_decay.erase(22); // remove Bremsstrahlung
			if (trueDecay == current_decay)
				total_counter++;
		}
	}
	
	// do the normal stuff
	massReader::eventProcessing(); // here the decay indices are modified in loadCandidateVariables()
	
	reco_counter += decay_indices.size();
	
} // eventProcessing()

int kpReader::loadCandidateVariables(TAnaCand *pCand)
{
	int type, first_mu = 1;
	TAnaTrack *track;
	TAnaCand *jpsiCand;
	int result,j;
	
	// default initialization
	fMassJPsi = -1.0;
	fPlabMu1 = TVector3();
	fPlabMu2 = TVector3();
	fPlabKp = TVector3();
	
	if (pCand->fType % 1000 != 521) return 0;
	
	result = massReader::loadCandidateVariables(pCand);
	
	// set the momenta
	for (j = pCand->fSig1; j <= pCand->fSig2 && j >= 0; j++) {
		
		track = fpEvt->getSigTrack(j);
		type = abs(track->fMCID);
		track = fpEvt->getRecTrack(track->fIndex);
		if (type == 13) {
			// muon
			if (first_mu)	fPlabMu1 = track->fPlab;
			else			fPlabMu2 = track->fPlab;
			first_mu = 0;
		} else if (type == 321)
			fPlabKp = track->fPlab;
	}
	
	// muon1 ist usually that one with the bigger pt
	if (fPlabMu1.Perp() < fPlabMu2.Perp())
		swap(fPlabMu1,fPlabMu2);
	
	// set the jpsi mass
	for (j = pCand->fDau1; j <= pCand->fDau2 && j >= 0; j++) {
		
		jpsiCand = fpEvt->getCand(j);
		if (jpsiCand->fType % 1000 == 443) {
			fMassJPsi = jpsiCand->fMass;
			break;
		}
	}
	
	return result;
} // loadCandidateVariables()

void kpReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("mass_jpsi",&fMassJPsi,"mass_jpsi/F");
	reduced_tree->Branch("plab_mu1","TVector3",&fPlabMu1Ptr);
	reduced_tree->Branch("plab_mu2","TVector3",&fPlabMu2Ptr);
	reduced_tree->Branch("plab_kp","TVector3",&fPlabKpPtr);
} // bookHist()

int kpReader::checkTruth(TAnaCand *pCand)
{
	int result;
	multiset<int> particles;
	TAnaTrack *track;
	TGenCand *gen;
	
	result = massReader::checkTruth(pCand);
	if (!result) goto bail;
	
	// check if they come from the right decay channel!
	track = fpEvt->getSigTrack(pCand->fSig1);
	track = fpEvt->getRecTrack(track->fIndex);
		
	gen = fpEvt->getGenCand(track->fGenIndex);
	while (abs(gen->fID) != 521) // works because else, massReader::checkTruth fails
		gen = fpEvt->getGenCand(gen->fMom1);
	
	buildDecay(gen,&particles);
	particles.erase(22); // remove Bremsstrahlung
	
	result = (particles == trueDecay);
	
	if (result) // save the decay
		decay_indices.insert(gen->fNumber);
bail:
	return result;
} // checkTruth()
