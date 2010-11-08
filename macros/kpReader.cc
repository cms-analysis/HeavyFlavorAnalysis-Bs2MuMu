#include "kpReader.hh"

#include <utility>

using namespace std;

const static float kMassBplus = 5.2792;

kpReader::kpReader(TChain *tree, TString evtClassName) :
	massReader(tree,evtClassName),
	fCutTrackQual_mu1(0),
	fCutTrackQual_mu2(0),
	fCutTrackQual_kp(0),
	fCutOppSign_mu(false),
	total_counter(0),
	reco_counter(0),
	reco_single(0),
	reco_double(0)
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
	cout << "\t\treco single: " << reco_single << endl;
	cout << "\t\treco double: " << reco_double << endl;
	cout << "\ttotal counter: " << total_counter << endl;
} // ~kpReader()

void kpReader::eventProcessing()
{
	unsigned j,nc;
	TGenCand *pGen;
	multiset<int> current_decay;
	map<int,int>::const_iterator it;
	
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
	massReader::eventProcessing(); // here the decay_indices are modified in loadCandidateVariables()
	
	reco_counter += decay_indices.size();
	
	for (it = decay_indices.begin(); it != decay_indices.end(); ++it) {
		
		if(it->second == 1) reco_single++;
		else if(it->second == 2) reco_double++;
		else
			cout << "Reconstructed candidate with " << it->second << " muons?" << endl;

	}
} // eventProcessing()

int kpReader::loadCandidateVariables(TAnaCand *pCand)
{
	int type, first_mu = 1;
	TAnaTrack *sigTrack,*recTrack;
	TAnaCand *jpsiCand;
	int result,j;
	
	// default initialization
	fMassJPsi = -1.0;
	fDeltaR = -1.0;
	fPlabMu1 = TVector3();
	fPlabMu2 = TVector3();
	fPlabKp = TVector3();
	fCtau = 0.0;
	
	if (pCand->fType % 1000 != 521) return 0;
	
	result = massReader::loadCandidateVariables(pCand);
	
	fCtau = kMassBplus / pCand->fPlab.Mag() * fD3;
	
	// set the momenta
	for (j = pCand->fSig1; j <= pCand->fSig2 && j >= 0; j++) {
		
		sigTrack = fpEvt->getSigTrack(j);
		type = abs(sigTrack->fMCID);
		recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
		if (type == 13) {
			// muon
			if (first_mu) {
				fPlabMu1 = sigTrack->fPlab;
				fTrackQual_mu1 = recTrack->fTrackQuality;
				fQ_mu1 = recTrack->fQ;
			} else {
				fPlabMu2 = sigTrack->fPlab;
				fTrackQual_mu2 = recTrack->fTrackQuality;
				fQ_mu2 = recTrack->fQ;
			}
			first_mu = 0;
		} else if (type == 321) {
			fPlabKp = sigTrack->fPlab;
			fQ_kp = recTrack->fQ;
			fTrackQual_kp = recTrack->fTrackQuality;
		}
	}
	
	// muon1 ist usually that one with the bigger pt
	if (fPlabMu1.Perp() < fPlabMu2.Perp()) {
		swap(fPlabMu1,fPlabMu2);
		swap(fTrackQual_mu1,fTrackQual_mu2);
		swap(fQ_mu1,fQ_mu2);
	}
	
	// set the jpsi mass
	for (j = pCand->fDau1; j <= pCand->fDau2 && j >= 0; j++) {
		
		jpsiCand = fpEvt->getCand(j);
		if (jpsiCand->fType % 1000 == 443) {
			fMassJPsi = jpsiCand->fMass;
			fD3_BpJpsi = jpsiCand->fVtx.fD3d;
			fD3e_BpJpsi = jpsiCand->fVtx.fD3dE;
			break;
		}
	}

	// set the deltaR of the J/Psi w.r.t. Kp
	if (fPlabKp.Perp() > 0 && (fPlabMu1 + fPlabMu2).Perp() > 0)
	  fDeltaR = (float)fPlabKp.DeltaR(fPlabMu1 + fPlabMu2);
	
	fPtMu1 = fPlabMu1.Perp();
	fPtMu2 = fPlabMu2.Perp();
	fPtKp = fPlabKp.Perp();

	return result;
} // loadCandidateVariables()

void kpReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("mass_jpsi",&fMassJPsi,"mass_jpsi/F");
	reduced_tree->Branch("deltaR",&fDeltaR,"deltaR/F");
	reduced_tree->Branch("plab_mu1","TVector3",&fPlabMu1Ptr);
	reduced_tree->Branch("plab_mu2","TVector3",&fPlabMu2Ptr);
	reduced_tree->Branch("plab_kp","TVector3",&fPlabKpPtr);
	reduced_tree->Branch("pt_mu1",&fPtMu1,"pt_mu1/F");
	reduced_tree->Branch("pt_mu2",&fPtMu2,"pt_mu2/F");
	reduced_tree->Branch("pt_kp",&fPtKp,"pt_kp/F");
	reduced_tree->Branch("track_qual_mu1",&fTrackQual_mu1,"track_qual_mu1/I");
	reduced_tree->Branch("track_qual_mu2",&fTrackQual_mu2,"track_qual_mu2/I");
	reduced_tree->Branch("track_qual_kp",&fTrackQual_kp,"track_qual_kp/I");
	reduced_tree->Branch("q_mu1",&fQ_mu1,"q_mu1/I");
	reduced_tree->Branch("q_mu2",&fQ_mu2,"q_mu2/I");
	reduced_tree->Branch("q_kp",&fQ_kp,"q_kp/I");
	reduced_tree->Branch("d3_bp_to_jpsi",&fD3_BpJpsi,"d3_bp_to_jpsi/F");
	reduced_tree->Branch("d3_bp_to_jpsi_e",&fD3e_BpJpsi,"d3_bp_to_jpsi_e/F");
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
	
	if (result && fCandidate == 300521) // save the decay
		decay_indices[gen->fNumber] = max((int)fNbrMuons, decay_indices[gen->fNumber]);
bail:
	return result;
} // checkTruth()

bool kpReader::parseCut(char *cutName, float cutValue, int dump)
{
	bool parsed;
	
	// parse the options..
	parsed = (strcmp(cutName,"TRACK_QUAL_MU1") == 0);
	if (parsed) {
		fCutTrackQual_mu1 = (int)cutValue;
		if (dump) cout << "TRACK_QUAL_MU1: " << fCutTrackQual_mu1 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRACK_QUAL_MU2") == 0);
	if (parsed) {
		fCutTrackQual_mu2 = (int)cutValue;
		if (dump) cout << "TRACK_QUAL_MU2: " << fCutTrackQual_mu2 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRACK_QUAL_KP") == 0);
	if (parsed) {
		fCutTrackQual_kp = (int)cutValue;
		if (dump) cout << "TRACK_QUAL_KP: " << fCutTrackQual_kp << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"OPPOSITE_SIGN_MU") == 0);
	if (parsed) {
		fCutOppSign_mu = (cutValue != 0.0);
		if (dump) cout << "OPPOSITE_SIGN_MU: " << fCutOppSign_mu << endl;
		goto bail;
	}
	
	// nothing parsed, then maybe massReader knows about it.
	parsed = massReader::parseCut(cutName, cutValue, dump);
	
bail:
	return parsed;
} // parseCut()

bool kpReader::applyCut()
{
	bool pass = true;
	
	// check track quality of mu1
	if (fCutTrackQual_mu1) {
		pass = (fTrackQual_mu1 & fCutTrackQual_mu1) != 0;
		if (!pass) goto bail;
	}
	
	// check track quality of mu2
	if (fCutTrackQual_mu2) {
		pass = (fTrackQual_mu2 & fCutTrackQual_mu2) != 0;
		if (!pass) goto bail;
	}
	
	// check track quality of kp
	if (fCutTrackQual_kp) {
		pass = (fTrackQual_kp & fCutTrackQual_kp) != 0;
		if (!pass) goto bail;
	}
	
	// check opposite sign of muons
	if (fCutOppSign_mu) {
		pass = fQ_mu1 != fQ_mu2;
		if (!pass) goto bail;
	}
	
	// check superclass cuts
	pass = massReader::applyCut();
	
bail:
	// return result
	return pass;
} // applyCut()
