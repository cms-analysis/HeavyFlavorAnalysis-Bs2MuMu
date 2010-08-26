#include "phiReader.hh"

#include <utility>
#include <TLorentzVector.h>

using namespace std;

phiReader::phiReader(TChain *tree, TString evtClassName) : massReader(tree,evtClassName), total_counter(0), reco_single(0), reco_double(0)
{
	cout << "Instantiating phiReader..." << endl;
	fPlabMu1Ptr = &fPlabMu1;
	fPlabMu2Ptr = &fPlabMu2;
	fPlabKp1Ptr = &fPlabKp1;
	fPlabKp2Ptr = &fPlabKp2;
	
	fTreeName = "phiReader reduced tree.";
	
	trueDecay.insert(13); // mu
	trueDecay.insert(13); // mu
	trueDecay.insert(321); // kp
	trueDecay.insert(321); // kp
	trueDecay.insert(443); // j/psi
	trueDecay.insert(333); // phi
	trueDecay.insert(531); // Bs
	
	cout << "phiReader instantiated..." << endl;
} // phiReader()

phiReader::~phiReader()
{
	// just dump some stuff
	cout << "~phiReader: " << endl;
	cout << "\treco counter: " << reco_single + reco_double << endl;
	cout << "\t\treco single: " << reco_single << endl;
	cout << "\t\treco double: " << reco_double << endl;
	cout << "\ttotal counter: " << total_counter << endl;
} // ~phiReader()

void phiReader::eventProcessing()
{
	multiset<int> current_decay;
	map<int,int>::const_iterator it;
	TGenCand *pGen;
	unsigned j,nc;
	
	decay_indices.clear();
	
	// count the number of real decays using the generator block
	nc = fpEvt->nGenCands();
	for (j = 0; j < nc; j++) {
		pGen = fpEvt->getGenCand(j);
		if (abs(pGen->fID) == 531) {
			current_decay.clear();
			buildDecay(pGen,&current_decay);
			current_decay.erase(22); // remove Bremsstrahlung
			if (trueDecay == current_decay)
				total_counter++;
		}
	}
	
	// do the normal stuff
	massReader::eventProcessing();
	
	// adjust the reco counters
	for (it = decay_indices.begin(); it != decay_indices.end(); ++it) {
		if (it->second == 1) reco_single++;
		else if(it->second == 2) reco_double++;
		else cout << "Reconstructed candidate with " << it->second << " muons??" << endl;
	}
} // eventProcessing()


int phiReader::loadCandidateVariables(TAnaCand *pCand)
{
	int result,type,j;
	int firstMu = 1, firstKp = 1;
	TAnaTrack *sigTrack,*recTrack;
	TAnaCand *dau;
	
	// default initialization
	fMassJPsi = fMassPhi = fDeltaR = fDeltaR_Kaons = -1.0f;
	fPtMu1 = fPtMu2 = fPtKp1 = fPtKp2 = 0.0f;
	
	fPlabMu1 = TVector3();
	fPlabMu2 = TVector3();
	fPlabKp1 = TVector3();
	fPlabKp2 = TVector3();
	
	fTrackQual_mu1 = fTrackQual_mu2 = fTrackQual_kp1 = fTrackQual_kp2 = 0;
	fQ_mu1 = fQ_mu2 = fQ_kp1 = fQ_kp2 = 0;
	fD3_BsJpsi = fD3e_BsJpsi = -1.0f;
	
	if (pCand->fType % 1000 != 531) return 0;
	
	result = massReader::loadCandidateVariables(pCand);
	
	// set the momenta
	for (j = pCand->fSig1; j <= pCand->fSig2; j++) {
		
		sigTrack = fpEvt->getSigTrack(j);
		type = abs(sigTrack->fMCID);
		recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
		if (type == 13) {
			if (firstMu) {
				fPlabMu1 = sigTrack->fPlab;
				fTrackQual_mu1 = recTrack->fTrackQuality;
				fQ_mu1 = recTrack->fQ;
			} else {
				fPlabMu2 = sigTrack->fPlab;
				fTrackQual_mu2 = recTrack->fTrackQuality;
				fQ_mu2 = recTrack->fQ;
			}
			firstMu = 0;
		} else if (type == 321) {
			if (firstKp) {
				fPlabKp1 = sigTrack->fPlab;
				fTrackQual_kp1 = recTrack->fTrackQuality;
				fQ_kp1 = recTrack->fQ;
			} else {
				fPlabKp2 = sigTrack->fPlab;
				fTrackQual_kp2 = recTrack->fTrackQuality;
				fQ_kp2 = recTrack->fQ;
			}
			firstKp = 0;
		}
	}
	
	// mu1 is with higher pt
	if (fPlabMu1.Perp() < fPlabMu2.Perp()) {
		swap(fPlabMu1,fPlabMu2);
		swap(fTrackQual_mu1,fTrackQual_mu2);
		swap(fQ_mu1,fQ_mu2);
	}
	
	// kp1 is with higher pt
	if (fPlabKp1.Perp() < fPlabKp2.Perp()) {
		swap(fPlabKp1,fPlabKp2);
		swap(fTrackQual_kp1,fTrackQual_kp2);
		swap(fQ_kp1,fQ_kp2);
	}
	
	// set the j/psi mass
	for (j = pCand->fDau1; j <= pCand->fDau2 && j>=0; j++) {
		
		dau = fpEvt->getCand(j);
		if (dau->fType % 1000 == 443)
			fMassJPsi = dau->fMass;
		else if(dau->fType % 1000 == 333)
			fMassPhi = dau->fMass;
	}
	
	if (fMassJPsi < 0.0f) { // no jpsi subcandidate found
		TLorentzVector m1,m2,jpsi;
		
		m1.SetXYZM(fPlabMu1.X(),fPlabMu1.Y(),fPlabMu1.Z(),MMUON);
		m2.SetXYZM(fPlabMu2.X(),fPlabMu2.Y(),fPlabMu2.Z(),MMUON);
		
		jpsi = m1 + m2;
		fMassJPsi = jpsi.M();
	}
	
	if (fMassPhi < 0.0f) { // no phi subcandidate found
		TLorentzVector k1,k2,phi;
		
		k1.SetXYZM(fPlabKp1.X(),fPlabKp1.Y(),fPlabKp1.Z(),MKAON);
		k2.SetXYZM(fPlabKp2.X(),fPlabKp2.Y(),fPlabKp2.Z(),MKAON);
		
		phi = k1 + k2;
		fMassPhi = phi.M();
	}
	
	// set the deltaR of the J/Psi w.r.t. Phi
	if ( (fPlabMu1 + fPlabMu2).Perp() > 0 && (fPlabKp1 + fPlabKp2).Perp() > 0 )
		fDeltaR = (float)(fPlabMu1 + fPlabMu2).DeltaR(fPlabKp1 + fPlabKp2);
	
	// set the deltaR of the Kaons
	if ( fPlabKp1.Perp() > 0 && fPlabKp2.Perp() > 0 )
		fDeltaR_Kaons = (float)fPlabKp1.DeltaR(fPlabKp2);
	
	// set the transveral momenta
	fPtMu1 = fPlabMu1.Perp();
	fPtMu2 = fPlabMu2.Perp();
	fPtKp1 = fPlabKp1.Perp();
	fPtKp2 = fPlabKp2.Perp();
	
	return result;
} // loadCandidateVariables()

void phiReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("mass_jpsi",&fMassJPsi,"mass_jpsi/F");
	reduced_tree->Branch("mass_phi",&fMassPhi,"mass_phi/F");
	reduced_tree->Branch("deltaR",&fDeltaR,"deltaR/F");
	reduced_tree->Branch("deltaR_kaons",&fDeltaR_Kaons,"deltaR_kaons/F");
	reduced_tree->Branch("plab_mu1","TVector3",&fPlabMu1Ptr);
	reduced_tree->Branch("plab_mu2","TVector3",&fPlabMu2Ptr);
	reduced_tree->Branch("plab_kp1","TVector3",&fPlabKp1Ptr);
	reduced_tree->Branch("plab_kp2","TVector3",&fPlabKp2Ptr);
	reduced_tree->Branch("pt_mu1",&fPtMu1,"pt_mu1/F");
	reduced_tree->Branch("pt_mu2",&fPtMu2,"pt_mu2/F");
	reduced_tree->Branch("pt_kp1",&fPtKp1,"pt_kp1/F");
	reduced_tree->Branch("pt_kp2",&fPtKp2,"pt_kp2/F");
	reduced_tree->Branch("track_qual_mu1",&fTrackQual_mu1,"track_qual_mu1/I");
	reduced_tree->Branch("track_qual_mu2",&fTrackQual_mu2,"track_qual_mu2/I");
	reduced_tree->Branch("track_qual_kp1",&fTrackQual_kp1,"track_qual_kp1/I");
	reduced_tree->Branch("track_qual_kp2",&fTrackQual_kp2,"track_qual_kp2/I");
	reduced_tree->Branch("q_mu1",&fQ_mu1,"q_mu1/I");
	reduced_tree->Branch("q_mu2",&fQ_mu2,"q_mu2/I");
	reduced_tree->Branch("q_kp1",&fQ_kp1,"q_kp1/I");
	reduced_tree->Branch("q_kp2",&fQ_kp2,"q_kp2/I");
	reduced_tree->Branch("d3_bs_to_jpsi",&fD3_BsJpsi,"d3_bs_to_jpsi/F");
	reduced_tree->Branch("d3e_bs_to_jpsi",&fD3e_BsJpsi,"d3e_bs_to_jpsi/F");
} // bookHist()

int phiReader::checkTruth(TAnaCand *pCand)
{
	int result;
	multiset<int> particles;
	TAnaTrack *track;
	TGenCand *gen;
	
	result = massReader::checkTruth(pCand);
	if(!result) goto bail;
	
	// check if the decay conincides
	track = fpEvt->getSigTrack(pCand->fSig1);
	track = fpEvt->getRecTrack(track->fIndex);
	
	gen = fpEvt->getGenCand(track->fGenIndex);
	while (abs(gen->fID) != 531)
		gen = fpEvt->getGenCand(gen->fMom1);
	
	buildDecay(gen,&particles);
	particles.erase(22); // remove Bremsstrahlung
	
	result = (particles == trueDecay);
	
	if (result && fCandidate == 300531)
		decay_indices[gen->fNumber] = max((int)fNbrMuons,decay_indices[gen->fNumber]);
bail:
	return result;
} // checkTruth()
