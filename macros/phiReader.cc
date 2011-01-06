#include "phiReader.hh"

#include <utility>
#include <TLorentzVector.h>

using namespace std;

const static float kMassBs = 5.3663;

phiReader::phiReader(TChain *tree, TString evtClassName) :
	massReader(tree,evtClassName),
	fCutTrackQual_mu1(0),
	fCutTrackQual_mu2(0),
	fCutTrackQual_kp1(0),
	fCutTrackQual_kp2(0),
	fCutMuID_mask(0),
	fCutMuID_reqall(false),
	fCutOppSign_mu(false),
	fCutOppSign_kp(false),
	fCutMass_JPsiLow(0.0),
	fCutMass_JPsiHigh(0.0),
	fCutMass_PhiLow(0.0),
	fCutMass_PhiHigh(0.0),
	fCutPt_Kp2(0.0),
	fCutDeltaR_Kaons(0.0),
	fCutPt_Mu2(0.0),
	total_counter(0),
	reco_single(0),
	reco_double(0)
{
	cout << "Instantiating phiReader..." << endl;	
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
	TVector3 plabMu1;
	TVector3 plabMu2;
	TVector3 plabKp1;
	TVector3 plabKp2;
	TLorentzVector mu1;
	TLorentzVector mu2;
	map<int,int> cand_tracks;
	
	// default initialization
	fMassJPsi = fMassPhi = fDeltaR = fDeltaR_Kaons = -1.0f;
	fPtMu1 = fPtMu2 = fPtKp1 = fPtKp2 = 0.0f;
	
	fTrackQual_mu1 = fTrackQual_mu2 = fTrackQual_kp1 = fTrackQual_kp2 = 0;
	fQ_mu1 = fQ_mu2 = fQ_kp1 = fQ_kp2 = 0;
	fD3_BsJpsi = fD3e_BsJpsi = -1.0f;
	fCtau = 0.0;
	fMuID1 = 0;
	fMuID2 = 0;
	fEtaMu1 = 0.0f;
	fEtaMu2 = 0.0f;
	fMuID_Kp1 = 0;
	fMuID_Kp2 = 0;
	
	if (pCand->fType % 1000 != 531) return 0;
	
	result = massReader::loadCandidateVariables(pCand);
	
	fCtau = kMassBs / pCand->fPlab.Mag() * fD3; // estimate the proper time!
	
	// set the constraint mass value
	findCandStructure(pCand,&cand_tracks);
	dau = findCandidate(400531,&cand_tracks);
	if (dau) fMassConstraint = dau->fMass;
	cand_tracks.clear();
	
	// set the momenta
	findAllTrackIndices(pCand,&cand_tracks);
	for (map<int,int>::const_iterator it = cand_tracks.begin(); it!=cand_tracks.end(); ++it) {
	
		sigTrack = fpEvt->getSigTrack(it->second);
		type = abs(sigTrack->fMCID);
		recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
		if (type == 13) {
			if (firstMu) {
				plabMu1 = sigTrack->fPlab;
				fTrackQual_mu1 = recTrack->fTrackQuality;
				fQ_mu1 = recTrack->fQ;
				mu1.SetXYZM(recTrack->fPlab.X(),recTrack->fPlab.Y(),recTrack->fPlab.Z(),MUMASS);
				fMuID1 = recTrack->fMuID;
				if (fMuID1 < 0) fMuID1 = 0;
				fEtaMu1 = sigTrack->fPlab.Eta();
			} else {
				plabMu2 = sigTrack->fPlab;
				fTrackQual_mu2 = recTrack->fTrackQuality;
				fQ_mu2 = recTrack->fQ;
				mu2.SetXYZM(recTrack->fPlab.X(),recTrack->fPlab.Y(),recTrack->fPlab.Z(),MUMASS);
				fMuID2 = recTrack->fMuID;
				if (fMuID2 < 0) fMuID2 = 0;
				fEtaMu2 = sigTrack->fPlab.Eta();
			}
			firstMu = 0;
		} else if (type == 321) {
			if (firstKp) {
				plabKp1 = sigTrack->fPlab;
				fTrackQual_kp1 = recTrack->fTrackQuality;
				fQ_kp1 = recTrack->fQ;
				fMuID_Kp1 = recTrack->fMuID;
				if (fMuID_Kp1 < 0) fMuID_Kp1 = 0;
			} else {
				plabKp2 = sigTrack->fPlab;
				fTrackQual_kp2 = recTrack->fTrackQuality;
				fQ_kp2 = recTrack->fQ;
				fMuID_Kp2 = recTrack->fMuID;
				if (fMuID_Kp2 < 0) fMuID_Kp2 = 0;
			}
			firstKp = 0;
		}
	}
	
	// mu1 is with higher pt
	if (plabMu1.Perp() < plabMu2.Perp()) {
		swap(plabMu1,plabMu2);
		swap(fTrackQual_mu1,fTrackQual_mu2);
		swap(fQ_mu1,fQ_mu2);
		swap(fMuID1,fMuID2);
		swap(fEtaMu1,fEtaMu2);
	}
	
	// kp1 is with higher pt
	if (plabKp1.Perp() < plabKp2.Perp()) {
		swap(plabKp1,plabKp2);
		swap(fTrackQual_kp1,fTrackQual_kp2);
		swap(fQ_kp1,fQ_kp2);
		swap(fMuID_Kp1,fMuID_Kp2);
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
		
		m1.SetXYZM(plabMu1.X(),plabMu1.Y(),plabMu1.Z(),MMUON);
		m2.SetXYZM(plabMu2.X(),plabMu2.Y(),plabMu2.Z(),MMUON);
		
		jpsi = m1 + m2;
		fMassJPsi = jpsi.M();
	}
	
	fMassJPsiRec = (mu1 + mu2).M();
	
	if (fMassPhi < 0.0f) { // no phi subcandidate found
		TLorentzVector k1,k2,phi;
		
		k1.SetXYZM(plabKp1.X(),plabKp1.Y(),plabKp1.Z(),MKAON);
		k2.SetXYZM(plabKp2.X(),plabKp2.Y(),plabKp2.Z(),MKAON);
		
		phi = k1 + k2;
		fMassPhi = phi.M();
	}
	
	// set the deltaR of the J/Psi w.r.t. Phi
	if ( (plabMu1 + plabMu2).Perp() > 0 && (plabKp1 + plabKp2).Perp() > 0 )
		fDeltaR = (float)(plabMu1 + plabMu2).DeltaR(plabKp1 + plabKp2);
	
	// set the deltaR of the Kaons
	if ( plabKp1.Perp() > 0 && plabKp2.Perp() > 0 )
		fDeltaR_Kaons = (float)plabKp1.DeltaR(plabKp2);
	
	// set the transveral momenta
	fPtMu1 = plabMu1.Perp();
	fPtMu2 = plabMu2.Perp();
	fPtKp1 = plabKp1.Perp();
	fPtKp2 = plabKp2.Perp();
	
	return result;
} // loadCandidateVariables()

void phiReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("mass_jpsi",&fMassJPsi,"mass_jpsi/F");
	reduced_tree->Branch("mass_jpsi_rec",&fMassJPsiRec,"mass_jpsi_rec/F");
	reduced_tree->Branch("mass_phi",&fMassPhi,"mass_phi/F");
	reduced_tree->Branch("deltaR",&fDeltaR,"deltaR/F");
	reduced_tree->Branch("deltaR_kaons",&fDeltaR_Kaons,"deltaR_kaons/F");
	reduced_tree->Branch("pt_mu1",&fPtMu1,"pt_mu1/F");
	reduced_tree->Branch("pt_mu2",&fPtMu2,"pt_mu2/F");
	reduced_tree->Branch("pt_kp1",&fPtKp1,"pt_kp1/F");
	reduced_tree->Branch("pt_kp2",&fPtKp2,"pt_kp2/F");
	reduced_tree->Branch("id_mu1",&fMuID1,"id_mu1/I");
	reduced_tree->Branch("id_mu2",&fMuID2,"id_mu2/I");
	reduced_tree->Branch("eta_mu1",&fEtaMu1,"eta_mu1/F");
	reduced_tree->Branch("eta_mu2",&fEtaMu2,"eta_mu2/F");
	reduced_tree->Branch("id_kp1",&fMuID_Kp1,"id_kp1/I");
	reduced_tree->Branch("id_kp2",&fMuID_Kp2,"id_kp2/I");
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

bool phiReader::parseCut(char *cutName, float cutLow, float cutHigh, int dump)
{
	bool parsed;
	
	// parse the options...
	parsed = (strcmp(cutName,"TRACK_QUAL_MU1") == 0);
	if (parsed) {
		fCutTrackQual_mu1 = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_MU1: " << fCutTrackQual_mu1 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRACK_QUAL_MU2") == 0);
	if (parsed) {
		fCutTrackQual_mu2 = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_MU2: " << fCutTrackQual_mu2 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRACK_QUAL_KP1") == 0);
	if (parsed) {
		fCutTrackQual_kp1 = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_KP1: " << fCutTrackQual_kp1 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"TRACK_QUAL_KP2") == 0);
	if (parsed) {
		fCutTrackQual_kp2 = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_KP2: " << fCutTrackQual_kp2 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"OPPOSITE_SIGN_MU") == 0);
	if (parsed) {
		fCutOppSign_mu = (cutLow != 0.0f);
		if (dump) cout << "OPPOSITE_SIGN_MU: " << fCutOppSign_mu << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"OPPOSITE_SIGN_KP") == 0);
	if (parsed) {
		fCutOppSign_kp = (cutLow != 0.0f);
		if (dump) cout << "OPPOSITE_SIGN_KP: " << fCutOppSign_kp << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MASS_JPSI") == 0);
	if (parsed) {
		fCutMass_JPsiLow = cutLow;
		fCutMass_JPsiHigh = cutHigh;
		if (dump) cout << "MASS_JPSI: (" << fCutMass_JPsiLow << ", " << fCutMass_JPsiHigh << ")" << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MASS_PHI") == 0);
	if (parsed) {
		fCutMass_PhiLow = cutLow;
		fCutMass_PhiHigh = cutHigh;
		if (dump) cout << "MASS_PHI: (" << fCutMass_PhiLow << ", " << fCutMass_PhiHigh << ")" << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"PT_KAON2") == 0);
	if (parsed) {
		fCutPt_Kp2 = cutLow;
		if (dump) cout << "PT_KAON2: " << fCutPt_Kp2 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"PT_MU2") == 0);
	if (parsed) {
		fCutPt_Mu2 = cutLow;
		if (dump) cout << "PT_MU2: " << fCutPt_Mu2 << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"DELTAR_KAONS") == 0);
	if (parsed) {
		fCutDeltaR_Kaons = cutLow;
		if (dump) cout << "DELTAR_KAONS: " << fCutDeltaR_Kaons << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MUID") == 0);
	if (parsed) {
		fCutMuID_mask = (int)cutLow;
		fCutMuID_reqall = (cutHigh != 0.0);
		if (dump) cout << "MUID: mask " << fCutMuID_mask << " and require all " << fCutMuID_reqall << endl;
		goto bail;
	}
	
	// ask superclass if nothing parsed yet
	parsed = massReader::parseCut(cutName, cutLow, cutHigh, dump);
bail:
	return parsed;
} // parseCut()


bool phiReader::applyCut()
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
	
	// check track quality of kp1
	if (fCutTrackQual_kp1) {
		pass = (fTrackQual_kp1 & fCutTrackQual_kp1) != 0;
		if (!pass) goto bail;
	}
	
	// check track quality of kp2
	if (fCutTrackQual_kp2) {
		pass = (fTrackQual_kp2 & fCutTrackQual_kp2) != 0;
		if (!pass) goto bail;
	}
	
	// check opposite sign of muons
	if (fCutOppSign_mu) {
		pass = fQ_mu1 != fQ_mu2;
		if (!pass) goto bail;
	}
	
	// check opposite sign of kaons
	if (fCutOppSign_kp) {
		pass = fQ_kp1 != fQ_kp2;
		if (!pass) goto bail;
	}
	
	if (fCutMass_JPsiLow > 0.0) {
		pass = fMassJPsi > fCutMass_JPsiLow;
		if (!pass) goto bail;
	}
	
	if (fCutMass_JPsiHigh > 0.0) {
		pass = fMassJPsi < fCutMass_JPsiHigh;
		if (!pass) goto bail;
	}
	
	if (fCutMass_PhiLow > 0.0) {
		pass = fMassPhi > fCutMass_PhiLow;
		if (!pass) goto bail;
	}
	
	if (fCutMass_PhiHigh > 0.0) {
		pass = fMassPhi < fCutMass_PhiHigh;
		if (!pass) goto bail;
	}
	
	if (fCutPt_Mu2 > 0.0) {
		pass = fPtMu2 > fCutPt_Mu2;
		if (!pass) goto bail;
	}
	
	if (fCutPt_Kp2 > 0.0) {
		pass = fPtKp2 > fCutPt_Kp2;
		if (!pass) goto bail;
	}
	
	if (fCutDeltaR_Kaons > 0.0) {
		pass = fDeltaR_Kaons < fCutDeltaR_Kaons;
		if (!pass) goto bail;
	}
	
	if (fCutMuID_mask != 0) {
		
		if (fCutMuID_reqall)
			pass = ((fMuID1 & fCutMuID_mask) == fCutMuID_mask) && ((fMuID2 & fCutMuID_mask) == fCutMuID_mask);
		else
			pass = ((fMuID1 & fCutMuID_mask) != 0) && ((fMuID2 & fCutMuID_mask) != 0);
		
		if (!pass) goto bail;
	}
	
	// check superclass cuts
	pass = massReader::applyCut();
	
bail:
	// return result
	return pass;
} // applyCut()
