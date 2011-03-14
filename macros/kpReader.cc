#include "kpReader.hh"

#include <TLorentzVector.h>
#include <utility>

using namespace std;

const static float kMassBplus = 5.2792;

kpReader::kpReader(TChain *tree, TString evtClassName) :
	massReader(tree,evtClassName),
	fCutTrackQual_mu1(0),
	fCutTrackQual_mu2(0),
	fCutTrackQual_kp(0),
	fCutMuID_mask(0),
	fCutMuID_reqall(false),
	fCutOppSign_mu(false),
	fCutMass_JPsiLow(0.0),
	fCutMass_JPsiHigh(0.0),
	fCutPt_Kaon(0.0),
	total_counter(0),
	reco_counter(0),
	reco_single(0),
	reco_double(0)
{
	fTreeName = "kpReader reduced tree";
	fTruthType = 521; // searching B+
	
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

void kpReader::clearVariables()
{
	massReader::clearVariables();
	
	fMassJPsi = 0.0;
	fMassJPsiRec = 0.0;
	fDeltaR = 0.0;
	
	fChi2Jpsi = 0.0;
	fMassJpsiKp = 0.0;
	
	fPtMu1 = 0.0;
	fPtMu2 = 0.0;
	fPtKp = 0.0;
	
	fMuID1 = 0; fMuID2 = 0;
	fEtaMu1 = 0.0; fEtaMu2 = 0.0;
	
	fPtKp_Gen = 0.0;
	fEtaKp_Gen = 0.0;
	
	fTrackQual_mu1 = 0;
	fTrackQual_mu2 = 0;
	fTrackQual_kp = 0;
	
	fQ_mu1 = 0.0;
	fQ_mu2 = 0.0;
	fQ_kp = 0.0;
	
	fD3_BpJpsi = 0.0;
	fD3e_BpJpsi = 0.0;
} // clearVariables()

int kpReader::loadCandidateVariables(TAnaCand *pCand)
{
	int type, first_mu = 1;
	TAnaTrack *sigTrack,*recTrack;
	TAnaCand *jpsiCand;
	TGenCand *kaonGen;
	TVector3 plabMu1;
	TVector3 plabMu2;
	TVector3 plabKp;
	int result,j,k;
	map<int,int> cand_tracks;
	
	TLorentzVector mu1;
	TLorentzVector mu2;
	
	// default initialization
	fMassJPsi = -1.0f;
	fMassJPsiRec = -1.0f;
	fDeltaR = -1.0f;
	fCtau = 0.0f;
	fChi2Jpsi = -1.0f;
	fMassJpsiKp = -1.0f;
	fMuID1 = 0;
	fMuID2 = 0;
	fEtaMu1 = 0.0f;
	fEtaMu2 = 0.0f;
		
	result = massReader::loadCandidateVariables(pCand);
	
	fCtau = kMassBplus / pCand->fPlab.Mag() * fD3;
	
	// set the constraint mass value
	findCandStructure(pCand,&cand_tracks);
	jpsiCand = findCandidate(400521,&cand_tracks);
	if (jpsiCand) fMassConstraint = jpsiCand->fMass;
	cand_tracks.clear();
	
	// set the momenta, iterate
	findAllTrackIndices(pCand,&cand_tracks);
	for (map<int,int>::const_iterator it = cand_tracks.begin(); it!=cand_tracks.end(); ++it) {
		
		sigTrack = fpEvt->getSigTrack(it->second);
		type = abs(sigTrack->fMCID);
		recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
		if (type == 13) {
			// muon
			if (first_mu) {
				plabMu1 = sigTrack->fPlab;
				fTrackQual_mu1 = recTrack->fTrackQuality;
				fQ_mu1 = recTrack->fQ;
				mu1.SetXYZM(recTrack->fPlab.X(),recTrack->fPlab.Y(),recTrack->fPlab.Z(),MUMASS);
				fMuID1 = recTrack->fMuID;
				fEtaMu1 = sigTrack->fPlab.Eta();
				
			} else {
				plabMu2 = sigTrack->fPlab;
				fTrackQual_mu2 = recTrack->fTrackQuality;
				fQ_mu2 = recTrack->fQ;
				mu2.SetXYZM(recTrack->fPlab.X(),recTrack->fPlab.Y(),recTrack->fPlab.Z(),MUMASS);
				fMuID2 = recTrack->fMuID;
				fEtaMu2 = sigTrack->fPlab.Eta();
			}
			first_mu = 0;
		} else if (type == 321) {
			plabKp = sigTrack->fPlab;
			fQ_kp = recTrack->fQ;
			fTrackQual_kp = recTrack->fTrackQuality;
		}
	}
	
	// muon1 ist usually that one with the bigger pt
	if (plabMu1.Perp() < plabMu2.Perp()) {
		swap(plabMu1,plabMu2);
		swap(fTrackQual_mu1,fTrackQual_mu2);
		swap(fQ_mu1,fQ_mu2);
		swap(fMuID1,fMuID2);
		swap(fEtaMu1,fEtaMu2);
	}
	
	// set the jpsi mass
	for (j = pCand->fDau1; j <= pCand->fDau2 && j >= 0; j++) {
		
		jpsiCand = fpEvt->getCand(j);
		if (jpsiCand->fType % 1000 == 443) {
			TLorentzVector *m1 = NULL,*m2 = NULL,*kp = NULL;
			
			fMassJPsi = jpsiCand->fMass;
			fChi2Jpsi = jpsiCand->fVtx.fChi2;
			fD3_BpJpsi = jpsiCand->fVtx.fD3d;
			fD3e_BpJpsi = jpsiCand->fVtx.fD3dE;
			
			// calculate the mass of the k+ by computing from j/psi sig tracks + kaon rectrack
			for (k = jpsiCand->fSig1; 0 <= k && k <= jpsiCand->fSig2; k++) {
				
				sigTrack = fpEvt->getSigTrack(k);
				if (abs(sigTrack->fMCID) == 13) {
					if (!m1) {
						m1 = new TLorentzVector;
						m1->SetXYZM(sigTrack->fPlab.X(),sigTrack->fPlab.Y(),sigTrack->fPlab.Z(),MMUON);
					} else if (!m2) {
						m2 = new TLorentzVector;
						m2->SetXYZM(sigTrack->fPlab.X(),sigTrack->fPlab.Y(),sigTrack->fPlab.Z(),MMUON);
					}
				}
			}
			
			for (k = pCand->fSig1; 0 <= k && k <= pCand->fSig2; k++) {
				sigTrack = fpEvt->getSigTrack(k);
				if (abs(sigTrack->fMCID) == 321 && !kp) {
					kp = new TLorentzVector;
					kp->SetXYZM(sigTrack->fPlab.X(),sigTrack->fPlab.Y(),sigTrack->fPlab.Z(),MKAON);
				}
			}
			
			fMassJpsiKp = ((*m1) + (*m2) + (*kp)).M();
			
			delete m1;
			delete m2;
			delete kp;
			break;
		}
	}
	
	fMassJPsiRec = (mu1 + mu2).M();

	// set the deltaR of the J/Psi w.r.t. Kp
	if (plabKp.Perp() > 0 && (plabMu1 + plabMu2).Perp() > 0)
	  fDeltaR = (float)plabKp.DeltaR(plabMu1 + plabMu2);
	
	fPtMu1 = plabMu1.Perp();
	fPtMu2 = plabMu2.Perp();
	fPtKp = plabKp.Perp();
	
	// set the generator variables
	cand_tracks.clear();
	findAllTrackIndices(pCand,&cand_tracks);
	for (map<int,int>::const_iterator it = cand_tracks.begin(); it != cand_tracks.end(); ++it) {
		recTrack = fpEvt->getRecTrack(it->first);
		if (recTrack->fGenIndex < 0 || recTrack->fGenIndex >= fpEvt->nGenCands())
			continue;
		
		kaonGen = fpEvt->getGenCand(recTrack->fGenIndex);
		if (abs(kaonGen->fID) == 321) {
			fPtKp_Gen = kaonGen->fP.Perp();
			fEtaKp_Gen = kaonGen->fP.Eta();
		}
	}	

	return result;
} // loadCandidateVariables()

void kpReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("mass_jpsi",&fMassJPsi,"mass_jpsi/F");
	reduced_tree->Branch("mass_jpsi_rec",&fMassJPsiRec,"mass_jpsi_rec/F");
	reduced_tree->Branch("deltaR",&fDeltaR,"deltaR/F");
	reduced_tree->Branch("pt_mu1",&fPtMu1,"pt_mu1/F");
	reduced_tree->Branch("pt_mu2",&fPtMu2,"pt_mu2/F");
	reduced_tree->Branch("pt_kp",&fPtKp,"pt_kp/F");
	reduced_tree->Branch("id_mu1",&fMuID1,"id_mu1/I");
	reduced_tree->Branch("id_mu2",&fMuID2,"id_mu2/I");
	reduced_tree->Branch("eta_mu1",&fEtaMu1,"eta_mu1/F");
	reduced_tree->Branch("eta_mu2",&fEtaMu2,"eta_mu2/F");
	reduced_tree->Branch("pt_kp_gen",&fPtKp_Gen,"pt_kp_gen/F");
	reduced_tree->Branch("eta_kp_gen",&fEtaKp_Gen,"eta_kp_gen/F");
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

bool kpReader::parseCut(char *cutName, float cutLow, float cutHigh, int dump)
{
	bool parsed;
	
	// parse the options..
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
	
	parsed = (strcmp(cutName,"TRACK_QUAL_KP") == 0);
	if (parsed) {
		fCutTrackQual_kp = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_KP: " << fCutTrackQual_kp << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"OPPOSITE_SIGN_MU") == 0);
	if (parsed) {
		fCutOppSign_mu = (cutLow != 0.0);
		if (dump) cout << "OPPOSITE_SIGN_MU: " << fCutOppSign_mu << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MASS_JPSI") == 0);
	if (parsed) {
		fCutMass_JPsiLow = cutLow;
		fCutMass_JPsiHigh = cutHigh;
		if (dump) cout << "MASS_JPSI: (" << fCutMass_JPsiLow << ", " << fCutMass_JPsiHigh << ")" << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"PT_KAON") == 0);
	if (parsed) {
		fCutPt_Kaon = cutLow;
		if (dump) cout << "PT_KAON: " << fCutPt_Kaon << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MUID") == 0);
	if (parsed) {
		fCutMuID_mask = (int)cutLow;
		fCutMuID_reqall = (cutHigh != 0.0);
		if (dump) cout << "MUID: mask " << fCutMuID_mask << " and require all " << fCutMuID_reqall << endl;
		goto bail;
	}
	
	// nothing parsed, then maybe massReader knows about it.
	parsed = massReader::parseCut(cutName, cutLow, cutHigh, dump);
	
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
	
	if (fCutMass_JPsiLow > 0.0) {
		pass = fMassJPsi > fCutMass_JPsiLow;
		if(!pass) goto bail;
	}
	
	if (fCutMass_JPsiHigh > 0.0) {
		pass = fMassJPsi < fCutMass_JPsiHigh;
		if(!pass) goto bail;
	}
	
	if (fCutPt_Kaon > 0.0) {
		pass = fPtKp > fCutPt_Kaon;
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
