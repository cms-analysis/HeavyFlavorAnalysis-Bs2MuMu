#include "kpReader.hh"

#include <TLorentzVector.h>
#include <utility>

using namespace std;

const static float kMassBplus = 5.2792;

kpReader::kpReader(TChain *tree, TString evtClassName) :
	massReader(tree,evtClassName),
	fCutTrackQual_kp(0),
	fCutMass_JPsiLow(0.0),
	fCutMass_JPsiHigh(0.0),
	fCutPt_Kaon(0.0)
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
} // ~kpReader()

void kpReader::clearVariables()
{
	massReader::clearVariables();
	
	fMassJPsi = 0.0;
	fMassJPsiRec = 0.0;
	fDeltaR = 0.0;
	
	fChi2Jpsi = 0.0;
	
	fPtKp = 0.0;
	fEtaKp = 0.0;
	
	fPtKp_Gen = 0.0;
	fEtaKp_Gen = 0.0;
	
	fTrackQual_kp = 0;
	
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
	int result,j;
	map<int,int> cand_tracks;
	
	TLorentzVector mu1;
	TLorentzVector mu2;
	
	// default initialization
	fMassJPsi = -1.0f;
	fMassJPsiRec = -1.0f;
	fChi2Jpsi = -1.0f;
			
	result = massReader::loadCandidateVariables(pCand);
	
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
				mu1.SetXYZM(recTrack->fPlab.X(),recTrack->fPlab.Y(),recTrack->fPlab.Z(),MUMASS);
			} else {
				plabMu2 = sigTrack->fPlab;
				mu2.SetXYZM(recTrack->fPlab.X(),recTrack->fPlab.Y(),recTrack->fPlab.Z(),MUMASS);
			}
			first_mu = 0;
		} else if (type == 321) {
			plabKp = sigTrack->fPlab;
			fEtaKp = sigTrack->fPlab.Eta();
			fQ_kp = recTrack->fQ;
			fTrackQual_kp = recTrack->fTrackQuality;
		}
	}
	
	// set the jpsi mass
	for (j = pCand->fDau1; j <= pCand->fDau2 && j >= 0; j++) {
		
		jpsiCand = fpEvt->getCand(j);
		if (jpsiCand->fType % 1000 == 443) {
			
			fMassJPsi = jpsiCand->fMass;
			fChi2Jpsi = jpsiCand->fVtx.fChi2;
			fD3_BpJpsi = jpsiCand->fVtx.fD3d;
			fD3e_BpJpsi = jpsiCand->fVtx.fD3dE;
			
			break;
		}
	}
	
	fMassJPsiRec = (mu1 + mu2).M();
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
	reduced_tree->Branch("pt_kp",&fPtKp,"pt_kp/F");
	reduced_tree->Branch("eta_kp",&fEtaKp,"eta_kp/F");
	reduced_tree->Branch("pt_kp_gen",&fPtKp_Gen,"pt_kp_gen/F");
	reduced_tree->Branch("eta_kp_gen",&fEtaKp_Gen,"eta_kp_gen/F");
	reduced_tree->Branch("track_qual_kp",&fTrackQual_kp,"track_qual_kp/I");
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
	
bail:
	return result;
} // checkTruth()

bool kpReader::parseCut(char *cutName, float cutLow, float cutHigh, int dump)
{
	bool parsed;
	
	parsed = (strcmp(cutName,"TRACK_QUAL_KP") == 0);
	if (parsed) {
		fCutTrackQual_kp = (int)cutLow;
		if (dump) cout << "TRACK_QUAL_KP: " << fCutTrackQual_kp << endl;
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
	
	// nothing parsed, then maybe massReader knows about it.
	parsed = massReader::parseCut(cutName, cutLow, cutHigh, dump);
	
bail:
	return parsed;
} // parseCut()

bool kpReader::applyCut()
{
	bool pass = true;
	
	// check track quality of kp
	if (fCutTrackQual_kp) {
		pass = (fTrackQual_kp & fCutTrackQual_kp) != 0;
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
	
	// check superclass cuts
	pass = massReader::applyCut();
	
bail:
	// return result
	return pass;
} // applyCut()
