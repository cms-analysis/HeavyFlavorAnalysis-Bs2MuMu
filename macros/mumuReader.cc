/*
 *  mumuReader.cpp
 *  macros
 *
 *  Created by Christoph on 09.01.11.
 *  Copyright 2011 PSI. All rights reserved.
 *
 */

#include "mumuReader.hh"

const static float kMassBs = 5.3663;

using namespace std;

mumuReader::mumuReader(TChain *tree, TString evtClassName) :
	massReader(tree, evtClassName),
	fCutTrackQual_mu1(0),
	fCutTrackQual_mu2(0),
	fCutMuID_mask(0),
	fCutMuID_reqall(false),
	fCutOppSign_mu(false)
{
	cout << "Instantiating mumuReader..." << flush;
	fTreeName = "mumuReader reduced tree.";
	fTruthType = 531;
	
	trueDecay.insert(13); // mu
	trueDecay.insert(13); // mu
	trueDecay.insert(fTruthType); // Bs
	
	cout << " ok" << endl;
} // mumuReader()

void mumuReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("pt_mu1",&fPtMu1,"pt_mu1/F");
	reduced_tree->Branch("pt_mu2",&fPtMu2,"pt_mu2/F");
	reduced_tree->Branch("id_mu1",&fMuID1,"id_mu1/I");
	reduced_tree->Branch("id_mu2",&fMuID2,"id_mu2/I");
	reduced_tree->Branch("deltaR",&fDeltaR,"deltaR/F");
	reduced_tree->Branch("eta_mu1",&fEtaMu1,"eta_mu1/F");
	reduced_tree->Branch("eta_mu2",&fEtaMu2,"eta_mu2/F");
	reduced_tree->Branch("track_qual_mu1",&fTrackQual_mu1,"track_qual_mu1/I");
	reduced_tree->Branch("track_qual_mu2",&fTrackQual_mu2,"track_qual_mu2/I");
	reduced_tree->Branch("q_mu1",&fQ_mu1,"q_mu1/I");
	reduced_tree->Branch("q_mu2",&fQ_mu2,"q_mu2/I");
} // bookHist()

void mumuReader::clearVariables()
{
	massReader::clearVariables();
	
	fPtMu1 = 0.0; fPtMu2 = 0.0;
	fMuID1 = 0; fMuID2 = 0;
	fDeltaR = 0.0;
	fEtaMu1 = 0.0; fEtaMu2 = 0.0;
	fTrackQual_mu1 = 0;
	fTrackQual_mu2 = 0;
	fQ_mu1 = 0;
	fQ_mu2 = 0;
} // clearVariables()

int mumuReader::loadCandidateVariables(TAnaCand *pCand)
{
	int result,type;
	bool firstMu = true;
	map<int,int> cand_tracks;
	TAnaTrack *sigTrack,*recTrack;
	TVector3 plabMu1,plabMu2;
	
	// default initialization
	fPtMu1 = fPtMu2 = 0.0f;
	fMuID1 = fMuID2 = 0;
	fDeltaR = 0.0f;
	fEtaMu1 = fEtaMu2 = 0.0f;
	fTrackQual_mu1 = fTrackQual_mu2 = 0;
	fQ_mu1 = fQ_mu2 = 0;
	fCtau = 0.0f;
	
	if (BLIND && 5.1 < pCand->fMass && pCand->fMass < 5.5) return 0;
	
	result = massReader::loadCandidateVariables(pCand);
	fCtau = kMassBs / pCand->fPlab.Mag() * fD3; // estimate the proper time!
	
	// set the momenta
	findAllTrackIndices(pCand,&cand_tracks);
	for (map<int,int>::const_iterator it = cand_tracks.begin(); it!=cand_tracks.end(); ++it) {
		
		sigTrack = fpEvt->getSigTrack(it->second);
		recTrack = fpEvt->getRecTrack(sigTrack->fIndex);
		type = abs(sigTrack->fMCID);
		
		switch (type) {
			case 13:
				// muons
				if (firstMu) {
					plabMu1 = sigTrack->fPlab;
					fPtMu1 = sigTrack->fPlab.Perp();
					fTrackQual_mu1 = recTrack->fTrackQuality;
					fQ_mu1 = recTrack->fQ;
					fMuID1 = recTrack->fMuID > 0 ? recTrack->fMuID : 0;
					fEtaMu1 = sigTrack->fPlab.Eta();
				} else {
					plabMu2 = sigTrack->fPlab;
					fPtMu2 = sigTrack->fPlab.Perp();
					fTrackQual_mu2 = recTrack->fTrackQuality;
					fQ_mu2 = recTrack->fQ;
					fMuID2 = recTrack->fMuID > 0 ? recTrack->fMuID : 0;
					fEtaMu2 = sigTrack->fPlab.Eta();
				}
				firstMu = false;
				break;
			default:
				break;
		}
	}
	
	// mu1 is with higher pt
	if (fPtMu1 < fPtMu2) {
		swap(fPtMu1,fPtMu2);
		swap(fTrackQual_mu1,fTrackQual_mu2);
		swap(fQ_mu1,fQ_mu2);
		swap(fMuID1,fMuID2);
		swap(fEtaMu1,fEtaMu2);
	}
	
	if (fPtMu2 > 0)
		fDeltaR = plabMu1.DeltaR(plabMu2);
	
	return result;
} // loadCandidateVariables()

int mumuReader::checkTruth(TAnaCand *pCand)
{
	int result;
	multiset<int> particles;
	TAnaTrack *track;
	TGenCand *gen;
	
	result = massReader::checkTruth(pCand);
	if (!result) goto bail;
	
	// check if the decay coincides...
	track = fpEvt->getSigTrack(pCand->fSig1);
	track = fpEvt->getRecTrack(track->fIndex);
	
	gen = fpEvt->getGenCand(track->fGenIndex);
	while (abs(gen->fID) != fTruthType) // truthtype set in constructor...
		gen = fpEvt->getGenCand(gen->fMom1);
	
	buildDecay(gen,&particles);
	particles.erase(22); // remove Bremsstrahlung
	
	result = (particles == trueDecay);
	
bail:
	return result;
} // checkTruth()

bool mumuReader::parseCut(char *cutName, float cutLow, float cutHigh, int dump)
{
	bool parsed;
	
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
	
	parsed = (strcmp(cutName,"OPPOSITE_SIGN_MU") == 0);
	if (parsed) {
		fCutOppSign_mu = (cutLow != 0.0f);
		if (dump) cout << "OPPOSITE_SIGN_MU: " << fCutOppSign_mu << endl;
		goto bail;
	}
	
	parsed = (strcmp(cutName,"MUID") == 0);
	if (parsed) {
		fCutMuID_mask = (int)cutLow;
		fCutMuID_reqall = (cutHigh != 0.0);
		if (dump) cout << "MUID: mask " << fCutMuID_mask << " and require all " << fCutMuID_reqall << endl;
		goto bail;
	}
	
	parsed = massReader::parseCut(cutName, cutLow, cutHigh, dump);
bail:
	return parsed;
} // parseCut()

bool mumuReader::applyCut()
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
	
	// check the Muon ID
	if (fCutMuID_mask != 0) {
		if (fCutMuID_reqall)
			pass = ((fMuID1 & fCutMuID_mask) == fCutMuID_mask) && ((fMuID2 & fCutMuID_mask) == fCutMuID_mask);
		else
			pass = ((fMuID1 & fCutMuID_mask) != 0) && ((fMuID2 & fCutMuID_mask) != 0);
		
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
	return pass;
} // applyCut()
