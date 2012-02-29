#include "ksReader.hh"

#include <utility>

using namespace std;

ksReader::ksReader(TChain *tree, TString evtClassName) : massReader(tree,evtClassName),total_counter(0),reco_counter(0)
{
	fPlabMu1Ptr = &fPlabMu1;
	fPlabMu2Ptr = &fPlabMu2;
	fPlabPi1Ptr = &fPlabPi1;
	fPlabPi2Ptr = &fPlabPi2;
	fPlabKsPtr = &fPlabKs;
	fPlabKsTruePtr = &fPlabKsTrue;
	
	fTreeName = "ksReader reduced tree";
		
	cout << "ksReader instantiated..." << endl;
} // ksReader()

ksReader::~ksReader()
{
	cout << "~ksReader:" << endl;
	cout << "\treco counter: " << reco_counter << endl;
	cout << "\ttotal counter: " << total_counter << endl;
} // ~ksReader()

int ksReader::loadCandidateVariables(TAnaCand *pCand)
{
	TAnaCand *subCand = NULL;
	TAnaTrack *track;
	TGenCand *ksTrue,*bdTrue;
	int result;
	int j;
	int firstMu = 1;
	int firstPi = 1;
	int type;
	
	// default initialization
	fMassJPsi = -1.0;
	fMassKs = -1.0;
	fAlphaKs = -1.0;
	fDxyKs = -1.0;
	fDxyeKs = -1.0;
	fChi2Ks = -1.0;
	fDeltaKs = -1.0;
	fDeltaKsTrue = -1.0;
	fMaxDocaKs = -1.0;

	fPlabMu1 = TVector3();
	fPlabMu2 = TVector3();
	fPlabPi1 = TVector3();
	fPlabPi2 = TVector3();
	fPlabKs = TVector3();
	fPlabKsTrue = TVector3();

	// we handle only the Bd decays.
	if (pCand->fType % 1000 != 511) return 0;
	
	result = massReader::loadCandidateVariables(pCand);
	
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
	// mu1 has higher pt, p1 has higher pt
	if (fPlabMu1.Perp() < fPlabMu2.Perp())
		swap(fPlabMu1,fPlabMu2);
	if (fPlabPi1.Perp() < fPlabPi2.Perp())
		swap(fPlabPi1,fPlabPi2);
	
	// Ks subvariables
	for (j = pCand->fDau1; j <= pCand->fDau2 && j>= 0; j++) {
		
		TAnaCand *tmpCand = fpEvt->getCand(j);
		if(tmpCand->fType % 1000 == 310) { // ks
			subCand = tmpCand;
			break;
		}
	}
	if (subCand) {
		fMassKs = subCand->fMass;
		fAlphaKs = subCand->fPlab.Angle(subCand->fVtx.fPoint - pCand->fVtx.fPoint);
		fDxyKs = subCand->fVtx.fDxy;
		fDxyeKs = subCand->fVtx.fDxyE;
		fChi2Ks = subCand->fVtx.fChi2;
		fPlabKs = subCand->fPlab;
		fMaxDocaKs = subCand->fMaxDoca;
		
		// check the truth distance...
		ksTrue = bdTrue = NULL;
		j = getGenIndex(subCand, 310);
		if (j >= 0) ksTrue = fpEvt->getGenCand(j);
		j = getGenIndex(pCand, 511);
		if (j >= 0) bdTrue = fpEvt->getGenCand(j);
		
		// do the difference...
		if (ksTrue && bdTrue) {
			
			TVector3 diff = subCand->fVtx.fPoint - pCand->fVtx.fPoint;
			fDeltaKs = diff.Mag();
			
			diff = ksTrue->fV - bdTrue->fV;
			fDeltaKsTrue = diff.Mag();
			
			fPlabKsTrue = ksTrue->fP.Vect();
		}
	}
	
	// jpsi subvariables
	subCand = NULL;
	for (j = pCand->fDau1; j >= 0 && j <= pCand->fDau2; j++) {
		
		TAnaCand *tmpCand = fpEvt->getCand(j);
		if (tmpCand->fType % 1000 == 443) { // J/Psi
			subCand = tmpCand;
			break;
		}
	}
	if (subCand) {
		fMassJPsi = subCand->fMass;
	}
	
	return result;
} // ksReader()

void ksReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("mass_jpsi",&fMassJPsi,"mass_jpsi/F");
	reduced_tree->Branch("alpha_ks",&fAlphaKs,"alpha_ks/F");
	reduced_tree->Branch("mass_ks",&fMassKs,"mass_ks/F");
	reduced_tree->Branch("dxy_ks",&fDxyKs,"dxy_ks/F");
	reduced_tree->Branch("dxye_ks",&fDxyeKs,"dxye_ks/F");
	reduced_tree->Branch("chi2_ks",&fChi2Ks,"chi2_ks/F");
	reduced_tree->Branch("max_doca_ks",&fMaxDocaKs,"max_doca_ks/F");
	reduced_tree->Branch("plab_mu1","TVector3",&fPlabMu1Ptr);
	reduced_tree->Branch("plab_mu2","TVector3",&fPlabMu2Ptr);
	reduced_tree->Branch("plab_pi1","TVector3",&fPlabPi1Ptr);
	reduced_tree->Branch("plab_pi2","TVector3",&fPlabPi2Ptr);
	reduced_tree->Branch("delta_ks",&fDeltaKs,"delta_ks/F");
	reduced_tree->Branch("delta_ks_true",&fDeltaKsTrue,"delta_ks_true/F");
	reduced_tree->Branch("plab_ks","TVector3",&fPlabKsPtr);
	reduced_tree->Branch("plab_ks_true","TVector3",&fPlabKsTruePtr);
} // bookHist()

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
