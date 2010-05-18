#include "ksReader.hh"

using namespace std;

ksReader::ksReader(TChain *tree, TString evtClassName) : massReader(tree,evtClassName)
{
	fPlabMu1Ptr = &fPlabMu1;
	fPlabMu2Ptr = &fPlabMu2;
	fPlabPi1Ptr = &fPlabPi1;
	fPlabPi2Ptr = &fPlabPi2;
} // ksReader()

ksReader::~ksReader() {}

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
	if (!ksCand) goto jpsi;
	
	fMassKs = ksCand->fMass;
	fAlphaKs = ksCand->fPlab.Angle(ksCand->fVtx.fPoint - pCand->fVtx.fPoint);
	fDxyKs = ksCand->fVtx.fDxy;
	fDxyeKs = ksCand->fVtx.fDxyE;
	fChi2Ks = ksCand->fVtx.fChi2;
	
jpsi:
	// jpsi subvariables
	ksCand = NULL;
	for (j = pCand->fDau1; j <= pCand->fDau2; j++) {
		
		TAnaCand *tmpCand = fpEvt->getCand(j);
		if (tmpCand->fType % 1000 == 443) { // J/Psi
			ksCand = tmpCand;
			break;
		}
	}
	if(!ksCand) goto bail;
	fMassJPsi = ksCand->fMass;
	
	result = 1;
	
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
