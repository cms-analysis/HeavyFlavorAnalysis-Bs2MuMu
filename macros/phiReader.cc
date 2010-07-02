#include "phiReader.hh"

#include <utility>

using namespace std;

phiReader::phiReader(TChain *tree, TString evtClassName) : massReader(tree,evtClassName)
{
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

int phiReader::loadCandidateVariables(TAnaCand *pCand)
{
	int result,type,j;
	int firstMu = 1, firstKp = 1;
	TAnaTrack *track;
	TAnaCand *dau;
	
	// default initialization
	fMassJPsi = -1.0;
	fMassPhi = -1.0;
	
	fPlabMu1 = TVector3();
	fPlabMu2 = TVector3();
	fPlabKp1 = TVector3();
	fPlabKp2 = TVector3();
	
	if (pCand->fType % 1000 != 531) return 0;
	
	result = massReader::loadCandidateVariables(pCand);
	
	// set the momenta
	for (j = pCand->fSig1; j < pCand->fSig2; j++) {
		
		track = fpEvt->getSigTrack(j);
		type = abs(track->fMCID);
		track = fpEvt->getRecTrack(track->fIndex);
		if (type == 13) {
			if (firstMu)	fPlabMu1 = track->fPlab;
			else			fPlabMu2 = track->fPlab;
		} else if (type == 321) {
			if (firstKp)	fPlabKp1 = track->fPlab;
			else			fPlabKp2 = track->fPlab;
		}
	}
	if (fPlabMu1.Perp() < fPlabMu2.Perp())
		swap(fPlabMu1,fPlabMu2);
	if (fPlabKp1.Perp() < fPlabKp2.Perp())
		swap(fPlabKp1,fPlabKp2);
	
	// set the masses
	for (j = pCand->fDau1; j <= pCand->fDau2 && j>=0; j++) {
		
		dau = fpEvt->getCand(j);
		if (dau->fType % 1000 == 443)
			fMassJPsi = dau->fMass;
		else if(dau->fType % 1000 == 333)
			fMassPhi = dau->fMass;
	}
	
	return result;
} // loadCandidateVariables()

void phiReader::bookHist()
{
	massReader::bookHist();
	reduced_tree->Branch("mass_jpsi",&fMassJPsi,"mass_jpsi/F");
	reduced_tree->Branch("mass_phi",&fMassPhi,"mass_phi/F");
	reduced_tree->Branch("plab_mu1","TVector3",&fPlabMu1Ptr);
	reduced_tree->Branch("plab_mu2","TVector3",&fPlabMu2Ptr);
	reduced_tree->Branch("plab_kp1","TVector3",&fPlabKp1Ptr);
	reduced_tree->Branch("plab_kp2","TVector3",&fPlabKp2Ptr);
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
	
bail:
	return result;
} // checkTruth()
