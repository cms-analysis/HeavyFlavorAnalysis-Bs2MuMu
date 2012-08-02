/*
 *  pixelReader.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 30.07.12.
 *
 */

#include "pixelReader.hh"

pixelReader::pixelReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName), reduced_tree(NULL)
{} // pixelReader()

pixelReader::~pixelReader()
{} // ~pixelReader()
	
void pixelReader::bookHist()
{
	reduced_tree = new TTree("T","Pixel2012 Reduced Tree");
	
	reduced_tree->Branch("mass",&fMass,"mass/F");
	reduced_tree->Branch("pt",&fPt,"pt/F");
	reduced_tree->Branch("doca",&fDoca,"doca/F");
} // bookHist()

void pixelReader::eventProcessing()
{
	int j,nc;
	
	nc = fpEvt->nCands();
	for (j = 0; j < nc; j++) {
		if (loadCandidateVariables(fpEvt->getCand(j)))
			reduced_tree->Fill();
	}
} // eventProcessing()

void pixelReader::closeHistFile()
{
	fpHistFile = reduced_tree->GetCurrentFile();
	treeReader01::closeHistFile();
} // closeHistFile()

bool pixelReader::loadCandidateVariables(TAnaCand *pCand)
{
	// set variables
	
	fMass = pCand->fMass;
	fPt = pCand->fPlab.Perp();
	fDoca = (float)calcMaxDoca(pCand);
	
	return true;
} // loadCandidateVariables()

double pixelReader::calcMaxDoca(TAnaCand *pCand)
{
	TGenCand *g1,*g2;
	TAnaTrack *trk;
	int j,k;
	double doca;
	double result = -99.0;
	TVector3 q;
	
	for (j = pCand->fSig1; j < pCand->fSig2; j++) {
		
		trk = fpEvt->getSigTrack(j);
		g1 = fpEvt->getGenCand(trk->fGenIndex);
		for (k = j+1; k <= pCand->fSig2; k++) {
			trk = fpEvt->getSigTrack(k);
			g2 = fpEvt->getGenCand(trk->fGenIndex);
			
			// compute the distance
			q = g1->fP.Vect().Cross(g2->fP.Vect());
			doca = TMath::Abs((g1->fV - g2->fV) * q);
			doca /= q.Mag();
			
			if (doca > result)
				result = doca;
		}
	}
	
	return result;
} // calcDoca()
