#include "massReader.hh"

massReader::massReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName),reduced_tree(NULL)
{
  fMomentumPtr = &fMomentum;
} // massReader()

massReader::~massReader()
{
	if (reduced_tree) {
		reduced_tree->Write();
		delete reduced_tree;
	}
} // ~massReader()

void massReader::eventProcessing()
{
	int j,nc;
	TAnaCand *pCand;
	
	// Fill a reduced tree with the following types
	//	- momentum
	//	- mass
	//	- type
	nc = fpEvt->nCands();
	for (j=0; j<nc; j++) {
		pCand = fpEvt->getCand(j);
		
		// Save in the tree
		fCandidate = pCand->fType;
		fMomentum = pCand->fPlab;
		fMass = pCand->fMass;
		reduced_tree->Fill();		
	}
} // massReader::eventProcessing()

void massReader::bookHist()
{
	// create the tree
	reduced_tree = new TTree("T","Candidate Mass / Eta / Charge");
	
	// and add the branches
	reduced_tree->Branch("candidate",&fCandidate,"candidate/I");
	reduced_tree->Branch("p","TVector3",&fMomentumPtr);
	reduced_tree->Branch("mass",&fMass,"mass");
} // massReader::bookHist()
