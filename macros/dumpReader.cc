#include "dumpReader.hh"

dumpReader::dumpReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName)
{
}

dumpReader::~dumpReader() {}

void dumpReader::eventProcessing()
{
	int j,ncand;
	TGenCand *cand;
	static int nbr = 0;
	
	cout << "Processing Event " << nbr++ << endl;
	ncand = fpEvt->nGenCands();
	for (j=0; j<ncand; j++) {
		cand = fpEvt->getGenCand(j);
		if(cand) {
			cout << "Cand " << j << " ID " << cand->fID << " Mom1 " << cand->fMom1 << " Mom2 " << cand->fMom2;
			cout << " Dau1 " << cand->fDau1 << " Dau2 " << cand->fDau2 << endl;
		} else {
			cerr << "Event returned invalid command" << endl;
		}
	}
}
