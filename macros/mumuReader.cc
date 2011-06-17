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
	massReader(tree, evtClassName)
{
	cout << "Instantiating mumuReader..." << flush;
	fTreeName = "mumuReader reduced tree.";
	fTruthType = 531;
	
	trueDecay.insert(13); // mu
	trueDecay.insert(13); // mu
	trueDecay.insert(fTruthType); // Bs
	
	cout << " ok" << endl;
} // mumuReader()

int mumuReader::loadCandidateVariables(TAnaCand *pCand)
{
	if (BLIND && 5.2 < pCand->fMass && pCand->fMass < 5.45) return 0;
	
	return massReader::loadCandidateVariables(pCand);
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
