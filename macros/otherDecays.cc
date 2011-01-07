/*
 *  otherDecays.cc
 *  macros
 *
 *  Created by Christoph on 06.01.11.
 *  Copyright 2011 PSI. All rights reserved.
 *
 */

#include "otherDecays.hh"

otherDecays::otherDecays(TChain* tree, TString evtClassName) : phiReader(tree, evtClassName)
{ } // otherDecays()

int otherDecays::loadCandidateVariables(TAnaCand *pCand)
{
	using std::set;
	using std::multiset;
	using std::map;
	using std::cout;
	using std::endl;
	
	int result = phiReader::loadCandidateVariables(pCand);
	int nGens;
	multiset<int> particles;
	multiset<int>::const_iterator partIt;
	set<int> bhadrons;
	map<int,int> candTracks;
	map<int,int>::const_iterator it;
	TGenCand *pGen = NULL;
	TAnaTrack *pTrack;
	
	if(!result || !applyCut() || (fTruthFlags == 0)) goto bail;
	
	// survived? dump the mass together with the decay
	bhadrons.insert(511); // B0
	bhadrons.insert(521); // B+
	bhadrons.insert(531); // Bs
	bhadrons.insert(5122); // Lambda_b	
	
	findAllTrackIndices(pCand,&candTracks);
	
	// iterate through all tracks and see if they have the same origin
	nGens = fpEvt->nGenCands();
	
	for (it = candTracks.begin(); it!=candTracks.end(); ++it) {
		
		pTrack = fpEvt->getRecTrack(it->first);
		pGen = fpEvt->getGenCand(pTrack->fGenIndex);
		
		// this works as truth_flags != 0
		while (bhadrons.count(abs(pGen->fID)) == 0)
			pGen = fpEvt->getGenCand(pGen->fMom1);
		
		if(bhadrons.count(abs(pGen->fID)) > 0) break;
	}
	
	buildDecay(pGen,&particles);
	
	// dump to stdout...
	cout << "mass " << pCand->fMass;
	for (partIt = particles.begin(); partIt != particles.end(); ++partIt) {
		cout << " " << *partIt;
	}
	cout << endl;
	
bail:
	return result;
} // loadCandidateVariables()
