/*
 *  decayCounter.cpp
 *  macros
 *
 *  Created by Christoph on 12.08.10.
 *  Copyright 2010 PSI. All rights reserved.
 *
 */

#include "decayCounter.hh"

decayCounter::decayCounter(TChain *tree, TString evtClassName) : treeReader01(tree,evtClassName),counter(0)
{
	trueDecay.insert(13); // mu
	trueDecay.insert(13); // mu
	trueDecay.insert(321); // kp
	trueDecay.insert(443); // J/Psi
	trueDecay.insert(521); // B+
} // decayCounter()

decayCounter::~decayCounter()
{
	using std::cout;
	using std::endl;
	
	cout << "===> decayCounter destructor!" << endl;
	cout << "Found counter = " << counter << endl;
} // ~decayCounter()

void decayCounter::eventProcessing()
{
	int nGens,j;
	TGenCand *gen;
	multiset<int> genDecay;
	
	// count the number of events with a trueDecay...
	nGens = fpEvt->nGenCands();
	for (j = 0; j < fpEvt->nGenCands(); j++) {
		gen = fpEvt->getGenCand(j);
		
		if (abs(gen->fID) == 521) { // B+
			genDecay.clear();
			buildDecay(gen,&genDecay);
			if (trueDecay == genDecay) {
				counter++;
				break;
			}
		}
	}
} // eventProcessing()

void decayCounter::buildDecay(TGenCand *gen, multiset<int> *particles)
{
	particles->insert(abs(gen->fID));
	for (int j = gen->fDau1; j <= gen->fDau2; j++)
		buildDecay(fpEvt->getGenCand(j),particles);
} // buildDecay()
