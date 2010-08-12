#ifndef DECAYCOUNTER_H
#define DECAYCOUNTER_H

#include "treeReader01.hh"
#include <set>

using std::multiset;

class decayCounter : public treeReader01 {
	
	public:
		explicit decayCounter(TChain *tree, TString evtClassName);
		virtual ~decayCounter();
		virtual void eventProcessing();
	
	private:
		unsigned long long counter;
		multiset<int> trueDecay;
		
		void buildDecay(TGenCand *gen, multiset<int> *particles);
};

#endif
