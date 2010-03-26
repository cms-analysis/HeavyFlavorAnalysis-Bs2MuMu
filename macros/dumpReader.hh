#ifndef DUMP_READER_H
#define DUMP_READER_H

#include "treeReader01.hh"

using namespace std;

class dumpReader : public treeReader01 {
	
	public:
		dumpReader(TChain *tree, TString evtClassName);
		~dumpReader();
		
		void eventProcessing();
};

#endif
