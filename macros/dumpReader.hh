#ifndef DUMP_READER_H
#define DUMP_READER_H

#include "treeReader01.hh"

#include <set>

using namespace std;

class dumpReader : public treeReader01 {
	
	public:
		dumpReader(TChain *tree, TString evtClassName);
		~dumpReader();
		
		void eventProcessing();
	private:
		void dumpCandidate(TAnaCand *pCand, unsigned indent = 0);
		void dumpGenerator(TGenCand *pGen);
		void dumpUnrecovered(TGenCand *pGen);
		
		void buildTracks(int genIx, set<int> *trackIndices, int ID = 0);
		int getGenIndex(TAnaCand *pCand);

		int checkReconstruction(int genIx, multiset<int> *particles, unsigned *nbrMuons);
		
		map<multiset<int>, int> decay_counter;
		multiset<int> true_channel;
		
		map<int,int> evt_nbr;
		map<int,int> rec_nbr;
		map<int,int> muon_nbr;
		
		unsigned unreco_counter;
};

#endif
