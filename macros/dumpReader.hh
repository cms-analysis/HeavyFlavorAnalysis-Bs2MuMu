#ifndef DUMP_READER_H
#define DUMP_READER_H

#include "treeReader01.hh"

#include <set>

class dumpReader : public treeReader01 {
	
	public:
		dumpReader(TChain *tree, TString evtClassName);
		~dumpReader();
		
		void eventProcessing();
	private:
		void dumpCandidate(TAnaCand *pCand, unsigned indent = 0);
		void dumpGenerator(TGenCand *pGen);
		void dumpUnrecovered(TGenCand *pGen);
		
		void buildTracks(int genIx, std::set<int> *trackIndices, int ID = 0);
		int getGenIndex(TAnaCand *pCand);

		int checkReconstruction(int genIx, std::multiset<int> *particles, unsigned *nbrMuons);
		
		std::map<std::multiset<int>, int> decay_counter;
		std::multiset<int> true_channel;
		
		std::map<int,int> evt_nbr;
		std::map<int,int> rec_nbr;
		std::map<int,int> muon_nbr;
		
		unsigned unreco_counter;
};

#endif
