/*
 *  runMergeTrees.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 25.08.10.
 *
 */

// C++ headers
#include <cstring>
#include <cstdio>
#include <iostream>

// ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TUnixSystem.h>

// HeavyFlavorObjects headers
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna01Event.hh"

using std::cout;
using std::endl;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: ./runMergeTrees -c <chain_file> -o <outputfile>
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

static void dump_help_message()
{
	std::cerr << "usage: runMergeTrees -c <chain_file> -o <outputfile>" << endl;
} // dump_help_message()

int main(int argc, const char *argv [])
{
	TAna01Event *evt = new TAna01Event;
	int processID = gSystem->GetPid();
	const char *chainfile = NULL;
	const char *outputfile = NULL;
	TChain *chain = NULL;
	TTree *tree = NULL;
	TFile* root_out = NULL;
	FILE *fchain = NULL;
	char buffer[1024];
	TString line;
	int j;
	
	cout << "Running under process ID  " << processID << endl;
	
	// split the tree into 4 GB chunks
	TTree::SetMaxTreeSize(4000000000ll); // 4 GB
	
	// parse the command line
	for (j = 1; j < argc; j++) {
		
		if (!strcmp(argv[j],"-c")) chainfile = argv[++j];
		else if (!strcmp(argv[j],"-o")) outputfile = argv[++j];
		else if (!strcmp(argv[j],"-h")) {
			dump_help_message();
			goto bail;
		}
		else {
			std::cerr << "Unknown argument '" << argv[j] << "'" << endl;
			dump_help_message();
			goto bail;
		}
	}
	
	if (!chainfile || !outputfile) {
		dump_help_message();
		goto bail;
	}
	
	cout << "chain file " << chainfile << endl;
	cout << "output file " << outputfile << endl;
	cout << "Setup the chain..." << endl;
	chain = new TChain("T1");
	fchain = fopen(chainfile,"r");
	while (fgets(buffer,sizeof(buffer),fchain)) {
		line.Append(buffer);
		if (line[line.Length()-1] != '\n')
			continue;
		
		line.Remove(line.Length()-1);
		
		cout << "Adding file: " << line.Data() << endl;
		chain->Add(line);
		
		line = TString(); // empty string
	}
	
	cout << "Chain has a total size of " << chain->GetEntries() << endl;
	chain->SetBranchAddress("TAna01Event",&evt); // we need to set this buffer (?), else crash
	
	// open the output file
	root_out = new TFile(outputfile,"recreate");
	root_out->cd();
	
	// copy the tree
	tree = chain->CloneTree();
	root_out->Write();
	
	cout << "Copied " << tree->GetEntries() << " number of entries" << endl;
	
bail:
	if (fchain) fclose(fchain);
	delete chain;
	delete root_out;
	delete evt;
	
	return 0;
} // main()
