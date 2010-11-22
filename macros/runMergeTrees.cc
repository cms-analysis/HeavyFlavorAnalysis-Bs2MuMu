/*
 *  runMergeTrees.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 25.08.10.
 *
 */

// C++ headers
#include <cstring>
#include <cstdio>
#include <cstdlib>
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
// %% Usage: ./runMergeTrees -c <chain_file> -o <outputfile> [-s <maxfilesize[GB]>]
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

static void dump_help_message()
{
	std::cerr << "usage: runMergeTrees -c <chain_file> -o <outputfile> [-s <maxfilesize[MB]> -t <treename>]" << endl;
} // dump_help_message()

int main(int argc, const char *argv [])
{
	TAna01Event *evt = new TAna01Event;
	int processID = gSystem->GetPid();
	const char *chainfile = NULL;
	const char *outputfile = NULL;
	const char *treename = "T1";
	TChain *chain = NULL;
	TTree *tree = NULL;
	TFile* root_out = NULL;
	FILE *fchain = NULL;
	char buffer[1024];
	TString line;
	int j;
	long long maxTreeSize = 4000000000ll;
	
	cout << "Running under process ID  " << processID << endl;
	
	// parse the command line
	for (j = 1; j < argc; j++) {
		
		if (!strcmp(argv[j],"-c")) chainfile = argv[++j];
		else if (!strcmp(argv[j],"-o")) outputfile = argv[++j];
		else if (!strcmp(argv[j],"-s")) maxTreeSize = 1000000ll * (long long)atoi(argv[++j]);
		else if (!strcmp(argv[j],"-t")) treename = argv[++j];
		else if (!strcmp(argv[j],"-h")) {
			dump_help_message();
			goto bail;
		} else {
			std::cerr << "Unknown argument '" << argv[j] << "'" << endl;
			dump_help_message();
			goto bail;
		}
	}
	
	if (!chainfile || !outputfile) {
		dump_help_message();
		goto bail;
	}
	
	cout << "chain file: " << chainfile << endl;
	cout << "output file: " << outputfile << endl;
	cout << "max tree size: " << maxTreeSize << endl;
	
	// set the max tree size
	TTree::SetMaxTreeSize(maxTreeSize);
	
	cout << "Setup the chain..." << endl;
	chain = new TChain(treename);
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
	root_out = tree->GetCurrentFile();
	root_out->Write();
bail:
	if (fchain) fclose(fchain);
	delete chain;
	delete root_out;
	delete evt;
	
	return 0;
} // main()
