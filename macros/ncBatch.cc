/*
 *  ncBatch.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 01.05.12.
 */

#include "ncMVA.hh"

// Standard headers
#include <string>
#include <iostream>

// ROOT headers
#include <TFile.h>
#include <TCut.h>

static void usage()
{
	std::cout << "usage: ncBatch <MLP options> <signalFile> <bkgFile>" << std::endl;
	abort();
} // usage()

// 1st argument MLP options
// 2st argument signal root file
// 3nd argument background root file
int main(int argc, char *argv [])
{
	TFile *sigFile = NULL,*bkgFile = NULL;
	TTree *sigTree, *bkgTree;
	double bkgWeight = 0.200000;
	double signalWeight = 0.000156;
	TCut barrelCut("TMath::Abs(eta_mu1) < 1.4 && TMath::Abs(eta_mu2) < 1.4");
	std::string methodOptions;
	
	if (argc < 4)
		usage();
	
	methodOptions = std::string(argv[1]);
	sigFile = new TFile(argv[2]);
	bkgFile = new TFile(argv[3]);
	
	sigTree = (TTree*)sigFile->Get("T");
	bkgTree = (TTree*)bkgFile->Get("T");
	
	std::cout << "Training with options: '" << methodOptions << "'" << std::endl;
	ncRunTraining(sigTree, signalWeight, bkgTree, bkgWeight, 0, barrelCut, &methodOptions);
	ncRunTraining(sigTree, signalWeight, bkgTree, bkgWeight, 1, !barrelCut, &methodOptions);
	
	delete sigFile;
	delete bkgFile;
	
	return 0;
} // main()
