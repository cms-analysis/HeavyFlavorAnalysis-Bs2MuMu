/*
 *  ncBatch.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 01.05.12.
 */

#include "ncMVA.h"

// Standard headers
#include <string>
#include <iostream>
#include <stdlib.h>

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
	ncMVA mva;
	std::string methodOptions;
	
	if (argc < 4)
		usage();
	
	methodOptions = std::string(argv[1]);
	mva.setMVAOpts(0,methodOptions);
	mva.setMVAOpts(1,methodOptions);
	mva.setMCPath(argv[2]);
	mva.setDataPath(argv[3]);
	
	std::cout << "Training with options: '" << methodOptions << "'" << std::endl;
	mva.runAllTrainings();
	
	return 0;
} // main()
