/*
 *  ncMVA.hh
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.04.12.
 */

#ifndef __NCMVA__
#define __NCMVA__

#include <string>
#include <TTree.h>
#include <TCut.h>

TTree* ncEvalAll(TTree *tree, bool verbose = true);
void ncEvalAll(const char *treeFileName, bool verbose = true);

void ncRunTraining(TTree *signalTree, double signalWeight, TTree *bkgTree, double bkgWeight, int channelIx, TCut channelCut, std::string methodOptions);
void ncRunDefaultTraining(const char *mcFile, const char *dataFile);

#endif
