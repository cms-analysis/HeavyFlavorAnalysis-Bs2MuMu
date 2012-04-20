/*
 *  ncMVA.hh
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.04.12.
 */

#ifndef __NCMVA__
#define __NCMVA__

#include <TTree.h>

TTree* ncEvalAll(TTree *tree, bool verbose = true);
void ncEvalAll(const char *treeFileName, bool verbose = true);

void ncRunTraining(TTree *signalTree, double signalWeight, TTree *bkgTree, double bkgWeight);
void ncRunDefaultTraining(const char *mcFile, const char *dataFile);

#endif
