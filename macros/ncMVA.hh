/*
 *  ncMVA.hh
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 11.04.12.
 */

#ifndef __NCMVA__
#define __NCMVA__

#include <TTree.h>

TTree* ncEvalAll(TTree *tree);
void ncEvalAll(const char *treeFileName);

void ncRunTraining(TTree *signalTree, double signalWeight, TTree *bkgTree, double bkgWeight);
void ncRunDefaultTraining(const char *mcFile, const char *dataFile);

#endif
