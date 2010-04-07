/*
 *  copyReader.cc
 *  massReader
 *
 *  Created by Christoph on 31.3.10.
 *  Copyright 2010 Christoph Nägeli. All rights reserved.
 *
 */

#include "copyReader.hh"

copyReader::copyReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName), copy_tree(NULL)
{ } // copyReader()

void copyReader::bookHist()
{
	// create the copy tree...
	copy_tree = new TTree("T","Events");
	
	copy_evt = new TAna01Event(0);
	copy_tree->Branch("TAna01Event","TAna01Event",&copy_evt);
} // bookHist()

void copyReader::eventProcessing()
{
	*copy_evt = *fpEvt;
	copy_tree->Fill();
} // eventProcessing()
