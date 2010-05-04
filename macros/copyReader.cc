/*
 *  copyReader.cc
 *  massReader
 *
 *  Created by Christoph on 31.3.10.
 *  Copyright 2010 Christoph Naegeli. All rights reserved.
 *
 */

#include "copyReader.hh"


// ----------------------------------------------------------------------
copyReader::copyReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName), copy_tree(NULL), fCopied(0) { 
  fTypeList.clear();
} 

// ----------------------------------------------------------------------
copyReader::~copyReader() { 
  cout << "Copied a total of " << fCopied << " events" << endl;
} 


// ----------------------------------------------------------------------
void copyReader::bookHist() {
  // create the copy tree...
  fpHistFile->cd();
  copy_tree = new TTree("T1","Events");
  
  copy_evt = new TAna01Event(0);
  copy_tree->Branch("TAna01Event","TAna01Event",&copy_evt, 256000/8, 1);
} 


// ----------------------------------------------------------------------
void copyReader::eventProcessing() {
  bool doCopy(false); 
  TAnaCand *pCand(0); 
  int type(-1); 
  for (int j = 0; j < fpEvt->nCands(); ++j) {
    pCand = fpEvt->getCand(j);
    for (unsigned int ic = 0; ic < fTypeList.size(); ++ic) {
      if (pCand->fType == fTypeList[ic]) {
	doCopy = true; 
	type = pCand->fType;
	break;
      }
    }
  }

  if (doCopy) {
    ++fCopied;
    cout << "Found a candidate with type " << type << ", copying, among ncands = " << fpEvt->nCands() << endl;
    *copy_evt = *fpEvt;
    copy_tree->Fill();
  }

}


// ----------------------------------------------------------------------
void copyReader::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> copyReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);
  
  TString fn(fCutFile.Data());
  
  if (dump) {
    cout << "====================================" << endl;
    cout << "==> copyReader: Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }
  
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
    
    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      fTypeList.push_back(CutValue); 
      if (dump) cout << " added TYPE:      " << TYPE << endl;
    }
  }

  cout << "Filtering on the following candidates" << endl;
  for (unsigned int i = 0; i < fTypeList.size(); ++i) {
    cout << Form("%3i %i", i, fTypeList[i]) << endl;
  }
}
