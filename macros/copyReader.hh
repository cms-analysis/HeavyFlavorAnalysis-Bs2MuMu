/*
 *  copyReader.h
 *
 *  Created by Christoph on 31.3.10.
 *  Copyright 2010 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef COPY_READER_H
#define COPY_READER_H


#include "treeReader01.hh"

class copyReader : public treeReader01 {
  
public:
  copyReader(TChain *tree, TString evtClassName);
  ~copyReader();
  
  void bookHist();
  void eventProcessing();
  void readCuts(TString filename, int dump = 1);
  
  
private:
  TTree *copy_tree;
  TAna01Event *copy_evt;
  std::vector<int> fTypeList; 

  int fCopied; 

  int FILTERONMASS;
  double MASSLO, MASSHI; 
};

#endif
