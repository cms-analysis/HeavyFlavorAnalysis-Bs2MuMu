#ifndef ANALYSISCUTS
#define ANALYSISCUTS

#include "TString.h"
#include "TObject.h"

#include <iostream>
#include <fstream>


#define MAXCUTS 100

using namespace::std;

// ----------------------------------------------------------------------
// AnalysisCuts
// ------------
// Class to organize cut sequences in an analysis. It allows to 
// ask the object whether a given cut is fulfilled, for cumulative cut 
// fulfilment (in the order of definition), and whether 
// "all other cuts" are fulfilled. 
//
// You can access cuts by index (good for loops) or name (good for 
// less confusion). 
//
// Don't forget to call update() before using it (i.e. after the cuts
// have been computed in your analysis)
//
// Bugs: Should probably be inlined to speed. But this seems not to work 
//       in CINT (???)
// ----------------------------------------------------------------------


class AnalysisCuts: public TObject {

public:

  AnalysisCuts(); 
  ~AnalysisCuts(); 

  int getIndex(const char *name); 
  const char* getName(int i); 
  const char* getDescription(int i); 

  void addCut(const char *name, int &location); 
  void addCut(const char *name, const char *description, int &location); 

  int singleCutTrue(int i); 
  int singleCutTrue(const char *name); 

  int cumulativeCutTrue(int i); 
  int cumulativeCutTrue(const char *name); 
  
  int allOtherCutsTrue(int i); 
  int allOtherCutsTrue(const char *name); 

  void update(); 
  
private: 
  
  int fNcuts; 
  int fUpdated; 
  TString fCutName[MAXCUTS]; 
  TString fDescription[MAXCUTS]; 
  int     fCutValue[MAXCUTS]; 
  int     fAocValue[MAXCUTS]; 
  int     fCumValue[MAXCUTS]; 
  
  int     *fCutLocation[MAXCUTS]; 

  ClassDef(AnalysisCuts,1) //Testing AnalysisCuts

}; 

#endif
