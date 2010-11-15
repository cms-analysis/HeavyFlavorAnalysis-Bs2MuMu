#ifndef OPTIMIZER
#define OPTIMIZER

#include "TString.h"
#include "TObject.h"
#include "TH1.h"

#include <iostream>
#include <fstream>

using namespace::std;


#define MAXCUTS 100

// ----------------------------------------------------------------------
// Optimizer
// ---------
// Class to optimize a cut given 
//  o 
//  o 
//  o 
//  o 
// ----------------------------------------------------------------------


class Optimizer: public TObject {

public:

  Optimizer(TH1D *hs, TH1D *hb); 
  ~Optimizer(); 

  TH1D* SoverB(int upperCut = 1);
  TH1D* S2overSplusB(int upperCut = 1);
  TH1D* EffS2overEffB(int upperCut = 1);

private: 
 
  TH1D *fpSignal, *fpBackground;


  ClassDef(Optimizer,1) //Testing Optimizer

}; 

#endif
