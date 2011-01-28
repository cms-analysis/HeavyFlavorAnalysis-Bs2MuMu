#ifndef INITFUNC
#define INITFUNC

#include "TString.h"
#include "TObject.h"
#include "TH1.h"
#include "TF1.h"

#include <iostream>


class initFunc: public TObject {

public:

  initFunc(); 
  ~initFunc(); 

  TF1* pol0(TH1 *h); 
  TF1* pol1(TH1 *h); 

  TF1* pol1gauss(TH1 *h); 
  TF1* pol1Gauss(TH1 *h); 


  TF1* pol0BsBlind(TH1 *h); 
  TF1* pol1BsBlind(TH1 *h); 


  void initPol1(double &p0, double &p1, TH1 *h);

  double fLo, fHi;

private: 

  ClassDef(initFunc,1) //Testing initFunc


}; 

#endif
