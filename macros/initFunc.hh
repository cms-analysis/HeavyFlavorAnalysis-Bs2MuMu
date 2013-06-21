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

  void resetLimits(); 
  void limitPar(int ipar, double lo, double hi); 

  // -- background functions
  TF1* pol0(TH1 *h); 
  TF1* pol0(double lo, double hi); 

  TF1* pol1(double lo, double hi); 
  TF1* pol1(TH1 *h); 
  TF1* pol1(TH1 *h, double lo, double hi); 

  TF1* pol2local(TH1 *h, double width = 0.1); 

  TF1* expo(double lo, double hi); 
  TF1* expo(TH1 *h); 
  TF1* expo(TH1 *h, double lo, double hi); 

  TF1* pol1Err(double lo, double hi); 
  TF1* expoErr(double lo, double hi); 

  // -- signal+background functions
  TF1* pol1gauss(TH1 *h, double peak = 5.3, double sigma = 0.04); 
  TF1* pol1Gauss(TH1 *h, double peak = 5.3, double sigma = 0.04); 

  TF1* pol1gauss2c(TH1 *h, double peak = 5.3, double sigma = 0.04); 

  TF1* expoGauss(TH1 *h, double peak = 5.3, double sigma = 0.04);

  TF1* expoErrGauss(TH1 *h, double peak = 5.3, double sigma = 0.04, double preco = 5.14); 
  TF1* expoErrgauss2c(TH1 *h, double peak = 5.3, double sigma1 = 0.04, double sigma2 = 0.1, double preco = 5.14); 
  TF1* expoErrgauss2(TH1 *h, double peak = 5.3, double sigma = 0.04, double peak2 = 5.3, double sigma2 = 0.1, double preco = 5.14); 
  TF1* pol1ErrGauss(TH1 *h, double peak = 5.3, double sigma = 0.04, double preco = 5.14); 

  TF1* expoErrgauss2f(TH1 *h, double peak = 5.3, double sigma1 = 0.04, double peak2 = 5.425, double sigma2 = 0.079, double fraction = -1.,
		      double preco = -1.); 

  TF1* crystalBall(TH1 *h, double peak = 5.3, double sigma = 0.04, double alpha = 1., double tailLength = 1.);
  TF1* pol1CrystalBall(TH1 *h, double peak = 5.3, double sigma = 0.04, double alpha = 1., double tailLength = 1.); 

  TF1* pol0BsBlind(TH1 *h); 
  TF1* pol1BsBlind(TH1 *h); 
  TF1* expoBsBlind(TH1 *h); 

  // Added to include landaus
  TF1* land(TH1 *h, double peak=5.3, double sigma=0.04); 
  TF1* landsimp(TH1 *h, double peak=5.3, double sigma=0.04); 
  TF1* expoErrgauss2Landau(TH1 *h, double peak = 5.3, double sigma = 0.04, double peak2 = 5.3, double sigma2 = 0.1, double preco = 5.14); 
  TF1* expoErrGaussLandau(TH1 *h, double peak = 5.3, double sigma = 0.04, double preco = 5.14); 


  void initPol1(double &p0, double &p1, TH1 *h);
  void initExpo(double &p0, double &p1, TH1 *h);

  double fLo, fHi;
  double fBgFractionLo, fBgFractionHi; 
  double fLimitLo[20], fLimitHi[20]; 
  bool fLimit[20];

private: 

  ClassDef(initFunc,1) //Testing initFunc


}; 

#endif
