#ifndef PLOTMISC
#define PLOTMISC

#include <TH1.h>
#include <TProfile.h>

#include "plotClass.hh"


class plotMisc: public plotClass {

public:

  plotMisc(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default", int mode = 0);
  ~plotMisc();
  
  virtual void loopFunction(int function, int mode = 0); 
  virtual void loopFunction1(int mode); 
  virtual void loopFunction2(int mode);
  int          etaBin(double eta); 
 
  void makeAll();

  void signalMass();
  void pidTableDisplay();
  void fakeRateOverlaysDK(std::string mode = "pions");
  void fakeRateOverlaysMM(std::string mode = "pions");
  void fakeRateOverlaysMG(std::string mode = "pions", std::string charge = "pos", std::string chan = "barrel");
  void massError();

  void effImpactTrkHit();

  void calcTauError();
  void calcCombSlopes();
  int bdtCat(double bdt, int chan);
  void bdtCatIdx(int cat, int chan, int &hmin, int &hmax);


  // -- histograms
  std::map<std::string, TH1D*> fh0M0, fh0M3, fh0M4;


  ClassDef(plotMisc,1) //Testing plotMisc

};


#endif

