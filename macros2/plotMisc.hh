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
  void fakeRateOverlays(std::string mode = "pions");
  void massError();



  // -- histograms
  std::map<std::string, TH1D*> fh0M0, fh0M3, fh0M4;


  ClassDef(plotMisc,1) //Testing plotMisc

};


#endif

