#ifndef PLOTOPTIMIZE
#define PLOTOPTIMIZE

#include "TTree.h"
#include "plotClass.hh"

class plotOptimize: public plotClass {

public:

  plotOptimize(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotOptimize();

  void makeAll(int nfiles = 50, int mode = 0);
  void optimizeULs(int nruns, int seed);
  void optimizeBdtULs(double minBdt, double maxBdt); 

  void bestUL(const char *fname, int mode);
  void readOptimize(int nfiles = 50, const char *fname="optimizeUL");
  void readFile(const char *fname, TTree *t);
  void recalcUL(int mode, double &ul0, double &ul1); 
  void displayCuts(const char *fname);

  int _chan, _file, _run, _closetrk, _cowboyVeto; 
  float _mlo, _mhi, _pt, _m1pt, _m2pt, _iso, _chi2dof, _alpha, _fls3d, _docatrk, _pvlip, _pvlips, _pvlip2, _pvlips2, _pvip, _pvips, 
    _maxdoca, _bdt; 
  float _ul, _ulC, _ulCP, _nobs, _nhlo, _nhhi, _nexp, _eff, _sig, _ssb, _ssb0, _ssb1, _ssb2; 
  
  ClassDef(plotOptimize,1) //Testing plotOptimize

};


#endif

