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

  void bestUL(const char *fname, int mode);
  void readOptimize(int nfiles = 50);
  void readFile(const char *fname, TTree *t);

  int _chan, _file, _run, _closetrk, _cowboyVeto; 
  float _mlo, _mhi, _pt, _m1pt, _m2pt, _iso, _chi2dof, _alpha, _fls3d, _docatrk, _pvlip, _pvlips; 
  float _ul, _nobs, _nexp, _eff, _sig, _ssb, _ssb0, _ssb1, _ssb2; 
  
  ClassDef(plotOptimize,1) //Testing plotOptimize

};


#endif

