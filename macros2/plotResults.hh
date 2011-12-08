#ifndef PLOTRESULTS
#define PLOTRESULTS

#include "plotClass.hh"

class plotResults: public plotClass {

public:

  plotResults(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotResults();

  void makeAll(int channels = 3);

  void computeNormUL();
  void computeCsBF();
  
  void rareBg();
  void printUlcalcNumbers();
  void printCsBFNumbers();
  double barlow(int nobs, double bg, double bgE, double sE);
  
  std::string fUlcalcFileName; 
  
  ClassDef(plotResults,1) //Testing plotResults

};


#endif

