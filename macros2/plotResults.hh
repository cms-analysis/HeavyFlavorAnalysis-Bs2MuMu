#ifndef PLOTRESULTS
#define PLOTRESULTS

#include "plotClass.hh"

class plotResults: public plotClass {

public:

  plotResults(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotResults();

  void makeAll(int channels = 15);

  void computeNormUL();
  void computeCsBF();

  void acceptancePerProcess();
  void rareBg();
  void scaledHist(int mode = 0);

  void allInvertedIso();
  void histInvertedIso(const char *var, int n, double lo, double hi);
  TH1D* invertedIso(int chan, const char *cuts);
  void invertedIsoPrediction();

  void printUlcalcNumbers(std::string fname);
  void createAllCfgFiles(std::string fname);
  void printCsBFNumbers();
  double barlow(int nobs, double bg, double bgE, double sE);
  double scaledYield(numbers *a, numbers *no, double chanbf, double fsfu); 

  
  bool fNormProcessed; 
  double fBlExp, fBlExpE, fBlObs, fBlObsE;
  
  ClassDef(plotResults,1) //Testing plotResults


};


#endif

