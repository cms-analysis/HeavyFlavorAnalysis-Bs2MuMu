#ifndef PLOTRESULTS
#define PLOTRESULTS

#include "plotClass.hh"

class plotResults: public plotClass {

public:

  plotResults(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default", int mode = 11);
  ~plotResults();

  void makeAll(int channels = 3);

  void computeNormUL();
  void computeCsBF();

  void setupNorm();

  void acceptancePerProcess();
  void rareBg(std::string mode = "nada");
  void scaledHist(int mode = 0);

  void fls3dEfficiency(std::string cuts, std::string pdfname);
  void fls3dVsX(std::string x, std::string cuts, std::string pdfname);

  // -- new:
  void invertedIsolationStudy();
  void determineInvertedIsolationYield(int print = 0);
  void plotInvertedIsolationScan(std::string pdfname, TH1D *h0, TH1D *h1, TH1D *x0, TH1D *x1);
  // -- old:
  void allInvertedIso();
  void histInvertedIso(const char *var, int n, double lo, double hi);
  TH1D* invertedIso(int chan, const char *cuts);
  void invertedIsoPrediction();
  std::pair<TH1D*, TH1D*> singleRelativeYield(std::string fstring);

  void computeErrors(std::vector<numbers*>); 
  void printUlcalcNumbers(std::string fname);
  void createAllCfgFiles(std::string fname);
  void printCsBFNumbers();
  double barlow(int nobs, double bg, double bgE, double sE);
  double scaledYield(numbers *a, numbers *no, double chanbf, double fsfu); 

  
  double fBlExp, fBlExpE, fBlObs, fBlObsE;

  double fBl0Exp, fBl0ExpE, fBl0Obs, fBl0ObsE;
  double fBl1Exp, fBl1ExpE, fBl1Obs, fBl1ObsE;

  
  ClassDef(plotResults,1) //Testing plotResults


};


#endif

