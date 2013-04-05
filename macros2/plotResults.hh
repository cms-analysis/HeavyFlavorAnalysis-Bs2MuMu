#ifndef PLOTRESULTS
#define PLOTRESULTS

#include "plotClass.hh"

class plotResults: public plotClass {

public:

  plotResults(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default", int mode = 0);
  ~plotResults();

  void makeAll(int channels = 1, int nevents = -1);

  void setupNorm();
  void fitHists(int chan = 0);

  void resetHistograms();
  void fillAndSaveHistograms(int nevents = -1); 
  void rareBgHists(std::string smode = "nada", int nevents = -1);
  void otherNumbers(std::string smode);
  void saveHists(std::string smode);
  
  void testAccEff(string smode);
  void numbersAfterLoopOverTree(int chan, int mode, numbers *aa, std::string directory);

  void play1(int mode); 
  void play2(int mode); 
  void play3(int mode); 
  void play4(int mode); 
  void play5(int mode); 

  void calculateNoNumbers(int chan, int mode = 1);
  void calculateCsNumbers(int chan, int mode = 1);
  void calculateNumbers(int mode); 
  void calculateRareBgNumbers(int chan);
  void calculateSgNumbers(int chan);
  void numbersFromHist(int chan, int mode, numbers *aa); 

  void saveLargerTree(std::string mode); 

  virtual void loopFunction(int function, int mode = 0); 
  virtual void loopFunction1(int mode); 
  virtual void loopFunction2(int mode); 

  void acceptancePerProcess();
  void scaledHist(int mode = 0);

  void fls3dEfficiency(std::string cuts, std::string pdfname);
  void fls3dVsX(std::string x, std::string cuts, std::string pdfname);

  void computeErrors(std::vector<numbers*>); 
  void printUlcalcNumbers(std::string fname);
  void createAllCfgFiles(std::string fname);
  void printCsBFNumbers();
  double scaledYield(numbers *a, numbers *no, std::string chan, double fsfu); 
  
  double fBlExp, fBlExpE, fBlObs, fBlObsE;

  double fBl0Exp, fBl0ExpE, fBl0Obs, fBl0ObsE;
  double fBl1Exp, fBl1ExpE, fBl1Obs, fBl1ObsE;

  
  ClassDef(plotResults,1) //Testing plotResults


};


#endif

