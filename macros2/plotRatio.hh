#ifndef PLOTRATIO
#define PLOTRATIO

#include "plotClass.hh"

class plotRatio: public plotClass {

public:

  plotRatio(const char *files="anaBmm.v11.files", const char *dir = "default", const char *cuts = "default", int mode = 11);
  ~plotRatio();

  void makeAll(int channels = 3);

  void computeCsNoRatio();
  void printCsBFNumbers();
  void computeCsNoRatioNew();
  
  virtual void loopFunction(int mode); 
  virtual void loopFunction(int mode, int imode); 
  virtual void loopFunction1(int imode); 
  virtual void loopFunction2(int imode);


//  std::vector<int> accNoNum, accCsNum, accNoDen, accCsDen;
//  std::vector<int> acc2NoNum, acc2CsNum, acc2NoDen, acc2CsDen;
  std::vector<int> effNoNum, effCsNum, effNoDen, effCsDen;
  std::vector<TH1D*> hNoSignal, hCsSignal;
  std::vector<TH1D*> hNoSignalC, hCsSignalC;
  std::vector<TH1D*> hNoPt, hCsPt;

  unsigned int fNPtbins;
  std::vector< std::vector<int> > accNoNumBin, accCsNumBin, accNoDenBin, accCsDenBin;
  std::vector< std::vector<int> > acc2NoNumBin, acc2CsNumBin, acc2NoDenBin, acc2CsDenBin;
  std::vector< std::vector<int> > effNoNumBin, effCsNumBin, effNoDenBin, effCsDenBin;
  std::vector< std::vector<TH1D*> > hNoSignalBin, hCsSignalBin;
  std::vector< std::vector<TH1D*> > hNoSignalCBin, hCsSignalCBin;

  std::vector<float> fitNoYield, fitNoYieldE, fitCsYield, fitCsYieldE;
  std::vector<float> fitNoYieldC, fitNoYieldCE, fitCsYieldC, fitCsYieldCE;
  std::vector< std::vector<float> > fitNoYieldBin, fitNoYieldEBin, fitCsYieldBin, fitCsYieldEBin;
  std::vector< std::vector<float> > fitNoYieldCBin, fitNoYieldCEBin, fitCsYieldCBin, fitCsYieldCEBin;


  ClassDef(plotRatio,1) //Testing plotRatio


};


#endif

