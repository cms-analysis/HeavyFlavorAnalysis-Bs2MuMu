#ifndef PLOTBDT
#define PLOTBDT

#include "plotClass.hh"
#include "preselection.hh"

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TProfile.h>
#include <TTree.h>
#include <TKey.h>


// ----------------------------------------------------------------------
class plotBDT: public plotClass {

public:

  plotBDT(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default", int mode = 11);
  ~plotBDT();
  
  void makeAll(int channels = 7); 

  void tmvaControlPlots();
  void dumpParameters();
  void tmvaPlots(std::string type);
  void variableRanking(); 
  void plotSSB();
  void ssb();
  void overlap();
  void overlayBdtOutput(); 

  void hackedMC(int chan = 0);
  void hackedMC1(std::string cuts = "1", double xmin = -0.4, double xmax = 0.5, std::string func = "pol0");
  void hackedMC2(double bdtCut1 = 0.0, double bdtCut2 = 0.2, std::string func = "pol0");
  void hackedMC3(int chan = 0);

  void allCorrelationPlots(double bdtcut, std::string fname = "TMVA-0");
  void correlationPlot(double bdtcut, std::string var, double xmin, double xmax, std::string fname);

  void validateAllDistributions(); 
  void validateDistributions(int channel, const char *type, int classID); 

  void bdtScan();
  void bdtDependencies(std::string mode = "SgData");
  void illustrateAntiMuonSample(const char *cuts = "hlt&&fls3d>5&&alpha<0.03&&chi2/dof<2&&iso>0.7"); 

  void testLoop(std::string mode);
  virtual void loopFunction(int function, int mode = 0); 
  virtual void loopFunction1(int mode); 
  virtual void loopFunction2(int mode);
  void  resetHistograms();
  void plotEffVsBg(int offset);

  std::string replaceLabelWithTex(string label);
  void SetFrameStyle( TH1* frame, Float_t scale = 1.0);
  void NormalizeHists( TH1* sig, TH1* bkg = 0);
  void SetSignalAndBackgroundStyle( TH1* sig, TH1* bkg, TH1* all = 0);
  void GetMethodTitle( TString & name, TKey * ikey);
  void GetMethodName( TString & name, TKey * mkey);
  void GetMethodTitle( TString & name, TDirectory * idir);
  int  GetNumberOfTargets( TDirectory *dir);
  int  GetNumberOfInputVariables( TDirectory *dir);

  std::string fXmlFile, fBdtString;
  TFile *fRootFile; 
  
  std::vector<string> fReaderVariables;
  
  // -- histograms/profiles  filled in loopFunction
  std::vector<TH1D*> fhMass, fhMassNoCuts;
  std::vector<TH1D*> fhBDT, fhLoBDT, fhInBDT, fhHiBDT; 

  std::vector<TProfile*> fpMassBDT, fpMassAcBDT, fpNpvBDT, fpNpvAcBDT; 

  // -- hacked MC plots
  std::vector<TH1D*> fhMcMass, fhMcBDT, fhMcRatio;
  std::vector<TH1D*> fhMcMass0, fhMcMass1, fhMcMass2, fhMcRatio1, fhMcRatio2;
  std::vector<TH1D*> fhMcBDT5;
  // -- other histograms
  std::vector<TH1D*> fhBgBDT, fhSgBDT;

  ClassDef(plotBDT,1) //Testing plotBDT

};


#endif

