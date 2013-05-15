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


struct dt {
  std::vector<int> vars;
  std::vector<int> cuts;
};


// ----------------------------------------------------------------------
class plotBDT: public plotClass {

public:

  plotBDT(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default", int mode = 0);
  ~plotBDT();
  
  void makeAll(int channels = 1); 

  void tmvaControlPlots();
  void dumpParameters();
  void tmvaPlots(std::string type);
  void variableRanking(); 

  void xmlParsing(); 
  void xmlParsingVariables(std::string xmlfile); 
  void xmlParsingReadTree(std::string xmlfile); 
  void xmlResetHists(); 

  void plotSSB();
  void ssb();
  void overlap();
  void overlayBdtOutput(); 

  void hackedMC();

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
  virtual void loopFunction3(int mode);
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
  std::vector<TH1D*> fhBDT, fhAnaBDT, fhLoBDT, fhInBDT, fhHiBDT; 

  std::vector<TProfile*> fpMassBDT, fpMassAcBDT, fpMassAdBDT, fpNpvBDT, fpNpvAcBDT, fpNpvAdBDT; 
  std::vector<TH1D*> fhNpvBDTchan0, fhNpvBDTchan1, fhNpvAcBDTchan0, fhNpvAcBDTchan1, fhNpvAdBDTchan0, fhNpvAdBDTchan1; 
  TH1D *fmeanNpvBDTchan0, *fmeanNpvBDTchan1, *fmeanNpvAcBDTchan0, *fmeanNpvAcBDTchan1, *fmeanNpvAdBDTchan0, *fmeanNpvAdBDTchan1;
  TH1D *frmsNpvBDTchan0, *frmsNpvBDTchan1, *frmsNpvAcBDTchan0, *frmsNpvAcBDTchan1, *frmsNpvAdBDTchan0, *frmsNpvAdBDTchan1;

  // -- hacked MC plots
  std::vector<TH1D*> fhMcMass, fhMcBDT, fhMcRatio;
  std::vector<TH1D*> fhMcMass0, fhMcMass1, fhMcMass2, fhMcRatio1, fhMcRatio2;
  std::vector<TH1D*> fhMcBDT5;
  // -- other histograms
  std::vector<TH1D*> fhBgBDT, fhSgBDT;

  // -- XML parsing
  std::map<int, std::string> fVariableMap; 
  TH1D *fhBdtVariables, *fhBdtVariablesW8, *fhBdtNodes, *fhBdtNodesW8;
  std::vector<TH1D*> fhBdtVariableCuts;
  std::vector<TH1D*> fhBdtVariableCutsW8;

  
  ClassDef(plotBDT,1) //Testing plotBDT

};


#endif

