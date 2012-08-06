#ifndef PLOTBDT
#define PLOTBDT

#include "plotClass.hh"

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TProfile.h>
#include <TTree.h>
#include <TKey.h>


// ----------------------------------------------------------------------
class plotBDT: public plotClass {

public:

  plotBDT(const char *files="anaBmm.v11.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotBDT();
  
  void makeAll(int channels = 7); 

  void tmvaControlPlots();
  void dumpParameters();
  void tmvaPlots();
  void variableRanking(); 
  void ssb();
  void overlap();

  void validateAllOddEven();
  void validateOddEven(const char *fnOdd, const char *fnEven, const char *type = "test", int classID = 0);
  
  void bdtScan();
  void bdtDependencies(std::string mode = "SgData");
  void illustrateAntiMuonSample(const char *cuts = "hlt&&fls3d>5&&alpha<0.03&&chi2/dof<2&&iso>0.7"); 

  void testLoop(std::string mode);
  virtual void loopFunction(int mode); 
  virtual void loopFunction1(); 
  virtual void loopFunction2();
  void  resetHistograms();
  void plotEffVsBg(int offset);

  std::string replaceLabelWithTex(string label);
  void SetFrameStyle( TH1* frame, Float_t scale = 1.0 );
  void NormalizeHists( TH1* sig, TH1* bkg = 0);
  void SetSignalAndBackgroundStyle( TH1* sig, TH1* bkg, TH1* all = 0 );
  void GetMethodTitle( TString & name, TKey * ikey );
  void GetMethodName( TString & name, TKey * mkey );
  void GetMethodTitle( TString & name, TDirectory * idir );

  std::string fXmlFile, fBdtString;
  TFile *fRootFile; 
  
  // -- histograms/profiles  filled in loopFunction
  std::vector<TH1D*> fhMass, fhMassNoCuts;
  std::vector<TH1D*> fhBDT, fhLoBDT, fhInBDT, fhHiBDT; 

  std::vector<TProfile*> fpMassBDT, fpMassAcBDT, fpNpvBDT, fpNpvAcBDT; 

  // -- other histograms
  std::vector<TH1D*> fhBgBDT, fhSgBDT;

  ClassDef(plotBDT,1) //Testing plotBDT

};


#endif

