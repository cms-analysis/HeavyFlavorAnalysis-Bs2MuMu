#ifndef PLOTEFFICIENCIES
#define PLOTEFFICIENCIES

#include "plotClass.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"

class plotEfficiencies: public plotClass {

public:

  plotEfficiencies(const char *files="anaBmm.default.files", const char *dir = "default", const char *cuts = "default", int mode = 0);
  ~plotEfficiencies();

  void makeAll(int channels = 1);
  
  void mcTriggerEffs();
  void tnpVsMC(double m1pt, double m2pt, std::string what = "default");
  void triggerSignal(std::string cuts = "fls3d>8&&chi2/dof<2");
  void triggerNormalization(std::string cuts = "fls3d>8&&chi2/dof<2&&alpha<0.05&&iso>0.5");
  void misid(int pdgid, double bdtCut, std::string sample);
  void pidtables(int pdgid, std::string sample); 

  // -- works only on specific large small trees
  void allOverlayStudy();
  void hltOverlayStudy(std::string mode, bool barrel);
  TH1D* hltHist(const char *var, double xmin, double xmax, const char *treename, bool barrel);

  virtual void loopFunction(int function, int mode = 0); 
  virtual void loopFunction1(int mode); 
  virtual void loopFunction2(int mode); 

  void resetHistograms();
  void saveHists(std::string smode, double m1pt, double m2pt, std::string what);
  void numbersFromHist(int chan, int mode, double m1pt, double m2pt, numbers *aa);
  void texNumbers(double m1pt, double m2pt, std::string what);

  void convertLucasHistograms();
  void readFile(const char *fname = "luca/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_data_all.root", 
				  const char *pfix = "DATA");
  void read2Files(PidTable &a, const char *f1name, const char *f2name, const char *hname);

  bool fSplitSeagullsFromCowboys; 

  std::vector<TH1D*> fhptPosAll, fhptPosPass, fhptNegAll, fhptNegPass;
  TH2D *fh2PosAll, *fh2PosPass, *fh2NegAll, *fh2NegPass;
  int fPdgId;
  double fBDTcut; 

  ClassDef(plotEfficiencies,1) //Testing plotEfficiencies

};


#endif

