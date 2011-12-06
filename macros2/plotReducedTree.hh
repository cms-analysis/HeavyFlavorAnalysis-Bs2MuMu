#ifndef PLOTREDUCEDTREE
#define PLOTREDUCEDTREE

#include "plotClass.hh"

struct numbers {
  std::string name;
  int index;
  double effGenFilter, effGenFilterE; 
  double fitYield, fitYieldE;
  double genFileYield, genYield, genChanYield, recoYield, chanYield, muidYield, trigYield, candYield;
  double genFileYieldE, genYieldE, genChanYieldE, recoYieldE, chanYieldE, muidYieldE, trigYieldE, candYieldE;
  double ana0Yield,  anaYield,  anaMuonYield,  anaTriggerYield,  anaWmcYield;
  double ana0YieldE, anaYieldE, anaMuonYieldE, anaTriggerYieldE, anaWmcYieldE;
  double cFrac, cFracE, acc, accE, accChan, accChanE;
  double effChan, effChanE; 
  double accMuidMC, accMuidMCE, accTrigMC, accTrigMCE;
  double effMuidMC, effMuidMCE, effTrigMC, effTrigMCE;
  double effMuidPid, effMuidPidE, effTrigPid, effTrigPidE;
  double effMuidPidMC, effMuidPidMCE, effTrigPidMC, effTrigPidMCE;
  double effCand, effCandE; 
  double effAna, effAnaE; 
  double effTot, effTotE, aEffProdMC, aEffProdMCE, effProdMC, effProdMCE, effProdPid, effProdPidE; 
  double effTotChan, effTotChanE; 
  double prodGenYield, combGenYield, chanGenYield; // eps*A corrected
  // -- signal stuff
  double expSignal;
  double bgObs, bgBsExp, bgBsExpE, bgBdExp, bgBdExpE; 
  double bsObs, bdObs; 
  double tauBs, tauBsE, tauBd, tauBdE; 
  double offRare, offRareE, bsRare, bsRareE, bdRare, bdRareE; 
  double pss, pssE, pdd, pddE;
  double psd, psdE, pds, pdsE;
  double mBdLo, mBdHi, mBsLo, mBsHi;
};


class plotReducedTree: public plotClass {

public:

  plotReducedTree(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotReducedTree();

  void makeAll(int channels = 3);
  
  void loopTree(int mode, int proc = -1);
  void tnpVsMC(double m1pt, double m2pt);
  void triggerSignal(std::string cuts = "fls3d>5&&alpha<0.05");
  void triggerNormalization(std::string cuts = "fls3d>5&&alpha<0.05");

  void accEffFromEffTree(std::string fname, std::string dname, numbers &a, cuts &b, int proc);
  void normYield(TH1 *h, int mode, double lo = 5.15, double hi=5.5);
  void csYield(TH1 *h, int mode, double lo = 5.25, double hi=5.6);
  void bgBlind(TH1 *h, int mode = 2, double lo = 4.5, double hi = 6.5); 

  void printNumbers(numbers &a, ostream &OUT);
  void initNumbers(numbers *a); 
  int  detChan(double m1eta, double m2eta);

  std::vector<TH1D*> fhMassAcc;
  std::vector<TH1D*> fhMassChan;
  std::vector<TH1D*> fhMassAbsNoCuts;
  std::vector<TH1D*> fhMassNoCuts;
  std::vector<TH1D*> fhMassNoCutsManyBins;
  std::vector<TH1D*> fhMassWithAnaCuts;
  std::vector<TH1D*> fhMassWithMuonCuts; 
  std::vector<TH1D*> fhMassWithTriggerCuts; 
  std::vector<TH1D*> fhMassWithAllCuts;
  std::vector<TH1D*> fhMassWithAllCutsBlind;
  std::vector<TH1D*> fhNorm;
  std::vector<TH1D*> fhMassWithMassCuts;
  std::vector<TH1D*> fhMassWithAnaCutsManyBins;
  std::vector<TH1D*> fhMassWithMuonCutsManyBins; 
  std::vector<TH1D*> fhMassWithTriggerCutsManyBins; 
  std::vector<TH1D*> fhMassWithAllCutsManyBins; 
  std::vector<TH1D*> fhMassWithMassCutsManyBins; 
  std::vector<TH1D*> fhMuId, fhMuIdMC;
  std::vector<TH1D*> fhMuTr, fhMuTrMC;
  std::vector<TH1D*> fh0PidTrigger, fh1PidTrigger, fh0PidMuID, fh1PidMuID;
  std::vector<TH1D*> fh0PidMCTrigger, fh1PidMCTrigger, fh0PidMCMuID, fh1PidMCMuID;
  std::vector<TH1D*> fh0MCTrigger, fh1MCTrigger, fh0MCMuID, fh1MCMuID;

  //  numbers fNumbersSig, fNumbersNorm, fNumbersCS; 
  std::vector<numbers*> fNumbersBs, fNumbersBd, fNumbersNo, fNumbersCs; 

  bool fDoApplyCowboyVeto; 
  double fBgExp, fBgExpE; 
  double fBgHist, fBgHistE; 
  double fBgHistExp, fBgHistExpE;

  double fNoSig, fNoSigE; 
  double fCsSig, fCsSigE; 

  double fBF, fu, fs, fsfu, fsfuE;

  bool fDoPrint;


  ClassDef(plotReducedTree,1) //Testing plotReducedTree

};


#endif

