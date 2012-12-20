#include "plotRatio.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "../macros/bayesianlimit.hh"
//#include "relativeYield.hh"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "THStack.h"
#include "TProfile.h"

using namespace std; 
using std::string; 

ClassImp(plotRatio)


// ----------------------------------------------------------------------
plotRatio::plotRatio(const char *files, const char *dir, const char *cuts, int mode) : plotClass(files, dir, cuts, mode) { 


  fDoPrint = true; 

  string hfname  = fDirectory + "/anaBmm.plotRatio." + fSuffix + ".root";
  cout << "open fHistFile: " << hfname << " RECREATE" << endl;
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

  fF["SgData"]->cd();
  TH1D *h1 = new TH1D("hasd", "", 100, 0., 100.); 
  h1->SetDirectory(fHistFile); 
  
  fNumbersFileName = fDirectory + "/anaBmm.plotRatio." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  printCuts(cout); 

// add stuff for my new method
  fDoUseBDT = true; 
  int NBINS = 100;

  fNPtbins = 4;
  accNoNumBin.resize(fNchan);  
  accCsNumBin.resize(fNchan);  
  accNoDenBin.resize(fNchan);  
  accCsDenBin.resize(fNchan);  
  acc2NoNumBin.resize(fNchan);  
  acc2CsNumBin.resize(fNchan);  
  acc2NoDenBin.resize(fNchan);  
  acc2CsDenBin.resize(fNchan);  
  effNoNumBin.resize(fNchan);  
  effCsNumBin.resize(fNchan);  
  effNoDenBin.resize(fNchan);  
  effCsDenBin.resize(fNchan);
  hNoSignalBin.resize(fNchan);
  hCsSignalBin.resize(fNchan);
  hNoSignalCBin.resize(fNchan);
  hCsSignalCBin.resize(fNchan);
  fitNoYieldBin.resize(fNchan);
  fitNoYieldEBin.resize(fNchan);
  fitCsYieldBin.resize(fNchan);
  fitCsYieldEBin.resize(fNchan);
  fitNoYieldCBin.resize(fNchan);
  fitNoYieldCEBin.resize(fNchan);
  fitCsYieldCBin.resize(fNchan);
  fitCsYieldCEBin.resize(fNchan);

  TH1D *h(0); 
  for (unsigned int i = 0; i < fNchan; ++i) {
    effNoNum.push_back(0);
    effCsNum.push_back(0);
    effNoDen.push_back(0);
    effCsDen.push_back(0);
    effNoNumBin[i].resize(fNPtbins);
    effCsNumBin[i].resize(fNPtbins);
    effNoDenBin[i].resize(fNPtbins);
    effCsDenBin[i].resize(fNPtbins);  
    hNoSignalBin[i].resize(fNPtbins);  
    hCsSignalBin[i].resize(fNPtbins);  
    hNoSignalCBin[i].resize(fNPtbins);  
    hCsSignalCBin[i].resize(fNPtbins);  
    
    for (unsigned int j = 0; j < fNPtbins; j++) {
      effNoNumBin[i][j] = 0;
      effCsNumBin[i][j] = 0;
      effNoDenBin[i][j] = 0;
      effNoDenBin[i][j] = 0;
    }
  }
  
  

  TH1D *h2 = new TH1D("hNoSignal0", "hNoSignal0", NBINS, 4.9, 5.9);
  h2->SetDirectory(fHistFile); 
  hNoSignal.push_back(h2);
  TH1D *h3 = new TH1D("hNoSignal1", "hNoSignal1", NBINS, 4.9, 5.9);
  h3->SetDirectory(fHistFile); 
  hNoSignal.push_back(h3);
  TH1D *h4 = new TH1D("hCsSignal0", "hCsSignal0", NBINS, 4.9, 5.9);
  h4->SetDirectory(fHistFile); 
  hCsSignal.push_back(h4);
  TH1D *h5 = new TH1D("hCsSignal1", "hCsSignal1", NBINS, 4.9, 5.9);
  h5->SetDirectory(fHistFile); 
  hCsSignal.push_back(h5);
  
  TH1D *hC2 = new TH1D("hNoSignalC0", "hNoSignalC0", NBINS, 4.9, 5.9);
  hC2->SetDirectory(fHistFile); 
  hNoSignalC.push_back(hC2);
  TH1D *hC3 = new TH1D("hNoSignalC1", "hNoSignalC1", NBINS, 4.9, 5.9);
  hC3->SetDirectory(fHistFile); 
  hNoSignalC.push_back(hC3);
  TH1D *hC4 = new TH1D("hCsSignalC0", "hCsSignalC0", NBINS, 4.9, 5.9);
  hC4->SetDirectory(fHistFile); 
  hCsSignalC.push_back(hC4);
  TH1D *hC5 = new TH1D("hCsSignalC1", "hCsSignalC1", NBINS, 4.9, 5.9);
  hC5->SetDirectory(fHistFile); 
  hCsSignalC.push_back(hC5);

  TH1D *hNo00 = new TH1D(Form("hNoSignalBin00"), Form("hNoSignalBin00"), NBINS, 4.9, 5.9);
  hNoSignalBin[0][0] = hNo00;
  TH1D *hNo01 = new TH1D(Form("hNoSignalBin01"), Form("hNoSignalBin01"), NBINS, 4.9, 5.9);
  hNoSignalBin[0][1] = hNo01;
  TH1D *hNo02 = new TH1D(Form("hNoSignalBin02"), Form("hNoSignalBin02"), NBINS, 4.9, 5.9);
  hNoSignalBin[0][2] = hNo02;
  TH1D *hNo03 = new TH1D(Form("hNoSignalBin03"), Form("hNoSignalBin03"), NBINS, 4.9, 5.9);
  hNoSignalBin[0][3] = hNo03;
  TH1D *hNo10 = new TH1D(Form("hNoSignalBin10"), Form("hNoSignalBin10"), NBINS, 4.9, 5.9);
  hNoSignalBin[1][0] = hNo10;
  TH1D *hNo11 = new TH1D(Form("hNoSignalBin11"), Form("hNoSignalBin11"), NBINS, 4.9, 5.9);
  hNoSignalBin[1][1] = hNo11;
  TH1D *hNo12 = new TH1D(Form("hNoSignalBin12"), Form("hNoSignalBin12"), NBINS, 4.9, 5.9);
  hNoSignalBin[1][2] = hNo12;
  TH1D *hNo13 = new TH1D(Form("hNoSignalBin13"), Form("hNoSignalBin13"), NBINS, 4.9, 5.9);
  hNoSignalBin[1][3] = hNo13;

  TH1D *hCs00 = new TH1D(Form("hCsSignalBin00"), Form("hCsSignalBin00"), NBINS, 4.9, 5.9);
  hCsSignalBin[0][0] = hCs00;
  TH1D *hCs01 = new TH1D(Form("hCsSignalBin01"), Form("hCsSignalBin01"), NBINS, 4.9, 5.9);
  hCsSignalBin[0][1] = hCs01;
  TH1D *hCs02 = new TH1D(Form("hCsSignalBin02"), Form("hCsSignalBin02"), NBINS, 4.9, 5.9);
  hCsSignalBin[0][2] = hCs02;
  TH1D *hCs03 = new TH1D(Form("hCsSignalBin03"), Form("hCsSignalBin03"), NBINS, 4.9, 5.9);
  hCsSignalBin[0][3] = hCs03;
  TH1D *hCs10 = new TH1D(Form("hCsSignalBin10"), Form("hCsSignalBin10"), NBINS, 4.9, 5.9);
  hCsSignalBin[1][0] = hCs10;
  TH1D *hCs11 = new TH1D(Form("hCsSignalBin11"), Form("hCsSignalBin11"), NBINS, 4.9, 5.9);
  hCsSignalBin[1][1] = hCs11;
  TH1D *hCs12 = new TH1D(Form("hCsSignalBin12"), Form("hCsSignalBin12"), NBINS, 4.9, 5.9);
  hCsSignalBin[1][2] = hCs12;
  TH1D *hCs13 = new TH1D(Form("hCsSignalBin13"), Form("hCsSignalBin13"), NBINS, 4.9, 5.9);
  hCsSignalBin[1][3] = hCs13;

  TH1D *hNoC00 = new TH1D(Form("hNoSignalCBin00"), Form("hNoSignalCBin00"), NBINS, 4.9, 5.9);
  hNoSignalCBin[0][0] = hNoC00;
  TH1D *hNoC01 = new TH1D(Form("hNoSignalCBin01"), Form("hNoSignalCBin01"), NBINS, 4.9, 5.9);
  hNoSignalCBin[0][1] = hNoC01;
  TH1D *hNoC02 = new TH1D(Form("hNoSignalCBin02"), Form("hNoSignalCBin02"), NBINS, 4.9, 5.9);
  hNoSignalCBin[0][2] = hNoC02;
  TH1D *hNoC03 = new TH1D(Form("hNoSignalCBin03"), Form("hNoSignalCBin03"), NBINS, 4.9, 5.9);
  hNoSignalCBin[0][3] = hNoC03;
  TH1D *hNoC10 = new TH1D(Form("hNoSignalCBin10"), Form("hNoSignalCBin10"), NBINS, 4.9, 5.9);
  hNoSignalCBin[1][0] = hNoC10;
  TH1D *hNoC11 = new TH1D(Form("hNoSignalCBin11"), Form("hNoSignalCBin11"), NBINS, 4.9, 5.9);
  hNoSignalCBin[1][1] = hNoC11;
  TH1D *hNoC12 = new TH1D(Form("hNoSignalCBin12"), Form("hNoSignalCBin12"), NBINS, 4.9, 5.9);
  hNoSignalCBin[1][2] = hNoC12;
  TH1D *hNoC13 = new TH1D(Form("hNoSignalCBin13"), Form("hNoSignalCBin13"), NBINS, 4.9, 5.9);
  hNoSignalCBin[1][3] = hNoC13;

  TH1D *hCsC00 = new TH1D(Form("hCsSignalCBin00"), Form("hCsSignalCBin00"), NBINS, 4.9, 5.9);
  hCsSignalCBin[0][0] = hCsC00;
  TH1D *hCsC01 = new TH1D(Form("hCsSignalCBin01"), Form("hCsSignalCBin01"), NBINS, 4.9, 5.9);
  hCsSignalCBin[0][1] = hCsC01;
  TH1D *hCsC02 = new TH1D(Form("hCsSignalCBin02"), Form("hCsSignalCBin02"), NBINS, 4.9, 5.9);
  hCsSignalCBin[0][2] = hCsC02;
  TH1D *hCsC03 = new TH1D(Form("hCsSignalCBin03"), Form("hCsSignalCBin03"), NBINS, 4.9, 5.9);
  hCsSignalCBin[0][3] = hCsC03;
  TH1D *hCsC10 = new TH1D(Form("hCsSignalCBin10"), Form("hCsSignalCBin10"), NBINS, 4.9, 5.9);
  hCsSignalCBin[1][0] = hCsC10;
  TH1D *hCsC11 = new TH1D(Form("hCsSignalCBin11"), Form("hCsSignalCBin11"), NBINS, 4.9, 5.9);
  hCsSignalCBin[1][1] = hCsC11;
  TH1D *hCsC12 = new TH1D(Form("hCsSignalCBin12"), Form("hCsSignalCBin12"), NBINS, 4.9, 5.9);
  hCsSignalCBin[1][2] = hCsC12;
  TH1D *hCsC13 = new TH1D(Form("hCsSignalCBin13"), Form("hCsSignalCBin13"), NBINS, 4.9, 5.9);
  hCsSignalCBin[1][3] = hCsC13;


  TH1D *h6 = new TH1D("hNoPt0", "hNoPt0", 50, 0, 50.);
  h6->SetDirectory(fHistFile); 
  hNoPt.push_back(h6);
  TH1D *h7 = new TH1D("hNoPt1", "hNoPt1", 50, 0, 50.);
  h7->SetDirectory(fHistFile); 
  hNoPt.push_back(h7);
  TH1D *h8 = new TH1D("hCsPt0", "hCsPt0", 50, 0, 50.);
  h8->SetDirectory(fHistFile); 
  hCsPt.push_back(h8);
  TH1D *h9 = new TH1D("hCsPt1", "hCsPt1", 50, 0, 50.);
  h9->SetDirectory(fHistFile); 
  hCsPt.push_back(h9);

  
  fDoUseBDT = true; 
  fInvertedIso = false; 
  fNormProcessed = false; 
}

// ----------------------------------------------------------------------
plotRatio::~plotRatio() {
  cout << "plotRatio dtor: " << fHistFile << endl;
  if (fHistFile) {
    fHistFile->Write();
    fHistFile->Close();
  }
}


// ----------------------------------------------------------------------
void plotRatio::makeAll(int channels) {

  zone(1);
  if (channels & 2) {
    fNormProcessed = false; 
    fDoUseBDT = true; 
//    computeCsNoRatio();
    computeCsNoRatioNew();
  }

}

// ----------------------------------------------------------------------
void plotRatio::computeCsNoRatio() {


  cout << "--> loopTree: CS MC" << endl;
  loopTree(20);  // CS signal eff
  c0->Modified(); c0->Update();
  
  cout << "--> loopTree: signal data" << endl;
  loopTree(25);  // control sample data 
  c0->Modified(); c0->Update();
  
  if (false == fNormProcessed) {
    fNormProcessed = true; 
    cout << "--> loopTree: norm MC" << endl;
    loopTree(10); // normalization eff
    c0->Modified(); c0->Update();
    
    cout << "--> loopTree: norm data" << endl;
    loopTree(15); // data normalization 
    c0->Modified(); c0->Update();
  }

  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;    
  fTEX << "% -- Control sample branching fraction" << endl;
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;    

  string cache = fSuffix; 
  if (fDoUseBDT) fSuffix = "bdt" + fSuffix; 
  
  double result, resultE; 
  for (int i = 0; i < 2; ++i) {
    result = (fNumbersCs[i]->fitYield/fNumbersNo[i]->fitYield)
      *(fu/fs)
      *(fNumbersNo[i]->acc/fNumbersCs[i]->acc)
      *(fNumbersNo[i]->effCand/fNumbersCs[i]->effCand)     
      *(fNumbersNo[i]->effAna/fNumbersCs[i]->effAna)
      *(fNumbersNo[i]->effMuidTNP/fNumbersCs[i]->effMuidTNP)
      *(fNumbersNo[i]->effTrigTNP/fNumbersCs[i]->effTrigTNP)
      * fBF["NoMc"];
cout << "chan " << i << ": Acceptance No = " << fNumbersNo[i]->acc << "  CS = " << fNumbersCs[i]->acc << endl;
    
    resultE = dRatio(fNumbersCs[i]->fitYield, fNumbersCs[i]->fitYieldE, fNumbersNo[i]->fitYield, fNumbersNo[i]->fitYieldE)
      *(fu/fs)
      *(fNumbersNo[i]->acc/fNumbersCs[i]->acc)
      *(fNumbersNo[i]->effCand/fNumbersCs[i]->effCand)     
      *(fNumbersNo[i]->effMuidTNP/fNumbersCs[i]->effMuidTNP)
      *(fNumbersNo[i]->effTrigTNP/fNumbersCs[i]->effTrigTNP)
      *(fNumbersNo[i]->effAna/fNumbersCs[i]->effAna)
      * fBF["NoMc"];

    cout << "chan " << i << ": PID fact branching fraction: " << result << "+/-" << resultE << endl;
    fTEX << formatTex(result, Form("%s:N-CSBF-TNP-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(resultE, Form("%s:N-CSBF-TNP-BS%i:err", fSuffix.c_str(), i), 6) << endl;

    result = (fNumbersCs[i]->fitYield/fNumbersNo[i]->fitYield)
      *(fu/fs)
      *(fNumbersNo[i]->acc/fNumbersCs[i]->acc)
      *(fNumbersNo[i]->effCand/fNumbersCs[i]->effCand)     
      *(fNumbersNo[i]->effAna/fNumbersCs[i]->effAna)
      *(fNumbersNo[i]->effMuidMC/fNumbersCs[i]->effMuidMC)
      *(fNumbersNo[i]->effTrigMC/fNumbersCs[i]->effTrigMC)
      * fBF["NoMc"];
cout << "chan " << i << ": Acceptance No = " << fNumbersNo[i]->acc << "  CS = " << fNumbersCs[i]->acc << endl;

    
    resultE = dRatio(fNumbersCs[i]->fitYield, fNumbersCs[i]->fitYieldE, fNumbersNo[i]->fitYield, fNumbersNo[i]->fitYieldE)
      *(fu/fs)
      *(fNumbersNo[i]->acc/fNumbersCs[i]->acc)
      *(fNumbersNo[i]->effCand/fNumbersCs[i]->effCand)     
      *(fNumbersNo[i]->effMuidMC/fNumbersCs[i]->effMuidMC)
      *(fNumbersNo[i]->effTrigMC/fNumbersCs[i]->effTrigMC)
      *(fNumbersNo[i]->effAna/fNumbersCs[i]->effAna)
      * fBF["NoMc"];

    cout << "chan " << i << ": MC fact branching fraction: " << result << "+/-" << resultE << endl;
    fTEX << formatTex(result, Form("%s:N-CSBF-MC-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(resultE, Form("%s:N-CSBF-MC-BS%i:err", fSuffix.c_str(), i), 6) << endl;
    
    result = (fNumbersCs[i]->fitYield/fNumbersNo[i]->fitYield)
      *(fu/fs)
      *(fNumbersNo[i]->effTot/fNumbersCs[i]->effTot)
      * fBF["NoMc"];
      
      cout << "  ********************************************************************** " << endl;
      cout << "yields: chan " << i << "  No = " <<  fNumbersNo[i]->fitYield << " +- " << fNumbersNo[i]->fitYieldE << 
                                      "  CS = " << fNumbersCs[i]->fitYield  << " +- " << fNumbersCs[i]->fitYieldE << endl;
      cout << "  ********************************************************************** " << endl;

      
      cout << "total eff, chan " << i << "  No = " <<  fNumbersNo[i]->effTot << "  Cs = " << fNumbersCs[i]->effTot << endl;

    cout << "chan " << i << ": branching fraction: " << result << "+/-" << resultE << endl;

    fTEX << formatTex(result, Form("%s:N-CSBF-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(resultE, Form("%s:N-CSBF-BS%i:err", fSuffix.c_str(), i), 6) << endl;

  }

  fSuffix = cache; 

  printCsBFNumbers();
}

// ----------------------------------------------------------------------
void plotRatio::printCsBFNumbers() {

  string cache = fSuffix; 
  if (fDoUseBDT) fSuffix = "bdt" + fSuffix; 
  
  fTEX << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  for (unsigned int i = 0; i < fNchan; ++i) {
    fTEX << "% -- CONTROL SAMPLE " << i << endl;
    fTEX << formatTex(fNumbersCs[i]->effTot, Form("%s:N-EFF-TOT-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(fNumbersCs[i]->effTotE, Form("%s:N-EFF-TOT-BS%i:err", fSuffix.c_str(), i), 6) << endl;

    fTEX << formatTex(fNumbersCs[i]->acc, Form("%s:N-ACC-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->accE, Form("%s:N-ACC-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->effMuidTNP, Form("%s:N-EFF-MU-PID-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->effMuidTNPE, Form("%s:N-EFF-MU-PID-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->effMuidTNPMCE, Form("%s:N-EFF-MU-PIDMC-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->effMuidMCE, Form("%s:N-EFF-MU-MC-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->effTrigTNP, Form("%s:N-EFF-TRIG-PID-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->effTrigTNPE, Form("%s:N-EFF-TRIG-PID-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->effTrigTNPMCE, Form("%s:N-EFF-TRIG-PIDMC-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->effTrigMCE, Form("%s:N-EFF-TRIG-MC-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->effCand, Form("%s:N-EFF-CAND-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->effCandE, Form("%s:N-EFF-CAND-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->effAna, Form("%s:N-EFF-ANA-BS%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersCs[i]->effAnaE, Form("%s:N-EFF-ANA-BS%i:err", fSuffix.c_str(), i), 4) << endl;

    fTEX << formatTex(fNumbersCs[i]->fitYield, Form("%s:N-OBS-BS%i:val", fSuffix.c_str(), i), 0) << endl;
    fTEX << formatTex(fNumbersCs[i]->fitYieldE, Form("%s:N-OBS-BS%i:err", fSuffix.c_str(), i), 0) << endl;
  }

  fSuffix = cache; 
}

void plotRatio::computeCsNoRatioNew() {

cout << "just to test. bdt barrel cut = " << fCuts[0]->bdt << "  and endcap = " << fCuts[1]->bdt << endl;

  int eventsToLoop = -1;

// I don't understand where these values come from
// they are supplied externally in the input file...
  float effFilterCs = 0.0162;
  float effFilterNo = 0.0181;
  float fracCsKst[2];
  float effTurnon[2];
  
  fIsMC=true;
 
  // genFileYield actually comes from the NoMc or CsMc file effTree->GetEntries(). Now sure what Acc trees are used for....

  string mode = "CsMc";
  TTree *t = (TTree*)fF[mode]->Get("candAnaBs2JpsiPhi/effTree");
  cout << "Cs acc tree entries = " << t->GetEntries() << endl;
  effCsDen[0] = t->GetEntries()/effFilterCs;
  effCsDen[1] = t->GetEntries()/effFilterCs;
  effCsDenBin[0][0] = t->GetEntries()/effFilterCs;
  effCsDenBin[0][1] = t->GetEntries()/effFilterCs;
  effCsDenBin[0][2] = t->GetEntries()/effFilterCs;
  effCsDenBin[0][3] = t->GetEntries()/effFilterCs;
  effCsDenBin[1][0] = t->GetEntries()/effFilterCs;
  effCsDenBin[1][1] = t->GetEntries()/effFilterCs;
  effCsDenBin[1][2] = t->GetEntries()/effFilterCs;
  effCsDenBin[1][3] = t->GetEntries()/effFilterCs;

  t = getTree(mode);
  setupTree(t, mode);
  loopOverTree(t, mode, 1, eventsToLoop);
  
  mode = "NoMc";
  t = (TTree*)fF[mode]->Get("candAnaBu2JpsiK/effTree");
  cout << "No acc tree entries = " << t->GetEntries() << endl;
  effNoDen[0] = t->GetEntries()/effFilterNo;
  effNoDen[1] = t->GetEntries()/effFilterNo;
  effNoDenBin[0][0] = t->GetEntries()/effFilterNo;
  effNoDenBin[0][1] = t->GetEntries()/effFilterNo;
  effNoDenBin[0][2] = t->GetEntries()/effFilterNo;
  effNoDenBin[0][3] = t->GetEntries()/effFilterNo;
  effNoDenBin[1][0] = t->GetEntries()/effFilterNo;
  effNoDenBin[1][1] = t->GetEntries()/effFilterNo;
  effNoDenBin[1][2] = t->GetEntries()/effFilterNo;
  effNoDenBin[1][3] = t->GetEntries()/effFilterNo;
  
  t = getTree(mode);
  setupTree(t, mode);
  loopOverTree(t, mode, 1, eventsToLoop);

  fIsMC=false;
  mode = "CsData";
  t = getTree(mode);
  setupTree(t, mode);
  loopOverTree(t, mode, 1, eventsToLoop);
  
  mode = "NoData";
  t = getTree(mode);
  setupTree(t, mode);
  loopOverTree(t, mode, 1, eventsToLoop);

  //now fit the mass histograms for yields
  cout << "before integrated fits" << endl;
  for (unsigned int i = 0; i < fNchan; ++i) {
//first do fits to mass constrained distributions
cout << "before normyield2 integrated fit i = " << i << endl;
    normYield2(hNoSignalC[i],i,5.0, 5.6);
    fitNoYieldC.push_back(fNoSig);
    fitNoYieldCE.push_back(fNoSigE);
              
    effTurnon[i] = fNoErrTurnon;
    cout << "turn on value for mode " << i << " " << fNoErrTurnon << endl;
	                                                           
cout << "before csyield2 integrated fit i = " << i << endl;
    csYield2(hCsSignalC[i],i,5.1, 5.60);
    fitCsYieldC.push_back(fCsSig);
    fitCsYieldCE.push_back(fCsSigE);
cout << "before normyield2 integrated unconstrainted fit i = " << i << endl;

//second do fits to unconstrainted distributions
    normYield2(hNoSignal[i],i,4.9, 5.6, effTurnon[0]);
    fitNoYield.push_back(fNoSig);
    fitNoYieldE.push_back(fNoSigE);
    
cout << "before csyield2 integrated unconstrained fit i = " << i << endl;
    csYield2(hCsSignal[i],i,5.1, 5.60); //defaults fraction = -1 (floating), preco = -1 (no step function)
    fitCsYield.push_back(fCsSig);
    fitCsYieldE.push_back(fCsSigE);
    //get fraction of K* from integrated fits here
    fracCsKst[i] = fCsKstFrac;
    cout << "K*0 fraction for mode " << i << " " << fCsKstFrac << endl;

    for (unsigned int j = 0; j < 4; j++){
      normYield2(hNoSignalBin[i][j],i,4.9,5.6, effTurnon[0]);
      fitNoYieldBin[i].push_back(fNoSig);
      fitNoYieldEBin[i].push_back(fNoSigE);

      csYield2(hCsSignalBin[i][j],i,5.1,5.65, fracCsKst[i]);
      fitCsYieldBin[i].push_back(fCsSig);
      fitCsYieldEBin[i].push_back(fCsSigE);

      normYield2(hNoSignalCBin[i][j],i,5.0,5.6);
      fitNoYieldCBin[i].push_back(fNoSig);
      fitNoYieldCEBin[i].push_back(fNoSigE);

      csYield2(hCsSignalCBin[i][j],i,5.1,5.65);
      fitCsYieldCBin[i].push_back(fCsSig);
      fitCsYieldCEBin[i].push_back(fCsSigE);
    }

  }

  fHistFile->cd();
  


  cout << "got some efficiency results: " << endl;
  cout << "first for barrel: No num = " << effNoNum[0] << " den = " << effNoDen[0] << " = " << (float)effNoNum[0]/(float)effNoDen[0] <<
   " Cs num = " << effCsNum[0] << " den = " << effCsDen[0] << " = " << (float)effCsNum[0]/(float)effCsDen[0] << endl;
  cout << "next for endcap:  No num = " << effNoNum[1] << " den = " << effNoDen[1] << " = " << (float)effNoNum[1]/(float)effNoDen[1] <<
   " Cs num = " << effCsNum[1] << " den = " << effCsDen[1] << " = " << (float)effCsNum[1]/(float)effCsDen[1] << endl;

  cout << "got some data yields: " << endl;
  cout << "yield normalization barrel = " << fitNoYield[0] << " +- " << fitNoYieldE[0] << "   endcap = " << fitNoYield[1] << " +- " << fitNoYieldE[1] << endl;
  cout << "yield control barrel       = " << fitCsYield[0] << " +- " << fitCsYieldE[0] << "   endcap = " << fitCsYield[1] << " +- " << fitCsYieldE[1] << endl;
  cout << "or yields with JPsi mass constraints:" << endl;
  cout << "yield normalization barrel = " << fitNoYieldC[0] << " +- " << fitNoYieldCE[0] << "   endcap = " << fitNoYieldC[1] << " +- " << fitNoYieldCE[1] << endl;
  cout << "yield control barrel       = " << fitCsYieldC[0] << " +- " << fitCsYieldCE[0] << "   endcap = " << fitCsYieldC[1] << " +- " << fitCsYieldCE[1] << endl;

cout << "    ================================== " << endl;

float barrelBsBF = (fitCsYield[0]/fitNoYield[0])*(fu/fs)*(((float)effNoNum[0]/(float)effNoDen[0])/((float)effCsNum[0]/(float)effCsDen[0]))*fBF["NoMc"];
float endcapBsBF = (fitCsYield[1]/fitNoYield[1])*(fu/fs)*(((float)effNoNum[1]/(float)effNoDen[1])/((float)effCsNum[1]/(float)effCsDen[1]))*fBF["NoMc"];

cout << "  Bs BF (to compare to old code): barrel = " << barrelBsBF << " endcap = " << endcapBsBF << endl;

float barrelfsfu = (fitCsYield[0]/fitNoYield[0])*(((float)effNoNum[0]/(float)effNoDen[0])/((float)effCsNum[0]/(float)effCsDen[0]))*((0.001016)/(0.00109*0.489)); //last factor if BF B->JPsiK/BF Bs->JPsiPhi(kk)
float endcapfsfu = (fitCsYield[1]/fitNoYield[1])*(((float)effNoNum[1]/(float)effNoDen[1])/((float)effCsNum[1]/(float)effCsDen[1]))*((0.001016)/(0.00109*0.489));

cout << "  Or our value of fs/fu: barrel = " << barrelfsfu << "  endcap = " << endcapfsfu << endl;

  cout << "==============================" << endl;
  cout << "Next let's get some binned results" << endl;
  
  //2011 values
  //float noEffSyst0 = 0.071;
  //float csEffSyst0 = 0.071;
  //float noEffSyst1 = 0.108;
  //float csEffSyst1 = 0.105;
  //2012 values
  float noEffSyst0 = 0.040;
  float csEffSyst0 = 0.040;
  float noEffSyst1 = 0.040;
  float csEffSyst1 = 0.040;

  float effNo[2];
  effNo[0] = ((float)effNoNum[0]/(float)effNoDen[0]);
  effNo[1] = ((float)effNoNum[1]/(float)effNoDen[1]);
  float effCs[2];
  effCs[0] = ((float)effCsNum[0]/(float)effCsDen[0]);
  effCs[1] = ((float)effCsNum[1]/(float)effCsDen[1]);
  float effNoE[2];
  effNoE[0] = sqrt((float)effNoNum[0])/(float)effNoDen[0];
  effNoE[1] = sqrt((float)effNoNum[1])/(float)effNoDen[1];
  float effCsE[2];
  effCsE[0] = sqrt((float)effCsNum[0])/(float)effCsDen[0];
  effCsE[1] = sqrt((float)effCsNum[1])/(float)effCsDen[1];
  float effNoBin[2][4];
  float effCsBin[2][4];
  float effNoEBin[2][4];
  float effCsEBin[2][4];
  float fsfu[2];
  float fsfuEStat[2];
  float fsfuESyst[2];
  fsfu[0] = barrelfsfu;
  fsfu[1] = endcapfsfu;
  fsfuEStat[0] = fsfu[0]*sqrt( (fitNoYieldE[0]/fitNoYield[0])*(fitNoYieldE[0]/fitNoYield[0]) + (fitCsYieldE[0]/fitCsYield[0])*(fitCsYieldE[0]/fitCsYield[0]) + (effNoE[0]/effNo[0])*(effNoE[0]/effNo[0]) + (effCsE[0]/effCs[0])*(effCsE[0]/effCs[0]) );
  fsfuEStat[1] = fsfu[1]*sqrt( (fitNoYieldE[1]/fitNoYield[1])*(fitNoYieldE[1]/fitNoYield[1]) + (fitCsYieldE[1]/fitCsYield[1])*(fitCsYieldE[1]/fitCsYield[1]) + (effNoE[1]/effNo[1])*(effNoE[1]/effNo[1]) + (effCsE[1]/effCs[1])*(effCsE[1]/effCs[1]) );
  fsfuESyst[0] = fsfu[0]*sqrt(0.05*0.05 + 0.05*0.05 + noEffSyst0*noEffSyst0 + csEffSyst0*csEffSyst0);
  fsfuESyst[1] = fsfu[1]*sqrt(0.05*0.05 + 0.05*0.05 + noEffSyst1*noEffSyst1 + csEffSyst1*csEffSyst1);
  float fsfuBin[2][4];
  float fsfuBinEStat[2][4];
  float fsfuBinESyst[2][4];
  float fsfuBinETot[2][4];

  
  for (int ib=0; ib<4; ib++){
    cout << "==========  starting bin " << ib << " ===========" << endl;
    cout << " Fit yields: No barrel = " << fitNoYieldBin[0][ib] << " +- " << fitNoYieldEBin[0][ib] << "  endcap = " << fitNoYieldBin[1][ib] << " +- " << fitNoYieldEBin[1][ib] <<
                   "      Cs barrel = " << fitCsYieldBin[0][ib] << " +- " << fitCsYieldEBin[0][ib] << "  endcap = " << fitCsYieldBin[1][ib] << " +- " << fitCsYieldEBin[1][ib] << endl;
    cout << " Fit yields with JPsi mass constraint: No barrel = " << fitNoYieldCBin[0][ib] << " +- " << fitNoYieldCEBin[0][ib] << "  endcap = " << fitNoYieldCBin[1][ib] << " +- " << fitNoYieldCEBin[1][ib] <<
                   "                                Cs barrel = " << fitCsYieldCBin[0][ib] << " +- " << fitCsYieldCEBin[0][ib] << "  endcap = " << fitCsYieldCBin[1][ib] << " +- " << fitCsYieldCEBin[1][ib] << endl;
    cout << " Eff nums: No barrel = " << effNoNumBin[0][ib] << "  endcap = " << effNoNumBin[1][ib] << 
                    "   Cs barrel = " << effCsNumBin[0][ib] << "  endcap = " << effCsNumBin[1][ib] <<  endl;
    cout << " Eff: No barrel = " << (float)(effNoNumBin[0][ib])/(float)(effNoDenBin[0][ib]) << "  endcap = " << (float)(effNoNumBin[1][ib])/(float)(effNoDenBin[1][ib]) << 
               "   Cs barrel = " << (float)(effCsNumBin[0][ib])/(float)(effCsDenBin[0][ib]) << "  endcap = " << (float)(effCsNumBin[1][ib])/(float)(effCsDenBin[1][ib]) << endl;
    cout << endl;
    cout << "or Bs BF" << endl;
    cout << "  Bs BF (to compare to old code): barrel = " <<    (fitCsYieldBin[0][ib]/fitNoYieldBin[0][ib])*(fu/fs)*(((float)(effNoNumBin[0][ib])/(float)(effNoDenBin[0][ib]))/((float)(effCsNumBin[0][ib])/(float)(effCsDenBin[0][ib])))*fBF["NoMc"] <<
                                         " endcap = " << (fitCsYieldBin[1][ib]/fitNoYieldBin[1][ib])*(fu/fs)*(((float)(effNoNumBin[1][ib])/(float)(effNoDenBin[1][ib]))/((float)(effCsNumBin[1][ib])/(float)(effCsDenBin[1][ib])))*fBF["NoMc"] << endl;
    cout << "or can quote it as ratio of Bs/B+" << endl;
    cout << "  Bs/Bp BF : barrel = " << (fitCsYieldBin[0][ib]/fitNoYieldBin[0][ib])*(fu/fs)*(((float)(effNoNumBin[0][ib])/(float)(effNoDenBin[0][ib]))/((float)(effCsNumBin[0][ib])/(float)(effCsDenBin[0][ib]))) <<
                        " endcap = " << (fitCsYieldBin[1][ib]/fitNoYieldBin[1][ib])*(fu/fs)*(((float)(effNoNumBin[1][ib])/(float)(effNoDenBin[1][ib]))/((float)(effCsNumBin[1][ib])/(float)(effCsDenBin[1][ib]))) << endl;
    cout << "or can quote it as calculated value of fs/fu" << endl;
    fsfuBin[0][ib] = (fitCsYieldBin[0][ib]/fitNoYieldBin[0][ib])*(((float)(effNoNumBin[0][ib])/(float)(effNoDenBin[0][ib]))/((float)(effCsNumBin[0][ib])/(float)(effCsDenBin[0][ib])))*((0.001016)/(0.00109*0.489));
    fsfuBin[1][ib] = (fitCsYieldBin[1][ib]/fitNoYieldBin[1][ib])*(((float)(effNoNumBin[1][ib])/(float)(effNoDenBin[1][ib]))/((float)(effCsNumBin[1][ib])/(float)(effCsDenBin[1][ib])))*((0.001016)/(0.00109*0.489));
    effNoBin[0][ib] = ((float)(effNoNumBin[0][ib])/(float)(effNoDenBin[0][ib]));
    effCsBin[0][ib] = ((float)(effCsNumBin[0][ib])/(float)(effCsDenBin[0][ib]));
    effNoEBin[0][ib] = sqrt((float)effNoNumBin[0][ib])/(float)effNoDenBin[0][ib];
    effCsEBin[0][ib] =  sqrt((float)effCsNumBin[0][ib])/(float)effCsDenBin[0][ib];
    effNoBin[1][ib] = ((float)(effNoNumBin[1][ib])/(float)(effNoDenBin[1][ib]));
    effCsBin[1][ib] = ((float)(effCsNumBin[1][ib])/(float)(effCsDenBin[1][ib]));
    effNoEBin[1][ib] = sqrt((float)effNoNumBin[1][ib])/(float)effNoDenBin[1][ib];
    effCsEBin[1][ib] =  sqrt((float)effCsNumBin[1][ib])/(float)effCsDenBin[1][ib];

    fsfuBinEStat[0][ib] = fsfuBin[0][ib]*sqrt( ((fitNoYieldEBin[0][ib]/fitNoYieldBin[0][ib])*(fitNoYieldEBin[0][ib]/fitNoYieldBin[0][ib])) + ((fitCsYieldEBin[0][ib]/fitCsYieldBin[0][ib])*(fitCsYieldEBin[0][ib]/fitCsYieldBin[0][ib])) + ((effNoEBin[0][ib]/effNoBin[0][ib])*(effNoEBin[0][ib]/effNoBin[0][ib])) + ((effCsEBin[0][ib]/effCsBin[0][ib])*(effCsEBin[0][ib]/effCsBin[0][ib])) );
cout << "the pieces of the stat error are: " << fsfuBin[0][ib] << "  " << (fitNoYieldEBin[0][ib]/fitNoYieldBin[0][ib]) << "   " << (fitCsYieldEBin[0][ib]/fitCsYieldBin[0][ib]) << "  " << (effNoEBin[0][ib]/effNoBin[0][ib]) << "  " <<  (effCsEBin[0][ib]/effCsBin[0][ib]) << endl;
    fsfuBinEStat[1][ib] = fsfuBin[1][ib]*sqrt( (fitNoYieldEBin[1][ib]/fitNoYieldBin[1][ib])*(fitNoYieldEBin[1][ib]/fitNoYieldBin[1][ib]) + (fitCsYieldEBin[1][ib]/fitCsYieldBin[1][ib])*(fitCsYieldEBin[1][ib]/fitCsYieldBin[1][ib]) + (effNoEBin[1][ib]/effNoBin[1][ib])*(effNoEBin[1][ib]/effNoBin[1][ib]) + (effCsEBin[1][ib]/effCsBin[1][ib])*(effCsEBin[1][ib]/effCsBin[1][ib]) );
    fsfuBinESyst[0][ib] = fsfuBin[0][ib]*sqrt(0.05*0.05 + 0.05*0.05 + noEffSyst0*noEffSyst0 + csEffSyst0*csEffSyst0);
    fsfuBinESyst[1][ib] = fsfuBin[1][ib]*sqrt(0.05*0.05 + 0.05*0.05 + noEffSyst1*noEffSyst1 + csEffSyst1*csEffSyst1);
    fsfuBinETot[0][ib] = sqrt( fsfuBinEStat[0][ib]*fsfuBinEStat[0][ib] + fsfuBinESyst[0][ib]*fsfuBinESyst[0][ib] );
    fsfuBinETot[1][ib] = sqrt( fsfuBinEStat[1][ib]*fsfuBinEStat[1][ib] + fsfuBinESyst[1][ib]*fsfuBinESyst[1][ib] );
    cout << "  fs/fu : barrel = " << fsfuBin[0][ib] << " +- " << fsfuBinEStat[0][ib] << " +- " << fsfuBinESyst[0][ib] << 
                     " endcap = " << fsfuBin[1][ib] << " +- " << fsfuBinEStat[1][ib] << " +- " << fsfuBinESyst[1][ib] <<endl;
    cout << endl;
    cout << endl;
  }


cout << endl;
cout << "putting it all in a nicer format" << endl;
cout << "Barrel\t\t Total\t\t\t\t Bin1\t\t\t\t Bin2\t\t\t\t Bin3\t\t\t\t Bin4" << endl;
cout << "No yield\t" << fitNoYield[0] << " +- " << fitNoYieldE[0] << " +- " << 0.05*fitNoYield[0]<< "\t\t" << fitNoYieldBin[0][0] << " +- " << fitNoYieldEBin[0][0] << " +- " << 0.05*fitNoYieldBin[0][0]<<  "\t\t" << fitNoYieldBin[0][1] << " +- " << fitNoYieldEBin[0][1] << " +- " << 0.05*fitNoYieldBin[0][1]<<  "\t\t" << fitNoYieldBin[0][2] << " +- " << fitNoYieldEBin[0][2] << " +- " << 0.05*fitNoYieldBin[0][2]<< "\t\t" << fitNoYieldBin[0][3] << " +- " << fitNoYieldEBin[0][3] << " +- " << 0.05*fitNoYieldBin[0][3]<< endl;
cout << "Cs yield\t" << fitCsYield[0] << " +- " << fitCsYieldE[0] << " +- " << 0.05*fitNoYield[0]<< "\t\t" << fitCsYieldBin[0][0] << " +- " << fitCsYieldEBin[0][0] << " +- " << 0.05*fitCsYieldBin[0][0]<<  "\t\t" << fitCsYieldBin[0][1] << " +- " << fitCsYieldEBin[0][1] << " +- " << 0.05*fitCsYieldBin[0][1]<<  "\t\t" << fitCsYieldBin[0][2] << " +- " << fitCsYieldEBin[0][2] << " +- " << 0.05*fitCsYieldBin[0][2]<< "\t\t" << fitCsYieldBin[0][3] << " +- " << fitCsYieldEBin[0][3] << " +- " << 0.05*fitCsYieldBin[0][3]<< endl;
cout << "No eff.\t\t" << effNo[0] << " +- " << effNoE[0] << " +- " << noEffSyst0*effNo[0] << "\t" << effNoBin[0][0] << " +- " << effNoEBin[0][0] << " +- " << noEffSyst0*effNoBin[0][0] <<  "\t" << effNoBin[0][1] << " +- " << effNoEBin[0][1] << " +- " << noEffSyst0*effNoBin[0][1] <<  "\t" << effNoBin[0][2] << " +- " << effNoEBin[0][2] << " +- " << noEffSyst0*effNoBin[0][2]<< "\t" << effNoBin[0][3] << " +- " << effNoEBin[0][3] << " +- " << noEffSyst0*effNoBin[0][3]<< endl;
cout << "Cs eff.\t\t" << effCs[0] << " +- " << effCsE[0] << " +- " << csEffSyst0*effCs[0] << "\t" << effCsBin[0][0] << " +- " << effCsEBin[0][0] << " +- " << csEffSyst0*effCsBin[0][0] <<  "\t" << effCsBin[0][1] << " +- " << effCsEBin[0][1] << " +- " << csEffSyst0*effCsBin[0][1] <<  "\t" << effCsBin[0][2] << " +- " << effCsEBin[0][2] << " +- " << csEffSyst0*effCsBin[0][2]<< "\t" << effCsBin[0][3] << " +- " << effCsEBin[0][3] << " +- " << csEffSyst0*effCsBin[0][3]<< endl;
cout << "fs/fu\t\t" << fsfu[0] << " +- " << fsfuEStat[0] << " +- " << fsfuESyst[0] << "\t\t\t\t" << fsfuBin[0][0] << " +- " << fsfuBinEStat[0][0] << " +- " << fsfuBinESyst[0][0]<< "\t\t\t" << fsfuBin[0][1]<< " +- " << fsfuBinEStat[0][1] << " +- " << fsfuBinESyst[0][1] << "\t\t\t" << fsfuBin[0][2]<< " +- " << fsfuBinEStat[0][2] << " +- " << fsfuBinESyst[0][2] << "\t\t\t" << fsfuBin[0][3]<< " +- " << fsfuBinEStat[0][3] << " +- " << fsfuBinESyst[0][3] << endl;
cout << endl;
cout << endl;
cout << "Endcap\t\t Total\t\t\t\t Bin1\t\t\t\t Bin2\t\t\t\t Bin3\t\t\t\t Bin4" << endl;
cout << "No yield\t" << fitNoYield[1] << " +- " << fitNoYieldE[1] << " +- " << 0.05*fitNoYield[1]<< "\t\t" << fitNoYieldBin[1][0] << " +- " << fitNoYieldEBin[1][0] << " +- " << 0.05*fitNoYieldBin[1][0]<<  "\t\t" << fitNoYieldBin[1][1] << " +- " << fitNoYieldEBin[1][1] << " +- " << 0.05*fitNoYieldBin[1][1]<<  "\t\t\t" << fitNoYieldBin[1][2] << " +- " << fitNoYieldEBin[1][2] << " +- " << 0.05*fitNoYieldBin[1][2]<< "\t\t" << fitNoYieldBin[1][3] << " +- " << fitNoYieldEBin[1][3] << " +- " << 0.05*fitNoYieldBin[1][3]<< endl;
cout << "Cs yield\t" << fitCsYield[1] << " +- " << fitCsYieldE[1] << " +- " << 0.05*fitNoYield[1]<< "\t\t" << fitCsYieldBin[1][0] << " +- " << fitCsYieldEBin[1][0] << " +- " << 0.05*fitCsYieldBin[1][0]<<  "\t\t" << fitCsYieldBin[1][1] << " +- " << fitCsYieldEBin[1][1] << " +- " << 0.05*fitCsYieldBin[1][1]<<  "\t\t"   << fitCsYieldBin[1][2] << " +- " << fitCsYieldEBin[1][2] << " +- " << 0.05*fitCsYieldBin[1][2]<< "\t\t" << fitCsYieldBin[1][3] << " +- " << fitCsYieldEBin[1][3] << " +- " << 0.05*fitCsYieldBin[1][3]<< endl;
cout << "No eff.\t\t" << effNo[1] << " +- " << effNoE[1] << " +- " << noEffSyst1*effNo[1] << "\t" << effNoBin[1][0] << " +- " << effNoEBin[1][0] << " +- " << noEffSyst1*effNoBin[1][0] <<  "\t\t" << effNoBin[1][1] << " +- " << effNoEBin[1][1] << " +- " << noEffSyst1*effNoBin[1][1] <<  "\t" << effNoBin[1][2] << " +- " << effNoEBin[1][2] << " +- " << noEffSyst1*effNoBin[1][2]<< "\t" << effNoBin[1][3] << " +- " << effNoEBin[1][3] << " +- " << noEffSyst1*effNoBin[1][3]<< endl;
cout << "Cs eff.\t\t" << effCs[1] << " +- " << effCsE[1] << " +- " << csEffSyst1*effCs[1] << "\t" << effCsBin[1][0] << " +- " << effCsEBin[1][0] << " +- " << csEffSyst1*effCsBin[1][0] <<  "\t"   << effCsBin[1][1] << " +- " << effCsEBin[1][1] << " +- " << csEffSyst1*effCsBin[1][1] <<  "\t" << effCsBin[1][2] << " +- " << effCsEBin[1][2] << " +- " << csEffSyst1*effCsBin[1][2]<< "\t" << effCsBin[1][3] << " +- " << effCsEBin[1][3] << " +- " << csEffSyst1*effCsBin[1][3]<< endl;
cout << "fs/fu\t\t" << fsfu[1] << " +- " << fsfuEStat[1] << " +- " << fsfuESyst[1]<< "\t\t\t\t" << fsfuBin[1][0] << " +- " << fsfuBinEStat[1][0] << " +- " << fsfuBinESyst[1][0]<< "\t\t\t" << fsfuBin[1][1] << " +- " << fsfuBinEStat[1][1] << " +- " << fsfuBinESyst[1][1]<< "\t\t\t" << fsfuBin[1][2] << " +- " << fsfuBinEStat[1][2] << " +- " << fsfuBinESyst[1][2]<< "\t\t\t" << fsfuBin[1][3] << " +- " << fsfuBinEStat[1][3] << " +- " << fsfuBinESyst[1][3]<< endl;
cout << endl;
cout << endl;

//then let's make some plots.
   TCanvas c1;
   const int nbins = 4;
   float xbins[5] = {12.,16.,20.,25.,30.};
    TH1F *ratio0 = new TH1F("ratio0","ratio0",nbins,xbins);
    ratio0->SetBinContent(1,fsfuBin[0][0]);
    ratio0->SetBinError(1,fsfuBinETot[0][0]);
    cout << " bin 1 error = " << fsfuBinETot[0][0] << endl;
    ratio0->SetBinContent(2,fsfuBin[0][1]);
    ratio0->SetBinError(2,fsfuBinETot[0][1]);
    cout << " bin 1 error = " << fsfuBinETot[0][1] << endl;
    ratio0->SetBinContent(3,fsfuBin[0][2]);
    ratio0->SetBinError(3,fsfuBinETot[0][2]);
    cout << " bin 1 error = " << fsfuBinETot[0][2] << endl;
    ratio0->SetBinContent(4,fsfuBin[0][3]);
    ratio0->SetBinError(4,fsfuBinETot[0][3]);
    cout << " bin 1 error = " << fsfuBinETot[0][3] << endl;
    ratio0->GetXaxis()->SetTitle("B p_{T} (GeV)");
    ratio0->GetYaxis()->SetTitle("fs/fu");
    TH1F *ratio0b = (TH1F*) ratio0->Clone("ratio0b");
    TF1 *fitpol0 = new TF1("fitpol0","[0]", 0, 100);
    ratio0->Fit("fitpol0");
    ratio0->Draw();
    //stamp(0.18, fStampString, 0.67, fStampCms); 
    cout << "Barrel fits" << endl;
    cout << "pol0 fit = " << fitpol0->GetParameter(0) << " +- " << fitpol0->GetParError(0) << endl;

    c1.SaveAs("default/ratio0_pol0.pdf");
    TCanvas c2;


    TF1 *fitpol1 = new TF1("fitpol1","[0]+(x-22.)*[1]", 0, 100);
    
    ratio0b->Fit("fitpol1");
    cout << "pol1 fit = " << fitpol1->GetParameter(0) << " +- " << fitpol1->GetParError(0) << "   " <<
       fitpol1->GetParameter(1) << " +- " << fitpol1->GetParError(1) <<endl;
    //stamp(0.18, fStampString, 0.67, fStampCms); 
    c2.SaveAs("default/ratio0_pol1.pdf");

    TCanvas c3;

    TH1F *ratio1 = new TH1F("ratio1","ratio1",nbins,xbins);
    ratio1->SetBinContent(1,fsfuBin[1][0]);
    ratio1->SetBinError(1,fsfuBinETot[1][0]);
    ratio1->SetBinContent(2,fsfuBin[1][1]);
    ratio1->SetBinError(2,fsfuBinETot[1][1]);
    ratio1->SetBinContent(3,fsfuBin[1][2]);
    ratio1->SetBinError(3,fsfuBinETot[1][2]);
    ratio1->SetBinContent(4,fsfuBin[1][3]);
    ratio1->SetBinError(4,fsfuBinETot[1][3]);
    ratio1->GetXaxis()->SetTitle("B p_{T} (GeV)");
    ratio1->GetYaxis()->SetTitle("fs/fu");
    TH1F *ratio1b = (TH1F*) ratio1->Clone("ratio1b");





    ratio1->Fit("fitpol0");
    ratio1->Draw();

    cout << "Endcap fits" << endl;
    cout << "pol0 fit = " << fitpol0->GetParameter(0) << " +- " << fitpol0->GetParError(0) << endl;
    //stamp(0.18, fStampString, 0.67, fStampCms); 
    c3.SaveAs("default/ratio1_pol0.pdf");

    TCanvas c4;
    ratio1b->Fit("fitpol1");
    cout << "pol1 fit = " << fitpol1->GetParameter(0) << " +- " << fitpol1->GetParError(0) << "   " <<
       fitpol1->GetParameter(1) << " +- " << fitpol1->GetParError(1) <<endl;
    //stamp(0.18, fStampString, 0.67, fStampCms); 
    c4.SaveAs("default/ratio1_pol1.pdf");

}

// ----------------------------------------------------------------------
void plotRatio::loopFunction(int mode) {
  if (1 == mode) loopFunction1(0);
  if (2 == mode) loopFunction2(0);
}
// ----------------------------------------------------------------------
void plotRatio::loopFunction(int mode, int imode) {
  if (1 == mode) loopFunction1(imode);
  if (2 == mode) loopFunction2(imode);
}


// ----------------------------------------------------------------------
void plotRatio::loopFunction1(int mode) {
 
// get efficiency for MC and mass distributions for data

  // first set up pt binning.
  int ptbin = 0;
  if (fb.pt>=16 && fb.pt<20) ptbin = 1;
  if (fb.pt>=20 && fb.pt<25) ptbin = 2;
  if (fb.pt>=25) ptbin = 3;

  bool bp2jpsikp(false), bs2jpsiphi(false); 
  if (10 == mode)  bp2jpsikp = true;
  if (20 == mode)  bs2jpsiphi = true;
  						
  // -- channel index									      
					      								      
  fChan = detChan(fb.m1eta, fb.m2eta);  	
  if (fChan < 0)  return;


							  
  if (!fIsMC) {										  
    if (!fb.json) return;							  
  }											  

	
//then apply analysis cuts
// first, if BDT, muon id and hlt
  if(fDoUseBDT) {
    if ( (fBDT > fCuts[fChan]->bdt) && fb.gmuid && fb.hlt && fb.gtqual) {    
    // count efficiency numerator					
      if(fIsMC) {							
    	if (bp2jpsikp) {
    	  effNoNum[fChan]++;				
    	  effNoNumBin[fChan][ptbin]++;
    	}
    	if (bs2jpsiphi) {
    	  effCsNum[fChan]++;
    	  effCsNumBin[fChan][ptbin]++;
    	}
    	//also fill a histo of B+ pt and Bs pt to check distributions
    	if (bp2jpsikp) hNoPt[fChan]->Fill(fb.pt);
    	if (bs2jpsiphi) hCsPt[fChan]->Fill(fb.pt);			      
      } else {  							
      // fill mass plots for data					
    	if(bs2jpsiphi) {
    	  hCsSignal[fChan]->Fill(fb.m); //fb.cm is constrained mass instead
    	  hCsSignalC[fChan]->Fill(fb.cm); //fb.cm is constrained mass instead
    	  hCsSignalBin[fChan][ptbin]->Fill(fb.m);
	  hCsSignalCBin[fChan][ptbin]->Fill(fb.cm);
    	}
    	if(bp2jpsikp) {
    	  hNoSignal[fChan]->Fill(fb.m); 		
    	  hNoSignalC[fChan]->Fill(fb.cm); 		
    	  hNoSignalBin[fChan][ptbin]->Fill(fb.m);
	  hNoSignalCBin[fChan][ptbin]->Fill(fb.cm);
    	}
      } 								
    } // end passing BDT cuts
  } else {
  // apply cut-based selection instead
  
  
  cuts *pCuts(0);
  pCuts = fCuts[fChan]; 
    //cuts copied directly from plotClass loopTree method  
      // -- gen-level acceptance cuts
    if (fIsMC) {
      if (TMath::Abs(fb.g1eta) > 2.5) return;
      if (TMath::Abs(fb.g2eta) > 2.5) return;
      if (fb.g1pt < 1.0) return;
      if (fb.g2pt < 1.0) return;
      if (bp2jpsikp) {
	// gen-level cuts for Bu2JpsiKp
	if (fb.g1pt < 3.5) return;
	if (fb.g2pt < 3.5) return;
	if (TMath::Abs(fb.g3eta) > 2.5) return;
	if (fb.g3pt < 0.4) return;
      }
      
      if (bs2jpsiphi) {
	if (TMath::Abs(fb.g3eta) > 2.5) return;
	if (TMath::Abs(fb.g4eta) > 2.5) return;
	// gen-level cuts for Bs2JpsiPhi
	if (fb.g1pt < 3.5) return;
	if (fb.g2pt < 3.5) return;
	if (fb.g3pt < 0.4) return;
	if (fb.g4pt < 0.4) return;
      }
    } else {
      if (!fb.json) return;
    }
    // -- immutable cuts: require basic muon and trackQual cuts
    if (false == fb.gtqual)  return;
    if (TMath::Abs(fb.m1eta) > 2.4) return;
    if (TMath::Abs(fb.m2eta) > 2.4) return;

    if (fb.m1pt < 1.0) return;
    if (fb.m2pt < 1.0) return;

    if (bp2jpsikp) {
      if (TMath::Abs(fb.k1eta) > 2.4) return;
      if (fb.k1pt < 0.5) return;
    }
    if (bs2jpsiphi) {
      if (TMath::Abs(fb.k1eta) > 2.4) return;
      if (TMath::Abs(fb.k2eta) > 2.4) return;
      if (fb.k1pt < 0.5) return;
      if (fb.k2pt < 0.5) return;
    }

    // -- analysis cuts 				  
    if (fb.m1q*fb.m2q > 0) return;			  
    if (fb.m1pt < pCuts->m1pt) return;  		  
    if (fb.m2pt < pCuts->m2pt) return;  		  
    							  
    if (fb.fl3d > 1.5) return;				  
    if (fb.pvw8 < 0.6) return;				  
    							  
    if (fb.pt < pCuts->pt) return; 			  
    if (fb.iso < pCuts->iso) return; 			  
    if (fb.chi2/fb.dof > pCuts->chi2dof) return;	  
    if (TMath::IsNaN(fb.fls3d)) return; 		  
    if (fb.fls3d < pCuts->fls3d) return;		  
    if (fb.alpha > pCuts->alpha) return;		  
    if (fb.docatrk < pCuts->docatrk) return;		  
    							  
    if (bs2jpsiphi && fb.dr > 0.25) return;		  
    if (bs2jpsiphi && fb.mkk < 0.995) return;		  
    if (bs2jpsiphi && fb.mkk > 1.045) return;		  
    							  
    // -- new cuts					  
    if (fb.closetrk >= pCuts->closetrk) return; 	  
    if (TMath::Abs(fb.pvlip) > pCuts->pvlip) return;	  
    if (TMath::Abs(fb.pvlips) > pCuts->pvlips) return;	  
    if (TMath::Abs(fb.pvlip2) < pCuts->pvlip2) return;	  
    if (TMath::Abs(fb.pvlips2) < pCuts->pvlips2) return;  
    if (fb.maxdoca > pCuts->maxdoca) return;		  
    if (fb.pvips > pCuts->pvips) return;		  
    if (fb.pvip > pCuts->pvip) return;			  
    							  
    if (bs2jpsiphi || bp2jpsikp) {			  
      if (fb.mpsi > 3.2) return;
      if (fb.mpsi < 3.0) return;
      // -- cowboy veto 
      if (fb.psipt < 7) return;
    } 			
    				  
    //also require good muon id and the trigger:
    if (fb.gmuid == false) return;
    if (fb.hlt == false) return;

//passed all cuts

    // count efficiency numerator					
      if(fIsMC) {							
    	if (bp2jpsikp) {
    	  effNoNum[fChan]++;				
    	  effNoNumBin[fChan][ptbin]++;
    	}
    	if (bs2jpsiphi) {
    	  effCsNum[fChan]++;
    	  effCsNumBin[fChan][ptbin]++;
    	}
    	//also fill a histo of B+ pt and Bs pt to check distributions
    	if (bp2jpsikp) hNoPt[fChan]->Fill(fb.pt);
    	if (bs2jpsiphi) hCsPt[fChan]->Fill(fb.pt);			      
      } else {  							
      // fill mass plots for data					
    	if(bs2jpsiphi) {
    	  hCsSignal[fChan]->Fill(fb.m); //or is it fb.cm?   constrained mass?
    	  hCsSignalBin[fChan][ptbin]->Fill(fb.m);
    	}
    	if(bp2jpsikp) {
    	  hNoSignal[fChan]->Fill(fb.m); 		
    	  hNoSignalBin[fChan][ptbin]->Fill(fb.m);
    	}
      } 								
  } // end using cut-based selection
}
// ----------------------------------------------------------------------
void plotRatio::loopFunction2(int imode) {
  //not using this currently...
  return;

}
