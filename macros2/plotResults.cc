#include "plotResults.hh"
#include "plotOverlays.hh"
#include "plotPU.hh"
#include "plotEfficiencies.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "../macros/bayesianlimit.hh"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "THStack.h"

using namespace std; 
using std::string; 

ClassImp(plotResults)


// ----------------------------------------------------------------------
plotResults::plotResults(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

  fDoPrint = true; 

  fNumbersFileName = fDirectory + "/anaBmm.plotResults." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  printCuts(cout); 

  fBF = 0.0593*1.014e-3;
  // PDG 2010:
  fu  = 0.401;
  fs  = 0.113;
  // -- CMS with PDG input
  fsfu = 0.282;
  fsfuE = 0.037/0.282;

  fDoApplyCowboyVeto = true;   
  fDoApplyCowboyVetoAlsoInSignal = false;   

  fNormProcessed = false; 
}

// ----------------------------------------------------------------------
plotResults::~plotResults() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotResults::makeAll(int channels) {

  fNumbersNo[0]->effTot = 0.00053;
  fNumbersNo[1]->effTot = 0.00019;
  fNumbersNo[0]->fitYield = 36208;
  fNumbersNo[1]->fitYield = 12758;

  rareBg(); 
  return;


  if (channels & 2) {
    plotOverlays a1; 
    a1.makeAll(); 
  }

  if (channels & 4) {
    plotPU a2; 
    a2.makeAll(); 
  }

  if (channels & 8) {
    plotEfficiencies a3; 
    a3.makeAll(); 
  }

  if (channels & 1) {
    fNormProcessed = false; 
    fDoUseBDT = false; 
    fDoApplyCowboyVeto = true;   
    fDoApplyCowboyVetoAlsoInSignal = false;   
    computeNormUL();
    computeCsBF();

    fNormProcessed = false; 
    fDoUseBDT = true; 
    fDoApplyCowboyVeto = true;   
    fDoApplyCowboyVetoAlsoInSignal = false;   
    computeNormUL();
    computeCsBF();
  }
}

// ----------------------------------------------------------------------
void plotResults::computeNormUL() {

  cout << "--> Normalization" << endl;
  if (false == fNormProcessed) {
    fNormProcessed = true; 
    cout << "--> loopTree: norm MC" << endl;
    loopTree(10); // normalization eff
    c0->Modified(); c0->Update();
    cout << "--> loopTree: norm data" << endl;
    loopTree(15); // data normalization 
    c0->Modified(); c0->Update();
  }

  cout << "--> rareBg" << endl;
  rareBg();

  cout << "--> loopTree: signal MC" << endl;
  loopTree(0);  // signal eff
  c0->Modified(); c0->Update();
  loopTree(1);  // Bd2MuMu eff
  c0->Modified(); c0->Update();
  cout << "--> loopTree: signal data" << endl;
  loopTree(5);  // data signal
  c0->Modified(); c0->Update();


  string bla = fDirectory + "/anaBmm.plotResults." + fSuffix;
  if (fDoUseBDT) {
    bla += ".bdt.ulc";
  } else {
    bla += ".cnc.ulc";
  }
  cout << "===> Storing ULC numbers in file " << bla << endl;
  system(Form("/bin/rm -f %s", bla.c_str()));
  printUlcalcNumbers(bla);

  double   fNul = 0.;

  cout << "printing fNumbersBs[0]" << endl;
  printNumbers(*fNumbersBs[0], cout); 
  printNumbers(*fNumbersBs[0], fOUT); 

  cout << "printing fNumbersBd[0]" << endl;
  printNumbers(*fNumbersBd[0], cout); 
  printNumbers(*fNumbersBd[0], fOUT); 

  cout << "printing fNumbersNorm[0]" << endl;
  printNumbers(*fNumbersNo[0], cout); 
  printNumbers(*fNumbersNo[0], fOUT); 

  double fUL = (fNumbersBs[0]->bgBsExp/fNumbersNo[0]->fitYield) //????FIXME
    *(fu/fs)
    *(fNumbersNo[0]->acc/fNumbersBs[0]->acc)
    *(fNumbersNo[0]->effCand/fNumbersBs[0]->effCand)     
    *(fNumbersNo[0]->effMuidMC/fNumbersBs[0]->effMuidMC)
    *(fNumbersNo[0]->effTrigMC/fNumbersBs[0]->effTrigMC)
    *(fNumbersNo[0]->effAna/fNumbersBs[0]->effAna)
    * fBF;

  cout << "prod(eff) expected UL: " << fUL << endl;

  fUL = (fNul/fNumbersNo[0]->fitYield)
    *(fu/fs)
    *(fNumbersNo[0]->effTot/fNumbersBs[0]->effTot)
    * fBF;

  cout << "effTot expected UL:    " << fUL << endl;

  //  system(Form("../ulcalc/bin/ulcalc %s", fUlcalcFileName.c_str())); 

}


// ----------------------------------------------------------------------
void plotResults::computeCsBF() {
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
      * fBF;

    
    resultE = dRatio(fNumbersCs[i]->fitYield, fNumbersCs[i]->fitYieldE, fNumbersNo[i]->fitYield, fNumbersNo[i]->fitYieldE)
      *(fu/fs)
      *(fNumbersNo[i]->acc/fNumbersCs[i]->acc)
      *(fNumbersNo[i]->effCand/fNumbersCs[i]->effCand)     
      *(fNumbersNo[i]->effMuidTNP/fNumbersCs[i]->effMuidTNP)
      *(fNumbersNo[i]->effTrigTNP/fNumbersCs[i]->effTrigTNP)
      *(fNumbersNo[i]->effAna/fNumbersCs[i]->effAna)
      * fBF;

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
      * fBF;

    
    resultE = dRatio(fNumbersCs[i]->fitYield, fNumbersCs[i]->fitYieldE, fNumbersNo[i]->fitYield, fNumbersNo[i]->fitYieldE)
      *(fu/fs)
      *(fNumbersNo[i]->acc/fNumbersCs[i]->acc)
      *(fNumbersNo[i]->effCand/fNumbersCs[i]->effCand)     
      *(fNumbersNo[i]->effMuidMC/fNumbersCs[i]->effMuidMC)
      *(fNumbersNo[i]->effTrigMC/fNumbersCs[i]->effTrigMC)
      *(fNumbersNo[i]->effAna/fNumbersCs[i]->effAna)
      * fBF;

    cout << "chan " << i << ": MC fact branching fraction: " << result << "+/-" << resultE << endl;
    fTEX << formatTex(result, Form("%s:N-CSBF-MC-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(resultE, Form("%s:N-CSBF-MC-BS%i:err", fSuffix.c_str(), i), 6) << endl;
    
    result = (fNumbersCs[i]->fitYield/fNumbersNo[i]->fitYield)
      *(fu/fs)
      *(fNumbersNo[i]->effTot/fNumbersCs[i]->effTot)
      * fBF;

    cout << "chan " << i << ": branching fraction: " << result << "+/-" << resultE << endl;

    fTEX << formatTex(result, Form("%s:N-CSBF-BS%i:val", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(resultE, Form("%s:N-CSBF-BS%i:err", fSuffix.c_str(), i), 6) << endl;

  }

  fSuffix = cache; 

  printCsBFNumbers();
}



// ----------------------------------------------------------------------
void plotResults::printUlcalcNumbers(string fname) {
  ofstream OUT(fname.c_str());

  OUT << "######################################################################" << endl;
  fTEX << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  //  double sysMu(0.05/TMath::Sqrt(2.)), sysTr(0.02/TMath::Sqrt(2.)), 
  double sysMu(0.05), sysTr(0.02), 
    sysAna(0.08), sysAnaNo(0.04), sysCand(0.01), sysAcc(0.04), sysNorm(0.05), 
    sysPSS(0.05), sysTot(-1.); 
  double totE, err1, err2; 

  string cache = fSuffix; 
  if (fDoUseBDT) fSuffix = "bdt" + fSuffix; 
  
  if (1) {
    sysTot = 0.05*0.05 + 0.02*0.02 + sysAna*sysAna + sysCand*sysCand + sysAcc*sysAcc;
    sysTot = TMath::Sqrt(sysTot); 

    for (unsigned int i = 0; i < fNchan; ++i) {
      OUT << "# -- NORMALIZATION " << i << endl;
      fTEX << "% -- NORMALIZATION " << i << endl;

      OUT << "#EFF_TOT_BPLUS\t" << i << "\t" << fNumbersNo[i]->effTot << endl;
      err1 = fNumbersNo[i]->effTotE; 
      err2 = sysTot*fNumbersNo[i]->effTot;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effTot, Form("%s:N-EFF-TOT-BPLUS%i:val", fSuffix.c_str(), i), 5) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-TOT-BPLUS%i:err", fSuffix.c_str(), i), 5) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-TOT-BPLUS%i:sys", fSuffix.c_str(), i), 5) << endl;
      fTEX << formatTex(totE, Form("%s:N-EFF-TOT-BPLUS%i:tot", fSuffix.c_str(), i), 5) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTot, totE, Form("%s:N-EFF-TOT-BPLUS%i:all", fSuffix.c_str(), i), 1e-3, 2) 
	   << endl;

      OUT << "ACC_BPLUS\t" << i << "\t" << fNumbersNo[i]->acc 
	  <<"\t" << sysAcc*fNumbersNo[i]->acc 
	  << endl;
      err1 = fNumbersNo[i]->accE; 
      err2 = sysAcc*fNumbersNo[i]->acc;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->acc, Form("%s:N-ACC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err1, Form("%s:N-ACC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err2, Form("%s:N-ACC-BPLUS%i:sys", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(totE, Form("%s:N-ACC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->acc, totE, Form("%s:N-ACC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      err1 = fNumbersNo[i]->effMuidTNPE; 
      err2 = sysMu*fNumbersNo[i]->effMuidTNP;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effMuidTNP, Form("%s:N-EFF-MU-PID-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-MU-PID-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-MU-PID-BPLUS%i:sys", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(totE, Form("%s:N-EFF-MU-PID-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effMuidTNP, totE, Form("%s:N-EFF-MU-PID-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "EFF_MU_BPLUS\t" << i << "\t"  << fNumbersNo[i]->effMuidMC 
	  << "\t"  << sysMu*fNumbersNo[i]->effMuidMC 
	  << endl;
      err1 = fNumbersNo[i]->effMuidMCE; 
      err2 = sysMu*fNumbersNo[i]->effMuidMC;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:sys", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(totE, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effMuidTNPMC, totE, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      err1 = fNumbersNo[i]->effMuidMCE; 
      err2 = sysNorm*fNumbersNo[i]->effMuidMC;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-MU-MC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-MU-MC-BPLUS%i:sys", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(totE, Form("%s:N-EFF-MU-MC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effMuidMC, totE, Form("%s:N-EFF-MU-MC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      //    OUT << "EFF_TRIG_BPLUS\t" << i << "\t" << fNumbersNo[i]->effTrigTNP << endl;
      err1 = fNumbersNo[i]->effTrigTNPE; 
      err2 = sysTr*fNumbersNo[i]->effTrigTNP;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effTrigTNP, Form("%s:N-EFF-TRIG-PID-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-PID-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-PID-BPLUS%i:sys", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(totE, Form("%s:N-EFF-TRIG-PID-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTrigTNP, totE, Form("%s:N-EFF-TRIG-PID-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      err1 = fNumbersNo[i]->effTrigTNPMCE; 
      err2 = sysTr*fNumbersNo[i]->effTrigTNPMC;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:sys", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTrigTNPMC, totE, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "EFF_TRIG_BPLUS\t" << i << "\t" << fNumbersNo[i]->effTrigMC 
	  << "\t" << sysTr*fNumbersNo[i]->effTrigMC 
	  << endl;
      err1 = fNumbersNo[i]->effTrigMCE; 
      err2 = sysTr*fNumbersNo[i]->effTrigMC;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-MC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-MC-BPLUS%i:sys", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(totE, Form("%s:N-EFF-TRIG-MC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTrigMC, totE, Form("%s:N-EFF-TRIG-MC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "EFF_CAND_BPLUS\t" << i << "\t" << fNumbersNo[i]->effCand
	  << "\t" << sysCand*fNumbersNo[i]->effCand
	  << endl;
      err1 = fNumbersNo[i]->effCandE; 
      err2 = sysCand*fNumbersNo[i]->effCand;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effCand, Form("%s:N-EFF-CAND-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-CAND-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-CAND-BPLUS%i:sys", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(totE, Form("%s:N-EFF-CAND-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effCand, totE, Form("%s:N-EFF-CAND-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      //    OUT << "EFF_ANA_BPLUS\t" << i << "\t" << fNumbersNo[i]->effAna << endl;
      OUT << "EFF_ANA_BPLUS\t" << i << "\t" << fNumbersNo[i]->effAna
	  << "\t" << sysAnaNo*fNumbersNo[i]->effAna 
	  << endl;
      err1 = fNumbersNo[i]->effAnaE; 
      err2 = sysAnaNo*fNumbersNo[i]->effAna;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->effAna, Form("%s:N-EFF-ANA-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(err1, Form("%s:N-EFF-ANA-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(err2, Form("%s:N-EFF-ANA-BPLUS%i:sys", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(totE, Form("%s:N-EFF-ANA-BPLUS%i:tot", fSuffix.c_str(), i), 4) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effAna, totE, Form("%s:N-EFF-ANA-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "OBS_BPLUS\t" << i << "\t" << fNumbersNo[i]->fitYield
	  << "\t" << sysNorm*fNumbersNo[i]->fitYield 
	  << endl;
      err1 = fNumbersNo[i]->fitYieldE; 
      err2 = sysNorm*fNumbersNo[i]->fitYield;
      totE = TMath::Sqrt(err1*err1 + err2*err2); 
      fTEX << formatTex(fNumbersNo[i]->fitYield, Form("%s:N-OBS-BPLUS%i:val", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(err1, Form("%s:N-OBS-BPLUS%i:err", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(err2, Form("%s:N-OBS-BPLUS%i:sys", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(totE, Form("%s:N-OBS-BPLUS%i:tot", fSuffix.c_str(), i), 0) << endl;
      fTEX << scientificTex(fNumbersNo[i]->fitYield, totE, Form("%s:N-OBS-BPLUS%i:all", fSuffix.c_str(), i), 1e3, 0) << endl;
    } 
  } else {
    //     OUT << "TOT_BPLUS\t" << "0\t" << (fDataLumi[fSgData]/39.4)*440000 << endl;
    //     OUT << "TOT_BPLUS\t" << "1\t" << (fDataLumi[fSgData]/39.4)*383000 << endl;
  }

  OUT << "######################################################################" << endl;
  fTEX << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  double scale(0.), scaledSig, scaledSigE(0.), scaledSigS(0.); 
  for (unsigned int i = 0; i < fNchan; ++i) {
    OUT << "# -- SIGNAL " << i << endl;
    fTEX << "% -- SIGNAL " << i << endl;
    scale      = fLumi["SgData"]/fLumi["SgMc"];
    scaledSig  = fNumbersBs[i]->anaWmcYield*scale;
    scaledSigE = scaledSig*TMath::Sqrt(fNumbersBs[i]->anaWmcYield)/fNumbersBs[i]->anaWmcYield;
    fTEX << formatTex(fNumbersBs[i]->anaWmcYield, Form("%s:N-EXP-SIG-BSMM%d:wmcyield", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(scale, Form("%s:N-EXP-SIG-BSMM%d:scale", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(scaledSig, Form("%s:N-EXP-SIG-BSMM%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(scaledSigE, Form("%s:N-EXP-SIG-BSMM%d:err", fSuffix.c_str(), i), 2) << endl;
    //    fTEX << formatTex(0.2*scaledSig, Form("%s:N-EXP-SIG-BSMM%d:sys", fSuffix.c_str(), i), 2) << endl;

    double yield = scaledYield(fNumbersBs[i], fNumbersNo[i], 3.2e-9, fsfu);

    //          BF(Bs -> mu mu)   fs epstot(Bs) 
    //   n_s = -----------------  -- ---------  N(B+) 
    //          BF(B+ -> mu muK)  fu epstot(B+) 

    //                   mass reso   signal eff  norm eff    kaon          norm fit    muid        trigger
    double commonError = TMath::Sqrt(0.03*0.03 + 0.08*0.08 + 0.04*0.04 + 0.039*0.039 + 0.05*0.05 + 0.05*0.05 + 0.03*0.03);
    cout << "****** commonError = " << commonError << endl;

    double bf = 3.2e-9/6.0e-5;
    scaledSig  = bf * (fsfu) * fNumbersBs[i]->pss * (fNumbersBs[i]->effTot/fNumbersNo[i]->effTot) * fNumbersNo[i]->fitYield;
    scaledSigE = TMath::Sqrt(fsfuE*fsfuE + commonError*commonError + (0.2/3.2)*(0.2/3.2))*scaledSig;
    cout << "****** scaledSig(orig) = " << scaledSig << endl;

    scaledSig  = fNumbersBs[i]->bsNoScaled;
    scaledSigE = fNumbersBs[i]->bsNoScaledE;
    cout << "****** scaledSig(Bs) =   " << scaledSig << endl;
    OUT << "#EXP_SIG_BSMM\t" << i << "\t" << fNumbersBs[i]->bsNoScaled << endl;

    fTEX << formatTex(scaledSig, Form("%s:N-EXP2-SIG-BSMM%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(scaledSigE, Form("%s:N-EXP2-SIG-BSMM%d:err", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(scaledSigS, Form("%s:N-EXP2-SIG-BSMM%d:sys", fSuffix.c_str(), i), 2) << endl;

    scaledSig  = fNumbersBd[i]->anaWmcYield*fLumi["SgData"]/fLumi["BdMc"];
    scaledSigE = scaledSig*TMath::Sqrt(fNumbersBd[i]->anaWmcYield)/fNumbersBd[i]->anaWmcYield;
    fTEX << formatTex(scaledSig, Form("%s:N-EXP-SIG-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(scaledSigE, Form("%s:N-EXP-SIG-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    
    bf = 1.0e-10/6.0e-5;
    scaledSig  = bf * fNumbersBd[i]->pdd * (fNumbersBd[i]->effTot/fNumbersNo[i]->effTot) * fNumbersNo[i]->fitYield;
    scaledSigE = TMath::Sqrt(commonError*commonError + (0.1/1.0)*(0.1/1.0))*scaledSig;
    cout << "****** scaledSig(orig) = " << scaledSig << endl;

    yield = scaledYield(fNumbersBd[i], fNumbersNo[i], 1.0e-10, 1.);
    scaledSig  = fNumbersBd[i]->bdNoScaled;
    scaledSigE = fNumbersBd[i]->bdNoScaledE;
    cout << "****** scaledSig(Bd) =   " << scaledSig << endl;
    OUT << "#EXP_SIG_BDMM\t" << i << "\t" << fNumbersBd[i]->bdNoScaled << endl;

    fTEX << formatTex(scaledSig, Form("%s:N-EXP2-SIG-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(scaledSigE, Form("%s:N-EXP2-SIG-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;

    OUT << "OBS_BKG\t" << i << "\t" << fNumbersBs[i]->bgObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bgObs, Form("%s:N-OBS-BKG%d:val", fSuffix.c_str(), i), 0) << endl;
    fTEX << formatTex(fNumbersBs[i]->bgBsExp, Form("%s:N-EXP-BSMM%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bgBsExpE, Form("%s:N-EXP-BSMM%d:err", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bgBdExp, Form("%s:N-EXP-BDMM%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bgBdExpE, Form("%s:N-EXP-BDMM%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "LOW_BD\t" << i << "\t" << fNumbersBs[i]->mBdLo << endl;
    fTEX << formatTex(fNumbersBs[i]->mBdLo, Form("%s:N-LOW-BD%d:val", fSuffix.c_str(), i), 3) << endl;
    
    OUT << "HIGH_BD\t" << i << "\t" << fNumbersBs[i]->mBdHi << endl;
    fTEX << formatTex(fNumbersBs[i]->mBdHi, Form("%s:N-HIGH-BD%d:val", fSuffix.c_str(), i), 3) << endl;

    OUT << "LOW_BS\t" << i << "\t" << fNumbersBs[i]->mBsLo << endl;
    fTEX << formatTex(fNumbersBs[i]->mBsLo, Form("%s:N-LOW-BS%d:val", fSuffix.c_str(), i), 3) << endl;

    OUT << "HIGH_BS\t" << i << "\t" << fNumbersBs[i]->mBsHi << endl;
    fTEX << formatTex(fNumbersBs[i]->mBsHi, Form("%s:N-HIGH-BS%d:val", fSuffix.c_str(), i), 3) << endl;

    OUT << "PSS\t" << i << "\t" << fNumbersBs[i]->pss
	<< "\t" << sysPSS*fNumbersBs[i]->pss 
	<< endl;
    fTEX << formatTex(fNumbersBs[i]->pss, Form("%s:N-PSS%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->pssE, Form("%s:N-PSS%d:err", fSuffix.c_str(), i), 3) << endl;
    //    fTEX << formatTex(sysPSS*fNumbersBs[i]->pss, Form("%s:N-PSS%d:sys", fSuffix.c_str(), i), 3) << endl;

    OUT << "PSD\t" << i << "\t" << fNumbersBd[i]->psd
	<< "\t" << sysPSS*fNumbersBd[i]->psd 
	<< endl;
    fTEX << formatTex(fNumbersBd[i]->psd, Form("%s:N-PSD%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->psdE, Form("%s:N-PSD%d:err", fSuffix.c_str(), i), 3) << endl;
    //    fTEX << formatTex(sysPSS*fNumbersBd[i]->psdE, Form("%s:N-PSD%d:sys", fSuffix.c_str(), i), 3) << endl;

    OUT << "PDS\t" << i << "\t" << fNumbersBs[i]->pds
	<< "\t" << sysPSS*fNumbersBs[i]->pds 
	<< endl;
    fTEX << formatTex(fNumbersBs[i]->pds, Form("%s:N-PDS%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->pdsE, Form("%s:N-PDS%d:err", fSuffix.c_str(), i), 3) << endl;
    //    fTEX << formatTex(sysPSS*fNumbersBs[i]->pds, Form("%s:N-PDS%d:sys", fSuffix.c_str(), i), 3) << endl;

    OUT << "PDD\t" << i << "\t" << fNumbersBd[i]->pdd
	<< "\t" << sysPSS*fNumbersBd[i]->pdd 
	<< endl;
    fTEX << formatTex(fNumbersBd[i]->pdd, Form("%s:N-PDD%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->pddE, Form("%s:N-PDD%d:err", fSuffix.c_str(), i), 3) << endl;
    //    fTEX << formatTex(sysPSS*fNumbersBd[i]->pdd, Form("%s:N-PDD%d:sys", fSuffix.c_str(), i), 3) << endl;

    OUT << "#EFF_TOT_BSMM\t" << i << "\t" << fNumbersBs[i]->effTot << endl;
    err1 = fNumbersBs[i]->effTotE; 
    err2 = sysTot*fNumbersBs[i]->effTot;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effTot, Form("%s:N-EFF-TOT-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-TOT-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-TOT-BSMM%d:sys", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-TOT-BSMM%d:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTot, totE, Form("%s:N-EFF-TOT-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "ACC_BSMM\t" << i << "\t" << fNumbersBs[i]->acc
	<< "\t" << sysAcc*fNumbersBs[i]->acc 
	<< endl;
    err1 = fNumbersBs[i]->accE; 
    err2 = sysAcc*fNumbersBs[i]->acc;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->acc, Form("%s:N-ACC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-ACC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-ACC-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-ACC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->acc, totE, Form("%s:N-ACC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    //    OUT << "EFF_MU_BSMM\t" << i << "\t" << fNumbersBs[i]->effMuidTNP << endl;
    err1 = fNumbersBs[i]->effMuidTNPE; 
    err2 = sysMu*fNumbersBs[i]->effMuidTNP;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effMuidTNP, Form("%s:N-EFF-MU-PID-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-MU-PID-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-MU-PID-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-MU-PID-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effMuidTNP, totE, Form("%s:N-EFF-MU-PID-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    err1 = fNumbersBs[i]->effMuidTNPMCE; 
    err2 = sysMu*fNumbersBs[i]->effMuidTNPMC;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-MU-PIDMC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-MU-PIDMC-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-MU-PIDMC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effMuidTNPMC, totE, Form("%s:N-EFF-MU-PIDMC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_MU_BSMM\t" << i << "\t" << fNumbersBs[i]->effMuidMC 
	<< "\t" << sysMu*fNumbersBs[i]->effMuidMC << endl;
    err1 = fNumbersBs[i]->effMuidMCE; 
    err2 = sysMu*fNumbersBs[i]->effMuidMC;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-MU-MC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-MU-MC-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-MU-MC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effMuidMC, totE, Form("%s:N-EFF-MU-MC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    //    OUT << "EFF_TRIG_BSMM\t" << i << "\t" << fNumbersBs[i]->effTrigTNP << endl;
    err1 = fNumbersBs[i]->effTrigTNPE; 
    err2 = sysTr*fNumbersBs[i]->effTrigTNP;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effTrigTNP, Form("%s:N-EFF-TRIG-PID-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-PID-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-PID-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-TRIG-PID-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTrigTNP, totE, Form("%s:N-EFF-TRIG-PID-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    err1 = fNumbersBs[i]->effTrigTNPMCE; 
    err2 = sysTr*fNumbersBs[i]->effTrigTNPMC;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTrigTNPMC, totE, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_TRIG_BSMM\t" << i << "\t" << fNumbersBs[i]->effTrigMC 
	<< "\t" << sysTr*fNumbersBs[i]->effTrigMC 
	<< endl;
    err1 = fNumbersBs[i]->effTrigMCE; 
    err2 = sysTr*fNumbersBs[i]->effTrigMC;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-MC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-MC-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-TRIG-MC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTrigMC, totE, Form("%s:N-EFF-TRIG-MC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_CAND_BSMM\t" << i << "\t" << fNumbersBs[i]->effCand
	<< "\t" << sysCand*fNumbersBs[i]->effCand << endl;
    err1 = fNumbersBs[i]->effCandE; 
    err2 = sysCand*fNumbersBs[i]->effCand;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effCand, Form("%s:N-EFF-CAND-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-CAND-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-CAND-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-CAND-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effCand, totE, Form("%s:N-EFF-CAND-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_ANA_BSMM\t" << i << "\t" << fNumbersBs[i]->effAna
	<< "\t" << sysAna*fNumbersBs[i]->effAna 
	<< endl;
    err1 = fNumbersBs[i]->effAnaE; 
    err2 = sysAna*fNumbersBs[i]->effAna;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBs[i]->effAna, Form("%s:N-EFF-ANA-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-ANA-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-ANA-BSMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-ANA-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effAna, totE, Form("%s:N-EFF-ANA-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "#EFF_TOT_BDMM\t" << i << "\t" << fNumbersBd[i]->effTot << endl;
    err1 = fNumbersBs[i]->effTotE; 
    err2 = sysTot*fNumbersBs[i]->effTot;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effTot, Form("%s:N-EFF-TOT-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-TOT-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-TOT-BDMM%d:sys", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-TOT-BDMM%d:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTot, totE, Form("%s:N-EFF-TOT-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "ACC_BDMM\t" << i << "\t" << fNumbersBd[i]->acc 
	<< "\t" << sysAcc*fNumbersBd[i]->acc
	<< endl;
    err1 = fNumbersBd[i]->accE; 
    err2 = sysAcc*fNumbersBd[i]->acc;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->acc, Form("%s:N-ACC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-ACC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-ACC-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-ACC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->acc, totE, Form("%s:N-ACC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    //    OUT << "EFF_MU_BDMM\t" << i << "\t" << fNumbersBd[i]->effMuidTNP << endl;
    err1 = fNumbersBd[i]->effMuidTNPE; 
    err2 = sysMu*fNumbersBd[i]->effMuidTNP;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effMuidTNP, Form("%s:N-EFF-MU-PID-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-MU-PID-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-MU-PID-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-MU-PID-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effMuidTNP, totE, Form("%s:N-EFF-MU-PID-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    err1 = fNumbersBd[i]->effMuidTNPMCE; 
    err2 = sysMu*fNumbersBd[i]->effMuidTNPMC;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-MU-PIDMC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-MU-PIDMC-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-MU-PIDMC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effMuidTNPMC, totE, Form("%s:N-EFF-MU-PIDMC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_MU_BDMM\t" << i	<< "\t" << fNumbersBd[i]->effMuidMC
	<< "\t" << sysMu*fNumbersBd[i]->effMuidMC
	<< endl;
    err1 = fNumbersBd[i]->effMuidMCE; 
    err2 = sysMu*fNumbersBd[i]->effMuidMC;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-MU-MC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-MU-MC-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-MU-MC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effMuidMC, totE, Form("%s:N-EFF-MU-MC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    //    OUT << "EFF_TRIG_BDMM\t" << i << "\t" << fNumbersBd[i]->effTrigTNP << endl;
    err1 = fNumbersBd[i]->effTrigTNPE; 
    err2 = sysTr*fNumbersBd[i]->effTrigTNP;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effTrigTNP, Form("%s:N-EFF-TRIG-PID-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-PID-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-PID-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-TRIG-PID-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effTrigTNP, totE, Form("%s:N-EFF-TRIG-PID-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    err1 = fNumbersBd[i]->effTrigTNPMCE; 
    err2 = sysTr*fNumbersBd[i]->effTrigTNPMC;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effTrigTNPMC, totE, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_TRIG_BDMM\t" << i << "\t" << fNumbersBd[i]->effTrigMC << endl;
    err1 = fNumbersBd[i]->effTrigMCE; 
    err2 = sysTr*fNumbersBd[i]->effTrigMC;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-TRIG-MC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-TRIG-MC-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-TRIG-MC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effTrigMC, totE, Form("%s:N-EFF-TRIG-MC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_CAND_BDMM\t" << i << "\t" << fNumbersBd[i]->effCand
	<< "\t" << sysCand*fNumbersBd[i]->effCand 
	<< endl;
    err1 = fNumbersBd[i]->effCandE; 
    err2 = sysCand*fNumbersBd[i]->effCand;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effCand, Form("%s:N-EFF-CAND-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-CAND-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-CAND-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-CAND-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(totE, fNumbersBd[i]->effCand, Form("%s:N-EFF-CAND-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_ANA_BDMM\t" << i << "\t" << fNumbersBd[i]->effAna 
	<< "\t" << sysAna*fNumbersBd[i]->effAna
	<< endl;
    err1 = fNumbersBd[i]->effAnaE; 
    err2 = sysAna*fNumbersBd[i]->effAna;
    totE = TMath::Sqrt(err1*err1 + err2*err2); 
    fTEX << formatTex(fNumbersBd[i]->effAna, Form("%s:N-EFF-ANA-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err1, Form("%s:N-EFF-ANA-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(err2, Form("%s:N-EFF-ANA-BDMM%d:sys", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(totE, Form("%s:N-EFF-ANA-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effAna, totE, Form("%s:N-EFF-ANA-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "# Observed in signal boxes" << endl;
    OUT << "OBS_BSMM\t" << i << "\t" << fNumbersBs[i]->bsObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bsObs, Form("%s:N-OBS-BSMM%d:val", fSuffix.c_str(), i), 0) << endl;

    OUT << "OBS_BDMM\t" << i << "\t" << fNumbersBs[i]->bdObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bdObs, Form("%s:N-OBS-BDMM%d:val", fSuffix.c_str(), i), 0) << endl;

    OUT << "PEAK_BKG_OFF\t" << i << "\t" << fNumbersBs[i]->offRare 
	<< "\t" << fNumbersBs[i]->offRareE
	<< endl;
    fTEX << formatTex(fNumbersBs[i]->offRare, Form("%s:N-OFF-RARE%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->offRareE,Form("%s:N-OFF-RARE%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "PEAK_BKG_BS\t" << i << "\t" << fNumbersBs[i]->bsRare
	<< "\t" << fNumbersBs[i]->bsRareE
	<< endl;
    fTEX << formatTex(fNumbersBs[i]->bsRare, Form("%s:N-PEAK-BKG-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bsRareE,Form("%s:N-PEAK-BKG-BS%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "PEAK_BKG_BD\t" << i << "\t"	<< fNumbersBs[i]->bdRare
	<< "\t" << fNumbersBs[i]->bdRareE
	<< endl;
    fTEX << formatTex(fNumbersBs[i]->bdRare, Form("%s:N-PEAK-BKG-BD%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bdRareE,Form("%s:N-PEAK-BKG-BD%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "TAU_BS\t" << i << "\t" << fNumbersBs[i]->tauBs 
	<< "\t" << fNumbersBs[i]->tauBsE
	<< endl;
    fTEX << formatTex(fNumbersBs[i]->tauBs, Form("%s:N-TAU-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->tauBsE, Form("%s:N-TAU-BS%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "TAU_BD\t" << i << "\t" << fNumbersBs[i]->tauBd << "\t" << fNumbersBs[i]->tauBdE << endl;
    fTEX << formatTex(fNumbersBs[i]->tauBd, Form("%s:N-TAU-BD%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->tauBdE, Form("%s:N-TAU-BD%d:err", fSuffix.c_str(), i), 2) << endl;

    
    double bsExpObs = fNumbersBs[i]->tauBs*fNumbersBs[i]->bgObs
      + fNumbersBs[i]->bsRare
      + fNumbersBs[i]->bsNoScaled;

    double bsExpObsE = TMath::Sqrt(fNumbersBs[i]->tauBs*fNumbersBs[i]->tauBs*fNumbersBs[i]->bgObs
				   + fNumbersBs[i]->bsRareE*fNumbersBs[i]->bsRareE
				   + fNumbersBs[i]->bsNoScaledE*fNumbersBs[i]->bsNoScaledE);

    OUT << "#EXP_OBS_BSMM\t" << i << "\t" << bsExpObs << "\t" << bsExpObsE << endl;
    fTEX << formatTex(bsExpObs, Form("%s:N-EXP-OBS-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(bsExpObsE, Form("%s:N-EXP-OBS-BS%d:err", fSuffix.c_str(), i), 2) << endl;

    double bdExpObs = fNumbersBs[i]->tauBd*fNumbersBs[i]->bgObs
      + fNumbersBs[i]->bdRare
      + fNumbersBd[i]->bdNoScaled;

    double bdExpObsE = TMath::Sqrt(fNumbersBs[i]->tauBd*fNumbersBs[i]->tauBd*fNumbersBs[i]->bgObs
				   + fNumbersBs[i]->bdRareE*fNumbersBs[i]->bdRareE
				   + fNumbersBd[i]->bdNoScaledE*fNumbersBd[i]->bdNoScaledE);

    OUT << "#EXP_OBS_BDMM\t" << i << "\t" << bdExpObs << "\t" << bdExpObsE << endl;
    fTEX << formatTex(bdExpObs, Form("%s:N-EXP-OBS-BD%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(bdExpObsE, Form("%s:N-EXP-OBS-BD%d:err", fSuffix.c_str(), i), 2) << endl;


  }

  OUT.close();

  fSuffix = cache; 
}


// ----------------------------------------------------------------------
void plotResults::printCsBFNumbers() {

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


// ----------------------------------------------------------------------
// returns the number of events expected 
// (corresponding to the efftot given, careful about mass cuts!)
double  plotResults::scaledYield(numbers *a, numbers *no, double chanbf, double fsfu) {

  bool verbose(true); 

    //          BF(Bs -> mu mu)   fs epstot(Bs) 
    //   n_s = -----------------  -- ---------  N(B+) 
    //          BF(B+ -> mu muK)  fu epstot(B+) 

  if (verbose) {
    cout << "** scale yields from " << a->name << " to " << no->name << endl;
    cout << "   pss:  " << a->pss << endl;
    cout << "   fsfu: " << fsfu << endl;
    cout << "   bf:   " << chanbf << endl;
    cout << "   eps:  " << a->effTot << "/" << no->effTot << endl;
    cout << "   N(B+):" << no->fitYield << endl;
  }

  double relError(0.10); 
  if (fsfu < 0.9) relError = 0.15; 
  double yield = (chanbf/6.0e-5) * (fsfu) * (a->effTot/no->effTot) * no->fitYield; 

  if (string::npos != a->name.find("signal Bs2MuMu")) {
    cout << "SCALING SIGNAL BS2MUMU" << endl;
    a->bsNoScaled  = a->pss * yield;
    a->bsNoScaledE = relError*a->bsNoScaled; 
  } else if (string::npos != a->name.find("signal Bd2MuMu")) {
    cout << "SCALING SIGNAL BD2MUMU" << endl;
    a->bdNoScaled  = a->pdd * yield;
    a->bdNoScaledE = relError*a->bdNoScaled; 
  } else if (string::npos != a->name.find("Bla")) {
    cout << "SCALING BACKGROUND" << endl;
    a->bsRare  = a->pss * yield;
    a->bsRareE = relError*a->bsRare; 
    a->bdRare  = a->pdd * yield;
    a->bdRareE = relError*a->bdRare; 
  }

  if (verbose) {
    cout << "      fs:      " << (chanbf/6.0e-5) * (fsfu) * a->pss * (a->effTot/no->effTot) << endl;
    cout << "   Ns(X):      " << a->bsRare << endl;
    cout << "   Nd(X):      " << a->bdRare << endl;
    cout << "******* Yield: " << yield << endl;
  }

  return yield; 
}


// ----------------------------------------------------------------------
void plotResults::rareBg() {

  string cache = fSuffix; 
  if (fDoUseBDT) fSuffix = "bdt" + fSuffix; 

  c0->Clear();
  gStyle->SetOptStat(0);

  TH1D *eRare = (TH1D*)(fhMassWithAllCuts[0]->Clone("eRare"));  eRare->SetLineColor(kBlack); 
  TH1D *bRare = (TH1D*)(fhMassWithAllCuts[0]->Clone("bRare"));  bRare->SetLineColor(kBlack); 

  THStack *hRareBg0 = new THStack("hRareBg0","");
  THStack *hRareBg1 = new THStack("hRareBg1","");

  std::map<string, int> colors, hatches;
  std::map<string, double> mscale;  
  std::map<string, double> err;  
  std::map<string, double> chanbf;  
  double epsPi(0.0015), errPi2(0.15*0.15); // relative errors on misid rates are statistical error from Danek's fits
  double epsKa(0.0017), errKa2(0.15*0.15); 
  double epsPr(0.0005), errPr2(0.15*0.15);
  colors.insert(make_pair("bgLb2KP", 46)); hatches.insert(make_pair("bgLb2KP", 3004)); mscale.insert(make_pair("bgLb2KP", epsPi*epsPr)); 
  chanbf.insert(make_pair("bgLb2KP", 5.6e-6)); 
  err.insert(make_pair("bgLb2KP", TMath::Sqrt(0.3*0.3 + errPi2 + errPr2))); 

  colors.insert(make_pair("bgLb2PiP", 49)); hatches.insert(make_pair("bgLb2PiP", 3005)); mscale.insert(make_pair("bgLb2PiP", epsKa*epsPr)); 
  chanbf.insert(make_pair("bgLb2PiP", 3.5e-6)); 
  err.insert(make_pair("bgLb2PiP", TMath::Sqrt(0.31*0.31 + errKa2 + errPr2))); 

  colors.insert(make_pair("bgLb2PMuNu", 48)); hatches.insert(make_pair("bgLb2PMuNu", 3006)); mscale.insert(make_pair("bgLb2PMuNu", epsPr)); 
  chanbf.insert(make_pair("bgLb2PMuNu", 1.3e-4)); 
  err.insert(make_pair("bgLb2PMuNu", TMath::Sqrt(0.31*0.31 + errPr2))); 


  colors.insert(make_pair("bgBs2KK", 30)); hatches.insert(make_pair("bgBs2KK", 3004)); mscale.insert(make_pair("bgBs2KK", epsKa*epsKa)); 
  chanbf.insert(make_pair("bgBs2KK", 2.7e-5)); 
  err.insert(make_pair("bgBs2KK", TMath::Sqrt(0.15*0.15 + errKa2 + errKa2))); 

  colors.insert(make_pair("bgBs2KPi", 32)); hatches.insert(make_pair("bgBs2KPi", 3005)); mscale.insert(make_pair("bgBs2KPi", epsPi*epsKa)); 
  chanbf.insert(make_pair("bgBs2KPi", 5.0e-6)); 
  err.insert(make_pair("bgBs2KPi", TMath::Sqrt(0.22*0.22 + errPi2 + errKa2))); 

  colors.insert(make_pair("bgBs2PiPi", 33)); hatches.insert(make_pair("bgBs2PiPi", 3007)); mscale.insert(make_pair("bgBs2PiPi", epsPi*epsPi)); 
  chanbf.insert(make_pair("bgBs2PiPi", 1.2e-6)); 
  err.insert(make_pair("bgBs2PiPi", TMath::Sqrt(1.0*1.0 + errPi2 + errPi2))); 

  colors.insert(make_pair("bgBs2KMuNu", 34)); hatches.insert(make_pair("bgBs2KMuNu", 3008)); mscale.insert(make_pair("bgBs2KMuNu", epsKa)); 
  chanbf.insert(make_pair("bgBs2KMuNu", 1.3e-4)); 
  err.insert(make_pair("bgBs2KMuNu", TMath::Sqrt(0.2*0.2 + errKa2))); 


  colors.insert(make_pair("bgBd2KK", 40)); hatches.insert(make_pair("bgBd2KK", 3004)); mscale.insert(make_pair("bgBd2KK", epsKa*epsKa)); 
  chanbf.insert(make_pair("bgBd2KK", 1.5e-7)); 
  err.insert(make_pair("bgBd2KK", TMath::Sqrt(0.73*0.73 + errKa2 + errKa2))); 

  colors.insert(make_pair("bgBd2KPi", 41)); hatches.insert(make_pair("bgBd2KPi", 3005)); mscale.insert(make_pair("bgBd2KPi", epsKa*epsPi)); 
  chanbf.insert(make_pair("bgBd2KPi", 1.9e-5)); 
  err.insert(make_pair("bgBd2KPi", TMath::Sqrt(0.05*0.05 + errKa2 + errPi2))); 

  colors.insert(make_pair("bgBd2PiPi", 42)); hatches.insert(make_pair("bgBd2PiPi", 3007)); mscale.insert(make_pair("bgBd2PiPi", epsPi*epsPi)); 
  chanbf.insert(make_pair("bgBd2PiPi", 5.2e-6)); 
  err.insert(make_pair("bgBd2PiPi", TMath::Sqrt(0.04*0.04 + errPi2 + errPi2))); 

  colors.insert(make_pair("bgBd2PiMuNu", 43)); hatches.insert(make_pair("bgBd2PiMuNu", 3008)); mscale.insert(make_pair("bgBd2PiMuNu", epsPi)); 
  chanbf.insert(make_pair("bgBd2PiMuNu", 1.3e-4)); 
  err.insert(make_pair("bgBd2PiMuNu", TMath::Sqrt(0.2*0.2 + errPi2))); 

  newLegend(0.55, 0.3, 0.80, 0.85); 

  double error(0.), errorInc(0.); 
  double teff[] = {0.85, 0.68};
  double rareBs[]  = {0., 0.};
  double rareBsE[] = {0., 0.};
  double rareBd[]  = {0., 0.};
  double rareBdE[] = {0., 0.};

  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << "% --- rare Background per-channel numbers" << endl;  

  c0->Divide(1,2);

  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    if (string::npos == imap->first.find("bg")) {
      continue;
    }

    if (0 == fF[imap->first]) {
      continue; 
    }

    //     if (string::npos == imap->first.find("bgBd2PiPi")) {
    //       continue;
    //     }

    double misid = mscale[imap->first];
    double ngenfile = ((TH1D*)fF[imap->first]->Get("monEvents"))->GetBinContent(1); 
 
    fF[imap->first]->cd("candAnaMuMu");

    TH1D *hRare[2]; 
    double tot, bd, bs, efftot, pss, pdd;
    
    loopTree(99); 

    for (int ichan = 0; ichan < 2; ++ichan) {
      
      hRare[ichan] = (TH1D*)(fhMassWithAllCuts[ichan]->Clone(Form("h1Rare%d", ichan)));  
      
      tot =  fhMassWithAllCutsManyBins[ichan]->GetSumOfWeights(); 
      bd = fhMassWithAllCutsManyBins[ichan]->Integral(fhMassWithAllCutsManyBins[ichan]->FindBin(fCuts[ichan]->mBdLo), 
						      fhMassWithAllCutsManyBins[ichan]->FindBin(fCuts[ichan]->mBdHi));
      bs = fhMassWithAllCutsManyBins[ichan]->Integral(fhMassWithAllCutsManyBins[ichan]->FindBin(fCuts[ichan]->mBsLo), 
						      fhMassWithAllCutsManyBins[ichan]->FindBin(fCuts[ichan]->mBsHi));

      efftot = static_cast<double>(hRare[ichan]->GetSumOfWeights())/static_cast<double>(ngenfile)*fFilterEff[imap->first];
      pss = bs/tot;
      pdd = bd/tot;
      
      fNumbersBla[ichan]->effTot = efftot;
      fNumbersBla[ichan]->pss = pss;
      fNumbersBla[ichan]->pdd = pdd;
      double pRatio(0.);
      if (string::npos != imap->first.find("Bs")) pRatio = fsfu;
      if (string::npos != imap->first.find("Lb")) pRatio = fsfu;
      if (string::npos != imap->first.find("Bd")) pRatio = 1.;
      if (string::npos != imap->first.find("Bu")) pRatio = 1.;
      
      double yield = scaledYield(fNumbersBla[ichan], fNumbersNo[ichan], chanbf[imap->first], pRatio);

      cout << "====> efftot: " << hRare[ichan]->GetSumOfWeights() << "/" << ngenfile << " = " << efftot << endl;
      cout << "   misid*trig: " << misid*teff[ichan] << endl;
      // -- NB: rareBxE contains the SQUARE of the error
      rareBs[ichan]  += fNumbersBla[ichan]->bsRare*misid*teff[ichan];
      error =  misid*teff[ichan]*err[imap->first]*fNumbersBla[ichan]->bsRare;
      errorInc =  misid*teff[ichan]*misid*teff[ichan]*fNumbersBla[ichan]->bsRareE*fNumbersBla[ichan]->bsRareE + error*error;
      cout << " ->Bs increment: " << fNumbersBla[ichan]->bsRare*misid*teff[ichan] 
	   << " error: " << error << " errorInc: " << errorInc << " err[]: " << err[imap->first]
	   << endl;
      rareBsE[ichan] += errorInc; 
      errorInc = TMath::Sqrt(errorInc); 
      fTEX <<  Form("\\vdef{%s:%s:bsRare%d}   {\\ensuremath{{%7.6f } } }", fSuffix.c_str(), imap->first.c_str(), ichan, rareBs[ichan]) << endl;
      fTEX <<  Form("\\vdef{%s:%s:bsRare%dE}  {\\ensuremath{{%7.6f } } }", fSuffix.c_str(), imap->first.c_str(), ichan, errorInc) << endl;

      rareBd[ichan]  += fNumbersBla[ichan]->bdRare*misid*teff[ichan];
      error =  misid*teff[ichan]*err[imap->first]*fNumbersBla[ichan]->bdRare;
      errorInc = misid*teff[ichan]*misid*teff[ichan]*fNumbersBla[ichan]->bdRareE*fNumbersBla[ichan]->bdRareE + error*error;
      cout << " ->Bd increment: " << fNumbersBla[ichan]->bdRare*misid*teff[ichan]
	   << " error: " << error << " errorInc: " << errorInc << " err[]: " << err[imap->first]
	   << endl;
      rareBdE[ichan] += errorInc; 
      errorInc = TMath::Sqrt(errorInc); 
      fTEX <<  Form("\\vdef{%s:%s:bdRare%d}   {\\ensuremath{{%7.6f } } }", fSuffix.c_str(), imap->first.c_str(), ichan, rareBd[ichan]) << endl;
      fTEX <<  Form("\\vdef{%s:%s:bdRare%dE}  {\\ensuremath{{%7.6f } } }", fSuffix.c_str(), imap->first.c_str(), ichan, errorInc) << endl;

      hRare[ichan]->SetFillColor(colors[imap->first]);
      hRare[ichan]->SetFillStyle(1000);
      hRare[ichan]->Scale(yield*misid*teff[ichan]/tot);
      cout << "histogram bins: " << hRare[ichan]->FindBin(fCuts[ichan]->mBsLo) << " " << hRare[ichan]->FindBin(fCuts[ichan]->mBsHi) << endl;
      cout << "histogram hRareX: " << hRare[ichan]->Integral(hRare[ichan]->FindBin(fCuts[ichan]->mBsLo), 
							     hRare[ichan]->FindBin(fCuts[ichan]->mBsHi))
	   << " total: " << hRare[ichan]->GetSumOfWeights()
	   << endl;
      
    }

    bRare->Add(hRare[0]); 
    eRare->Add(hRare[1]); 

    legg->AddEntry(hRare[0], fName[imap->first].c_str(), "f"); 
    c0->cd(1); 
    hRare[0]->Draw();
    tl->DrawLatex(0.5, 0.92, imap->first.c_str());
    c0->cd(2); 
    hRare[1]->Draw();
    c0->Modified();
    c0->Update();

    hRareBg0->Add(hRare[0]); 
    hRareBg1->Add(hRare[1]); 
    
  }

  c0->Clear();

  rareBsE[0] = TMath::Sqrt(rareBsE[0]);
  rareBsE[1] = TMath::Sqrt(rareBsE[1]);
  rareBdE[0] = TMath::Sqrt(rareBdE[0]);
  rareBdE[1] = TMath::Sqrt(rareBdE[1]);
							 
  gStyle->SetOptStat(0);

  shrinkPad(0.12, 0.18); 

  hRareBg0->SetMaximum(1.0); 
  hRareBg0->Draw();
  TH1D *hhRareBg0 = (TH1D*)hRareBg0->GetHistogram(); 
  hhRareBg0->SetAxisRange(4.9, 5.9, "X"); 
  setTitles(hhRareBg0, "m_{#mu #mu} [GeV]", Form("Candidates/%4.3f GeV", hhRareBg0->GetBinWidth(1)), 0.06, 0.9, 1.5);
  legg->SetHeader("CMS simulation"); 
  legg->Draw(); 
  hhRareBg0->Draw("same");
  string pdfname = Form("%s/%s_rare0.pdf", fDirectory.c_str(), fSuffix.c_str());
  //  stamp(0.2, "CMS, 4.9 fb^{-1}", 0.65, "#sqrt{s} = 7 TeV"); 
  if (fDoPrint) {
    if (c0) c0->SaveAs(pdfname.c_str());
  }

  hRareBg1->SetMaximum(1.0); 
  hRareBg1->Draw();
  TH1D *hhRareBg1 = (TH1D*)hRareBg1->GetHistogram(); 
  hhRareBg1->SetAxisRange(4.9, 5.9, "X"); 
  setTitles(hhRareBg1, "m_{#mu #mu} [GeV]", Form("Candidates/%4.3f GeV", hhRareBg0->GetBinWidth(1)), 0.06, 0.9, 1.5);
  legg->SetHeader("CMS simulation"); 
  legg->Draw(); 
  hhRareBg1->Draw("same");
  pdfname = Form("%s/%s_rare1.pdf", fDirectory.c_str(), fSuffix.c_str());
  //  stamp(0.2, "CMS, 4.9 fb^{-1}", 0.65, "#sqrt{s} = 7 TeV"); 
  if (fDoPrint){
    if (c0) c0->SaveAs(pdfname.c_str());
  }

  fNumbersBs[0]->bsRare = rareBs[0]; 
  fNumbersBs[0]->bsRareE= rareBsE[0]; 
  fNumbersBs[0]->bdRare = rareBd[0]; 
  fNumbersBs[0]->bdRareE= rareBdE[0]; 

  fNumbersBs[0]->offRare = 0.; 
  fNumbersBs[0]->offRareE= 0.; 

  fNumbersBs[1]->bsRare = rareBs[1]; 
  fNumbersBs[1]->bsRareE= rareBsE[1]; 
  fNumbersBs[1]->bdRare = rareBd[1]; 
  fNumbersBs[1]->bdRareE= rareBdE[1]; 

  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << "% --- rare Background summary numbers" << endl;  
  fTEX <<  Form("\\vdef{%s:bsRare0}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), rareBs[0]) << endl;
  fTEX <<  Form("\\vdef{%s:bsRare0E}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), rareBsE[0]) << endl;
  fTEX <<  Form("\\vdef{%s:bsRare1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), rareBs[1]) << endl;
  fTEX <<  Form("\\vdef{%s:bsRare1E}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), rareBsE[1]) << endl;
  fTEX <<  Form("\\vdef{%s:bdRare0}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), rareBd[0]) << endl;
  fTEX <<  Form("\\vdef{%s:bdRare0E}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), rareBdE[0]) << endl;
  fTEX <<  Form("\\vdef{%s:bdRare1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), rareBd[1]) << endl;
  fTEX <<  Form("\\vdef{%s:bdRare1E}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), rareBdE[1]) << endl;

  TDirectory *pD = gFile; 
  fHistFile->cd();
  bRare->SetDirectory(gDirectory);
  bRare->Write(); 
  eRare->SetDirectory(gDirectory);
  eRare->Write(); 
  pD->cd();

  fSuffix = cache; 
}


// ----------------------------------------------------------------------
double plotResults::barlow(int nobs, double bg, double bgE, double sE) {
  double ul(-99.); 

  if (nobs > 10)  {
    cout << "barlow: dunno" << endl;
    return -99;
  }

  // -- "unified confidence limits" from PDG as intervals to probe 
  double uclLo[] = {0.00, 0.11, 0.53, 1.10, 1.47, 1.84, 2.21, 3.56, 3.96, 4.36, 5.50};
  double uclHi[] = {2.44, 4.36, 5.91, 7.42, 8.60, 9.99,11.47,12.53,13.99,15.30,16.50};

  double cl(0.1); // this is 90% UL!

  // -- Set up Gaussian pdf for expected background 
  double bgSig= (bgE > 0. ? bgE : bg/3.); 

  // -- Set up Gaussian pdf for sensitivity
  double sSig= (sE > 0. ? sE : 0.1); 

  // -- Set up histogram for counting
  TH1D *h = (TH1D*)gROOT->FindObject("hbarlow"); 
  if (h) {
    h->Reset();
  } else {
    //    h = new TH1D("hbarlow", "", 10000, 0., 100.); 
    h = new TH1D("hbarlow", "", 50, 0., 50.); 
  }

  int bin = h->FindBin(nobs); 

  // -- Run toy
  double tSignal(0.), tBackground(0.), tSensitivity(1.), tCounts(0.);
  double tFraction(0.); 
  double minDeviation(99.), bestAlpha(99.); 

  int NSTEPS(100); 
  double step = (uclHi[nobs] - uclLo[nobs])/NSTEPS; 

  if (nobs <= 2) {
    step   = 0.02; 
  } else {
    step   = 0.1; 
  }
  NSTEPS = int((uclHi[nobs] - uclLo[nobs])/step + 1);

  // -- This loop is quite rough ...?
  for (int k = 0; k < NSTEPS; ++k) {
    h->Reset();
    tSignal = uclLo[nobs] + k*step; 
    
    for (int i = 0; i < 100000; ++i) {
      tBackground  = gRandom->Gaus(bg, bgSig); 
      tSensitivity = gRandom->Gaus(1., sSig); 
    
      tCounts = gRandom->Poisson((tSignal+tBackground)*tSensitivity);
      h->Fill(tCounts); 
    }
    tFraction= h->Integral(1, bin) / h->Integral();
    if (TMath::Abs(tFraction - cl) < minDeviation) {
      minDeviation = TMath::Abs(tFraction - cl); 
      bestAlpha    = tFraction; 
      ul           = tSignal;
    }

    if (tFraction < (cl - minDeviation)) {
      //      cout << "break " << tFraction << " bestAlpha: " << bestAlpha << " (should be " << cl << ")" << endl;
      break;
    }
  }

  return ul; 
}
