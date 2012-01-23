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

//   fNumbersNo[0]->effTot = 0.00101;
//   fNumbersNo[1]->effTot = 0.00030;
//   fNumbersNo[0]->fitYield = 69177;
//   fNumbersNo[1]->fitYield = 19991;
  
//   rareBg(); 
//   return;


  if (channels & 16) {
    plotEfficiencies a3(fFiles.c_str()); 
    a3.makeAll(); 
  }

  if (channels & 8) {
    plotPU a2(fFiles.c_str()); 
    a2.makeAll(); 
  }

  if (channels & 4) {
    plotOverlays a1(fFiles.c_str()); 
    a1.makeAll(); 
  }

  if (channels & 1) {
    fNormProcessed = false; 
    fDoUseBDT = false; 
    fDoApplyCowboyVeto = false;   
    fDoApplyCowboyVetoAlsoInSignal = false;   
    computeNormUL();
    computeCsBF();
    acceptancePerProcess();
  }

  if (channels & 2) {
    fNormProcessed = false; 
    fDoUseBDT = true; 
    fDoApplyCowboyVeto = false;   
    fDoApplyCowboyVetoAlsoInSignal = false;   
    computeNormUL();
    computeCsBF();
    //    acceptancePerProcess();
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
  scaledHist(0);
  loopTree(1);  // Bd2MuMu eff
  c0->Modified(); c0->Update();
  scaledHist(1);
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
  createAllCfgFiles(bla); 

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

//   double fUL = (fNumbersBs[0]->bgBsExp/fNumbersNo[0]->fitYield) //????FIXME
//     *(fu/fs)
//     *(fNumbersNo[0]->acc/fNumbersBs[0]->acc)
//     *(fNumbersNo[0]->effCand/fNumbersBs[0]->effCand)     
//     *(fNumbersNo[0]->effMuidMC/fNumbersBs[0]->effMuidMC)
//     *(fNumbersNo[0]->effTrigMC/fNumbersBs[0]->effTrigMC)
//     *(fNumbersNo[0]->effAna/fNumbersBs[0]->effAna)
//     * fBF;

//   cout << "prod(eff) expected UL: " << fUL << endl;

//   fUL = (fNul/fNumbersNo[0]->fitYield)
//     *(fu/fs)
//     *(fNumbersNo[0]->effTot/fNumbersBs[0]->effTot)
//     * fBF;

//   cout << "effTot expected UL:    " << fUL << endl;

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
void plotResults::scaledHist(int mode) {

  double yield(0.), tot(0.);

  string cache("cnc"); 
  if (fDoUseBDT) cache = "bdt";  

  TH1D *h(0);
  for (int i = 0; i < 2; ++i) {
    if (0 == mode) {
      tot  = fhMassWithAllCuts[i]->Integral(0, fhMassWithAllCuts[i]->GetNbinsX()+1); 
      yield = scaledYield(fNumbersBs[i], fNumbersNo[i], 3.2e-9, fsfu);
      h = (TH1D*)(fhMassWithAllCuts[i]->Clone(Form("Bs_%s_chan%d", cache.c_str(), i)));  
      h->Scale(yield/tot);
    } 
    else if (1 == mode) {
      tot  = fhMassWithAllCuts[i]->Integral(0, fhMassWithAllCuts[i]->GetNbinsX()+1); 
      yield = scaledYield(fNumbersBd[i], fNumbersNo[i], 1.0e-10, 1.);
      h = (TH1D*)(fhMassWithAllCuts[i]->Clone(Form("Bd_%s_chan%d", cache.c_str(), i)));  
      h->Scale(yield/tot);
    }
    h->SetDirectory(fHistFile);
    h->Write();
  }
}


// ----------------------------------------------------------------------
void plotResults::createAllCfgFiles(string fname) { 

  vector<string> lines; 
  char  buffer[200];
  ifstream is(fname.c_str());
  while (is.getline(buffer, 200, '\n')) {
    lines.push_back(string(buffer));
  }
  is.close();

  float bsExp0(0.), bsExp1(0.), bdExp0(0.), bdExp1(0.), err(0.); 
  for (int i = 0; i < lines.size(); ++i) {
    cout << lines[i] << endl;   
    if (string::npos != lines[i].find("#EXP_OBS_BSMM\t0")) sscanf(lines[i].c_str(), "#EXP_OBS_BSMM\t0\t%f\t%f", &bsExp0, &err);
    if (string::npos != lines[i].find("#EXP_OBS_BDMM\t0")) sscanf(lines[i].c_str(), "#EXP_OBS_BDMM\t0\t%f\t%f", &bdExp0, &err);
    if (string::npos != lines[i].find("#EXP_OBS_BSMM\t1")) sscanf(lines[i].c_str(), "#EXP_OBS_BSMM\t1\t%f\t%f", &bsExp1, &err);
    if (string::npos != lines[i].find("#EXP_OBS_BDMM\t1")) sscanf(lines[i].c_str(), "#EXP_OBS_BDMM\t1\t%f\t%f", &bdExp1, &err);
  }
  cout << "bsExp0: " << bsExp0 << endl;
  cout << "bdExp0: " << bdExp0 << endl;
  cout << "bsExp1: " << bsExp1 << endl;
  cout << "bdExp1: " << bdExp1 << endl;


//   double nexp = cbg0 + _sig; 
//   int    nobs = static_cast<int>(nexp+0.5); 
//   double w8(0.), w8cum(0.); 
//   double nulbayes(0.), nulw8(0.); 
// //   cout << "Lo: " << _nhlo << " Hi: " << _nhhi << " -> comb. BG = " << cbg0 << " sig: " << _sig 
// //        << " -> nexp: " << nexp << " -> nobs<exp> = " << nobs
// //        << endl;
//   vector<int> bd0, bs0, bd1, bs1; 
//   for (int ibd0 = 0; ibd0 < 5*bdExp0; ++ibd0) {
//     w8bd0 = TMath::PoissonI(ibd0, bdExp0); 
//     w8bd0Cum += w8bd0; 
//     if (ibd0 < bdExp0 && w8bd0 < 0.01) continue;
//     if (ibd0 > bdExp0 && w8db0cum > 0.99) break;
//   }



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

    OUT << "# Expected in signal boxes" << endl;
    double bsExpObs = fNumbersBs[i]->bgBsExp
      + fNumbersBs[i]->bsRare
      + fNumbersBs[i]->bsNoScaled;

    double bsExpObsE = TMath::Sqrt(fNumbersBs[i]->bgBsExpE
				   + fNumbersBs[i]->bsRareE*fNumbersBs[i]->bsRareE
				   + fNumbersBs[i]->bsNoScaledE*fNumbersBs[i]->bsNoScaledE);

    OUT << "#EXP_OBS_BSMM\t" << i << "\t" << bsExpObs << "\t" << bsExpObsE << endl;
    fTEX << formatTex(bsExpObs, Form("%s:N-EXP-OBS-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(bsExpObsE, Form("%s:N-EXP-OBS-BS%d:err", fSuffix.c_str(), i), 2) << endl;

    double bdExpObs = fNumbersBs[i]->bgBdExp
      + fNumbersBs[i]->bdRare
      + fNumbersBd[i]->bdNoScaled;

    double bdExpObsE = TMath::Sqrt(fNumbersBs[i]->bgBdExpE
				   + fNumbersBs[i]->bdRareE*fNumbersBs[i]->bdRareE
				   + fNumbersBd[i]->bdNoScaledE*fNumbersBd[i]->bdNoScaledE);

    OUT << "#EXP_OBS_BDMM\t" << i << "\t" << bdExpObs << "\t" << bdExpObsE << endl;
    fTEX << formatTex(bdExpObs, Form("%s:N-EXP-OBS-BD%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(bdExpObsE, Form("%s:N-EXP-OBS-BD%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "# Observed in signal boxes" << endl;
    OUT << "OBS_BSMM\t" << i << "\t" << fNumbersBs[i]->bsObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bsObs, Form("%s:N-OBS-BSMM%d:val", fSuffix.c_str(), i), 0) << endl;

    OUT << "OBS_BDMM\t" << i << "\t" << fNumbersBs[i]->bdObs << endl;
    fTEX << formatTex(fNumbersBs[i]->bdObs, Form("%s:N-OBS-BDMM%d:val", fSuffix.c_str(), i), 0) << endl;

    OUT << "PEAK_BKG_OFF\t" << i << "\t" << fNumbersBs[i]->offLoRare+fNumbersBs[i]->offHiRare 
	<< "\t" << TMath::Sqrt(fNumbersBs[i]->offLoRareE*fNumbersBs[i]->offLoRareE + fNumbersBs[i]->offHiRareE*fNumbersBs[i]->offHiRareE)
	<< endl;
    fTEX << formatTex(fNumbersBs[i]->offLoRare, Form("%s:N-OFFLO-RARE%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->offLoRareE,Form("%s:N-OFFLO-RARE%d:err", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->offHiRare, Form("%s:N-OFFHI-RARE%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->offHiRareE,Form("%s:N-OFFHI-RARE%d:err", fSuffix.c_str(), i), 2) << endl;

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
// returns the number of events expected (total integral over entire histogram!)
// (corresponding to the efftot given, careful about mass cuts!)
double  plotResults::scaledYield(numbers *a, numbers *no, double chanbf, double lfsfu) {

  bool verbose(true); 

    //          BF(Bs -> mu mu)   fs epstot(Bs) 
    //   n_s = -----------------  -- ---------  N(B+) 
    //          BF(B+ -> mu muK)  fu epstot(B+) 

  if (verbose) {
    cout << "** scale yields from " << a->name << " to " << no->name << endl;
    cout << "   pss:  " << a->pss << endl;
    cout << "   pdd:  " << a->pdd << endl;
    cout << "   fsfu: " << lfsfu << endl;
    cout << "   bf:   " << chanbf << endl;
    cout << "   eps:  " << a->effTot << "/" << no->effTot << endl;
    cout << "   N(B+):" << no->fitYield << endl;
  }

  double relError(0.10); 
  if (lfsfu < 0.9) relError = 0.15; 
  double yield = (chanbf/6.0e-5) * (lfsfu) * (a->effTot/no->effTot) * no->fitYield; 

  if (string::npos != a->name.find("signal Bs2MuMu")) {
    cout << "SCALING SIGNAL BS2MUMU" << endl;
    a->bsNoScaled  = a->pss * yield;
    a->bsNoScaledE = relError*a->bsNoScaled; 
  } else if (string::npos != a->name.find("signal Bd2MuMu")) {
    cout << "SCALING SIGNAL BD2MUMU" << endl;
    a->bdNoScaled  = a->pdd * yield;
    a->bdNoScaledE = relError*a->bdNoScaled; 
  } else if (string::npos != a->name.find("Bla")) {
    cout << "SCALING RARE BACKGROUND" << endl;
    a->bsRare  = a->pss * yield;
    a->bsRareE = relError*a->bsRare; 
    a->bdRare  = a->pdd * yield;
    a->bdRareE = relError*a->bdRare; 
  }

  if (verbose) {
    cout << "******* Yield: " << yield << endl;
    cout << "      sf:      " << (chanbf/6.0e-5) * (lfsfu) * (a->effTot/no->effTot) << endl;
    cout << "   sf(s):      " << (chanbf/6.0e-5) * (lfsfu) * (a->effTot/no->effTot) * a->pss << endl;
    cout << "   Ns(X):      " << a->bsRare << endl;
    cout << "   sf(d):      " << (chanbf/6.0e-5) * (lfsfu) * (a->effTot/no->effTot) * a->pdd << endl;
    cout << "   Nd(X):      " << a->bdRare << endl;
  }

  return yield; 
}


// ----------------------------------------------------------------------
void plotResults::rareBg() {

  c0->Clear();
  gStyle->SetOptStat(0);

  string cache("cnc"); 
  if (fDoUseBDT) cache = "bdt";  
  
  TH1D *eRare = (TH1D*)(fhMassWithAllCuts[0]->Clone(Form("eRare_%s", cache.c_str())));  eRare->SetLineColor(kBlack); eRare->Reset();
  TH1D *bRare = (TH1D*)(fhMassWithAllCuts[0]->Clone(Form("bRare_%s", cache.c_str())));  bRare->SetLineColor(kBlack); bRare->Reset();

  cache = fSuffix; 
  if (fDoUseBDT) fSuffix = "bdt" + fSuffix; 

  THStack *hRareBg0 = new THStack("hRareBg0","");
  THStack *hRareBg1 = new THStack("hRareBg1","");

  std::map<string, int> colors, hatches;
  std::map<string, double> mscale;  
  std::map<string, double> err;  
  std::map<string, double> chanbf;  

  // -- tight muon values
  double epsMu(0.8); // FIXME: this should be different for barrel and endcap!
  double epsPi(0.0015), errPi2(0.15*0.15); // relative errors on misid rates are statistical error from Danek's fits
  double epsKa(0.0017), errKa2(0.15*0.15); 
  double epsPr(0.0005), errPr2(0.15*0.15);

//   // -- ATM values
//   double epsMu(0.9); // FIXME: this should be different for barrel and endcap!
//   double epsPi(0.009), errPi2(0.15*0.15); // relative errors on misid rates are statistical error from Danek's fits
//   double epsKa(0.01),  errKa2(0.15*0.15); 
//   double epsPr(0.003), errPr2(0.15*0.15);

  colors.insert(make_pair("bgLb2KP", 46)); hatches.insert(make_pair("bgLb2KP", 3004)); mscale.insert(make_pair("bgLb2KP", epsPi*epsPr)); 
  chanbf.insert(make_pair("bgLb2KP", 5.6e-6)); 
  err.insert(make_pair("bgLb2KP", TMath::Sqrt(0.3*0.3 + errPi2 + errPr2))); 

  colors.insert(make_pair("bgLb2PiP", 49)); hatches.insert(make_pair("bgLb2PiP", 3005)); mscale.insert(make_pair("bgLb2PiP", epsKa*epsPr)); 
  chanbf.insert(make_pair("bgLb2PiP", 3.5e-6)); 
  err.insert(make_pair("bgLb2PiP", TMath::Sqrt(0.31*0.31 + errKa2 + errPr2))); 

  colors.insert(make_pair("bgLb2PMuNu", 48)); hatches.insert(make_pair("bgLb2PMuNu", 3006)); mscale.insert(make_pair("bgLb2PMuNu", epsPr*epsMu)); 
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

  colors.insert(make_pair("bgBs2KMuNu", 34)); hatches.insert(make_pair("bgBs2KMuNu", 3008)); mscale.insert(make_pair("bgBs2KMuNu", epsKa*epsMu)); 
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

  colors.insert(make_pair("bgBd2PiMuNu", 43)); hatches.insert(make_pair("bgBd2PiMuNu", 3008)); mscale.insert(make_pair("bgBd2PiMuNu", epsPi*epsMu)); 
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

  //  gStyle->SetOptStat(1111111);
  c0->Divide(1,2);

  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    if (string::npos == imap->first.find("bg")) {
      continue;
    }

//          if (string::npos == imap->first.find("bgBd2PiMuNu")) {
//            continue;
//          }

    if (0 == fF[imap->first]) {
      continue; 
    }

    double misid = mscale[imap->first];
    double ngenfile = ((TH1D*)fF[imap->first]->Get("monEvents"))->GetBinContent(1); 
    cout << "======> " << imap->first << " with misid: " << misid << endl;

    fF[imap->first]->cd("candAnaMuMu");

    TH1D *hRare[2]; 
    double tot0, tot, bd, bs, efftot, pss, pdd;
    
    loopTree(99); 

    for (int ichan = 0; ichan < 2; ++ichan) {
      
      if (fDoUseBDT) 
	hRare[ichan] = (TH1D*)(fhMassWithAllCuts[ichan]->Clone(Form("h1Rare_BDT_%d", ichan)));  
      else 
	hRare[ichan] = (TH1D*)(fhMassWithAllCuts[ichan]->Clone(Form("h1Rare_CNC_%d", ichan)));  

      tot  = fhMassWithAllCutsManyBins[ichan]->Integral(0, fhMassWithAllCutsManyBins[ichan]->GetNbinsX()+1); 
      tot0 = fhMassWithAllCutsManyBins[ichan]->GetSumOfWeights(); 
      bd   = fhMassWithAllCutsManyBins[ichan]->Integral(fhMassWithAllCutsManyBins[ichan]->FindBin(fCuts[ichan]->mBdLo), 
							fhMassWithAllCutsManyBins[ichan]->FindBin(fCuts[ichan]->mBdHi));
      bs   = fhMassWithAllCutsManyBins[ichan]->Integral(fhMassWithAllCutsManyBins[ichan]->FindBin(fCuts[ichan]->mBsLo), 
							fhMassWithAllCutsManyBins[ichan]->FindBin(fCuts[ichan]->mBsHi));

      efftot = tot/static_cast<double>(ngenfile)*fFilterEff[imap->first];
      
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

      cout << "====> efftot: " << tot << "/" << ngenfile << "*" << fFilterEff[imap->first] << " = " << efftot << endl;
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
      cout << "histogram hRareX"
	   << " total (b/w8): " << tot
	   << " total (a/w8): " << hRare[ichan]->Integral(0, hRare[ichan]->GetNbinsX()+1)
	   << " 4.9 - 6.0:    " << hRare[ichan]->Integral(hRare[ichan]->FindBin(4.9001), hRare[ichan]->GetNbinsX()+1)
	   << endl;
      cout << "yield:      " << yield << endl;
      cout << "misid*teff: " << misid*teff[ichan] << endl;
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
							 
  //  gStyle->SetOptStat(0);

  shrinkPad(0.12, 0.18); 

  hRareBg0->SetMaximum(1.0); 
  hRareBg0->Draw();
  TH1D *hhRareBg0 = (TH1D*)hRareBg0->GetHistogram(); 
  hhRareBg0->SetAxisRange(4.9, 5.9, "X"); 
  setTitles(hhRareBg0, "m_{#mu #mu} [GeV]", Form("Candidates/%4.3f GeV", hhRareBg0->GetBinWidth(1)), 0.06, 0.9, 1.5);
  legg->SetHeader("CMS simulation"); 
  legg->Draw(); 
  hhRareBg0->Draw("same");
  string pdfname;
  if (fDoUseBDT) pdfname = Form("%s/%s_bdt_rare0.pdf", fDirectory.c_str(), fSuffix.c_str());
  else  pdfname = Form("%s/%s_cnc_rare0.pdf", fDirectory.c_str(), fSuffix.c_str());
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
  if (fDoUseBDT) pdfname = Form("%s/%s_bdt_rare1.pdf", fDirectory.c_str(), fSuffix.c_str());
  else  pdfname = Form("%s/%s_cnc_rare1.pdf", fDirectory.c_str(), fSuffix.c_str());
  //  stamp(0.2, "CMS, 4.9 fb^{-1}", 0.65, "#sqrt{s} = 7 TeV"); 
  if (fDoPrint){
    if (c0) c0->SaveAs(pdfname.c_str());
  }

  fNumbersBs[0]->bsRare = rareBs[0]; 
  fNumbersBs[0]->bsRareE= rareBsE[0]; 
  fNumbersBs[0]->bdRare = rareBd[0]; 
  fNumbersBs[0]->bdRareE= rareBdE[0]; 
  double relErr = (fNumbersBs[0]->bdRareE/fNumbersBs[0]->bdRare); 

  fNumbersBs[0]->offLoRare = bRare->Integral(bRare->FindBin(4.9001), bRare->FindBin(5.199));
  fNumbersBs[0]->offLoRareE= relErr*fNumbersBs[0]->offLoRare; 
  fNumbersBs[0]->offHiRare = bRare->Integral(bRare->FindBin(5.451), bRare->FindBin(5.899));
  fNumbersBs[0]->offHiRareE= relErr*fNumbersBs[0]->offHiRare; 

  fNumbersBs[1]->bsRare = rareBs[1]; 
  fNumbersBs[1]->bsRareE= rareBsE[1]; 
  fNumbersBs[1]->bdRare = rareBd[1]; 
  fNumbersBs[1]->bdRareE= rareBdE[1]; 
  relErr = (fNumbersBs[1]->bdRareE/fNumbersBs[1]->bdRare); 

  fNumbersBs[1]->offLoRare = eRare->Integral(eRare->FindBin(4.9001), eRare->FindBin(5.199));
  fNumbersBs[1]->offLoRareE= relErr*fNumbersBs[1]->offLoRare; 
  fNumbersBs[1]->offHiRare = eRare->Integral(eRare->FindBin(5.451), eRare->FindBin(5.899));
  fNumbersBs[1]->offHiRareE= relErr*fNumbersBs[1]->offHiRare; 

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


// ----------------------------------------------------------------------
void plotResults::acceptancePerProcess() {

  vector<int> processes;
  processes.push_back(0); 
  processes.push_back(10); 

  double sgV[2][3];
  double sgE[2][3];
  double noV[2][3];
  double noE[2][3];

  for (int i = 0; i < processes.size(); ++i) {
    int mode = processes[i];

    if (0 == mode) fF["SgMcAcc"]->cd();
    else if (10 == mode) fF["NoMcAcc"]->cd();
    else break;

    for (int chan = 0; chan < 2; ++chan) {

      fTEX << "% -- acceptancePerProcess: " << mode << " " << chan << endl;
      
      string cuts; 
      
      if (0 == mode) {
	if (0 == chan) {
	  cuts = "g1pt>1&&g2pt>1&&abs(g1eta)<1.4&&abs(g2eta)<1.4";
	  cuts += "&&m1pt>1&&m2pt>1&&abs(m1eta)<1.4&&abs(m2eta)<1.4&&m1gt&&m2gt";
	} else {
	  cuts = "g1pt>1&&g2pt>1&&(abs(g1eta)>1.4||abs(g2eta)>1.4)&&abs(g1eta)<2.5&&abs(g2eta)<2.5";
	  cuts += "&&m1pt>1&&m2pt>1&&(abs(m1eta)>1.4||abs(m2eta)>1.4)&&abs(m1eta)<2.4&&abs(m2eta)<2.4&&m1gt&&m2gt";
	}
      } else if (10 == mode) {
	if (0 == chan) {
	  cuts = "g1pt>1&&g2pt>1&&g3pt>0.4&&abs(g1eta)<1.4&&abs(g2eta)<1.4&&abs(g3eta)<2.5";
	  cuts += "&&m1pt>1&&m2pt>1&&k1pt>0.5&&abs(m1eta)<1.4&&abs(m2eta)<1.4&&abs(k1eta)<2.4&&m1gt&&m2gt&&k1gt";
	} else {
	  cuts = "g1pt>1&&g2pt>1&&g3pt>0.4&&(abs(g1eta)>1.4||abs(g2eta)>1.4)&&abs(g1eta)<2.5&&abs(g2eta)<2.5&&abs(g3eta)<2.5";
	  cuts += "&&m1pt>1&&m2pt>1&&k1pt>0.5&&(abs(m1eta)>1.4||abs(m2eta)>1.4)&&abs(m1eta)<2.4&&abs(m2eta)<2.4&&abs(k1eta)<2.4&&m1gt&&m2gt&&k1gt";
	}
      }
      
      TTree *t; 
      if (0 == mode) t = (TTree*)gFile->Get("candAnaMuMu/effTree"); 
      else if (10 == mode) t = (TTree*)gFile->Get("candAnaBu2JpsiK/effTree"); 
      if (0 == t) {
	cout << "no tree effTree found " << (mode==0?"candAnaMuMu/effTree":"candAnaBu2JpsiK/effTree") << endl;
	return;
      }
      double c40 = t->Draw("g1pt", Form("%s&&procid==40", cuts.c_str()));
      double n40 = t->Draw("g1pt", "procid==40");
      
      double c41 = t->Draw("g1pt", Form("%s&&procid==41", cuts.c_str()));
      double n41 = t->Draw("g1pt", "procid==41");
      
      double c42 = t->Draw("g1pt", Form("%s&&procid==42", cuts.c_str()));
      double n42 = t->Draw("g1pt", "procid==42");
      
      cout << Form("GGF: %4.3f +/- %4.3f", c40/n40, dEff(static_cast<int>(c40), static_cast<int>(n40))) << endl;
      cout << Form("FEX: %4.3f +/- %4.3f", c41/n41, dEff(static_cast<int>(c41), static_cast<int>(n41))) << endl;
      cout << Form("GSP: %4.3f +/- %4.3f", c42/n42, dEff(static_cast<int>(c42), static_cast<int>(n42))) << endl;
      
      double r = c40/n40; 
      double rE = dEff(static_cast<int>(c40), static_cast<int>(n40));
      fTEX << formatTex(r, Form("%s:ggf%i-%i:val", fSuffix.c_str(), mode, chan), 3) << endl;
      fTEX << formatTex(rE, Form("%s:ggf%i-%i:err", fSuffix.c_str(), mode, chan), 3) << endl;
      if (0 == mode) {
	sgV[chan][0] = r; 
	sgE[chan][0] = rE; 
      } else {
	noV[chan][0] = r; 
	noE[chan][0] = rE; 
      }

      r = c41/n41; 
      rE = dEff(static_cast<int>(c41), static_cast<int>(n41));
      fTEX << formatTex(r, Form("%s:fex%i-%i:val", fSuffix.c_str(), mode, chan), 3) << endl;
      fTEX << formatTex(rE, Form("%s:fex%i-%i:err", fSuffix.c_str(), mode, chan), 3) << endl;
      if (0 == mode) {
	sgV[chan][1] = r; 
	sgE[chan][1] = rE; 
      } else {
	noV[chan][1] = r; 
	noE[chan][1] = rE; 
      }

      
      r = c42/n42; 
      rE = dEff(static_cast<int>(c42), static_cast<int>(n42));
      fTEX << formatTex(r, Form("%s:gsp%i-%i:val", fSuffix.c_str(), mode, chan), 3) << endl;
      fTEX << formatTex(rE, Form("%s:gsp%i-%i:err", fSuffix.c_str(), mode, chan), 3) << endl;
      if (0 == mode) {
	sgV[chan][2] = r; 
	sgE[chan][2] = rE; 
      } else {
	noV[chan][2] = r; 
	noE[chan][2] = rE; 
      }
    }
  }

  double r, rE; 
  for (int i = 0; i < 2; ++i) {
    r = sgV[i][0]/noV[i][0];
    rE= dRatio(sgV[i][0], sgE[i][0], noV[i][0], noE[i][0]); 
    fTEX << formatTex(r, Form("%s:ggfRatio%i:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(rE, Form("%s:ggfRatio%i:err", fSuffix.c_str(), i), 3) << endl;

    r = sgV[i][1]/noV[i][1];
    rE= dRatio(sgV[i][1], sgE[i][1], noV[i][1], noE[i][1]); 
    fTEX << formatTex(r, Form("%s:fexRatio%i:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(rE, Form("%s:fexRatio%i:err", fSuffix.c_str(), i), 3) << endl;

    r = sgV[i][2]/noV[i][2];
    rE= dRatio(sgV[i][2], sgE[i][2], noV[i][2], noE[i][2]); 
    fTEX << formatTex(r, Form("%s:gspRatio%i:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(rE, Form("%s:gspRatio%i:err", fSuffix.c_str(), i), 3) << endl;
  }

}



// ----------------------------------------------------------------------
void plotResults::allInvertedIso() {
  histInvertedIso("fls3d>", 10, 8., 28.); 
  histInvertedIso("chi2/dof<", 10, 1., 2.); 
  histInvertedIso("alpha<", 10, 0.01, 0.06); 
  histInvertedIso("m1pt>", 10, 4., 9.0); 
  histInvertedIso("m2pt>", 10, 4., 6.0); 
  histInvertedIso("pt>", 10, 6.0, 24.0); 
  histInvertedIso("docatrk>", 10, 0.00, 0.05); 
  histInvertedIso("pvips<", 10, 0.5, 4.); 
}

// ----------------------------------------------------------------------  
void plotResults::histInvertedIso(const char *var, int n, double lo, double hi) {

  gStyle->SetOptStat(0);
  
  vector<string> cuts;
  vector<double> cutv;
  vector<string> cutt;
  vector<string> cutf;
  string filename; 

  TH1D *BE = new TH1D("BE", "", n, lo, hi); BE->Sumw2(); 
  TH1D *BO = new TH1D("BO", "", n, lo, hi); BO->Sumw2(); 
  TH1D *EE = new TH1D("EE", "", n, lo, hi); EE->Sumw2(); 
  TH1D *EO = new TH1D("EO", "", n, lo, hi); EO->Sumw2(); 

  for (int ichan = 0; ichan < 2; ++ichan) {
    cuts.clear(); 
    cutv.clear(); 
    cutf.clear();
    
    cuts.push_back("fls3d>");    cutv.push_back(fCuts[ichan]->fls3d);   
    cutt.push_back("l_{3d}/#sigma(l_{3d}) >");
    cutf.push_back("fls3d");

    cuts.push_back("pvips<");    cutv.push_back(fCuts[ichan]->pvips);   
    cutt.push_back("#delta_{3D}/#sigma(#delta_{3d}) <");
    cutf.push_back("pvips");

    cuts.push_back("chi2/dof<"); cutv.push_back(fCuts[ichan]->chi2dof);
    cutt.push_back("#chi^{2}/dof < ");
    cutf.push_back("chi2dof");

    cuts.push_back("alpha<");    cutv.push_back(fCuts[ichan]->alpha);
    cutt.push_back("#alpha < ");
    cutf.push_back("alpha");

    cuts.push_back("m1pt>");     cutv.push_back(fCuts[ichan]->m1pt);
    cutt.push_back("p_{T,#mu1} > ");
    cutf.push_back("m1pt");

    cuts.push_back("m2pt>");     cutv.push_back(fCuts[ichan]->m2pt);
    cutt.push_back("p_{T,#mu2} > ");
    cutf.push_back("m2pt");

    cuts.push_back("docatrk>");  cutv.push_back(fCuts[ichan]->docatrk);
    cutt.push_back("d^{0}_{trk} > ");
    cutf.push_back("docatrk");

    cuts.push_back("pt>");  cutv.push_back(fCuts[ichan]->pt);
    cutt.push_back("p_{T,B} > ");
    cutf.push_back("pt");
    
    string ocuts, tcuts, rcuts;
    for (unsigned int i = 0; i < cuts.size(); ++i) {
      if (!strcmp(cuts[i].c_str(), var)) {
	setTitles(BE, cutt[i].c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 
	setTitles(EE, cutt[i].c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 
	filename = cutf[i];
	continue;
      }
      tcuts += Form("%s %f && ", cuts[i].c_str(), cutv[i]);
    }
    
    string::size_type n1 = tcuts.find_last_of("&&"); 
    ocuts = tcuts.substr(0, n1-1); 
    
    double steps = (hi-lo)/n;
    double cut; 
    for (int i = 0; i < n; ++i) {
      cut = lo + i*steps; 
      rcuts = ocuts + Form(" && %s %f", var, cut); 
      //      cout << "-->" << rcuts << endl;
      
      invertedIso(ichan, rcuts.c_str()); 
      if (0 == ichan) {
	BE->SetBinContent(i+1, fBlExp); BE->SetBinError(i+1, fBlExpE); 
	BO->SetBinContent(i+1, fBlObs); BO->SetBinError(i+1, fBlObsE); 
      } else {
	EE->SetBinContent(i+1, fBlExp); EE->SetBinError(i+1, fBlExpE); 
	EO->SetBinContent(i+1, fBlObs); EO->SetBinError(i+1, fBlObsE); 
      }	
    }
  }
  
  c0->Clear();
  c0->Divide(1,2);
  
  c0->cd(1);
  shrinkPad(0.2, 0.15); 
  double maxi = BE->GetMaximum(); 
  if (BO->GetMaximum() > maxi) {
    maxi = 1.2*BO->GetMaximum();
  } else {
    maxi = 1.2*BE->GetMaximum();
  }
  setFilledHist(BE);
  BE->SetMaximum(maxi); 
  BE->SetMinimum(0.); 
  BE->Draw("hist");
  BO->Draw("esame");
  tl->DrawLatex(0.15, 0.92, "Barrel");

  maxi = EE->GetMaximum(); 
  if (EO->GetMaximum() > maxi) {
    maxi = 1.2*EO->GetMaximum();  
  } else {
    maxi = 1.2*EE->GetMaximum();  
  }
  c0->cd(2);
  shrinkPad(0.2, 0.15); 
  setFilledHist(EE);
  EE->SetMaximum(maxi); 
  EE->SetMinimum(0.); 
  EE->Draw("hist");
  EO->Draw("esame");
  tl->DrawLatex(0.15, 0.92, "Endcap");

  c0->SaveAs(Form("%s/invertedIso-%s.pdf", fDirectory.c_str(), filename.c_str())); 

}


// ----------------------------------------------------------------------
void plotResults::invertedIsoPrediction() {
  
  gStyle->SetOptStat(0); 
  c0->Clear();
  
  string cuts; 
  TH1D *h1;
  for (int i = 0; i < 2; ++i) {
    cuts = string(Form("fls3d>%f", fCuts[i]->fls3d))
      + string(Form("&&chi2/dof<%f", fCuts[i]->chi2dof))
      + string(Form("&&alpha<%f", fCuts[i]->alpha))
      + string(Form("&&pt>%f", fCuts[i]->pt))
      + string(Form("&&m1pt>%f", fCuts[i]->m1pt))
      + string(Form("&&m2pt>%f", fCuts[i]->m2pt))
      //NO!      + string(Form("&&iso5>%f", fCuts[i]->iso1))
      + string(Form("&&docatrk>%f", fCuts[i]->docatrk))
      + string(Form("&&!TMath::IsNaN(fls3d)"))
      ;
    h1 = invertedIso(i, cuts.c_str()); 
    setTitles(h1, "m [GeV]", "Entries/bin"); 
    h1->DrawCopy();
    tl->DrawLatex(0.2, 0.92, (i==0?"Barrel":"Endcap"));
    double lo = h1->Integral(h1->FindBin(fBgLo), h1->FindBin(5.2));
    double hi = h1->Integral(h1->FindBin(5.45), h1->FindBin(fBgHi));
    double bs = h1->Integral(h1->FindBin(fCuts[i]->mBsLo), h1->FindBin(fCuts[i]->mBsHi)); 
    double bd = h1->Integral(h1->FindBin(fCuts[i]->mBdLo), h1->FindBin(fCuts[i]->mBdHi)); 
    double taus= (fCuts[i]->mBsHi - fCuts[i]->mBsLo)/(fBgHi - fBgLo - 0.25);
    double taud= (fCuts[i]->mBdHi - fCuts[i]->mBdLo)/(fBgHi - fBgLo - 0.25);
    double preds = (lo+hi)*taus;
    double relE  = TMath::Sqrt(lo+hi)/(lo+hi);
    double predd = (lo+hi)*taud;
    cout << "channel " << i << endl;
    cout << "taus = " << taus << " taud = " << taud << endl;
    cout << "lo: " << lo << " hi: " << hi << endl;
    cout << "predS = " << preds << " obs = " << bs << endl;
    cout << "predD = " << predd << " obs = " << bd << endl;

    fTEX << formatTex(preds,      Form("%s:invIsoPredBs%i:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(preds*relE, Form("%s:invIsoPredBs%i:err", fSuffix.c_str(), i), 2) << endl;

    fTEX << formatTex(predd,      Form("%s:invIsoPredBd%i:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(predd*relE, Form("%s:invIsoPredBd%i:err", fSuffix.c_str(), i), 2) << endl;

    fTEX << formatTex(bs,              Form("%s:invIsoObsBs%i:val", fSuffix.c_str(), i), 0) << endl;
    fTEX << formatTex(TMath::Sqrt(bs), Form("%s:invIsoObsBs%i:err", fSuffix.c_str(), i), 2) << endl;

    fTEX << formatTex(bd,              Form("%s:invIsoObsBd%i:val", fSuffix.c_str(), i), 0) << endl;
    fTEX << formatTex(TMath::Sqrt(bd), Form("%s:invIsoObsBd%i:err", fSuffix.c_str(), i), 2) << endl;

    c0->SaveAs(Form("%s/%s_invertedIsoPrediction%d.pdf", fDirectory.c_str(), fSuffix.c_str(), i)); 

  }


}



// ----------------------------------------------------------------------
TH1D* plotResults::invertedIso(int chan, const char *cuts) {
  string baseCuts = "gmuid&&hlt&&gmupt&&gmueta&&iso<0.7&&abs(pvlip)<0.05&&abs(pvlips)<2&&closetrk<3";
  string chanCut; 
  if (0 == chan) {
    chanCut = "(abs(m1eta)<1.4&&abs(m2eta)<1.4)";
  } else {
    chanCut = "(abs(m1eta)>1.4||abs(m2eta)>1.4)";
  }
  
  string allCuts  = baseCuts + "&&" + chanCut + "&&" + cuts; 
  cout << allCuts << endl;
  
  fF["SgData"]->cd();
  TH1D *h1(0); 
  if (h1) delete h1; 
  h1 = new TH1D("h1", "", 1000, 4.9, 5.9); 
  TTree *t = (TTree*)gFile->Get("candAnaMuMu/events"); 
  t->Draw("m>>h1", allCuts.c_str(), "goff");
  h1->Draw();

  int blo = h1->FindBin(5.2+0.001); 
  int bhi = h1->FindBin(5.45-0.001); 
  //  cout << "blo: " << blo << " bhi: " << bhi << endl;

  cout << "fBgLo: " << fBgLo << " fBgHi: " << fBgHi << endl;

  double lo = h1->Integral(1, blo);
  double bl = h1->Integral(blo, bhi);
  double blE= TMath::Sqrt(bl); 
  double hi = h1->Integral(bhi, h1->GetNbinsX());
  double ex = (lo+hi)*0.25/(fBgHi-fBgLo-0.25);
  double exE= TMath::Sqrt(ex); 
  double df = bl - ex; 
  double dfE=TMath::Sqrt(blE*blE + exE*exE);
  
  //  cout << "bl: " << bl << " lo: " << lo << " hi: " << hi << endl;
  cout << Form("seen: %4.0f+/-%3.1f", bl, blE) 
       << Form(" expected: %4.1f+/-%3.1f", ex, exE) 
       << Form(" difference: %4.1f+/-%3.1f", df, dfE) 
       << endl;

  fBlExp  = ex; 
  fBlExpE = exE; 
  fBlObs  = (bl>0.?bl:0.01);
  fBlObsE = (bl>0.5?blE:1.);

  return h1; 
}

