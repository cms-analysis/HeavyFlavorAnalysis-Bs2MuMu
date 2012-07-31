#include "plotResults.hh"

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

ClassImp(plotResults)


// ----------------------------------------------------------------------
plotResults::plotResults(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

  fDoPrint = true; 

  string hfname  = fDirectory + "/anaBmm.plotResults." + fSuffix + ".root";
  //  if (fHistFile) fHistFile->Close();
  cout << "open fHistFile: " << hfname << " RECREATE" << endl;
  fHistFile = new TFile(hfname.c_str(), "RECREATE");

  fF["SgData"]->cd();
  TH1D *h1 = new TH1D("hasd", "", 100, 0., 100.); 
  h1->SetDirectory(fHistFile); 
  
  fNumbersFileName = fDirectory + "/anaBmm.plotResults." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  printCuts(cout); 

  fDoUseBDT = true; 
  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 
  fInvertedIso = false; 
  fNormProcessed = false; 
}

// ----------------------------------------------------------------------
plotResults::~plotResults() {
  cout << "plotResults dtor: " << fHistFile << endl;
  if (fHistFile) {
    fHistFile->cd();
    fHistFile->Write();
    fHistFile->Close();
  }
}


// ----------------------------------------------------------------------
void plotResults::makeAll(int channels) {

  //  return;

//   fls3dVsX("pt", "m1pt>4.5&&m2pt>4.2&&hlt", "fls3d_prof_pt.pdf");
//   fls3dVsX("pt", "m1pt>4.5&&m2pt>4.2&&hlt&&docatrk>0.015", "fls3d_prof_docatrk_pt.pdf");
  
//   fls3dVsX("abs(eta)", "m1pt>4.5&&m2pt>4.2&&hlt&&docatrk>0.015", "fls3d_prof_docatrk_eta.pdf");
//   fls3dVsX("abs(eta)", "m1pt>4.5&&m2pt>4.2&&hlt", "fls3d_prof_eta.pdf");

//   fls3dEfficiency("1", "fls3dStudy_muoncuts.pdf");
//   fls3dEfficiency("alpha<0.05", "fls3dStudy_alpha_muoncuts.pdf");
//   fls3dEfficiency("alpha<0.05&&docatrk>0.015", "fls3dStudy_docatrk_alpha_muoncuts.pdf");

//   if (channels & 4) invertedIsolationStudy();

//   zone(1);
//   if (channels & 1) {
//     fNormProcessed = false; 
//     fDoUseBDT = false; 
//     fDoApplyCowboyVeto = false;   
//     fDoApplyCowboyVetoAlsoInSignal = false;   
//     computeNormUL();
//     computeCsBF();
//     acceptancePerProcess();
//   }

  zone(1);
  if (channels & 2) {
    fNormProcessed = false; 
    fDoUseBDT = true; 
    fDoApplyCowboyVeto = false;   
    fDoApplyCowboyVetoAlsoInSignal = false;   
    computeNormUL();
    computeCsBF();
    //acceptancePerProcess();
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

  double yield(0.); 

  for (int i = 0; i < 2; ++i) {
    yield = scaledYield(fNumbersBs[i], fNumbersNo[i], fBF["SgMc"], fsfu);
    yield = scaledYield(fNumbersBd[i], fNumbersNo[i], fBF["BdMc"], 1.);

    fNumbersBs[i]->bsExpObs = fNumbersBs[i]->bgBsExp + fNumbersBs[i]->bsRare + fNumbersBs[i]->bsNoScaled;
  
    fNumbersBs[i]->bsExpObsE = TMath::Sqrt(fNumbersBs[i]->bgBsExpE*fNumbersBs[i]->bgBsExpE //??
					   + fNumbersBs[i]->bsRareE*fNumbersBs[i]->bsRareE
					   + fNumbersBs[i]->bsNoScaledE*fNumbersBs[i]->bsNoScaledE);
    
    fNumbersBs[i]->bdExpObs = fNumbersBs[i]->bgBdExp + fNumbersBs[i]->bdRare + fNumbersBd[i]->bdNoScaled;
    
    fNumbersBs[i]->bdExpObsE = TMath::Sqrt(fNumbersBs[i]->bgBdExpE*fNumbersBs[i]->bgBdExpE //??
					   + fNumbersBs[i]->bdRareE*fNumbersBs[i]->bdRareE
					   + fNumbersBd[i]->bdNoScaledE*fNumbersBd[i]->bdNoScaledE);
  }

  computeErrors(fNumbersNo); 
  computeErrors(fNumbersBs); 
  computeErrors(fNumbersBd); 

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
      * fBF["NoMc"];

    
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
    //    h->Write();
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
void plotResults::computeErrors(std::vector<numbers*> a) {

  double sysMuid[]  = {0.04, 0.08};
  double sysTrig[]  = {0.03, 0.06};
  double sysTrk[]   = {0.04, 0.04};
  double sysAnaSg[] = {0.03, 0.03};
  double sysAnaCs[] = {0.03, 0.03};
  double sysAnaNo[] = {0.04, 0.04};
  double sysCand[]  = {0.01, 0.01};
  double sysAcc[]   = {0.035, 0.05};
  double sysNorm[]  = {0.05, 0.05};
  double sysPSS[]   = {0.05, 0.05};
  double sysTau[]   = {0.04, 0.04};
  
  double sysAna[2]; 
  double syst(0.), total(0.);
  
  string mode("nada"); 
  if (string::npos != a[0]->name.find("signal")) {
    mode = "sg"; 
    sysAna[0] = sysAnaSg[0];
    sysAna[1] = sysAnaSg[1];
  } 

  if (string::npos != a[0]->name.find("normalization")) {
    mode = "no"; 
    sysAna[0] = TMath::Sqrt(sysAnaNo[0]*sysAnaNo[0] + sysTrk[0]*sysTrk[0]);
    sysAna[1] = TMath::Sqrt(sysAnaNo[1]*sysAnaNo[1] + sysTrk[1]*sysTrk[1]);
  } 

  if (string::npos != a[0]->name.find("control")) {
    mode = "cs"; 
    sysAna[0] = TMath::Sqrt(sysAnaCs[0]*sysAnaCs[0] + 2*sysTrk[0]*sysTrk[0]);
    sysAna[1] = TMath::Sqrt(sysAnaCs[1]*sysAnaCs[1] + 2*sysTrk[1]*sysTrk[1]);
  } 

  for (int i = 0; i < 2; ++i) {
    cout << "Setting errors for " << a[i]->name << endl;

    // -- acceptance and selection
    syst = sysAcc[i]*a[i]->acc; 
    a[i]->accTE = TMath::Sqrt(a[i]->accE*a[i]->accE + syst*syst);

    syst = sysCand[i]*a[i]->effCand;
    a[i]->effCandTE = TMath::Sqrt(a[i]->effCandE*a[i]->effCandE + syst*syst);

    syst = sysAna[i]*a[i]->effAna;
    a[i]->effAnaTE = TMath::Sqrt(a[i]->effAnaE*a[i]->effAnaE + syst*syst); 

    // -- total efficiency
    total= TMath::Sqrt(sysAcc[i]*sysAcc[i] + sysAna[i]*sysAna[i] + sysCand[i]*sysCand[i] 
		       + sysMuid[i]*sysMuid[i] + sysTrig[i]*sysTrig[i]); 
    syst = total*a[i]->effTot;
    a[i]->effTotTE = TMath::Sqrt(a[i]->effTotE*a[i]->effTotE + syst*syst); 

    // -- muon id 
    syst = sysMuid[i]*a[i]->effMuidMC;
    a[i]->effMuidTE = TMath::Sqrt(a[i]->effMuidMCE*a[i]->effMuidMCE + syst*syst); 

    syst = sysMuid[i]*a[i]->effMuidMC;
    a[i]->effMuidMCTE = TMath::Sqrt(a[i]->effMuidMCE*a[i]->effMuidMCE + syst*syst); 

    syst = sysMuid[i]*a[i]->effMuidTNP;
    a[i]->effMuidTNPTE = TMath::Sqrt(a[i]->effMuidTNPE*a[i]->effMuidTNPE + syst*syst); 

    syst = sysMuid[i]*a[i]->effMuidTNPMC;
    a[i]->effMuidTNPMCTE = TMath::Sqrt(a[i]->effMuidTNPMCE*a[i]->effMuidTNPMCE + syst*syst); 

    // -- trigger
    syst = sysTrig[i]*a[i]->effTrigMC;
    a[i]->effTrigTE = TMath::Sqrt(a[i]->effTrigMCE*a[i]->effTrigMCE + syst*syst); 

    syst = sysTrig[i]*a[i]->effTrigMC;
    a[i]->effTrigMCTE = TMath::Sqrt(a[i]->effTrigMCE*a[i]->effTrigMCE + syst*syst); 

    syst = sysTrig[i]*a[i]->effTrigTNP;
    a[i]->effTrigTNPTE = TMath::Sqrt(a[i]->effTrigTNPE*a[i]->effTrigTNPE + syst*syst); 

    syst = sysTrig[i]*a[i]->effTrigTNPMC;
    a[i]->effTrigTNPMCTE = TMath::Sqrt(a[i]->effTrigTNPMCE*a[i]->effTrigTNPMCE + syst*syst); 

    // -- yields
    syst = sysNorm[i]*a[i]->fitYield;
    a[i]->fitYieldTE = TMath::Sqrt(a[i]->fitYieldE*a[i]->fitYieldE + syst*syst);
    
    // -- bg estimates
    if (mode == "sg") {
      syst = sysTau[i]*a[i]->bgBsExp; 
      a[i]->bgBsExpTE = TMath::Sqrt(a[i]->bgBsExpE*a[i]->bgBsExpE + syst*syst); 
      
      syst = sysTau[i]*a[i]->bgBdExp; 
      a[i]->bgBdExpTE = TMath::Sqrt(a[i]->bgBdExpE*a[i]->bgBdExpE + syst*syst); 

      syst = sysTau[i]*a[i]->bgBsExp; // the other errors already include the (dominant) systematic contribution
      a[i]->bsExpObsTE = TMath::Sqrt(a[i]->bsExpObsE*a[i]->bsExpObsE + syst*syst);

      syst = sysTau[i]*a[i]->bgBdExp; // the other errors already include the (dominant) systematic contribution
      a[i]->bdExpObsTE = TMath::Sqrt(a[i]->bsExpObsE*a[i]->bsExpObsE + syst*syst);

      syst = sysPSS[i]*a[i]->pss; 
      a[i]->pssTE = TMath::Sqrt(a[i]->pssE*a[i]->pssE + syst*syst);

      syst = sysPSS[i]*a[i]->pds; 
      a[i]->pdsTE = TMath::Sqrt(a[i]->pdsE*a[i]->pdsE + syst*syst);

      syst = sysPSS[i]*a[i]->pdd; 
      a[i]->pddTE = TMath::Sqrt(a[i]->pddE*a[i]->pddE + syst*syst);

      syst = sysPSS[i]*a[i]->psd; 
      a[i]->psdTE = TMath::Sqrt(a[i]->psdE*a[i]->psdE + syst*syst);
    }
  }

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
      fTEX << formatTex(fNumbersNo[i]->effTot, Form("%s:N-EFF-TOT-BPLUS%i:val", fSuffix.c_str(), i), 5) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTotE, Form("%s:N-EFF-TOT-BPLUS%i:err", fSuffix.c_str(), i), 6) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTotTE, Form("%s:N-EFF-TOT-BPLUS%i:tot", fSuffix.c_str(), i), 5) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTot, fNumbersNo[i]->effTotTE, 
			    Form("%s:N-EFF-TOT-BPLUS%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;

      OUT << "ACC_BPLUS\t" << i << "\t" << fNumbersNo[i]->acc <<"\t" << fNumbersNo[i]->accTE  << endl;
      fTEX << formatTex(fNumbersNo[i]->acc, Form("%s:N-ACC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->accE, Form("%s:N-ACC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->accTE, Form("%s:N-ACC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->acc, fNumbersNo[i]->accTE, 
			    Form("%s:N-ACC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      fTEX << formatTex(fNumbersNo[i]->effMuidTNP, Form("%s:N-EFF-MU-PID-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPE, Form("%s:N-EFF-MU-PID-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPTE, Form("%s:N-EFF-MU-PID-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effMuidTNP, fNumbersNo[i]->effMuidTNPTE, 
			    Form("%s:N-EFF-MU-PID-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "EFF_MU_BPLUS\t" << i << "\t"  << fNumbersNo[i]->effMuidMC << "\t"  << "0." << endl;
      //	  << "\t"  << sysMu*fNumbersNo[i]->effMuidMC 
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPMCE, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPMCTE, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effMuidTNPMC, fNumbersNo[i]->effMuidTNPMCTE, 
			    Form("%s:N-EFF-MU-PIDMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      fTEX << formatTex(fNumbersNo[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidMCE, Form("%s:N-EFF-MU-MC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidMCTE, Form("%s:N-EFF-MU-MC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effMuidMC, fNumbersNo[i]->effMuidMCTE, 
			    Form("%s:N-EFF-MU-MC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      fTEX << formatTex(fNumbersNo[i]->effTrigTNP, Form("%s:N-EFF-TRIG-PID-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigTNPE, Form("%s:N-EFF-TRIG-PID-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigTNPTE, Form("%s:N-EFF-TRIG-PID-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTrigTNP, fNumbersNo[i]->effTrigTNPTE, 
			    Form("%s:N-EFF-TRIG-PID-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      fTEX << formatTex(fNumbersNo[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigTNPMCE, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigTNPMCTE, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTrigTNPMC, fNumbersNo[i]->effTrigTNPMCTE, 
			    Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "EFF_TRIG_BPLUS\t" << i << "\t" << fNumbersNo[i]->effTrigMC << "\t" << "0." << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigMCE, Form("%s:N-EFF-TRIG-MC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigMCTE, Form("%s:N-EFF-TRIG-MC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTrigMC, fNumbersNo[i]->effTrigMCTE, 
			    Form("%s:N-EFF-TRIG-MC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "EFF_CAND_BPLUS\t" << i << "\t" << fNumbersNo[i]->effCand << "\t" << fNumbersNo[i]->effCandTE << endl;
      fTEX << formatTex(fNumbersNo[i]->effCand, Form("%s:N-EFF-CAND-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effCandE, Form("%s:N-EFF-CAND-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effCandTE, Form("%s:N-EFF-CAND-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effCand, fNumbersNo[i]->effCandTE, 
			    Form("%s:N-EFF-CAND-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "EFF_ANA_BPLUS\t" << i << "\t" << fNumbersNo[i]->effAna  << "\t" << fNumbersNo[i]->effAnaTE << endl;
      fTEX << formatTex(fNumbersNo[i]->effAna, Form("%s:N-EFF-ANA-BPLUS%i:val", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNo[i]->effAnaE, Form("%s:N-EFF-ANA-BPLUS%i:err", fSuffix.c_str(), i), 4) << endl;
      fTEX << formatTex(fNumbersNo[i]->effAnaTE, Form("%s:N-EFF-ANA-BPLUS%i:tot", fSuffix.c_str(), i), 4) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effAna, fNumbersNo[i]->effAnaTE, 
			    Form("%s:N-EFF-ANA-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      OUT << "OBS_BPLUS\t" << i << "\t" << fNumbersNo[i]->fitYield << "\t" << fNumbersNo[i]->fitYieldTE << endl;
      fTEX << formatTex(fNumbersNo[i]->fitYield, Form("%s:N-OBS-BPLUS%i:val", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(fNumbersNo[i]->fitYieldE, Form("%s:N-OBS-BPLUS%i:err", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(fNumbersNo[i]->fitYieldTE, Form("%s:N-OBS-BPLUS%i:tot", fSuffix.c_str(), i), 0) << endl;
      fTEX << scientificTex(fNumbersNo[i]->fitYield, fNumbersNo[i]->fitYieldTE, 
			    Form("%s:N-OBS-BPLUS%i:all", fSuffix.c_str(), i), 1e3, 0) << endl;

      fTEX << formatTex(fNumbersNo[i]->fitYieldC, Form("%s:N-OBS-CBPLUS%i:val", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(fNumbersNo[i]->fitYieldCE, Form("%s:N-OBS-CBPLUS%i:err", fSuffix.c_str(), i), 0) << endl;
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

    scaledSig  = fNumbersBs[i]->bsNoScaled;
    scaledSigE = fNumbersBs[i]->bsNoScaledE;
    cout << "****** scaledSig(Bs) =   " << scaledSig << endl;
    OUT << "#EXP_SIG_BSMM\t" << i << "\t" << fNumbersBs[i]->bsNoScaled << endl;

    fTEX << formatTex(scaledSig, Form("%s:N-EXP2-SIG-BSMM%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(scaledSigE, Form("%s:N-EXP2-SIG-BSMM%d:err", fSuffix.c_str(), i), 2) << endl;

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

    OUT << "PSS\t" << i << "\t" << fNumbersBs[i]->pss << "\t" << fNumbersBs[i]->pssTE << endl;
    fTEX << formatTex(fNumbersBs[i]->pss, Form("%s:N-PSS%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->pssE, Form("%s:N-PSS%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->pssTE, Form("%s:N-PSS%d:tot", fSuffix.c_str(), i), 3) << endl;

    OUT << "PSD\t" << i << "\t" << fNumbersBd[i]->psd << "\t" << fNumbersBd[i]->psdTE << endl;
    fTEX << formatTex(fNumbersBd[i]->psd, Form("%s:N-PSD%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->psdE, Form("%s:N-PSD%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->psdTE, Form("%s:N-PSD%d:tot", fSuffix.c_str(), i), 3) << endl;

    OUT << "PDS\t" << i << "\t" << fNumbersBs[i]->pds << "\t" << fNumbersBs[i]->pdsTE << endl;
    fTEX << formatTex(fNumbersBs[i]->pds, Form("%s:N-PDS%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->pdsE, Form("%s:N-PDS%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->pdsTE, Form("%s:N-PDS%d:tot", fSuffix.c_str(), i), 3) << endl;

    OUT << "PDD\t" << i << "\t" << fNumbersBd[i]->pdd << "\t" << fNumbersBd[i]->pddTE << endl;
    fTEX << formatTex(fNumbersBd[i]->pdd, Form("%s:N-PDD%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->pddE, Form("%s:N-PDD%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->pddTE, Form("%s:N-PDD%d:tot", fSuffix.c_str(), i), 3) << endl;

    OUT << "#EFF_TOT_BSMM\t" << i << "\t" << fNumbersBs[i]->effTot << endl;
    fTEX << formatTex(fNumbersBs[i]->effTot, Form("%s:N-EFF-TOT-BSMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTotE, Form("%s:N-EFF-TOT-BSMM%d:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTotTE, Form("%s:N-EFF-TOT-BSMM%d:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTot, fNumbersBs[i]->effTotTE, 
			  Form("%s:N-EFF-TOT-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "ACC_BSMM\t" << i << "\t" << fNumbersBs[i]->acc << "\t" << fNumbersBs[i]->accTE << endl;
    fTEX << formatTex(fNumbersBs[i]->acc, Form("%s:N-ACC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->accE, Form("%s:N-ACC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->accTE, Form("%s:N-ACC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->acc, fNumbersBs[i]->accTE, 
			  Form("%s:N-ACC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBs[i]->effMuidTNP, Form("%s:N-EFF-MU-PID-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidTNPE, Form("%s:N-EFF-MU-PID-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidTNPTE, Form("%s:N-EFF-MU-PID-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effMuidTNP, fNumbersBs[i]->effMuidTNPTE, 
			  Form("%s:N-EFF-MU-PID-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBs[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidTNPMCE, Form("%s:N-EFF-MU-PIDMC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidTNPMCTE, Form("%s:N-EFF-MU-PIDMC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effMuidTNPMC, fNumbersBs[i]->effMuidTNPMCTE, 
			  Form("%s:N-EFF-MU-PIDMC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_MU_BSMM\t" << i << "\t" << fNumbersBs[i]->effMuidMC << "\t" << fNumbersBs[i]->effMuidMCTE << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidMCE, Form("%s:N-EFF-MU-MC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidMCTE, Form("%s:N-EFF-MU-MC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effMuidMC, fNumbersBs[i]->effMuidMCTE, 
			  Form("%s:N-EFF-MU-MC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBs[i]->effTrigTNP, Form("%s:N-EFF-TRIG-PID-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigTNPE, Form("%s:N-EFF-TRIG-PID-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigTNPTE, Form("%s:N-EFF-TRIG-PID-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTrigTNP, fNumbersBs[i]->effTrigTNPTE, 
			  Form("%s:N-EFF-TRIG-PID-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBs[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigTNPMCE, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigTNPMCTE, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTrigTNPMC, fNumbersBs[i]->effTrigTNPMCTE, 
			  Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_TRIG_BSMM\t" << i << "\t" << fNumbersBs[i]->effTrigMC << "\t" << fNumbersBs[i]->effTrigMCTE << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigMCE, Form("%s:N-EFF-TRIG-MC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigMCTE, Form("%s:N-EFF-TRIG-MC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTrigMC, fNumbersBs[i]->effTrigMCTE, 
			  Form("%s:N-EFF-TRIG-MC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_CAND_BSMM\t" << i << "\t" << fNumbersBs[i]->effCand << "\t" << fNumbersBs[i]->effCandTE << endl;
    fTEX << formatTex(fNumbersBs[i]->effCand, Form("%s:N-EFF-CAND-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effCandE, Form("%s:N-EFF-CAND-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effCandTE, Form("%s:N-EFF-CAND-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effCand, fNumbersBs[i]->effCandTE, 
			  Form("%s:N-EFF-CAND-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_ANA_BSMM\t" << i << "\t" << fNumbersBs[i]->effAna << "\t" << fNumbersBs[i]->effAnaTE << endl;
    fTEX << formatTex(fNumbersBs[i]->effAna, Form("%s:N-EFF-ANA-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effAnaE, Form("%s:N-EFF-ANA-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effAnaTE, Form("%s:N-EFF-ANA-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effAna, fNumbersBs[i]->effAnaTE, 
			  Form("%s:N-EFF-ANA-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "#EFF_TOT_BDMM\t" << i << "\t" << fNumbersBd[i]->effTot << endl;
    fTEX << formatTex(fNumbersBd[i]->effTot, Form("%s:N-EFF-TOT-BDMM%d:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTotE, Form("%s:N-EFF-TOT-BDMM%d:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTotTE, Form("%s:N-EFF-TOT-BDMM%d:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTot, fNumbersBs[i]->effTotTE, 
			  Form("%s:N-EFF-TOT-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "ACC_BDMM\t" << i << "\t" << fNumbersBd[i]->acc << "\t" << fNumbersBd[i]->accTE << endl;
    fTEX << formatTex(fNumbersBd[i]->acc, Form("%s:N-ACC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->accE, Form("%s:N-ACC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->accTE, Form("%s:N-ACC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->acc, fNumbersBd[i]->accTE, 
			  Form("%s:N-ACC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBd[i]->effMuidTNP, Form("%s:N-EFF-MU-PID-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidTNPE, Form("%s:N-EFF-MU-PID-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidTNPTE, Form("%s:N-EFF-MU-PID-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effMuidTNP, fNumbersBd[i]->effMuidTNPTE, 
			  Form("%s:N-EFF-MU-PID-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBd[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidTNPMCE, Form("%s:N-EFF-MU-PIDMC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidTNPMCTE, Form("%s:N-EFF-MU-PIDMC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effMuidTNPMC, fNumbersBd[i]->effMuidTNPMCTE, 
			  Form("%s:N-EFF-MU-PIDMC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_MU_BDMM\t" << i	<< "\t" << fNumbersBd[i]->effMuidMC << "\t" << fNumbersBd[i]->effMuidMCTE << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidMC, Form("%s:N-EFF-MU-MC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidMCE, Form("%s:N-EFF-MU-MC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidMCTE, Form("%s:N-EFF-MU-MC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effMuidMC, fNumbersBd[i]->effMuidMCTE, 
			  Form("%s:N-EFF-MU-MC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBd[i]->effTrigTNP, Form("%s:N-EFF-TRIG-PID-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigTNPE, Form("%s:N-EFF-TRIG-PID-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigTNPTE, Form("%s:N-EFF-TRIG-PID-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effTrigTNP, fNumbersBd[i]->effTrigTNPTE, 
			  Form("%s:N-EFF-TRIG-PID-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBd[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigTNPMCE, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigTNPMCTE, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effTrigTNPMC, fNumbersBd[i]->effTrigTNPMCTE, 
			  Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_TRIG_BDMM\t" << i << "\t" << fNumbersBd[i]->effTrigMC << "\t" << fNumbersBd[i]->effTrigMCTE << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigMC, Form("%s:N-EFF-TRIG-MC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigMCE, Form("%s:N-EFF-TRIG-MC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigMCTE, Form("%s:N-EFF-TRIG-MC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effTrigMC, fNumbersBd[i]->effTrigMCTE, 
			  Form("%s:N-EFF-TRIG-MC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_CAND_BDMM\t" << i << "\t" << fNumbersBd[i]->effCand << "\t" << fNumbersBd[i]->effCandTE << endl;
    fTEX << formatTex(fNumbersBd[i]->effCand, Form("%s:N-EFF-CAND-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effCandE, Form("%s:N-EFF-CAND-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effCandTE, Form("%s:N-EFF-CAND-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effCand, fNumbersBd[i]->effCandTE,
			  Form("%s:N-EFF-CAND-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "EFF_ANA_BDMM\t" << i << "\t" << fNumbersBd[i]->effAna << "\t" << fNumbersBd[i]->effAnaTE << endl;
    fTEX << formatTex(fNumbersBd[i]->effAna, Form("%s:N-EFF-ANA-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effAnaE, Form("%s:N-EFF-ANA-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effAnaTE, Form("%s:N-EFF-ANA-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effAna, fNumbersBd[i]->effAnaTE, 
			  Form("%s:N-EFF-ANA-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "# Expected in signal boxes" << endl;
    OUT << "#EXP_OBS_BSMM\t" << i << "\t" << fNumbersBs[i]->bsExpObs << "\t" << fNumbersBs[i]->bsExpObsTE << endl;
    fTEX << formatTex(fNumbersBs[i]->bsExpObs, Form("%s:N-EXP-OBS-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bsExpObsE, Form("%s:N-EXP-OBS-BS%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "#EXP_OBS_BDMM\t" << i << "\t" << fNumbersBs[i]->bdExpObs << "\t" << fNumbersBs[i]->bdExpObsTE << endl;
    fTEX << formatTex(fNumbersBs[i]->bdExpObs, Form("%s:N-EXP-OBS-BD%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bdExpObsE, Form("%s:N-EXP-OBS-BD%d:err", fSuffix.c_str(), i), 2) << endl;

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

    OUT << "PEAK_BKG_BS\t" << i << "\t" << fNumbersBs[i]->bsRare << "\t" << fNumbersBs[i]->bsRareE << endl;
    fTEX << formatTex(fNumbersBs[i]->bsRare, Form("%s:N-PEAK-BKG-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bsRareE,Form("%s:N-PEAK-BKG-BS%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "PEAK_BKG_BD\t" << i << "\t"	<< fNumbersBs[i]->bdRare << "\t" << fNumbersBs[i]->bdRareE << endl;
    fTEX << formatTex(fNumbersBs[i]->bdRare, Form("%s:N-PEAK-BKG-BD%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->bdRareE,Form("%s:N-PEAK-BKG-BD%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "TAU_BS\t" << i << "\t" << fNumbersBs[i]->tauBs << "\t" << fNumbersBs[i]->tauBsE << endl;
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
void plotResults::rareBg(std::string mode) {

  if (!fNormProcessed) {
    setupNorm();
  }


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

  // -- 'tight muon' values
  double epsMu(0.75); // FIXME: this should be different for barrel and endcap!
  double epsPi(0.001),  errPi2(0.15*0.15); // relative errors on misid rates are statistical error from Danek's fits
  double epsKa(0.001),  errKa2(0.15*0.15); 
  double epsPr(0.0005), errPr2(0.15*0.15);

//   // -- not quite tight muon values
//   double epsMu(0.8); // FIXME: this should be different for barrel and endcap!
//   double epsPi(0.0015), errPi2(0.15*0.15); // relative errors on misid rates are statistical error from Danek's fits
//   double epsKa(0.0017), errKa2(0.15*0.15); 
//   double epsPr(0.0005), errPr2(0.15*0.15);

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

  colors.insert(make_pair("bgBd2PiMuNu", 43));hatches.insert(make_pair("bgBd2PiMuNu", 3008));mscale.insert(make_pair("bgBd2PiMuNu", epsPi*epsMu));
  chanbf.insert(make_pair("bgBd2PiMuNu", 1.3e-4)); 
  err.insert(make_pair("bgBd2PiMuNu", TMath::Sqrt(0.2*0.2 + errPi2))); 


  colors.insert(make_pair("bgBu2PiMuMu", 50));hatches.insert(make_pair("bgBu2PiMuMu", 3004));mscale.insert(make_pair("bgBu2PiMuMu", epsMu*epsMu));
  chanbf.insert(make_pair("bgBu2PiMuMu", 1.0e-8)); 
  err.insert(make_pair("bgBu2PiMuMu", TMath::Sqrt(0.2*0.2))); 

  colors.insert(make_pair("bgBu2KMuMu", 51));hatches.insert(make_pair("bgBu2KMuMu", 3005));mscale.insert(make_pair("bgBu2KMuMu", epsMu*epsMu));
  chanbf.insert(make_pair("bgBu2KMuMu", 4.3e-7)); 
  err.insert(make_pair("bgBu2KMuMu", TMath::Sqrt(0.2*0.2))); 

  colors.insert(make_pair("bgBu2KMuMu", 52));hatches.insert(make_pair("bgBu2KMuMu", 3005));mscale.insert(make_pair("bgBu2KMuMu", epsMu*epsMu));
  chanbf.insert(make_pair("bgBu2KMuMu", 5.0e-7)); 
  err.insert(make_pair("bgBu2KMuMu", TMath::Sqrt(0.1*0.1))); 

  colors.insert(make_pair("bgBs2PhiMuMu",53));hatches.insert(make_pair("bgBs2PhiMuMu",3006));mscale.insert(make_pair("bgBs2PhiMuMu",epsMu*epsMu));
  chanbf.insert(make_pair("bgBs2PhiMuMu", 1.4e-6)); 
  err.insert(make_pair("bgBs2PhiMuMu", TMath::Sqrt(0.4*0.4))); 

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

    if (mode != "nada" && string::npos == imap->first.find(mode)) {
      continue;
    }

    if (0 == fF[imap->first]) {
      continue; 
    }

    double misid = mscale[imap->first];
    double ngenfile = ((TH1D*)fF[imap->first]->Get("monEvents"))->GetBinContent(1); 
    cout << "======> " << imap->first << " with misid: " << misid << endl;

    fF[imap->first]->cd("candAnaMuMu");

    TH1D *hRare[2]; 
    double tot0, tot, bd, bs, efftot, pss, pdd;
    
    fRareName = imap->first; 
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

      double loSb = hRare[ichan]->Integral(hRare[ichan]->FindBin(4.9001), hRare[ichan]->FindBin(5.199));
      
      fTEX <<  Form("\\vdef{%s:%s:loSideband%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), imap->first.c_str(), ichan, 
		    loSb) << endl;
      fTEX <<  Form("\\vdef{%s:%s:loSideband%d:err}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), imap->first.c_str(), ichan, 
		    loSb*0.3) << endl;
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
  if (mode == "nada") hhRareBg0->SetAxisRange(4.9, 5.9, "X"); 
  setTitles(hhRareBg0, "m_{#mu #mu} [GeV]", Form("Candidates/%4.3f GeV", hhRareBg0->GetBinWidth(1)), 0.06, 0.9, 1.5);
  legg->SetHeader("CMS simulation"); 
  legg->Draw(); 
  hhRareBg0->Draw("same");
  stamp(0.18, fStampString, 0.67, fStampCms); 
  double size = tl->GetTextSize();
  tl->SetTextSize(0.07); 
  tl->DrawLatex(0.25, 0.8, "Barrel");   
  tl->SetTextSize(size); 

  string pdfname;
  if (fDoUseBDT) pdfname = Form("%s/%s_bdt_rare0.pdf", fDirectory.c_str(), fSuffix.c_str());
  else  pdfname = Form("%s/%s_cnc_rare0.pdf", fDirectory.c_str(), fSuffix.c_str());
  if (mode != "nada") pdfname = Form("%s/%s_cnc_%s_0.pdf", fDirectory.c_str(), mode.c_str(), fSuffix.c_str());
  if (fDoPrint) {
    if (c0) c0->SaveAs(pdfname.c_str());
  }

  hRareBg1->SetMaximum(1.0); 
  hRareBg1->Draw();
  TH1D *hhRareBg1 = (TH1D*)hRareBg1->GetHistogram(); 
  if (mode == "nada") hhRareBg1->SetAxisRange(4.9, 5.9, "X"); 
  setTitles(hhRareBg1, "m_{#mu #mu} [GeV]", Form("Candidates/%4.3f GeV", hhRareBg0->GetBinWidth(1)), 0.06, 0.9, 1.5);
  legg->SetHeader("CMS simulation"); 
  legg->Draw(); 
  hhRareBg1->Draw("same");
  if (fDoUseBDT) pdfname = Form("%s/%s_bdt_rare1.pdf", fDirectory.c_str(), fSuffix.c_str());
  else  pdfname = Form("%s/%s_cnc_rare1.pdf", fDirectory.c_str(), fSuffix.c_str());
  if (mode != "nada") pdfname = Form("%s/%s_cnc_%s_1.pdf", fDirectory.c_str(), mode.c_str(), fSuffix.c_str());
  stamp(0.18, fStampString, 0.67, fStampCms); 
  tl->SetTextSize(0.07); 
  tl->DrawLatex(0.25, 0.8, "Endcap");   
  tl->SetTextSize(size); 
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
  //  bRare->Write(); 
  eRare->SetDirectory(gDirectory);
  //  eRare->Write(); 
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
void plotResults::invertedIsolationStudy() {

  fInvertedIso = true; 

  fTEX.close();
  fTEX.open(Form("%s/anaBmm.invertedIsolation.default-11.tex", fDirectory.c_str()));

  // -- determine the inverted isolation prediction and observation without cut variation
  readCuts(fCutsFileName.c_str()); 
  determineInvertedIsolationYield(1);

  // -- scan a few cuts
  TH1D *h0, *h1, *x0, *x1;
  string axis(""), pdfname(""); 
  double lo(5.0), hi(15.0), cut(0.); 
  int npoints(5); 

  // -- m1pt
  readCuts(fCutsFileName.c_str()); 
  npoints = 5; 
  lo = 4.0; 
  hi = 5.0;
  h0 = new TH1D("h0", "", npoints, lo, hi); 
  h1 = new TH1D("h1", "", npoints, lo, hi); 
  x0 = new TH1D("x0", "", npoints, lo, hi); x0->Sumw2();
  x1 = new TH1D("x1", "", npoints, lo, hi); x1->Sumw2();
  axis = "p_{T,#mu1} > ";
  pdfname = "muon1pt";
  setTitles(x0, axis.c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 
  setTitles(x1, axis.c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 

  for (int j = 0; j < npoints; ++j) {
    cut = lo + j*(hi-lo)/npoints;
    fCuts[0]->m1pt = cut;
    fCuts[1]->m1pt = cut;
    printCuts(cout);
    determineInvertedIsolationYield();
    x0->SetBinContent(j+1, fBl0Exp); x0->SetBinError(j+1, fBl0ExpE); 
    x1->SetBinContent(j+1, fBl1Exp); x1->SetBinError(j+1, fBl1ExpE); 
    h0->SetBinContent(j+1, fBl0Obs); 
    h1->SetBinContent(j+1, fBl1Obs); 
  }
 
  plotInvertedIsolationScan(pdfname, h0, h1, x0, x1);


  // -- fls3d
  readCuts(fCutsFileName.c_str()); 
  npoints = 5; 
  lo = 5.0; 
  hi = 15.0;
  h0 = new TH1D("h0", "", npoints, lo, hi); 
  h1 = new TH1D("h1", "", npoints, lo, hi); 
  x0 = new TH1D("x0", "", npoints, lo, hi); x0->Sumw2();
  x1 = new TH1D("x1", "", npoints, lo, hi); x1->Sumw2();
  axis = "l/#sigma(l_{3d}) > ";
  pdfname = "fls3d";
  setTitles(x0, axis.c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 
  setTitles(x1, axis.c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 

  for (int j = 0; j < npoints; ++j) {
    cut = lo + j*(hi-lo)/npoints;
    fCuts[0]->fls3d = cut;
    fCuts[1]->fls3d = cut;
    printCuts(cout);
    determineInvertedIsolationYield();
    x0->SetBinContent(j+1, fBl0Exp); x0->SetBinError(j+1, fBl0ExpE); 
    x1->SetBinContent(j+1, fBl1Exp); x1->SetBinError(j+1, fBl1ExpE); 
    h0->SetBinContent(j+1, fBl0Obs); 
    h1->SetBinContent(j+1, fBl1Obs); 
  }
 
  plotInvertedIsolationScan(pdfname, h0, h1, x0, x1);


  // -- alpha
  readCuts(fCutsFileName.c_str()); 
  npoints = 5; 
  lo = 0.05; 
  hi = 0.15;
  h0 = new TH1D("h0", "", npoints, lo, hi); 
  h1 = new TH1D("h1", "", npoints, lo, hi); 
  x0 = new TH1D("x0", "", npoints, lo, hi); x0->Sumw2();
  x1 = new TH1D("x1", "", npoints, lo, hi); x1->Sumw2();
  axis = "#alpha < ";
  pdfname = "alpha";
  setTitles(x0, axis.c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 
  setTitles(x1, axis.c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 

  for (int j = 0; j < npoints; ++j) {
    cut = lo + j*(hi-lo)/npoints;
    fCuts[0]->alpha = cut;
    fCuts[1]->alpha = cut;
    printCuts(cout);
    determineInvertedIsolationYield();
    x0->SetBinContent(j+1, fBl0Exp); x0->SetBinError(j+1, fBl0ExpE); 
    x1->SetBinContent(j+1, fBl1Exp); x1->SetBinError(j+1, fBl1ExpE); 
    h0->SetBinContent(j+1, fBl0Obs); 
    h1->SetBinContent(j+1, fBl1Obs); 
  }
 
  plotInvertedIsolationScan(pdfname, h0, h1, x0, x1);


  // -- chi2/dof
  readCuts(fCutsFileName.c_str()); 
  npoints = 5; 
  lo = 2.0; 
  hi = 4.0;
  h0 = new TH1D("h0", "", npoints, lo, hi); 
  h1 = new TH1D("h1", "", npoints, lo, hi); 
  x0 = new TH1D("x0", "", npoints, lo, hi); x0->Sumw2();
  x1 = new TH1D("x1", "", npoints, lo, hi); x1->Sumw2();
  axis = "#chi^{2}/dof < ";
  pdfname = "chi2dof";
  setTitles(x0, axis.c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 
  setTitles(x1, axis.c_str(), "Candidates", 0.08, 1., 0.9, 0.07); 

  for (int j = 0; j < npoints; ++j) {
    cut = lo + j*(hi-lo)/npoints;
    fCuts[0]->chi2dof = cut;
    fCuts[1]->chi2dof = cut;
    printCuts(cout);
    determineInvertedIsolationYield();
    x0->SetBinContent(j+1, fBl0Exp); x0->SetBinError(j+1, fBl0ExpE); 
    x1->SetBinContent(j+1, fBl1Exp); x1->SetBinError(j+1, fBl1ExpE); 
    h0->SetBinContent(j+1, fBl0Obs); 
    h1->SetBinContent(j+1, fBl1Obs); 
  }
 
  plotInvertedIsolationScan(pdfname, h0, h1, x0, x1);

  // -- reset cuts to default
  readCuts(fCutsFileName.c_str()); 

  fTEX.close();
  fTEX.open(fNumbersFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotResults::plotInvertedIsolationScan(string pdfname, TH1D *h0, TH1D *h1, TH1D *x0, TH1D *x1) {
  zone(1,2);
  c0->cd(1);
  shrinkPad(0.2, 0.15); 
  double maxi = h0->GetMaximum(); 
  if (x0->GetMaximum() > maxi) {
    maxi = 1.2*x0->GetMaximum();
  } else {
    maxi = 1.2*h0->GetMaximum();
  }
  //  setFilledHist(x0);
  x0->SetMaximum(maxi); 
  x0->SetMinimum(0.); 
  x0->SetMarkerSize(0.); 
  x0->Draw("ehist");
  h0->Draw("e1same");
  tl->SetTextColor(kBlack);
  tl->DrawLatex(0.15, 0.92, "Barrel");
  
  maxi = x1->GetMaximum(); 
  if (h1->GetMaximum() > maxi) {
    maxi = 1.2*h1->GetMaximum();  
  } else {
    maxi = 1.2*x1->GetMaximum();  
  }
  c0->cd(2);
  shrinkPad(0.2, 0.15); 
  //  setFilledHist(x1);
  x1->SetMarkerSize(0.); 
  x1->SetMaximum(maxi); 
  x1->SetMinimum(0.); 
  x1->Draw("ehist");
  h1->Draw("e1same");
  tl->DrawLatex(0.15, 0.92, "Endcap");
  
  c0->SaveAs(Form("%s/invertedIso-%s.pdf", fDirectory.c_str(), pdfname.c_str())); 
  c0->Clear();
}


// ----------------------------------------------------------------------
void plotResults::determineInvertedIsolationYield(int print) {
  
  zone(2,2);
  loopTree(5); 
  TH1D *hd[] = {(TH1D*)(fhMassWithAllCuts[0]->Clone("inviso_cnc_0")),
		(TH1D*)(fhMassWithAllCuts[1]->Clone("inviso_cnc_1"))
  };

  // -- must run on norm sample to get correct scaling factor (which depends on the cuts!)
  loopTree(15); 
  loopTree(10); 

  c0->cd(2);
  pair<TH1D*, TH1D*> hrare = singleRelativeYield("bgBd2PiMuNu");
  TH1D *hr[] = {hrare.first, hrare.second};
  // -- correct for missing rest
  hr[0]->Scale(2.0);
  hr[1]->Scale(2.0);
  
  double eps(0.0001);
  for (int i = 0; i < 2; ++i) {
    double dlo = hd[i]->Integral(hd[i]->FindBin(fBgLo), hd[i]->FindBin(5.2-eps));
    double dhi = hd[i]->Integral(hd[i]->FindBin(5.45+eps), hd[i]->FindBin(fBgHi));
    double dbs = hd[i]->Integral(hd[i]->FindBin(fCuts[i]->mBsLo+eps), hd[i]->FindBin(fCuts[i]->mBsHi-eps)); 
    double dbd = hd[i]->Integral(hd[i]->FindBin(fCuts[i]->mBdLo+eps), hd[i]->FindBin(fCuts[i]->mBdHi-eps)); 

    double rlo = hr[i]->Integral(hr[i]->FindBin(fBgLo), hr[i]->FindBin(5.2-eps));
    double rhi = hr[i]->Integral(hr[i]->FindBin(5.45+eps), hr[i]->FindBin(fBgHi));
    double rbs = hr[i]->Integral(hr[i]->FindBin(fCuts[i]->mBsLo+eps), hr[i]->FindBin(fCuts[i]->mBsHi-eps)); 
    double rbd = hr[i]->Integral(hr[i]->FindBin(fCuts[i]->mBdLo+eps), hr[i]->FindBin(fCuts[i]->mBdHi-eps)); 

    // -- per decay mode in signal windows
    double taus= (fCuts[i]->mBsHi - fCuts[i]->mBsLo)/(fBgHi - fBgLo - 0.25);
    double taud= (fCuts[i]->mBdHi - fCuts[i]->mBdLo)/(fBgHi - fBgLo - 0.25);
    double preds = (dlo+dhi-rlo-rhi)*taus + rbs;
    double relE  = TMath::Sqrt(dlo+dhi)/(dlo+dhi);
    double predd = (dlo+dhi-rlo-rhi)*taud + rbd;
    cout << "channel " << i << endl;
    cout << "taus = " << taus << " taud = " << taud << endl;
    cout << "dlo: " << dlo << " dhi: " << dhi << " rlo: " << rlo << " rhi: " << rhi << endl;
    cout << "predS = " << preds << " obs = " << dbs << endl;
    cout << "predD = " << predd << " obs = " << dbd << endl;

    if (1 == print) {
      fTEX << formatTex(preds,      Form("%s:invIsoPredBs%i:val", fSuffix.c_str(), i), 2) << endl;
      fTEX << formatTex(preds*relE, Form("%s:invIsoPredBs%i:err", fSuffix.c_str(), i), 2) << endl;
      
      fTEX << formatTex(predd,      Form("%s:invIsoPredBd%i:val", fSuffix.c_str(), i), 2) << endl;
      fTEX << formatTex(predd*relE, Form("%s:invIsoPredBd%i:err", fSuffix.c_str(), i), 2) << endl;
      
      fTEX << formatTex(dbs,              Form("%s:invIsoObsBs%i:val", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(TMath::Sqrt(dbs), Form("%s:invIsoObsBs%i:err", fSuffix.c_str(), i), 2) << endl;
      
      fTEX << formatTex(dbd,              Form("%s:invIsoObsBd%i:val", fSuffix.c_str(), i), 0) << endl;
      fTEX << formatTex(TMath::Sqrt(dbd), Form("%s:invIsoObsBd%i:err", fSuffix.c_str(), i), 2) << endl;
    
      zone(1);
      hd[i]->Draw();
      tl->SetTextColor(kBlack);
      if (0 == i) tl->DrawLatex(0.15, 0.92, "Barrel");
      if (1 == i) tl->DrawLatex(0.15, 0.92, "Endcap");
      c0->SaveAs(Form("%s/%s_invertedIsoPrediction%d.pdf", fDirectory.c_str(), fSuffix.c_str(), i)); 
    }

    // -- in blinding window
    int blo = hd[i]->FindBin(5.2+eps); 
    int bhi = hd[i]->FindBin(5.45-eps); 
    double dbds = hd[i]->Integral(blo, bhi);

    taus = 0.25/(fBgHi-fBgLo-0.25);
    double predds = (dlo+dhi-rlo-rhi)*taus + rbs;

    if (0 == i) {
      fBl0Exp  = predds; 
      fBl0ExpE = predds*relE; 
      fBl0Obs  = (dbds>0.?dbds:0.01);
      fBl0ObsE = (dbds>0.5?dbds*relE:1.);
    } else {
      fBl1Exp  = predds;      
      fBl1ExpE = predds*relE; 
      fBl1Obs  = (dbds>0.?dbds:0.01);
      fBl1ObsE = (dbds>0.5?dbds*relE:1.);
    }


  }

  //  c0->Clear();
  zone(2,2);
  hd[0]->Draw();

  c0->cd(2); 
  hd[1]->Draw();

  c0->cd(3); 
  hr[0]->Draw();

  c0->cd(4); 
  hr[1]->Draw();

  c0->Clear();
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
  tl->SetTextColor(kBlack);
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
    tl->SetTextColor(kBlack);
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



// ----------------------------------------------------------------------
pair<TH1D*, TH1D*> plotResults::singleRelativeYield(std::string fstring) {

  fF[fstring]->cd("candAnaMuMu");
  gDirectory->pwd();
  loopTree(99); 

  TH1D *h[2];
  cout << "fstring:   " << fstring << endl;
  double ngenfile = ((TH1D*)fF[fstring.c_str()]->Get("monEvents"))->GetBinContent(1); 
  cout << "ngenfile:  " << ngenfile << endl;
  double filtereff = fFilterEff[fstring.c_str()];
  if (filtereff < 0.001) filtereff = 1; 
  cout << "filtereff: " << filtereff << endl;

  // -- definition of (mis)id and trigger efficiency numbers
  double mutrig[2] = {0.821, 0.726};
  double mumuid[2] = {0.77, 0.895};
  double pimuid[2] = {1.0e-3, 1.0e-3};
  double kamuid[2] = {1.0e-3, 1.0e-3};
  double prmuid[2] = {0.5e-3, 0.5e-3};
  
  string mode(""); 

  if (string::npos != fstring.find("KP")) mode = "KP";
  if (string::npos != fstring.find("KPi")) mode = "KPi";
  if (string::npos != fstring.find("KK")) mode = "KK";

  if (string::npos != fstring.find("PiPi")) mode = "PiPi";

  if (string::npos != fstring.find("KMu")) mode = "KMu";
  if (string::npos != fstring.find("PiMu")) mode = "PiMu";
  if (string::npos != fstring.find("PMu")) mode = "PMu";

  double bfn = fBF["NoMc"]; 
  double bf1 = fBF[fstring]; 
  cout << " ==> BF: " << bf1 << " mode = " << mode << endl;

  double total, scaledYield, nx;
  double eff, muScale, brScale, prScale, scale; 
  for (int ichan = 0; ichan < 2; ++ichan) {

    if (mode == "KP")   muScale = kamuid[ichan]*prmuid[ichan];
    if (mode == "KPi")  muScale = kamuid[ichan]*pimuid[ichan];
    if (mode == "KK")   muScale = kamuid[ichan]*kamuid[ichan];
    if (mode == "PiPi") muScale = pimuid[ichan]*pimuid[ichan];

    if (mode == "KMu")  muScale = kamuid[ichan]*mumuid[ichan];
    if (mode == "PiMu") muScale = pimuid[ichan]*mumuid[ichan];
    if (mode == "PMu")  muScale = prmuid[ichan]*mumuid[ichan];
   
    muScale *= mutrig[ichan];
 
    if (fDoUseBDT) {
      h[ichan] = (TH1D*)(fhMassWithAllCuts[ichan]->Clone(Form("h1Rare_BDT_%d", ichan)));  
      h[ichan]->SetTitle(Form("%s (BDT, %s)", fstring.c_str(), (0 == ichan?"barrel":"endcap")));
    } else {
      h[ichan] = (TH1D*)(fhMassWithAllCuts[ichan]->Clone(Form("h1Rare_CNC_%d", ichan)));  
      h[ichan]->SetTitle(Form("%s (baseline, %s)", fstring.c_str(), (0 == ichan?"barrel":"endcap")));
    }
    
    double orig0 = h[ichan]->GetSumOfWeights();
    double orig1 = h[ichan]->Integral(h[ichan]->FindBin(4.9001), h[ichan]->FindBin(5.8999)); 
    total   = h[ichan]->Integral(0, h[ichan]->GetNbinsX()+1); 
    eff     = total/static_cast<double>(ngenfile)*filtereff;
    brScale = fBF[fstring]/fBF["NoMc"];
    prScale = fProdR[fstring];

    scale = muScale * fProdR[fstring] * (eff/fNumbersNo[ichan]->effTot) * brScale;

    scaledYield = scale * fNumbersNo[ichan]->fitYield; 
    nx          = h[ichan]->Integral(0, h[ichan]->GetNbinsX()+1);
    h[ichan]->Scale(scaledYield/nx);

    cout << "channel " << ichan << ": " << total << " -> eff = " << eff << endl;
    cout << " ==> mutrig * muid * misid: " << muScale << endl;
    cout << " ==> eff:         " << eff << endl;
    cout << " ==> brScale:     " << brScale << " (" << fBF[fstring] << "/" << fBF["NoMc"] << ")" << endl;
    cout << " ==> prScale:     " << prScale << endl;
    cout << " ==> scale:       " << scale << endl;
    cout << " ==> n(B+):       " << fNumbersNo[ichan]->fitYield << endl;
    cout << " ==> nx:          " << nx 
	 << " in histogram: " << orig0
	 << " in range: " << orig1
	 << endl;
    cout << " ==> scaledYield: " << scaledYield
	 << " in histogram: " << h[ichan]->GetSumOfWeights() 
	 << " in range: " << h[ichan]->Integral(h[ichan]->FindBin(4.9001), h[ichan]->FindBin(5.8999))
	 << endl;

  }

  return(make_pair(h[0], h[1])); 
  
}


// ----------------------------------------------------------------------
void plotResults::fls3dEfficiency(string cuts, string pdfname) {

  gStyle->SetOptStat(0);

  zone(1,2);
  fF["SgMc"]->cd("candAnaMuMu");
  TH1D *h1 = new TH1D("h1", "", 100, 0., 5); h1->Sumw2(); setHist(h1, kBlack);
  TH1D *h2 = new TH1D("h2", "", 100, 0., 5); h2->Sumw2(); setHist(h2, kBlue);
  TH1D *h3 = new TH1D("h3", "", 100, 0., 5); h3->Sumw2(); setHist(h3, kRed);
  TH1D *hr = new TH1D("hr", "", 100, 0., 5); hr->Sumw2();
  TH1D *hs = new TH1D("hs", "", 100, 0., 5); hs->Sumw2();
  TTree *t = (TTree*)gDirectory->Get("events"); 
  
  string baseCuts = "m1pt>4.5&&m2pt>4.2&&hlt&&" + cuts; 
  t->Draw("1e12*tau>>h1", baseCuts.c_str(), "goff");
  setTitles(h1, "#tau [ps]", "a.u.");
  h1->Draw("hist");
  
  string allCuts = baseCuts + " && fls3d>10";
  t->Draw("1e12*tau>>h2", allCuts.c_str(), "goff");
  h2->Draw("histsame");

  allCuts = baseCuts + " && fls3d>15";
  t->Draw("1e12*tau>>h3", allCuts.c_str(), "goff");
  h3->Draw("histsame");

  tl->SetTextSize(0.03);
  tl->SetTextColor(kBlack);
  tl->DrawLatex(0.2, 0.95, baseCuts.c_str());

  tl->SetTextSize(0.07);
  tl->SetTextColor(kBlue);
  tl->DrawLatex(0.6, 0.8, "l/#sigma(l_{3D}) > 10");

  tl->SetTextColor(kRed);
  tl->DrawLatex(0.6, 0.70, "l/#sigma(l_{3D}) > 15");

  hr->Divide(h2, h1, 1., 1., "b"); setHist(hr, kBlue); 
  hs->Divide(h3, h1, 1., 1., "b"); setHist(hs, kRed); 

  c0->cd(2);
  setTitles(hr, "#tau [ps]", "#epsilon'");
  hr->Draw();
  hs->Draw("same");

  c0->SaveAs(Form("%s/%s", fDirectory.c_str(), pdfname.c_str())); 
}


// ----------------------------------------------------------------------
void plotResults::fls3dVsX(std::string x, std::string cuts, std::string pdfname) {

  gStyle->SetOptStat(0);

  double xmax(10);
  string xt;
  if (string::npos != x.find("pt")) {
    xmax = 40;
    xt = "pT(B) [GeV]";
  }
  if (string::npos != x.find("eta")) {
    xmax = 2.5;
    xt = "|#eta(B)| ";
  }    

  zone(1);
  fF["SgMc"]->cd("candAnaMuMu");
  TProfile *h1 = new TProfile("h1", "", 100, 0., xmax); 
  h1->SetMinimum(0.);
  h1->SetMaximum(80.);

  TTree *t = (TTree*)gDirectory->Get("events"); 
  
  string baseCuts = "m1pt>4.5&&m2pt>4.2&&hlt&&" + cuts; 
  t->Draw(Form("fls3d:%s>>h1", x.c_str()), baseCuts.c_str(), "prof");
  h1->SetXTitle(xt.c_str());
  h1->SetYTitle("l_{3D}/#sigma(l_{3D})"); 
  h1->Draw();

  c0->SaveAs(Form("%s/%s", fDirectory.c_str(), pdfname.c_str())); 
}


// ----------------------------------------------------------------------
void plotResults::setupNorm() {

  fNumbersNo[0]->effTot = 0.0011;
  fNumbersNo[1]->effTot = 0.00032;
  fNumbersNo[0]->fitYield = 82712;
  fNumbersNo[1]->fitYield = 23809;

}
