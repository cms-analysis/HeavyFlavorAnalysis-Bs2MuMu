#include "plotEfficiencies.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "TMath.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"

using namespace std; 
using std::string; 

ClassImp(plotEfficiencies)


// ----------------------------------------------------------------------
plotEfficiencies::plotEfficiencies(const char *files, const char *dir, const char *cuts, int mode) : plotClass(files, dir, cuts, mode) { 

  fDoPrint = true; 

  fDoUseBDT = true; 
  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 

  fNumbersFileName = fDirectory + "/anaBmm.plotEfficiencies." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

}

// ----------------------------------------------------------------------
plotEfficiencies::~plotEfficiencies() {
  fHistFile->Write();
  fHistFile->Close();
}

// ----------------------------------------------------------------------
void plotEfficiencies::makeAll(int channel) {
  if (channel &1) {

    // -- dump histograms
    string hfname  = fDirectory + "/anaBmm.plotEfficiencies." + fSuffix + ".root";
    //  string hfname  = "anaBmm.plotResults." + fSuffix + ".root";
    cout << "fHistFile: " << hfname;
    fHistFile = TFile::Open(hfname.c_str(), "RECREATE");
    cout << " opened " << endl;

    fDoApplyMuonPtCuts = false; 
    tnpVsMC(fCuts[0]->m1pt, fCuts[0]->m2pt, "analysis");

    fDoApplyMuonPtCuts = true; 
    
    for (int i = 0; i < 6; ++i) {
      double pt = 4.0 + i*1; 
      cout << "tnpVsMC(" << pt << ", " << pt << ")" << endl;
      tnpVsMC(pt, pt, "study");
    }

  }

  if (channel &2) {

    // -- read histograms
    string hfname  = fDirectory + "/anaBmm.plotEfficiencies." + fSuffix + ".root";
    //  string hfname  = "anaBmm.plotResults." + fSuffix + ".root";
    cout << "fHistFile: " << hfname;
    fHistFile = TFile::Open(hfname.c_str());
    cout << " opened " << endl;

    texNumbers(fCuts[0]->m1pt, fCuts[0]->m2pt, "analysis"); 

    for (int i = 0; i < 6; ++i) {
      double pt = 4.0 + i*1; 
      cout << "tnpVsMC(" << pt << ", " << pt << ")" << endl;
      texNumbers(pt, pt, "woCowboyVeto");
    }
  }

  if (channel &4) {
    // -- MC trigger efficiencies
    mcTriggerEffs();
  }
}

// ----------------------------------------------------------------------
void plotEfficiencies::mcTriggerEffs() {

  fDoApplyCowboyVeto = true;   
  fDoApplyCowboyVetoAlsoInSignal = false; 
  
  fF["SgMc3e33"]->cd("candAnaMuMu");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:SgMc3e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:SgMc3e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:SgMc3e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:SgMc3e33:trigEff1:err", fSuffix.c_str()), 3) << endl;

  fF["SgMc2e33"]->cd("candAnaMuMu");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:SgMc2e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:SgMc2e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:SgMc2e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:SgMc2e33:trigEff1:err", fSuffix.c_str()), 3) << endl;

  fF["SgMc1e33"]->cd("candAnaMuMu");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:SgMc1e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:SgMc1e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:SgMc1e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:SgMc1e33:trigEff1:err", fSuffix.c_str()), 3) << endl;


  // -- Bu -> J/psi phi
  fF["NoMc3e33"]->cd("candAnaBu2JpsiK");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:NoMc3e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:NoMc3e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:NoMc3e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:NoMc3e33:trigEff1:err", fSuffix.c_str()), 3) << endl;

  fF["NoMc2e33"]->cd("candAnaBu2JpsiK");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:NoMc2e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:NoMc2e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:NoMc2e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:NoMc2e33:trigEff1:err", fSuffix.c_str()), 3) << endl;


  fF["NoMc1e33"]->cd("candAnaBu2JpsiK");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:NoMc1e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:NoMc1e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:NoMc1e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:NoMc1e33:trigEff1:err", fSuffix.c_str()), 3) << endl;


  // -- Bs -> J/psi phi
  fF["CsMc3e33"]->cd("candAnaBs2JpsiPhi");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:CsMc3e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:CsMc3e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:CsMc3e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:CsMc3e33:trigEff1:err", fSuffix.c_str()), 3) << endl;

  fF["CsMc2e33"]->cd("candAnaBs2JpsiPhi");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:CsMc2e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:CsMc2e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:CsMc2e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:CsMc2e33:trigEff1:err", fSuffix.c_str()), 3) << endl;


  fF["CsMc1e33"]->cd("candAnaBs2JpsiPhi");
  loopTree(98); 
  cout << "trigger eff 0: " << fNumbersBla[0]->effTrigMC << endl;
  cout << "trigger eff 1: " << fNumbersBla[1]->effTrigMC << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMC, Form("%s:CsMc1e33:trigEff0:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[0]->effTrigMCE, Form("%s:CsMc1e33:trigEff0:err", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMC, Form("%s:CsMc1e33:trigEff1:val", fSuffix.c_str()), 3) << endl;
  fTEX << formatTex(fNumbersBla[1]->effTrigMCE, Form("%s:CsMc1e33:trigEff1:err", fSuffix.c_str()), 3) << endl;

}


// ----------------------------------------------------------------------
void plotEfficiencies::resetHistograms() {
  
  for (unsigned int i = 0; i < fNchan; ++i) {

    fhMuId[i]->Reset();
    fhMuTr[i]->Reset();
    fhMuIdMC[i]->Reset();
    fhMuTrMC[i]->Reset();
    
    fhMassAbsNoCuts[i]->Reset();
    fhMassNoCuts[i]->Reset();
    fhMassNoCutsManyBins[i]->Reset();

    fhMassWithAnaCuts[i]->Reset();
    fhMassWithAnaCutsManyBins[i]->Reset();

    fhMassWithMuonCuts[i]->Reset();
    fhMassWithMuonCutsManyBins[i]->Reset();

    fhMassWithTriggerCuts[i]->Reset();
    fhMassWithTriggerCutsManyBins[i]->Reset();

    fhMassWithAllCuts[i]->Reset();
    fhMassWithAllCutsBlind[i]->Reset();
    fhMassWithAllCutsManyBins[i]->Reset();

    fhNorm[i]->Reset();
    fhNormC[i]->Reset();

    fhMassWithMassCutsManyBins[i]->Reset();
    fhMassWithMassCuts[i]->Reset();
  }

}


// ----------------------------------------------------------------------
void plotEfficiencies::loopFunction(int function, int mode) {
  if (1 == function) loopFunction1(mode);
}



// ----------------------------------------------------------------------
// FIXME: This function (together with the identical version in plotClass) should migrate into plotClass
void plotEfficiencies::loopFunction1(int mode) {

  if (fChan < 0) return;
  
  double mass = fb.m; 

  bool bp2jpsikp(false), bs2jpsiphi(false); 
  if (10 == mode) bp2jpsikp = true; 
  if (15 == mode) bp2jpsikp = true; 

  if (20 == mode) bs2jpsiphi = true; 
  if (25 == mode) bs2jpsiphi = true; 
  
  fhMassAbsNoCuts[fChan]->Fill(mass);

  if (!fGoodAcceptance) return;

  // -- this is the base, after the raw acceptance cuts
  fhMassNoCuts[fChan]->Fill(mass);
  fhMassNoCutsManyBins[fChan]->Fill(mass); 

  if (fDoUseBDT) {
    if (fBDT < fCuts[fChan]->bdt) return;
  } else {
    if (!fGoodQ) return;
    if (!fGoodJpsiCuts) return;
    if (!fGoodPvAveW8) return;
    if (!fGoodMaxDoca) return;
    if (!fGoodLip) return;
    if (!fGoodLipS) return;
    if (!fGoodIp) return;
    if (!fGoodIpS) return;
    if (!fGoodPt) return;
    if (!fGoodEta) return;
    if (!fGoodAlpha) return;
    if (!fGoodChi2) return;
    if (!fGoodFLS) return;
    if (!fGoodCloseTrack) return;
    if (!fGoodIso) return;
    if (!fGoodDocaTrk) return;
  }

  fhMassWithAnaCuts[fChan]->Fill(mass); 
  fhMassWithAnaCutsManyBins[fChan]->Fill(mass); 

  double tr1w8(0.), tr2w8(0.), trw8(0.), m1w8(0.), m2w8(0.), mw8(0.0);
  // -- muon ID: Data PidTables
  m1w8 = fptM->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.);
  m2w8 = fptM->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.);
  mw8  = m1w8*m2w8; 
  fhMuId[fChan]->Fill(mw8); 
  
  // -- muon ID: MC PidTables
  m1w8 = fptMMC->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.);
  m2w8 = fptMMC->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.);
  mw8  = m1w8*m2w8; 
  if (mw8 > 0.) {
    fhMuIdMC[fChan]->Fill(mw8, 1./mw8); 
  }

  // -- muon trigger: Data PidTables
  tr1w8 = fptT1->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.)*fptT2->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.);
  tr2w8 = fptT1->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.)*fptT2->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.);
  trw8  = tr1w8*tr2w8; 
  if (bs2jpsiphi || bp2jpsikp) {
    if (fb.rr >= 3) {
      if (TMath::Abs(fb.m1eta) > 2.2) trw8 = 0; 
      if (TMath::Abs(fb.m2eta) > 2.2) trw8 = 0; 
    }
  }
  fhMuTr[fChan]->Fill(trw8); 
  
  // -- muon trigger: MC PidTables
  tr1w8 = fptT1MC->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.)*fptT2->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.);
  tr2w8 = fptT1MC->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.)*fptT2->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.);
  trw8  = tr1w8*tr2w8; 
  if (bs2jpsiphi || bp2jpsikp) {
    if (fb.rr >= 3) {
      if (TMath::Abs(fb.m1eta) > 2.2) trw8 = 0; 
      if (TMath::Abs(fb.m2eta) > 2.2) trw8 = 0; 
    }
  }
  fhMuTrMC[fChan]->Fill(trw8); 


  // -- MUON ID
  if (false == fb.gmuid) return;
  fhMassWithMuonCuts[fChan]->Fill(mass); 
  fhMassWithMuonCutsManyBins[fChan]->Fill(mass); 
  
  // -- TRIGGER
  if (false == fb.hlt) return;
  fhMassWithTriggerCuts[fChan]->Fill(mass); 
  fhMassWithTriggerCutsManyBins[fChan]->Fill(mass); 
  
  fhMassWithAllCuts[fChan]->Fill(mass); 
  if (5 == mode && !(5.2 < mass && mass < 5.45)) {
    fhMassWithAllCutsBlind[fChan]->Fill(mass); 
  }
  
  fhMassWithAllCutsManyBins[fChan]->Fill(mass); 
  
  fhNorm[fChan]->Fill(mass);
  fhNormC[fChan]->Fill(fb.cm);
  fhDstarPi[fChan]->Fill(mass);
  
  //FIXME  if (!fDoUseBDT) elist->Enter(jentry); 
  
  if (0 == mode && mass < fCuts[fChan]->mBsLo) return;
  if (0 == mode && mass > fCuts[fChan]->mBsHi) return;
  if (1 == mode && mass < fCuts[fChan]->mBdLo) return;
  if (1 == mode && mass > fCuts[fChan]->mBdHi) return;
  if (10 == mode && mass < fNoLo) return;
  if (10 == mode && mass > fNoHi) return;
  if (20 == mode && mass < fCsLo) return;
  if (20 == mode && mass > fCsHi) return;
  
  fhMassWithMassCuts[fChan]->Fill(mass);
  fhMassWithMassCutsManyBins[fChan]->Fill(mass); 
  
  fhBdtMass[fChan]->Fill(mass, fBDT); 
  //  if (fBDT > 0.3) cout << "mass = " << mass << " bdt = " << fBDT << endl;
}



// ----------------------------------------------------------------------
void plotEfficiencies::tnpVsMC(double m1pt, double m2pt, string what) {

  if (string::npos != what.find("study")) {
    for (unsigned int i = 0; i < fNchan; ++i) {
      fCuts[i]->m1pt = m1pt; 
      fCuts[i]->m2pt = m2pt; 
    }
  }

  TTree *t(0); 
  resetHistograms();
  fSetup = "SgMc"; 
  t = getTree(fSetup); 
  setupTree(t, fSetup); 
  loopOverTree(t, fSetup, 1);
  saveHists(fSetup, m1pt, m2pt, what);


  resetHistograms();
  fSetup = "NoMc"; 
  t = getTree(fSetup); 
  setupTree(t, fSetup); 
  loopOverTree(t, fSetup, 1);
  saveHists(fSetup, m1pt, m2pt, what);

}


// ----------------------------------------------------------------------
void plotEfficiencies::texNumbers(double m1pt, double m2pt, string what) {

  fTEX << "% -- tnpVsMC for pt = " << m1pt << " and " << m2pt << " and what = " << what << endl;

  string Suffix = fSuffix + what; 
  // -- reset to have constant name for the default analysis setup
  if (what == "analysis") {
    m1pt = 1.1; 
    m2pt = 1.1; 
  }

  double r(0.), rs(1.), rp(1.); 
  for (int i = 0; i < fNchan; ++i) {
    numbersFromHist(i,  0, m1pt, m2pt, fNumbersBs[i]);
    numbersFromHist(i, 10, m1pt, m2pt, fNumbersNo[i]);

    // -------
    // -- muid
    // -------
    r = fNumbersBs[i]->effMuidMC;
    fTEX << formatTex(r, Form("%s:MuidSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effMuidMC;
    fTEX << formatTex(r, Form("%s:MuidNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:TNPMCMuidSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:TNPMCMuidNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effMuidTNP;
    fTEX << formatTex(r, Form("%s:TNPMuidSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effMuidTNP;
    fTEX << formatTex(r, Form("%s:TNPMuidNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    // -- rho
    r = fNumbersBs[i]->effMuidMC/fNumbersBs[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:rhoMuidSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    rs = r; 

    r = fNumbersNo[i]->effMuidMC/fNumbersNo[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:rhoMuidNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    rp = r; 

    // -- ratios
    r = fNumbersBs[i]->effMuidMC/fNumbersNo[i]->effMuidMC;
    fTEX << formatTex(r, Form("%s:rMcMuid%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effMuidMC, fNumbersBs[i]->effMuidMCE, fNumbersNo[i]->effMuidMC, fNumbersNo[i]->effMuidMCE); 
    fTEX << formatTex(r, Form("%s:rMcMuid%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effMuidTNPMC/fNumbersNo[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:rTNPMCMuid%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effMuidTNPMC, fNumbersBs[i]->effMuidTNPMCE, fNumbersNo[i]->effMuidTNPMC, fNumbersNo[i]->effMuidTNPMCE); 
    fTEX << formatTex(r, Form("%s:rTNPMCMuid%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effMuidTNP/fNumbersNo[i]->effMuidTNP;
    fTEX << formatTex(r, Form("%s:rTNPMuid%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effMuidTNP, fNumbersBs[i]->effMuidTNPE, fNumbersNo[i]->effMuidTNPMC, fNumbersNo[i]->effMuidTNPMCE); 
    fTEX << formatTex(r, Form("%s:rTNPMuid%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    // -- TNP: ratios*rho
    r = (fNumbersBs[i]->effMuidTNPMC/fNumbersNo[i]->effMuidTNPMC)*(rs/rp);
    fTEX << formatTex(r, Form("%s:rRhoTNPMCMuid%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effMuidTNPMC*rs, fNumbersBs[i]->effMuidTNPMCE, fNumbersNo[i]->effMuidTNPMC*rp, fNumbersNo[i]->effMuidTNPMCE); 
    fTEX << formatTex(r, Form("%s:rRhoTNPMCMuid%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = (fNumbersBs[i]->effMuidTNP/fNumbersNo[i]->effMuidTNP)*(rs/rp);
    fTEX << formatTex(r, Form("%s:rRhoTNPMuid%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effMuidTNP*rs, fNumbersBs[i]->effMuidTNPE, fNumbersNo[i]->effMuidTNPMC*rp, fNumbersNo[i]->effMuidTNPMCE); 
    fTEX << formatTex(r, Form("%s:rRhoTNPMuid%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;



    // -------
    // -- trig
    // -------
    r = fNumbersBs[i]->effTrigMC;
    fTEX << formatTex(r, Form("%s:TrigSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effTrigMC;
    fTEX << formatTex(r, Form("%s:TrigNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:TNPMCTrigSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:TNPMCTrigNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effTrigTNP;
    fTEX << formatTex(r, Form("%s:TNPTrigSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effTrigTNP;
    fTEX << formatTex(r, Form("%s:TNPTrigNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    // -- rho
    r = fNumbersBs[i]->effTrigMC/fNumbersBs[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:rhoTrigSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    rs = r; 

    r = fNumbersNo[i]->effTrigMC/fNumbersNo[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:rhoTrigNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    rp = r; 

    // -- ratios
    r = fNumbersBs[i]->effTrigMC/fNumbersNo[i]->effTrigMC;
    fTEX << formatTex(r, Form("%s:rMcTrig%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effTrigMC, fNumbersBs[i]->effTrigMCE, fNumbersNo[i]->effTrigMC, fNumbersNo[i]->effTrigMCE); 
    fTEX << formatTex(r, Form("%s:rMcTrig%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effTrigTNPMC/fNumbersNo[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:rTNPMCTrig%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effTrigTNPMC, fNumbersBs[i]->effTrigTNPMCE, fNumbersNo[i]->effTrigTNPMC, fNumbersNo[i]->effTrigTNPMCE); 
    fTEX << formatTex(r, Form("%s:rTNPMCTrig%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effTrigTNP/fNumbersNo[i]->effTrigTNP;
    fTEX << formatTex(r, Form("%s:rTNPTrig%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effTrigTNP, fNumbersBs[i]->effTrigTNPE, fNumbersNo[i]->effTrigTNPMC, fNumbersNo[i]->effTrigTNPMCE); 
    fTEX << formatTex(r, Form("%s:rTNPTrig%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    // -- TNP ratios*rho
    r = (fNumbersBs[i]->effTrigTNPMC/fNumbersNo[i]->effTrigTNPMC)*(rs/rp);
    fTEX << formatTex(r, Form("%s:rRhoTNPMCTrig%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effTrigTNPMC*rs, fNumbersBs[i]->effTrigTNPMCE, fNumbersNo[i]->effTrigTNPMC*rp, fNumbersNo[i]->effTrigTNPMCE); 
    fTEX << formatTex(r, Form("%s:rRhoTNPMCTrig%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = (fNumbersBs[i]->effTrigTNP/fNumbersNo[i]->effTrigTNP)*(rs/rp);
    fTEX << formatTex(r, Form("%s:rRhoTNPTrig%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effTrigTNP*rs, fNumbersBs[i]->effTrigTNPE, fNumbersNo[i]->effTrigTNPMC*rp, fNumbersNo[i]->effTrigTNPMCE); 
    fTEX << formatTex(r, Form("%s:rRhoTNPTrig%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

  }

  // -- muid 
  if (1) {
    cout << "            SG 0 Eff  TNP  1 Eff  TNP  NO 0 Eff  TNP  1 Eff  TNP  0 rEff rTNP 1 rEff rTNP  SG rho0 rho1  NO rho0 rho1 " << endl;
    cout << Form("pT > %3.1f muid:   ", m2pt) 
	 << Form("%3.2f %3.2f   %3.2f %3.2f", 
		 fNumbersBs[0]->effMuidMC, fNumbersBs[0]->effMuidTNPMC, fNumbersBs[1]->effMuidMC, fNumbersBs[1]->effMuidTNPMC)
	 << Form("      %3.2f %3.2f   %3.2f %3.2f", 
		 fNumbersNo[0]->effMuidMC, fNumbersNo[0]->effMuidTNPMC, fNumbersNo[1]->effMuidMC, fNumbersNo[1]->effMuidTNPMC)
	 << Form("   %3.2f %3.2f   %3.2f %3.2f", 
		 fNumbersBs[0]->effMuidMC/fNumbersNo[0]->effMuidMC, fNumbersBs[0]->effMuidTNPMC/fNumbersNo[0]->effMuidTNPMC, 
		 fNumbersBs[1]->effMuidMC/fNumbersNo[1]->effMuidMC, fNumbersBs[1]->effMuidTNPMC/fNumbersNo[1]->effMuidTNPMC)
	 << Form("     %3.2f %3.2f", 
		 fNumbersBs[0]->effMuidMC/fNumbersBs[0]->effMuidTNPMC, fNumbersBs[1]->effMuidMC/fNumbersBs[1]->effMuidTNPMC)
	 << Form("     %3.2f %3.2f", 
		 fNumbersNo[0]->effMuidMC/fNumbersNo[0]->effMuidTNPMC, fNumbersNo[1]->effMuidMC/fNumbersNo[1]->effMuidTNPMC)
	 << endl;
  }

  // -- trig
  if (1) {
    cout << Form("pT > %3.1f trig: ", m2pt) 
	 << Form("%3.2f %3.2f %3.2f %3.2f", 
		 fNumbersBs[0]->effTrigMC, fNumbersBs[0]->effTrigTNPMC, 
		 fNumbersBs[1]->effTrigMC, fNumbersBs[1]->effTrigTNPMC)
	 << " " 
	 << Form("%3.2f %3.2f %3.2f %3.2f", 
		 fNumbersNo[0]->effTrigMC, fNumbersNo[0]->effTrigTNPMC, 
		 fNumbersNo[1]->effTrigMC, fNumbersNo[1]->effTrigTNPMC)
	 << " -> " 
	 << Form(" rho_B = %3.2f rho_E = %3.2f", 
		 fNumbersNo[0]->effTrigMC/fNumbersNo[0]->effTrigTNPMC, 
		 fNumbersNo[1]->effTrigMC/fNumbersNo[1]->effTrigTNPMC)
	 << Form(" r_MC = %3.2f %3.2f r_TNP(MC) = %3.2f %3.2f  r_TNP(D) = %3.2f %3.2f", 
		 fNumbersBs[0]->effTrigMC/fNumbersNo[0]->effTrigMC, 
		 fNumbersBs[1]->effTrigMC/fNumbersNo[1]->effTrigMC, 
		 fNumbersBs[0]->effTrigTNPMC/fNumbersNo[0]->effTrigTNPMC,
		 fNumbersBs[1]->effTrigTNPMC/fNumbersNo[1]->effTrigTNPMC,
		 fNumbersBs[0]->effTrigTNP/fNumbersNo[0]->effTrigTNP,
		 fNumbersBs[1]->effTrigTNP/fNumbersNo[1]->effTrigTNP
		 )
	 << endl;
  }

}


// ----------------------------------------------------------------------
void plotEfficiencies::saveHists(string smode, double m1pt, double m2pt, string what) {


  if (what != "study") {
    m1pt = 1.1; 
    m2pt = 1.1; 
  }

  fHistFile->cd();

  int mode(0); 
  if ("SgMc" == smode)   mode = 0; 
  if ("BdMc" == smode)   mode = 1; 
  if ("SgData" == smode) mode = 5; 

  if ("NoMc" == smode)   mode = 10; 
  if ("NoData" == smode) mode = 15; 

  if ("CsMc" == smode)   mode = 20; 
  if ("CsData" == smode) mode = 25; 
 
  TH1D *h1(0); 
  TH2D *h2(0); 

  for (unsigned int i = 0; i < fNchan; ++i) {
    string modifier = (fDoUseBDT?"bdt":"cnc"); 
    modifier = Form("%s:%2.1f:%2.1f", modifier.c_str(), m1pt, m2pt); 
    h1 = (TH1D*)(fhGenAndAccNumbers[i]->Clone(Form("hGenAndAccNumbers_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMuId[i]->Clone(Form("hMuId_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMuIdMC[i]->Clone(Form("hMuIdMC_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMuTr[i]->Clone(Form("hMuTr_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMuTrMC[i]->Clone(Form("hMuTrMC_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassAbsNoCuts[i]->Clone(Form("hMassAbsNoCuts_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassAbsNoCuts_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassNoCuts[i]->Clone(Form("hMassNoCuts_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassNoCuts_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassNoCuts[i]->Clone(Form("hMassNoCutsManyBins_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassNoCutsManyBins_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassWithAnaCuts[i]->Clone(Form("hMassWithAnaCuts_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassWithAnaCuts_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassWithMuonCuts[i]->Clone(Form("hMassWithMuonCuts_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassWithMuonCuts_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassWithTriggerCuts[i]->Clone(Form("hMassWithTriggerCuts_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassWithTriggerCuts_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassWithAllCuts[i]->Clone(Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassWithAllCuts_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassWithAllCutsManyBins[i]->Clone(Form("hMassWithAllCutsManyBins_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassWithAllCutsManyBins_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassWithMassCuts[i]->Clone(Form("hMassWithMassCuts_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassWithMassCuts_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhMassWithMassCutsManyBins[i]->Clone(Form("hMassWithMassCutsManyBins_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetTitle(Form("hMassWithMassCutsManyBins_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    if (string::npos != fSetup.find("No") || string::npos != fSetup.find("Cs")) {
      h1 = (TH1D*)(fhNorm[i]->Clone(Form("hNorm_%s_%d_chan%d", modifier.c_str(), mode, i)));       
      h1->SetTitle(Form("hNorm_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
      h1->SetDirectory(fHistFile); 
      h1->Write();

      h1 = (TH1D*)(fhNormC[i]->Clone(Form("hNormC_%s_%d_chan%d", modifier.c_str(), mode, i)));     
      h1->SetTitle(Form("hNormC_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
      h1->SetDirectory(fHistFile); 
      h1->Write();
    }

    if (string::npos != fSetup.find("DstarPi")) {
      h1 = (TH1D*)(fhDstarPi[i]->Clone(Form("hDstarPi_%s_%d_chan%d", modifier.c_str(), mode, i))); 
      h1->SetTitle(Form("hDstarPi_%s_%d_chan%d %s", modifier.c_str(), mode, i, smode.c_str())); 
      h1->SetDirectory(fHistFile); 
      h1->Write();
    }
    
    if (fDoUseBDT) {
      h2 = (TH2D*)(fhBdtMass[i]->Clone(Form("hBdtMass_%s_%d_chan%d", modifier.c_str(), mode, i))); 
      h2->SetDirectory(fHistFile); 
      h2->Write();
    }

  }

}


// ----------------------------------------------------------------------
void plotEfficiencies::numbersFromHist(int chan, int mode, double m1pt, double m2pt, numbers *aa) {
  TH1D *h1(0); 

  // -- efficiency and acceptance
  string modifier = (fDoUseBDT?"bdt":"cnc");
  modifier = Form("%s:%2.1f:%2.1f", modifier.c_str(), m1pt, m2pt); 

  cout << Form("hMassNoCuts_%s_%d_chan%d", modifier.c_str(), mode, chan) << endl;
  TH1D *hMassNoCuts              = (TH1D*)fHistFile->Get(Form("hMassNoCuts_%s_%d_chan%d", modifier.c_str(), mode, chan));
  TH1D *hMassWithAnaCuts         = (TH1D*)fHistFile->Get(Form("hMassWithAnaCuts_%s_%d_chan%d", modifier.c_str(), mode, chan));
  TH1D *hMassWithMuonCuts        = (TH1D*)fHistFile->Get(Form("hMassWithMuonCuts_%s_%d_chan%d", modifier.c_str(), mode, chan));
  TH1D *hMassWithTriggerCuts     = (TH1D*)fHistFile->Get(Form("hMassWithTriggerCuts_%s_%d_chan%d", modifier.c_str(), mode, chan));
  TH1D *hMassWithAllCuts         = (TH1D*)fHistFile->Get(Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), mode, chan));
  TH1D *hMassWithMassCuts        = (TH1D*)fHistFile->Get(Form("hMassWithMassCuts_%s_%d_chan%d", modifier.c_str(), mode, chan));
  TH1D *hMassWithAllCutsManyBins = (TH1D*)fHistFile->Get(Form("hMassWithAllCutsManyBins_%s_%d_chan%d", modifier.c_str(), mode, chan));

  TH1D *hMuId                    = (TH1D*)fHistFile->Get(Form("hMuId_%s_%d_chan%d", modifier.c_str(), mode, chan));
  TH1D *hMuIdMC                  = (TH1D*)fHistFile->Get(Form("hMuIdMC_%s_%d_chan%d", modifier.c_str(), mode, chan));

  TH1D *hMuTr                    = (TH1D*)fHistFile->Get(Form("hMuTr_%s_%d_chan%d", modifier.c_str(), mode, chan));
  TH1D *hMuTrMC                  = (TH1D*)fHistFile->Get(Form("hMuTrMC_%s_%d_chan%d", modifier.c_str(), mode, chan));

  double a = hMassNoCuts->GetSumOfWeights(); 
  double b = hMassWithAnaCuts->GetSumOfWeights();
  double c = hMassWithMuonCuts->GetSumOfWeights();
  double d = hMassWithTriggerCuts->GetSumOfWeights();
  double e = hMassWithAllCuts->GetSumOfWeights();
  double f = hMassWithMassCuts->GetSumOfWeights();
  aa->ana0Yield        = a;
  aa->ana0YieldE       = TMath::Sqrt(a);
  aa->anaYield         = b; 
  aa->anaYieldE        = TMath::Sqrt(b); 
  aa->anaMuonYield     = c; 
  aa->anaMuonYieldE    = TMath::Sqrt(c); 
  aa->anaTriggerYield  = d; 
  aa->anaTriggerYieldE = TMath::Sqrt(d); 
  aa->anaWmcYield      = f; 
  aa->anaWmcYieldE     = TMath::Sqrt(f);
  aa->effAna           = b/a*aa->effPtReco;
  aa->effAnaE          = dEff(static_cast<int>(b), static_cast<int>(a)); // FIXME add error from effPtReco
  aa->effMuidMC        = c/b;
  aa->effMuidMCE       = dEff(static_cast<int>(c), static_cast<int>(b));
  aa->effMuidTNP       = hMuId->GetMean();
  aa->effMuidTNPE      = hMuId->GetMeanError();
  aa->effMuidTNPMC     = hMuIdMC->GetMean();
  aa->effMuidTNPMCE    = hMuIdMC->GetMeanError();
  aa->effTrigMC        = d/c;
  aa->effTrigMCE       = dEff(static_cast<int>(d), static_cast<int>(c));
  aa->effTrigTNP       = hMuTr->GetMean();
  aa->effTrigTNPE      = hMuTr->GetMeanError();
  aa->effTrigTNPMC     = hMuTrMC->GetMean();
  aa->effTrigTNPMCE    = hMuTrMC->GetMeanError();

  printNumbers(*aa, cout); 
    
}

// ----------------------------------------------------------------------
void plotEfficiencies::triggerSignal(string cuts) { 
  static int version(0); 

  string defaultCuts = "4.9<m&&m<5.9&&gmuid&&gtqual&&m1pt>4&&m2pt>4&&m1q*m2q<0&&abs(eta)<2.4&&abs(m2eta)<2.4&&pt>5&&abs(m1ip)<2&&abs(m2ip)<2";

  string allCuts = defaultCuts + "&&" + cuts;
  string hltCuts = defaultCuts + "&&" + cuts + "&&hlt";

  int NBINS(4); 
  double HLO(0.), HHI(40.); 
  TH1D *M0 = new TH1D("M0", "", 40, 4.8, 6.0); 
  TH1D *H0 = new TH1D("H0", "", NBINS, HLO, HHI); 
  TH1D *H1 = new TH1D("H1", "", NBINS, HLO, HHI); 
  TH1D *EF = new TH1D("EF", "", NBINS, HLO, HHI); EF->Sumw2();
  TH1D *RUNS = new TH1D("RUNS", "", 210, 160000, 181000); 
  
  int yNBINS(4); 
  double  yHLO(0), yHHI(2.4);
  TH1D *Y0 = new TH1D("Y0", "", yNBINS, yHLO, yHHI); 
  TH1D *Y1 = new TH1D("Y1", "", yNBINS, yHLO, yHHI); 
  TH1D *YR = new TH1D("YR", "", yNBINS, yHLO, yHHI); YR->Sumw2(); 

  TTree *t; 

  //  loopTree(0); 

  cout << "allCuts: " << allCuts << endl;
  cout << "hltCuts: " << hltCuts << endl;

  c0->Clear(); 
  c0->Divide(2,2); 
  gStyle->SetOptStat(0); 

  int i(1); 
  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    if (string::npos == imap->first.find("Bmt")) continue;
    cout << "===> " << imap->first;
    cout << ": " << fF[imap->first]->GetName();
    cout << " -> " << imap->second << endl;

    fF[imap->first]->cd("candAnaMuMu1313");
    TH1D *m0 = new TH1D("m0", "", 40, 4.8, 6.0); 
    TH1D *h0 = new TH1D("h0", "", NBINS, HLO, HHI); 
    TH1D *h1 = new TH1D("h1", "", NBINS, HLO, HHI); 
    TH1D *y0 = new TH1D("y0", "", yNBINS, yHLO, yHHI); 
    TH1D *y1 = new TH1D("y1", "", yNBINS, yHLO, yHHI); 
    TH1D *runs = new TH1D("runs", "", 210, 160000, 181000); 
    t = (TTree*)gDirectory->Get("events"); 
    
    t->Draw("m>>m0", allCuts.c_str(), "goff");
    double n0   = t->Draw("pt>>h0", allCuts.c_str(), "goff");
    double n1   = t->Draw("pt>>h1", hltCuts.c_str(), "goff");
    double n2   = t->Draw("run>>runs", hltCuts.c_str(), "goff");
    t->Draw("eta>>y0", allCuts.c_str(), "goff");
    t->Draw("eta>>y1", hltCuts.c_str(), "goff");
    n0 = h0->Integral(1, h0->GetNbinsX()); 
    n1 = h1->Integral(1, h1->GetNbinsX()); 
    double eff  = n1/n0; 
    cout << "n0 = " << n0 << " n1 = " << n1 << " eff = " << eff << endl;
    double deff = dEff(static_cast<int>(n1), static_cast<int>(n0)); 

    M0->Add(m0); 
    H0->Add(h0); 
    H1->Add(h1); 
    RUNS->Add(runs); 

    Y0->Add(y0); 
    Y1->Add(y1); 
    
    c0->cd(i++); 
    setFilledHist(h0, kBlack, kRed, 3004); 
    setTitles(h0, "p_{T} [GeV]", "Entries/bin"); 
    h0->SetMinimum(0.01);
    h0->Draw();

    setFilledHist(h1, kBlack, kBlue, 3005); 
    h1->Draw("same");

    tl->SetTextSize(0.04); tl->DrawLatex(0.2, 0.92, imap->first.c_str()); 
    tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.92, Form("eff = %4.3f #pm %4.3f", eff, deff)); 
    
  }

  c0->cd(i++); 
  setFilledHist(H0, kBlack, kRed, 3004); 
  setTitles(H0, "p_{T} [GeV]", "Entries/bin"); 
  H0->SetMinimum(0.01);
  H0->Draw();
  
  setFilledHist(H1, kBlack, kBlue, 3005); 
  H1->Draw("same");

  //   double n0   = H0->GetSumOfWeights(); 
  //   double n1   = H1->GetSumOfWeights(); 
  double n0   = H0->Integral(1, H0->GetNbinsX()); 
  double n1   = H1->Integral(1, H1->GetNbinsX()); 
  cout << "n0 = " << n0 << " n1 = " << n1 << endl;
  double eff  = n1/n0; 
  double deff = dEff(static_cast<int>(n1), static_cast<int>(n0)); 

  tl->SetTextSize(0.04); tl->DrawLatex(0.2, 0.92, "Combined:"); 
  tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.92, Form("eff = %4.3f #pm %4.3f", eff, deff)); 
  
  //   double ave(0.), aveE(0.); 
  //   average(ave, aveE, 3, vEff, vEffE); 

  if (fDoPrint) c0->SaveAs(Form("%s/triggerSignal-overlay-pt-%d.pdf", fDirectory.c_str(), version));     

  c0->Clear(); 
  RUNS->Draw();
  if (fDoPrint) c0->SaveAs(Form("%s/triggerSignal-runs-%d.pdf", fDirectory.c_str(), version));     

  c0->Clear();
  EF->Divide(H1, H0, 1., 1., "B");
  setTitles(EF, "p_{T} [GeV]", "efficiency"); 
  EF->SetMinimum(0.); EF->SetMaximum(1.05); 
  EF->Draw();
  tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.92, Form("#epsilon = %4.3f #pm %4.3f", eff, deff)); 
  if (fDoPrint) c0->SaveAs(Form("%s/triggerSignal-efficiency-pt-%d.pdf", fDirectory.c_str(), version));     

  c0->Clear();
  YR->Divide(Y1, Y0, 1., 1., "B");
  setTitles(YR, "|#eta|", "efficiency"); 
  YR->SetMinimum(0.); YR->SetMaximum(1.05);
  YR->Draw();
  tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.92, Form("#epsilon = %4.3f #pm %4.3f", eff, deff)); 
  if (fDoPrint) c0->SaveAs(Form("%s/triggerSignal-efficiency-eta-%d.pdf", fDirectory.c_str(), version));     


  ++version;

}



// ----------------------------------------------------------------------
void plotEfficiencies::triggerNormalization(string cuts) { 
  static int version(0); 

  string defaultCuts = "2.9<mpsi&&mpsi<3.3&&4.8<m&&m<6&&gmuid&&gtqual&&m1pt>4&&m2pt>4&&pt>7&&m1q*m2q<0&&abs(eta)<2.4&&abs(m2eta)<2.4&&pt>5&&abs(m1ip)<2&&abs(m2ip)<2";

  string allCuts = defaultCuts + "&&" + cuts;
  string hltCuts = defaultCuts + "&&" + cuts + "&&hlt";

  int NBINS(8); 
  double HLO(0.), HHI(40.); 
  TH1D *M0 = new TH1D("M0", "", 40, 4.8, 6.0); 
  TH1D *H0 = new TH1D("H0", "", NBINS, HLO, HHI); 
  TH1D *H1 = new TH1D("H1", "", NBINS, HLO, HHI); 
  TH1D *HR = new TH1D("HR", "", NBINS, HLO, HHI); HR->Sumw2();  
  TH1D *RUNS = new TH1D("RUNS", "", 210, 160000, 181000); 

  int yNBINS(6); 
  double  yHLO(0), yHHI(2.4);
  TH1D *Y0 = new TH1D("Y0", "", yNBINS, yHLO, yHHI); 
  TH1D *Y1 = new TH1D("Y1", "", yNBINS, yHLO, yHHI); 
  TH1D *YR = new TH1D("YR", "", yNBINS, yHLO, yHHI); YR->Sumw2(); 
  
  TTree *t; 

  //  loopTree(0); 

  c0->Clear(); 
  c0->Divide(2,2); 
  gStyle->SetOptStat(0); 

  cout << "allCuts: " << allCuts << endl;
  cout << "hltCuts: " << hltCuts << endl;

  int i(1); 
  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    if (string::npos == imap->first.find("Bmt")) continue;
    cout << "===> " << imap->first;
    cout << ": " << fF[imap->first]->GetName();
    cout << " -> " << imap->second << endl;

    fF[imap->first]->cd("candAnaBu2JpsiK");
    TH1D *m0 = new TH1D("m0", "", 40, 4.8, 6.0); 
    TH1D *h0 = new TH1D("h0", "", NBINS, HLO, HHI); 
    TH1D *h1 = new TH1D("h1", "", NBINS, HLO, HHI); 
    TH1D *y0 = new TH1D("y0", "", yNBINS, yHLO, yHHI); 
    TH1D *y1 = new TH1D("y1", "", yNBINS, yHLO, yHHI); 
    TH1D *runs = new TH1D("runs", "", 210, 160000, 181000); 
    t = (TTree*)gDirectory->Get("events"); 
    
    //    t->Draw("m>>m0", allCuts.c_str(), "goff");
    double n0   = t->Draw("pt>>h0", allCuts.c_str(), "goff");
    double n1   = t->Draw("pt>>h1", hltCuts.c_str(), "goff");
    double n2   = t->Draw("run>>runs", hltCuts.c_str(), "goff");
    n0 = h0->Integral(1, h0->GetNbinsX()); 
    n1 = h1->Integral(1, h1->GetNbinsX()); 
    t->Draw("eta>>y0", allCuts.c_str(), "goff");
    t->Draw("eta>>y1", hltCuts.c_str(), "goff");
    double eff  = n1/n0; 
    double deff = dEff(static_cast<int>(n1), static_cast<int>(n0)); 

    M0->Add(m0); 
    H0->Add(h0); 
    H1->Add(h1); 
    RUNS->Add(runs); 

    Y0->Add(y0); 
    Y1->Add(y1); 
    
    c0->cd(i++); 
    setFilledHist(h0, kBlack, kRed, 3004); 
    setTitles(h0, "p_{T} [GeV]", "Entries/bin"); 
    h0->Draw();

    setFilledHist(h1, kBlack, kBlue, 3005); 
    h1->Draw("same");

    tl->SetTextSize(0.04); tl->DrawLatex(0.2, 0.92, imap->first.c_str()); 
    tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.92, Form("eff = %4.3f #pm %4.3f", eff, deff)); 
    
  }

  c0->cd(i++); 
  setFilledHist(H0, kBlack, kRed, 3004); 
  setTitles(H0, "p_{T} [GeV]", "Entries/bin"); 
  H0->Draw();
  
  setFilledHist(H1, kBlack, kBlue, 3005); 
  H1->Draw("same");

//   double n0   = H0->GetSumOfWeights(); 
//   double n1   = H1->GetSumOfWeights(); 
  double n0 = H0->Integral(1, H0->GetNbinsX()); 
  double n1 = H1->Integral(1, H1->GetNbinsX()); 
  double eff  = n1/n0; 
  cout << "n0 = " << n0 << " n1 = " << n1 << endl;
  double deff = dEff(static_cast<int>(n1), static_cast<int>(n0)); 

  tl->SetTextSize(0.04); tl->DrawLatex(0.2, 0.92, "Combined:"); 
  tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.92, Form("eff = %4.3f #pm %4.3f", eff, deff)); 
  
//   double ave(0.), aveE(0.); 
//   average(ave, aveE, 3, vEff, vEffE); 

  if (fDoPrint) c0->SaveAs(Form("%s/triggerNormalization-overlay-pt-%d.pdf", fDirectory.c_str(), version));     

  c0->Clear(); 
  RUNS->Draw();
  if (fDoPrint) c0->SaveAs(Form("%s/triggerNormalization-runs-%d.pdf", fDirectory.c_str(), version));     

  c0->Clear();
  HR->Divide(H1, H0, 1., 1., "B");
  setTitles(HR, "p_{T} [GeV]", "efficiency"); 
  HR->SetMinimum(0.); HR->SetMaximum(1.05); 
  HR->Draw();
  tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.92, Form("#epsilon = %4.3f #pm %4.3f", eff, deff)); 
  if (fDoPrint) c0->SaveAs(Form("%s/triggerNormalization-efficiency-pt-%d.pdf", fDirectory.c_str(), version));     

  c0->Clear();
  YR->Divide(Y1, Y0, 1., 1., "B");
  setTitles(YR, "|#eta|", "efficiency"); 
  YR->SetMinimum(0.); YR->SetMaximum(1.05); 
  YR->Draw();
  tl->SetTextSize(0.04); tl->DrawLatex(0.5, 0.92, Form("#epsilon = %4.3f #pm %4.3f", eff, deff)); 
  if (fDoPrint) c0->SaveAs(Form("%s/triggerNormalization-efficiency-eta-%d.pdf", fDirectory.c_str(), version));     
  

  

  ++version;

}



                                                   
// L1L2Efficiency_VBTF_Track2_data_all_histo.root     
// L1L2Efficiency_VBTF_Track2_datalike_mc_histo.root  
// L1L2Efficiency_VBTF_Track7_data_all_histo.root     
// L1L2Efficiency_VBTF_Track7_datalike_mc_histo.root  
// L3Efficiency_VBTF_data_all_histo.root              
// L3Efficiency_VBTF_datalike_mc_histo.root 


// MuonID_VBTF_Track2_data_all_histo.root    
// MuonID_VBTF_Track2_datalike_mc_histo.root 
// MuonID_VBTF_Track7_data_all_histo.root    
// MuonID_VBTF_Track7_datalike_mc_histo.root 


// ----------------------------------------------------------------------
void plotEfficiencies::convertLucasHistograms() {

  PidTable aTemplate, bTemplate, cTemplate;

  // -- no seagull/cowboy difference here
  readFile("../macros/pidtables/111210/L3Efficiency_VBTF_data_all_histo.root", "L3"); 
  readFile("../macros/pidtables/111210/L3Efficiency_VBTF_datalike_mc_histo.root", "L3"); 

  vector<string> top; 
  top.push_back("sg"); 
  top.push_back("cb"); 

 
  for (unsigned int i = 0; i < top.size(); ++i) {
    aTemplate.clear(); 
    read2Files(aTemplate, 
	       "../macros/pidtables/111210/L1L2Efficiency_VBTF_Track2_data_all_histo.root", 
	       "../macros/pidtables/111210/L1L2Efficiency_VBTF_Track7_data_all_histo.root", 
	       Form("L1L2_%s", top[i].c_str())); 
    aTemplate.shiftPmax(24.9, 999.); 
    aTemplate.dumpToFile(Form("../macros/pidtables/111210/L1L2Efficiency_VBTF_data_all_histo_%s.dat", top[i].c_str())); 


    aTemplate.clear();
    read2Files(aTemplate, 
	       "../macros/pidtables/111210/L1L2Efficiency_VBTF_Track2_datalike_mc_histo.root", 
	       "../macros/pidtables/111210/L1L2Efficiency_VBTF_Track7_datalike_mc_histo.root", 
	       Form("L1L2_%s", top[i].c_str())); 
    aTemplate.shiftPmax(24.9, 999.); 
    aTemplate.dumpToFile(Form("../macros/pidtables/111210/L1L2Efficiency_VBTF_datalike_mc_histo_%s.dat", top[i].c_str())); 


    aTemplate.clear();
    read2Files(aTemplate, 
	       "../macros/pidtables/111210/MuonID_VBTF_Track2_data_all_histo.root", 
	       "../macros/pidtables/111210/MuonID_VBTF_Track7_data_all_histo.root", 
	       Form("tight_%s", top[i].c_str())); 
    aTemplate.shiftPmax(24.9, 999.); 
    aTemplate.dumpToFile(Form("../macros/pidtables/111210/MuonID_VBTF_data_all_histo_%s.dat", top[i].c_str())); 


    aTemplate.clear();
    read2Files(aTemplate, 
	       "../macros/pidtables/111210/MuonID_VBTF_Track2_datalike_mc_histo.root", 
	       "../macros/pidtables/111210/MuonID_VBTF_Track7_datalike_mc_histo.root", 
	       Form("tight_%s", top[i].c_str())); 
    aTemplate.shiftPmax(24.9, 999.); 
    aTemplate.dumpToFile(Form("../macros/pidtables/111210/MuonID_VBTF_datalike_mc_histo_%s.dat", top[i].c_str())); 
  }

  // -- create combined sg/cb PidTables
  vector<string> tables; 
  tables.push_back("../macros/pidtables/111210/MuonID_VBTF_data_all_histo");
  tables.push_back("../macros/pidtables/111210/MuonID_VBTF_datalike_mc_histo");
  tables.push_back("../macros/pidtables/111210/L1L2Efficiency_VBTF_data_all_histo");
  tables.push_back("../macros/pidtables/111210/L1L2Efficiency_VBTF_datalike_mc_histo");
  for (int i = 0; i < tables.size(); ++i) {
    aTemplate.clear();
    bTemplate.clear();
    aTemplate.readFromFile(Form("%s_cb.dat", tables[i].c_str()));
    bTemplate.readFromFile(Form("%s_sg.dat", tables[i].c_str()));
    aTemplate.combine(bTemplate); 
    aTemplate.dumpToFile(Form("%s.dat", tables[i].c_str()));
  }

}


// ----------------------------------------------------------------------
void plotEfficiencies::readFile(const char *fname, const char *hname) {

  TFile *f = TFile::Open(fname); 
  f->ls(); 

  TH2F *hsg = (TH2F*)(f->Get(Form("%s", hname)))->Clone("hsg");
  cout << "hsg: " << hsg << endl;

  TH2F *h1 = (TH2F*)(f->Get(Form("%s", hname))); 
  cout << "h1: " << h1 << endl;
  h1->Clear();
  h1->SetTitle("h1");
  
  cout << hsg->GetEntries() << endl;
  
  zone(2,2);
  hsg->DrawCopy("colz");

  PidTable a;
  a.readFromEffHist(gFile, "hsg"); 
    
  a.eff2d(h1);
  c0->cd(2);
  h1->DrawCopy("colz");
  
  h1->Add(hsg, -1.); 

  c0->cd(3);
  h1->DrawCopy("colz");
  
  a.printAll();

  TString pname(fname); 
  pname.ReplaceAll(".root", ".dat"); 
  a.dumpToFile(pname.Data()); 

}


// ----------------------------------------------------------------------
void plotEfficiencies::read2Files(PidTable &a, const char *f1name, const char *f2name, const char *hname) {

  TFile *f1 = TFile::Open(f1name); 
  f1->ls(); 
  TH2F *hsg1 = (TH2F*)(f1->Get(hname))->Clone("hsg1");
  cout << "hsg1: " << hsg1 << endl;

  TFile *f2 = TFile::Open(f2name); 
  f2->ls(); 
  TH2F *hsg2 = (TH2F*)(f2->Get(hname))->Clone("hsg2");
  cout << "hsg2: " << hsg2 << endl;
  
  TH2F *h1 = (TH2F*)(f1->Get(hname))->Clone("h1"); 
  cout << "h1: " << h1 << endl;
  h1->Clear();
  h1->SetTitle("h1");
  
  PidTable a1;
  a1.readFromEffHist(f1, "hsg1"); 
  a1.printAll();

  PidTable a2;
  a2.readFromEffHist(f2, "hsg2"); 
  a2.printAll();
  
  a.fillTable(a1);
  a.fillTable(a2);
  
}


