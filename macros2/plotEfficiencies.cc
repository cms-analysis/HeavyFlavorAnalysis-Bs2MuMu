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
plotEfficiencies::plotEfficiencies(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

  fDoPrint = true; 

  fNumbersFileName = fDirectory + "/anaBmm.plotEfficiencies." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  int NBINS = (fMassHi - fMassLo)/0.025;

  int HBINS(15); 
  double HLO(0.), HHI(45.); 

  TH1D *h; 
  numbers *a(0); 
  cout << "----> fNchan = " << fNchan << " fCuts.size() = " << fCuts.size() << endl;
  for (unsigned int i = 0; i < fNchan; ++i) {
    cout << "--> set up structures for channel #" << i << endl;
    cout << "    etaMin: " << fCuts[i]->etaMin << " etaMax: " << fCuts[i]->etaMax << endl;
    // -- signal Bs2MuMu
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("signal Bs2MuMu %i", i); 
    fNumbersBs.push_back(a); 

    // -- signal Bd2MuMu
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("signal Bd2MuMu %i", i); 
    fNumbersBd.push_back(a); 

    // --  normalization
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name = Form("normalization %i", i); 
    fNumbersNo.push_back(a); 

    // --  control sample
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("control sample %i", i); 
    fNumbersCs.push_back(a); 

    h = new TH1D(Form("hMassWithMassCuts%d", i), Form("hMassWithMassCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMassCuts.push_back(h); 
    h = new TH1D(Form("fhMassWithMassCutsManyBins%d", i), Form("fhMassWithMassCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithMassCutsManyBins.push_back(h); 
  
    h = new TH1D(Form("hMassWithAllCuts%d", i), Form("hMassWithAllCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCuts.push_back(h); 
    h = new TH1D(Form("hMassWithAllCutsBlind%d", i), Form("hMassWithAllCutsBlind%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCutsBlind.push_back(h); 
    h = new TH1D(Form("hMassWithAllCutsManyBins%d", i), Form("hMassWithAllCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithAllCutsManyBins.push_back(h); 

    h = new TH1D(Form("hNorm%d", i), Form("hNorm%d", i), 60, 5.0, 5.6);
    fhNorm.push_back(h); 


    h = new TH1D(Form("hMassWithTriggerCuts%d", i), Form("hMassWithTriggerCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithTriggerCuts.push_back(h); 
    h = new TH1D(Form("hMassWithTriggerCutsManyBins%d", i), Form("hMassWithTriggerCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithTriggerCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassWithMuonCuts%d", i), Form("hMassWithMuonCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMuonCuts.push_back(h); 
    h = new TH1D(Form("hMassWithMuonCutsManyBins%d", i), Form("hMassWithMuonCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithMuonCutsManyBins.push_back(h); 

    h = new TH1D(Form("fhMassWithAnaCuts%d", i), Form("hMassChan%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAnaCuts.push_back(h); 
    h = new TH1D(Form("fhMassWithAnaCutsManyBins%d", i), Form("hMassChan%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithAnaCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassNoCuts%d", i), Form("hMassNoCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassNoCuts.push_back(h); 
    h = new TH1D(Form("fhMassNoCutsManyBins%d", i), Form("hMassChan%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassNoCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassAbsNoCuts%d", i), Form("hMassAbsNoCuts%d", i), 100, 0, 10);
    fhMassAbsNoCuts.push_back(h); 

    h = new TH1D(Form("hMuId%d", i), Form("hMuId%d", i), 100, 0., 1.);
    fhMuId.push_back(h); 
    h = new TH1D(Form("hMuTr%d", i), Form("hMuTr%d", i), 100, 0., 1.);
    fhMuTr.push_back(h); 

    h = new TH1D(Form("hMuIdMC%d", i), Form("hMuIdMC%d", i), 100, 0., 1.);
    fhMuIdMC.push_back(h); 
    h = new TH1D(Form("hMuTrMC%d", i), Form("hMuTrMC%d", i), 100, 0., 1.);
    fhMuTrMC.push_back(h); 

    h = new TH1D(Form("h0TNPTrigger%d", i), Form("hTNPTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0TNPTrigger.push_back(h); 
    h = new TH1D(Form("h1TNPTrigger%d", i), Form("hTNPTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1TNPTrigger.push_back(h); 
    h = new TH1D(Form("h0TNPMuID%d", i), Form("hTNPMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0TNPMuID.push_back(h); 
    h = new TH1D(Form("h1TNPMuID%d", i), Form("hTNPMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1TNPMuID.push_back(h); 

    h = new TH1D(Form("h0TNPMCTrigger%d", i), Form("hTNPMCTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0TNPMCTrigger.push_back(h); 
    h = new TH1D(Form("h1TNPMCTrigger%d", i), Form("hTNPMCTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1TNPMCTrigger.push_back(h); 
    h = new TH1D(Form("h0TNPMCMuID%d", i), Form("hTNPMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0TNPMCMuID.push_back(h); 
    h = new TH1D(Form("h1TNPMCMuID%d", i), Form("hTNPMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1TNPMCMuID.push_back(h); 


    h = new TH1D(Form("h0MCTrigger%d", i), Form("hMCTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0MCTrigger.push_back(h); 
    h = new TH1D(Form("h1MCTrigger%d", i), Form("hMCTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1MCTrigger.push_back(h); 
    h = new TH1D(Form("h0MCMuID%d", i), Form("hMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0MCMuID.push_back(h); 
    h = new TH1D(Form("h1MCMuID%d", i), Form("hMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1MCMuID.push_back(h); 
  }


}

// ----------------------------------------------------------------------
plotEfficiencies::~plotEfficiencies() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotEfficiencies::makeAll(int channels) {
  // -- BMT results
  triggerSignal(); 
  triggerNormalization(); 

  // -- efciciency ratios
  fDoApplyCowboyVeto = true;   
  fDoApplyCowboyVetoAlsoInSignal = false;   
  for (int i = 0; i < 6; ++i) {
    double pt = 4.0 + i*1.0; 
    cout << "tnpVsMC(" << pt << ", " << pt << ")" << endl;
    tnpVsMC(pt, pt, 1, "wCowboyVeto");
  }

  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 
  for (int i = 0; i < 6; ++i) {
    double pt = 4.0 + i*1; 
    cout << "tnpVsMC(" << pt << ", " << pt << ")" << endl;
    tnpVsMC(pt, pt, 1, "woCowboyVeto");
  }

  // -- MC trigger efficiencies
  mcTriggerEffs();
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
void plotEfficiencies::tnpVsMC(double m1pt, double m2pt, int chan, string what) {
  
  for (unsigned int i = 0; i < fNchan; ++i) {
    fCuts[i]->m1pt = m1pt; 
    fCuts[i]->m2pt = m2pt; 
  }
  
  loopTree(0); 
  loopTree(10); 

  fTEX << "% -- tnpVsMC for pt = " << m1pt << " and " << m2pt << " and what = " << what << endl;

  string Suffix = fSuffix + what; 

  double r(0.); 
  for (int i = 0; i < fNchan; ++i) {
    // -- muid
    r = fNumbersBs[i]->effMuidMC;
    fTEX << formatTex(r, Form("%s:MuidSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effMuidMC;
    fTEX << formatTex(r, Form("%s:MuidNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:TNPMuidSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:TNPMuidNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;


    r = fNumbersBs[i]->effMuidMC/fNumbersBs[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:rhoMuidSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersNo[i]->effMuidMC/fNumbersNo[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:rhoMuidNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effMuidMC/fNumbersNo[i]->effMuidMC;
    fTEX << formatTex(r, Form("%s:rMcMuid%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effMuidMC, fNumbersBs[i]->effMuidMCE, fNumbersNo[i]->effMuidMC, fNumbersNo[i]->effMuidMCE); 
    fTEX << formatTex(r, Form("%s:rMcMuid%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effMuidTNPMC/fNumbersNo[i]->effMuidTNPMC;
    fTEX << formatTex(r, Form("%s:rTNPMuid%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effMuidTNPMC, fNumbersBs[i]->effMuidTNPMCE, fNumbersNo[i]->effMuidTNPMC, fNumbersNo[i]->effMuidTNPMCE); 
    fTEX << formatTex(r, Form("%s:rTNPMuid%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    // -- trig
    r = fNumbersBs[i]->effTrigMC;
    fTEX << formatTex(r, Form("%s:TrigSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effTrigMC;
    fTEX << formatTex(r, Form("%s:TrigNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:TNPTrigSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
    r = fNumbersNo[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:TNPTrigNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effTrigMC/fNumbersBs[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:rhoTrigSg%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersNo[i]->effTrigMC/fNumbersNo[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:rhoTrigNo%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effTrigMC/fNumbersNo[i]->effTrigMC;
    fTEX << formatTex(r, Form("%s:rMcTrig%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effTrigMC, fNumbersBs[i]->effTrigMCE, fNumbersNo[i]->effTrigMC, fNumbersNo[i]->effTrigMCE); 
    fTEX << formatTex(r, Form("%s:rMcTrig%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = fNumbersBs[i]->effTrigTNPMC/fNumbersNo[i]->effTrigTNPMC;
    fTEX << formatTex(r, Form("%s:rTNPTrig%i-pT%2.0fpT%2.0f:val", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;

    r = dRatio(fNumbersBs[i]->effTrigTNPMC, fNumbersBs[i]->effTrigTNPMCE, fNumbersNo[i]->effTrigTNPMC, fNumbersNo[i]->effTrigTNPMCE); 
    fTEX << formatTex(r, Form("%s:rTNPTrig%i-pT%2.0fpT%2.0f:err", Suffix.c_str(), i, 10*m1pt, 10*m2pt), 3) << endl;
  }

  // -- muid 
  if (chan&1) {
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
  if (chan&2) {
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


