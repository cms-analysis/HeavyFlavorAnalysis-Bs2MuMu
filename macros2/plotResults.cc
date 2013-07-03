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
plotResults::plotResults(const char *files, const char *dir, const char *cuts, int mode) : plotClass(files, dir, cuts, mode) { 

  fDoPrint = true; 

  cout << "==> plotResults files: " << files << " dir: " << dir << " cuts: " << cuts << endl;

  fNumbersFileName = fDirectory + "/anaBmm.plotResults." + fSuffix + ".tex";
  //  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  printCuts(cout); 

  fDoUseBDT = true; 
  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 
  fInvertedIso = false; 
  fNormProcessed = false; 
  fSaveSmallTree = false;

  //  fStampString = "CMS, 5fb^{-1}"; 

}

// ----------------------------------------------------------------------
plotResults::~plotResults() {
  cout << "plotResults dtor: " << fHistFile << endl;
}


// ----------------------------------------------------------------------
void plotResults::makeAll(int channels, int nevents) {

  if (fDoUseBDT) {
    fStampString = "BDT preliminary"; 
  } else {
    fStampString = "CNC preliminary"; 
  }
 
  zone(1);
  //DBX  fSaveLargerTree = true; //FIXME

  if (channels & 1) {
    fillAndSaveHistograms(nevents); 
  }

  //DBX  return; // FIXME

  if (channels & 2) {
    calculateNumbers(2);
    if (fDoUseBDT) acceptanceMaps(1); 
  }


  // -- dimuons unblinded
  if (channels & 4) {
    TTree *t(0);
    string modifier("bdt"); 
    resetHistograms();
    fSaveSmallTree = true; 
    fSetup = "SgData"; 
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1);

    zone(1,2);
    gStyle->SetOptStat(0); 
    gStyle->SetOptTitle(0); 
    TH1D *h1(0); 
    for (unsigned int i = 0; i < fNchan; ++i) {
      c0->cd(i+1); 
      h1 = (TH1D*)(fhMassWithAllCutsManyBins[i]->Clone(Form("hMassWithAllCutsManyBins_%s_5_chan%d", modifier.c_str(), i))); 
      h1->SetAxisRange(4.9, 5.9, "X"); 
      h1->Draw();
      double ymax0 = 3; 
      double yleg = ymax0 - 0.3;
      h1->SetMaximum(ymax0); 
      double b0  = h1->Integral(h1->FindBin(5.2), h1->FindBin(5.3)); 
      double bs  = h1->Integral(h1->FindBin(5.3), h1->FindBin(5.45)); 
    
      drawArrow(0.6, 1, yleg); 
      drawArrow(0.4, 2, yleg-0.35); 
      tl->DrawLatex(0.4, yleg, Form("obs: %3.0f", bs));
      tl->DrawLatex(0.4, yleg-0.35, Form("obs: %3.0f", b0));
      tl->DrawLatex(0.2, 0.92, (0 == i?"Barrel":"Endcap"));
    }
  }
  c0->SaveAs(Form("%s/%s-unblinding.pdf", fDirectory.c_str(), fSuffix.c_str())); 
}



// ----------------------------------------------------------------------
void plotResults::acceptanceMaps(int mode) {

  string hfname  = fDirectory + "/anaBmm.plotResults." + fSuffix + ".root";
  TFile *f = TFile::Open(hfname.c_str()); 

  TH1D *h1a(0), *h1p(0), *h1r(0); 
  TH2D *h2a(0), *h2p(0), *h2r(0); 

  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 

  for (int i = 0; i < fNchan; ++i) {
    h1a = (TH1D*)f->Get(Form("hAccEtaAll_bdt_2_chan%d", i)); 
    h1p = (TH1D*)f->Get(Form("hAccEtaPass_bdt_2_chan%d", i)); 
    h1r = (TH1D*)h1a->Clone("eta_ratio"); h1r->Sumw2(); h1r->Reset();
    h1r->Divide(h1p, h1a, 1., 1., "b"); 
    zone(); 
    setTitles(h1r, "|#eta(B)| ", "acceptance", 0.05, 1.1, 1.5); 
    h1r->Draw(); 
    c0->SaveAs(Form("%s/%s-acc-eta-%d.pdf", fDirectory.c_str(), fSuffix.c_str(), i)); 

    h1a = (TH1D*)f->Get(Form("hAccPtAll_bdt_2_chan%d", i)); 
    h1p = (TH1D*)f->Get(Form("hAccPtPass_bdt_2_chan%d", i)); 
    h1r = (TH1D*)h1a->Clone("pt_ratio"); h1r->Sumw2(); h1r->Reset();
    h1r->Divide(h1p, h1a, 1., 1., "b"); 
    zone(); 
    setTitles(h1r, "pT(B) [GeV]", "acceptance", 0.05, 1.1, 1.5); 
    h1r->Draw(); 
    c0->SaveAs(Form("%s/%s-acc-pt-%d.pdf", fDirectory.c_str(), fSuffix.c_str(), i)); 

    h2a = (TH2D*)f->Get(Form("hAccAll_bdt_2_chan%d", i)); 
    h2p = (TH2D*)f->Get(Form("hAccPass_bdt_2_chan%d", i)); 
    h2r = (TH2D*)h2a->Clone("etapt_ratio"); h2r->Sumw2(); h2r->Reset();
    h2r->Divide(h2p, h2a, 1., 1., "b"); 
    zone(); 
    setTitles(h2r, "|#eta(B)|", "pT(B) [GeV]", 0.05, 1.1, 1.5); 
    h2r->Draw("colz"); 
    c0->SaveAs(Form("%s/%s-acc-etapt-%d.pdf", fDirectory.c_str(), fSuffix.c_str(), i)); 

  }

}


// ----------------------------------------------------------------------
void plotResults::play1(int mode) {

  double bdtmin = 0.; 

  int NBINS = 40;
  
  TFile *f1 = TFile::Open("2012/small-SgData.root"); 
  TFile *f2 = TFile::Open("2012/small-BdMc.root"); 

  TTree *t1 = (TTree*)f1->Get("SgData_bdt"); 
  TTree *t2 = (TTree*)f2->Get("BdMc_bdt"); 

  TH2D *h1 = new TH2D("h1", "", NBINS, 4.9, 5.9, 40, 0., 0.4);   h1->SetLineWidth(2);
  TH2D *h2 = new TH2D("h2", "", NBINS, 4.9, 5.9, 40, 0., 0.4); 
  
  string cuts; 

  makeCanvas(1);
  zone(2, 1, c1);

  for (int i = 0; i < 2; ++i) {
    c1->cd(i+1); 
    if (0 == i) {
      cuts = "(abs(m1eta)<1.4&&abs(m2eta)<1.4)"; 
    } else {
      cuts = "!(abs(m1eta)<1.4&&abs(m2eta)<1.4)"; 
    }
    
    t1->Draw("bdt:m>>h1", Form("%s&&bdt>%f", cuts.c_str(), bdtmin), "goff");
    t2->Draw("bdt:m>>h2", Form("%s&&bdt>%f", cuts.c_str(), bdtmin), "goff");
    
    h2->DrawCopy("col");
    h1->DrawCopy("boxsame");

    tl->DrawLatex(0.2, 0.92, (0==i? "Barrel": "Endcap")); 
  }

  
  newLegend(0.60, 0.6, 0.80, 0.85); 
  TLegend *legg = new TLegend(0.55, 0.6, 0.80, 0.85); 
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 
  legg->AddEntry(h1, "Data", "f"); 
  legg->AddEntry(h2, "B0 MC", "f"); 
  legg->Draw();


  c1->SaveAs(Form("%s/play1.pdf", fDirectory.c_str())); 
}



// ----------------------------------------------------------------------
void plotResults::play2(int chan) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  string hfname  = fDirectory + "/anaBmm.plotResults." + fSuffix + ".root";
  TFile *f = TFile::Open(hfname.c_str()); 
  
  // -- data
  TH1D *hd = (TH1D*)f->Get(Form("hMassWithAllCuts_bdt_5_chan%d", chan)); 

  // -- peaking backgrounds
  TH1D *hm1 = (TH1D*)f->Get(Form("hMassWithAllCuts_bdt_bgBd2KPi_chan%d", chan));  hm1->SetLineColor(41);  hm1->SetLineWidth(2);  
  TH1D *hm2 = (TH1D*)f->Get(Form("hMassWithAllCuts_bdt_bgBd2PiPi_chan%d", chan)); hm2->SetLineColor(42);  hm2->SetLineWidth(2); 
  TH1D *hm3 = (TH1D*)f->Get(Form("hMassWithAllCuts_bdt_bgBs2KK_chan%d", chan));   hm3->SetLineColor(30);  hm3->SetLineWidth(2); 

  TH1D *hs = (TH1D*)f->Get(Form("hMassWithAllCuts_bdt_1_chan%d", chan));   hs->SetLineColor(kBlue); 

  hd->SetLineWidth(2); 
  hd->SetAxisRange(4.9, 5.9, "X"); 
  hd->Draw();

  double peak(20.), base(1.36);
  if (1 == chan) {
    base = 0.77; 
    peak = 15.;
  }

  hm1->Scale(peak/hm1->GetSumOfWeights());
  hm2->Scale(peak/hm2->GetSumOfWeights()); 
  hm3->Scale(peak/hm3->GetSumOfWeights()); 

  hs->SetLineWidth(2); 
  hs->Scale(peak/hs->GetSumOfWeights()); 

  for (int i = 1; i < hm1->GetNbinsX(); ++i) {
    double x = hm1->GetBinCenter(i); 
    hm1->Fill(x, base);
    hm2->Fill(x, base);
    hm3->Fill(x, base);
    hs->Fill(x, base);
  }
  
//   hm1->Draw("same"); 
//   hm2->Draw("same"); 
//   hm3->Draw("same"); 
//   hs->Draw("same"); 



  zone(1,4);
  hd->Draw("e");
  hs->Draw("same"); 
  tl->SetTextSize(0.13); 
  tl->DrawLatex(0.6, 0.7, "B^{0} #rightarrow #mu #mu");

  c0->cd(2);
  hd->Draw("e");
  hm1->Draw("same"); 
  tl->DrawLatex(0.6, 0.7, "B^{0} #rightarrow K #pi");

  c0->cd(3);
  hd->Draw("e");
  hm2->Draw("same"); 
  tl->DrawLatex(0.6, 0.7, "B^{0} #rightarrow #pi #pi");

  c0->cd(4);
  hd->Draw("e");
  hm3->Draw("same"); 
  tl->DrawLatex(0.6, 0.7, "B^{0}_{s} #rightarrow K K");

  c0->SaveAs(Form("%s/%s-data-withBgOverlayed-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan)); 
  
}



// ----------------------------------------------------------------------
void plotResults::play3(int mode) {


}



// ----------------------------------------------------------------------
void plotResults::play4(int mode) {


}



// ----------------------------------------------------------------------
void plotResults::play5(int mode) {


}

// ----------------------------------------------------------------------
void plotResults::calculateNumbers(int mode) {

  // -- dump histograms
  string hfname  = fDirectory + "/anaBmm.plotResults." + fSuffix + ".root";
  //  string hfname  = "anaBmm.plotResults." + fSuffix + ".root";
  cout << "fHistFile: " << hfname;
  fHistFile = TFile::Open(hfname.c_str());
  cout << " opened " << endl;
  
  for (unsigned int chan = 0; chan < fNchan; ++chan) {
    fChan = chan; 
    initNumbers(fNumbersCs[chan], false); 
    initNumbers(fNumbersNo[chan], false); 
    initNumbers(fNumbersBs[chan], false); 
    initNumbers(fNumbersBd[chan], false); 
    calculateNoNumbers(chan, mode);
    calculateCsNumbers(chan, mode);
    calculateRareBgNumbers(chan);
    calculateSgNumbers(chan);
  }

  string bla = fDirectory + "/anaBmm.plotResults." + fSuffix + ".ulc";

  cout << "===> Storing ULC numbers in file " << bla << endl;
  system(Form("/bin/rm -f %s", bla.c_str()));
  printUlcalcNumbers(bla);
  printCsBFNumbers();

  //??  createAllCfgFiles(bla); 


  cout << "printing fNumbersBs" << endl;
  printNumbers(*fNumbersBs[0], cout); 
  printNumbers(*fNumbersBs[1], cout); 
  printNumbers(*fNumbersBs[0], fOUT); 
  printNumbers(*fNumbersBs[1], fOUT); 

  cout << "printing fNumbersBd" << endl;
  printNumbers(*fNumbersBd[0], cout); 
  printNumbers(*fNumbersBd[1], cout); 
  printNumbers(*fNumbersBd[0], fOUT); 
  printNumbers(*fNumbersBd[1], fOUT); 

  cout << "printing fNumbersNorm" << endl;
  printNumbers(*fNumbersNo[0], cout); 
  printNumbers(*fNumbersNo[1], cout); 
  printNumbers(*fNumbersNo[0], fOUT); 
  printNumbers(*fNumbersNo[1], fOUT); 

  fHistFile->Close();

}
    
// ----------------------------------------------------------------------
void plotResults::calculateNoNumbers(int chan, int mode) {

  c0->Clear();

  // -- efficiency and acceptance
  fSetup = "NoMc"; 
  numbersFromHist(chan, 10, fNumbersNo[chan]); 

  TH1D *h1(0); 
  // -- default mass
  string modifier = fDoUseBDT?"bdt":"cnc"; 
  string  name = Form("hNorm_%s_15_chan%d", modifier.c_str(), chan);
  h1 =  (TH1D*)fHistFile->Get(name.c_str()); 
  h1->DrawCopy(); 
  c0->cd(2);
  if (1 == mode) {
    normYield(h1, chan, 5.0); 
  } else if (2 == mode) {
    normYield2(h1, chan, 5.0, 5.6); 
  }
  fNumbersNo[chan]->fitYield  = fNoSig; 
  fNumbersNo[chan]->fitYieldE = fNoSigE; 

  // -- new numbers
  fNumbersNo[chan]->fFitNo = fNoSig; 
  fNumbersNo[chan]->fFitNoE1 = fNoSigE; 
  fNumbersNo[chan]->fFitNoE2 = 0.05*fNoSig; 

  // -- J/psi mass constrained
  name = Form("hNormC_%s_15_chan%d", modifier.c_str(), chan);
  h1 =  (TH1D*)fHistFile->Get(name.c_str()); 
  c0->cd(3);
  h1->DrawCopy(); 
  c0->cd(4);
  if (1 == mode) {
    normYield(h1, chan, 5.0); 
  } else if (2 == mode) {
    normYield2(h1, chan, 5.0); 
  }
  fNumbersNo[chan]->fitYieldC  = fNoSig; 
  fNumbersNo[chan]->fitYieldCE = fNoSigE; 

  // -- new numbers
  fNumbersNo[chan]->fFitNoC = fNoSig; 
  fNumbersNo[chan]->fFitNoCE1 = fNoSigE; 
  fNumbersNo[chan]->fFitNoCE2 = 0.05*fNoSig; 

  computeErrors(fNumbersNo); 

   
}


// ----------------------------------------------------------------------
void plotResults::calculateCsNumbers(int i, int mode) {

  // -- efficiency and acceptance
  fSetup = "CsMc"; 
  numbersFromHist(i, 20, fNumbersCs[i]); 

  TH1D *h1(0); 
  // -- default mass
  string modifier = fDoUseBDT?"bdt":"cnc"; 
  string  name = Form("hNorm_%s_25_chan%d", modifier.c_str(), i);
  h1 =  (TH1D*)fHistFile->Get(name.c_str()); 
  h1->DrawCopy(); 
  c0->cd(2);
  if (1 == mode) {
    csYield(h1, i, 5.0); 
  } else if (2 == mode) {
    csYield2(h1, i, 5.0); 
  }
  fNumbersCs[i]->fitYield  = fCsSig; 
  fNumbersCs[i]->fitYieldE = fCsSigE; 

  // -- new numbers
  fNumbersCs[i]->fFitCs   = fCsSig; 
  fNumbersCs[i]->fFitCsE1 = fCsSigE; 
  fNumbersCs[i]->fFitCsE2 = 0.07*fCsSig; 

  // -- J/psi mass constrained
  name = Form("hNormC_%s_25_chan%d", modifier.c_str(), i);
  h1 =  (TH1D*)fHistFile->Get(name.c_str()); 
  c0->cd(3);
  h1->DrawCopy(); 
  c0->cd(4);
  if (1 == mode) {
    csYield(h1, i, 5.0); 
  } else if (2 == mode) {
    csYield2(h1, i, 5.0); 
  }
  fNumbersCs[i]->fitYieldC  = fCsSig; 
  fNumbersCs[i]->fitYieldCE = fCsSigE; 

  // -- new numbers
  fNumbersCs[i]->fFitCsC   = fCsSig; 
  fNumbersCs[i]->fFitCsCE1 = fCsSigE; 
  fNumbersCs[i]->fFitCsCE2 = 0.07*fCsSig; 

  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;    
  fTEX << "% -- Control sample branching fraction " << endl;
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;    

  string cache = fSuffix; 
  fSuffix = (fDoUseBDT? "bdt":"cnc") + fSuffix; 
  
  double result, resultE; 
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
  
  // -- for adding some systematical errors
  double resultE2 = resultE*resultE; 
  resultE2 += 0.05*0.05*result*result + 0.05*0.05*result*result;
  double resultEsys = TMath::Sqrt(resultE2); 
  
  cout << "chan " << i << ": MC fact branching fraction: " << result << "+/-" << resultE << endl;
  fTEX << formatTex(result, Form("%s:N-CSBF-MC-BS%i:val", fSuffix.c_str(), i), 6) << endl;
  fTEX << formatTex(resultE, Form("%s:N-CSBF-MC-BS%i:err", fSuffix.c_str(), i), 6) << endl;
  fTEX << formatTex(resultEsys, Form("%s:N-CSBF-MC-BS%i:syst", fSuffix.c_str(), i), 6) << endl;
  
  result = (fNumbersCs[i]->fitYield/fNumbersNo[i]->fitYield)
    *(fu/fs)
    *(fNumbersNo[i]->effTot/fNumbersCs[i]->effTot)
    * fBF["NoMc"];
  
  cout << "chan " << i << ": branching fraction: " << result << "+/-" << resultE << endl;
  
  fTEX << formatTex(result, Form("%s:N-CSBF-BS%i:val", fSuffix.c_str(), i), 6) << endl;
  fTEX << formatTex(resultE, Form("%s:N-CSBF-BS%i:err", fSuffix.c_str(), i), 6) << endl;
  fTEX << formatTex(resultEsys, Form("%s:N-CSBF-BS%i:syst", fSuffix.c_str(), i), 6) << endl;
  

  fSuffix = cache; 


}



// ----------------------------------------------------------------------
void plotResults::calculateRareBgNumbers(int chan) {
  
  fSetup = "SgMc"; 
  numbersFromHist(chan, 0, fNumbersBs[chan]); 
  numbersFromHist(chan, 1, fNumbersBd[chan]); 

  string cache = fSuffix; 
  fSuffix = (fDoUseBDT? "bdt":"cnc") + fSuffix; 
  string modifier = (fDoUseBDT?"bdt":"cnc"); 

  TH1D *h1(0), *h2(0); 
  string name(""); 

//   // -- 'tight muon' values
//   double epsMu[]  = {0.83, 0.90}; // this is the square of the muon ID efficiency
//   double epsPi[]  = {0.001, 0.001}; 
//   double errPi2[] = {0.40*0.40, 0.40*0.40}; // changed to 40% uncertainty to be consistent with Michael's central value
//   double epsKa[]  = {0.001, 0.001}; 
//   double errKa2[] = {0.40*0.40, 0.40*0.40}; 
//   double epsPr[]  = {0.00015, 0.00015};
//   double errPr2[] = {1.00*1.00, 1.00*1.00};

  // -- mean 'MVA' values
  double epsMu[]  = {fNumbersBs[0]->effMuidMC, fNumbersBs[1]->effMuidMC}; // this is the square of the muon ID efficiency
  double epsPi[]  = {0.0009, 0.0009}; 
  double epsPi2[] = {epsPi[0]*epsPi[0], epsPi[1]*epsPi[1]}; 
  double errPi[]  = {0.50, 0.50};          
  double errPi2[] = {errPi[0]*errPi[0], errPi[1]*errPi[1]};
  double epsKa[]  = {0.0015, 0.0015}; 
  double epsKa2[] = {epsKa[0]*epsKa[0], epsKa[1]*epsKa[1]}; 
  double errKa[]  = {0.50, 0.50}; 
  double errKa2[] = {errKa[0]*errKa[0], errKa[1]*errKa[1]};
  double epsPr[]  = {0.00061, 0.00061};
  double epsPr2[] = {epsPr[0]*epsPr[0], epsPr[1]*epsPr[1]}; 
  double errPr[]  = {0.50, 0.50};
  double errPr2[] = {errKa[0]*errKa[0], errKa[1]*errKa[1]};

  double teff[] = {fNumbersBs[0]->effTrigMC, fNumbersBs[1]->effTrigMC};

  c0->Clear();
  gStyle->SetOptStat(0);

  std::map<string, int> colors, hatches;
  std::map<string, double> mscale;  
  std::map<string, double> err;  

  string cname; 
  
  cname = "bgLb2KP";
  colors.insert(make_pair(cname, 46)); hatches.insert(make_pair(cname, 3004)); mscale.insert(make_pair(cname, epsKa[chan]*epsPr[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + errPi2[chan] + errPr2[chan]))); 

  cname = "bgLb2PiP";
  colors.insert(make_pair(cname, 49)); hatches.insert(make_pair(cname, 3005)); mscale.insert(make_pair(cname, epsPi[chan]*epsPr[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + errKa2[chan] + errPr2[chan]))); 

  cname = "bgLb2PMuNu";
  colors.insert(make_pair(cname, 48)); hatches.insert(make_pair(cname, 3006)); mscale.insert(make_pair(cname, epsPr[chan]*epsMu[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + errPr2[chan]))); 


  cname = "bgBs2KK";
  colors.insert(make_pair(cname, 30)); hatches.insert(make_pair(cname, 3004)); mscale.insert(make_pair(cname, epsKa[chan]*epsKa[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + 1.))); //  + errKa2[chan] + errKa2[chan]))); 

  cname = "bgBs2KPi";
  colors.insert(make_pair(cname, 32)); hatches.insert(make_pair(cname, 3005)); mscale.insert(make_pair(cname, epsPi[chan]*epsKa[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + errPi2[chan] + errKa2[chan]))); 

  cname = "bgBs2PiPi"; 
  colors.insert(make_pair(cname, 33)); hatches.insert(make_pair(cname, 3007)); mscale.insert(make_pair(cname, epsPi[chan]*epsPi[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + 1.))); //  + errPi2[chan] + errPi2[chan]))); 

  cname = "bgBs2KMuNu";
  colors.insert(make_pair(cname, 34)); hatches.insert(make_pair(cname, 3008)); mscale.insert(make_pair(cname, epsKa[chan]*epsMu[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + errKa2[chan]))); 

  cname = "bgBs2MuMuGamma";
  colors.insert(make_pair(cname, 36));hatches.insert(make_pair(cname, 3009));mscale.insert(make_pair(cname, epsMu[chan]*epsMu[chan]));
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname]))); 


  cname = "bgBd2KK";
  colors.insert(make_pair(cname, 40)); hatches.insert(make_pair(cname, 3004)); mscale.insert(make_pair(cname, epsKa[chan]*epsKa[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + 1.))); //  + errKa2[chan] + errKa2[chan]))); 

  cname = "bgBd2KPi";
  colors.insert(make_pair(cname, 41)); hatches.insert(make_pair(cname, 3005)); mscale.insert(make_pair(cname, epsKa[chan]*epsPi[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + errKa2[chan] + errPi2[chan]))); 

  cname = "bgBd2PiPi";
  colors.insert(make_pair(cname, 42)); hatches.insert(make_pair(cname, 3007)); mscale.insert(make_pair(cname, epsPi[chan]*epsPi[chan])); 
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + 1.))); //  + errPi2[chan] + errPi2[chan]))); 

  cname = "bgBd2PiMuNu";
  colors.insert(make_pair(cname, 43));hatches.insert(make_pair(cname, 3008));mscale.insert(make_pair(cname, epsPi[chan]*epsMu[chan]));
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname] + errPi2[chan]))); 

  cname = "bgBd2Pi0MuMu";
  colors.insert(make_pair(cname, 44));hatches.insert(make_pair(cname, 3003));mscale.insert(make_pair(cname, epsMu[chan]*epsMu[chan]));
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname]))); 

  cname = "bgBd2MuMuGamma";
  colors.insert(make_pair(cname, 45));hatches.insert(make_pair(cname, 3002));mscale.insert(make_pair(cname, epsMu[chan]*epsMu[chan]));
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname]))); 


  cname = "bgBu2PiMuMu";
  colors.insert(make_pair(cname, 50));hatches.insert(make_pair(cname, 3004));mscale.insert(make_pair(cname, epsMu[chan]*epsMu[chan]));
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname]))); 

  cname = "bgBu2KMuMu";
  colors.insert(make_pair(cname, 52));hatches.insert(make_pair(cname, 3005));mscale.insert(make_pair(cname, epsMu[chan]*epsMu[chan]));
  err.insert(make_pair(cname, TMath::Sqrt(fBFE[cname]*fBFE[cname]))); 


  newLegend(0.55, 0.3, 0.80, 0.85); 
  TLegend *leggsl = new TLegend(0.55, 0.3, 0.80, 0.85); 
  leggsl->SetFillStyle(0); 
  leggsl->SetBorderSize(0); 
  leggsl->SetTextSize(0.04);  
  leggsl->SetFillColor(0); 
  leggsl->SetTextFont(42); 

  double valInc(0.), error(0.); 
  double vRare[4], vRareSl[4], vRareE[4], vRareSlE[4]; 
  double vRareHh[4], vRarePiPi[4], vRareKPi[4], vRareKK[4], vRarePiP[4], vRareKP[4]; 
  double vRareHhE[4], vRarePiPiE[4], vRareKPiE[4], vRareKKE[4], vRarePiPE[4], vRareKPE[4]; 
  for (int i = 0; i < 4; ++i) {
    vRare[i] = vRareSl[i] = vRareE[i] = vRareSlE[i] = 0.;

    vRareHh[i] = vRarePiPi[i] = vRareKPi[i] = vRareKK[i] = vRarePiP[i] = vRareKP[i] 
      = vRareHhE[i] = vRarePiPiE[i] = vRareKPiE[i] = vRareKKE[i] = vRarePiPE[i] = vRareKPE[i] = 0.;
  }

  double rareBs    = 0.;
  double rareBsE   = 0.;
  double rareBd    = 0.;
  double rareBdE   = 0.;

  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << "% --- rare Background per-channel numbers" << endl;  

  zone(1,2);

  TH1D *hMassWithAllCuts = (TH1D*)fHistFile->Get(Form("hMassWithAllCuts_%s_5_chan%d", modifier.c_str(),  chan)); 
  TH1D *bRare   = (TH1D*)(hMassWithAllCuts->Clone(Form("bRare_%s", modifier.c_str())));  bRare->SetLineColor(kBlack); bRare->Reset();
  TH1D *bslRare = (TH1D*)(hMassWithAllCuts->Clone(Form("bslRare_%s", modifier.c_str()))); bslRare->SetLineColor(kBlack); bslRare->Reset();

  THStack *hRareBg   = new THStack("hRareBg","");
  THStack *hslRareBg = new THStack("hslRareBg","");


  TH1D *hRare, *hRareNw8, *hRareW8; 
  double tot, bd, bs, efftot, pss, pdd;
  double eps(0.00001);
  int mbin(0); 
  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    fRareName = imap->first; 
    if (string::npos == fRareName.find("bg")) continue;
    if (0 == fF[fRareName])   continue; 


    double misid = mscale[fRareName];
    double ngenfile = ((TH1D*)fF[fRareName]->Get("monEvents"))->GetBinContent(1); 
    cout << "======> " << fRareName << " with ngenfile = " << ngenfile
	 << " filterEff = " << fFilterEff[fRareName] 
	 << " with old constant misid: " << misid << endl;

    name = Form("hMassWithAllCuts_%s_%s_chan%d", modifier.c_str(), fRareName.c_str(), chan);
    h1 = (TH1D*)fHistFile->Get(name.c_str()); 
    name = Form("hMassWithAllCutsManyBins_%s_%s_chan%d", modifier.c_str(), fRareName.c_str(), chan);
    h2 = (TH1D*)fHistFile->Get(name.c_str()); 
    hRareNw8 = h1; 

    tot  = h2->Integral(0, h2->GetNbinsX()+1); 
    bd   = h2->Integral(h2->FindBin(fCuts[chan]->mBdLo), h2->FindBin(fCuts[chan]->mBdHi));
    bs   = h2->Integral(h2->FindBin(fCuts[chan]->mBsLo), h2->FindBin(fCuts[chan]->mBsHi));

    efftot = tot/static_cast<double>(ngenfile)*fFilterEff[fRareName];

    cout << "NNNNNN%%%%%NNNNNNNNNN input: No = " << fNumbersNo[chan]->fitYield 
      //	 << " h1 : " << h1->GetName() << " " << h1->GetSumOfWeights()  
	 << " rare : " << hRareNw8->GetName() << " " << hRareNw8->GetSumOfWeights() 
	 << " tot: " << tot << " in B0 box: " << bd << " in Bs box: " << bs
	 << endl;

   
    if (tot>0) pss = bs/tot; 
    else pss = 0.;
    if (tot>0) pdd = bd/tot;
    else pdd = 0.;
    
    fNumbersBla[chan]->effTot = efftot;
    fNumbersBla[chan]->pss = pss;
    fNumbersBla[chan]->pdd = pdd;
    double pRatio(0.);
    if (string::npos != fRareName.find("Bs")) pRatio = fsfu;
    if (string::npos != fRareName.find("Lb")) pRatio = fsfu;
    if (string::npos != fRareName.find("Bd")) pRatio = 1.;
    if (string::npos != fRareName.find("Bu")) pRatio = 1.;
    
    // -- yield is the total expectation, based on eff/eff' and N(B+)
    //    it is used to normalize the rare bg histogram
    //    the relevant numbers are extracted from the integrals over specific regions
    double yield = scaledYield(fNumbersBla[chan], fNumbersNo[chan], fRareName, pRatio);
    hRareNw8->Scale(yield*misid*teff[chan]/tot);

    cout << "NNNNNN%%%%%NNNNNNNNNN scaled:        " << hRareNw8->GetSumOfWeights() << endl;
    // -- until here everything was unweighted!
    
    name = Form("hW8MassWithAllCuts_%s_%s_chan%d", modifier.c_str(), fRareName.c_str(), chan);
    hRareW8 = (TH1D*)fHistFile->Get(name.c_str());
    cout << "NNNNNN%%%%%NNNNNNNNNN raw weighted:  " << hRareW8->GetSumOfWeights() << endl;
    cout << "NNNNNN%%%%%NNNNNNNNNN ave fake rate: " << hRareW8->GetSumOfWeights()/hRareW8->GetEntries() << endl;
    hRareW8->Scale(yield*teff[chan]/tot); 
    cout << "NNNNNN%%%%%NNNNNNNNNN scale factor:  " << yield << "*" << teff[chan] << "/" << tot << " = " << yield*teff[chan]/tot << endl;
    cout << "NNNNNN%%%%%NNNNNNNNNN B+  weighted:  " << hRareW8->GetSumOfWeights() << endl;

    hRare = hRareW8;    

    mbin = 0; 
    valInc = hRare->Integral(hRare->FindBin(4.9), hRare->FindBin(fCuts[chan]->mBdLo-fEpsilon));
    error  = valInc*err[fRareName];
    fTEX <<  Form("\\vdef{%s:%s:loSideband%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fRareName.c_str(), chan, valInc) << endl;
    fTEX <<  Form("\\vdef{%s:%s:loSideband%d:err}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fRareName.c_str(), chan, error) << endl;
    // FIXME: add also MuMu to sl!
    //    if (string::npos == fRareName.find("Nu")) {
    if ((string::npos != fRareName.find("Nu")) || (string::npos != fRareName.find("MuMu"))) {
      vRareSl[mbin]  += valInc; 
      vRareSlE[mbin] += error*error; 
    } else {
      vRare[mbin]    += valInc; 
      vRareE[mbin]   += error*error; 
    }
    
    if ((string::npos != fRareName.find("PiPi")))     { vRarePiPi[mbin]  += valInc; }
    if ((string::npos != fRareName.find("KK")))       { vRareKK[mbin]  += valInc; }
    if ((string::npos != fRareName.find("KPi")))      { vRareKPi[mbin]  += valInc; }
    if ((string::npos != fRareName.find("bgLb2KP")))  { vRareKP[mbin]  += valInc; }
    if ((string::npos != fRareName.find("bgLb2PiP"))) { vRarePiP[mbin]  += valInc; }
    

    mbin = 1; 
    valInc = hRare->Integral(hRare->FindBin(fCuts[chan]->mBdLo+fEpsilon), hRare->FindBin(fCuts[chan]->mBsLo-fEpsilon));
    error  = valInc*err[fRareName];
    fTEX <<  Form("\\vdef{%s:%s:bdRare%d}   {\\ensuremath{{%7.6f } } }", fSuffix.c_str(), fRareName.c_str(), chan, valInc) << endl;
    fTEX <<  Form("\\vdef{%s:%s:bdRare%dE}  {\\ensuremath{{%7.6f } } }", fSuffix.c_str(), fRareName.c_str(), chan, error) << endl;
    //    if (string::npos == fRareName.find("Nu")) {
    if ((string::npos != fRareName.find("Nu")) || (string::npos != fRareName.find("MuMu"))) {
      vRareSl[mbin]  += valInc; 
      vRareSlE[mbin] += error*error; 
    } else {
      vRare[mbin]    += valInc; 
      vRareE[mbin]   += error*error; 
    }
    if ((string::npos != fRareName.find("PiPi")))     { vRarePiPi[mbin]  += valInc; }
    if ((string::npos != fRareName.find("KK")))       { vRareKK[mbin]  += valInc; }
    if ((string::npos != fRareName.find("KPi")))      { vRareKPi[mbin]  += valInc; }
    if ((string::npos != fRareName.find("bgLb2KP")))  { vRareKP[mbin]  += valInc; }
    if ((string::npos != fRareName.find("bgLb2PiP"))) { vRarePiP[mbin]  += valInc; }

    mbin = 2; 
    valInc = hRare->Integral(hRare->FindBin(fCuts[chan]->mBsLo+fEpsilon), hRare->FindBin(fCuts[chan]->mBsHi-fEpsilon));
    error  = valInc*err[fRareName];
    cout << "NNNNNN%%%%%NNNNNNNNNN inside Bs box: " << valInc << endl;
    fTEX <<  Form("\\vdef{%s:%s:bsRare%d}   {\\ensuremath{{%7.6f } } }", fSuffix.c_str(), fRareName.c_str(), chan, valInc) << endl;
    fTEX <<  Form("\\vdef{%s:%s:bsRare%dE}  {\\ensuremath{{%7.6f } } }", fSuffix.c_str(), fRareName.c_str(), chan, error) << endl;
    //    if (string::npos == fRareName.find("Nu")) {
    if ((string::npos != fRareName.find("Nu")) || (string::npos != fRareName.find("MuMu"))) {
      vRareSl[mbin]  += valInc; 
      vRareSlE[mbin] += error*error; 
    } else {
      vRare[mbin]    += valInc; 
      vRareE[mbin]   += error*error; 
    }
    if ((string::npos != fRareName.find("PiPi")))     { vRarePiPi[mbin]  += valInc; }
    if ((string::npos != fRareName.find("KK")))       { vRareKK[mbin]  += valInc; }
    if ((string::npos != fRareName.find("KPi")))      { vRareKPi[mbin]  += valInc; }
    if ((string::npos != fRareName.find("bgLb2KP")))  { vRareKP[mbin]  += valInc; }
    if ((string::npos != fRareName.find("bgLb2PiP"))) { vRarePiP[mbin]  += valInc; }

    mbin = 3; 
    valInc = hRare->Integral(hRare->FindBin(fCuts[chan]->mBsHi+fEpsilon), hRare->FindBin(5.9));
    error  = valInc*err[fRareName];
    fTEX <<  Form("\\vdef{%s:%s:hiSideband%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fRareName.c_str(), chan, valInc) << endl;
    fTEX <<  Form("\\vdef{%s:%s:hiSideband%d:err}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fRareName.c_str(), chan, error) << endl;
    //    if (string::npos == fRareName.find("Nu")) {
    if ((string::npos != fRareName.find("Nu")) || (string::npos != fRareName.find("MuMu"))) {
      vRareSl[mbin]  += valInc; 
      vRareSlE[mbin] += error*error; 
    } else {
      vRare[mbin]    += valInc; 
      vRareE[mbin]   += error*error; 
    }
    if ((string::npos != fRareName.find("PiPi")))     { vRarePiPi[mbin]  += valInc; }
    if ((string::npos != fRareName.find("KK")))       { vRareKK[mbin]  += valInc; }
    if ((string::npos != fRareName.find("KPi")))      { vRareKPi[mbin]  += valInc; }
    if ((string::npos != fRareName.find("bgLb2KP")))  { vRareKP[mbin]  += valInc; }
    if ((string::npos != fRareName.find("bgLb2PiP"))) { vRarePiP[mbin]  += valInc; }

    hRare->SetFillColor(colors[fRareName]);
    hRare->SetFillStyle(1000);
    
    if ((string::npos != fRareName.find("Nu")) || (string::npos != fRareName.find("MuMu"))) {
      cout << "&&&&&&&&&&&&&&&&& adding to bslRare: " << fRareName << endl;
      bslRare->Add(hRare); 
      hslRareBg->Add(hRare); 
      leggsl->AddEntry(hRare, fName[fRareName].c_str(), "f"); 
    } else {
      cout << "&&&&&&&&&&&&&&&&& adding to bRare: " << fRareName << endl;
      bRare->Add(hRare); 
      hRareBg->Add(hRare); 
      legg->AddEntry(hRare, fName[fRareName].c_str(), "f"); 
    }

    c0->cd(1); 
    hRare->Draw();
    tl->DrawLatex(0.5, 0.92, fRareName.c_str());
    c0->cd(2); 
    hRare->Draw();
    c0->Modified();
    c0->Update();
    c0->SaveAs(Form("%s/%s-%s_%s_chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), (fDoUseBDT?"bdt":"cnc"), fRareName.c_str(), chan));
  }



  c0->Clear();

  rareBsE = TMath::Sqrt(rareBsE);
  rareBdE = TMath::Sqrt(rareBdE);

  for (int i = 0; i < 4; ++i) {
    vRareE[i]   = TMath::Sqrt(vRareE[i]); 
    vRareSlE[i] = TMath::Sqrt(vRareSlE[i]); 

    vRareHh[i] = vRarePiPi[i] + vRareKK[i] + vRareKPi[i] + vRareKP[i] + vRarePiP[i];
    vRareHhE[i]= quadraticSum(3, 
			      (vRareKPi[i]/epsKa[0] + 2.*vRareKK[i]/epsKa[0] + vRareKP[i]/epsKa[0])*errKa[0]*epsKa[0], 
			      (vRareKPi[i]/epsPi[0] + 2.*vRarePiPi[i]/epsPi[0] + vRareKP[i]/epsPi[0])*errPi[0]*epsPi[0], 
			      (vRareKP[i]/epsPr[0] + 2.*vRarePiP[i]/epsPr[0])*errPr[0]*epsPr[0]
			      );
  }
	 
  shrinkPad(0.12, 0.18); 

  hRareBg->SetMaximum(1.3*hRareBg->GetMaximum()); 
  hRareBg->Draw();
  TH1D *hhRareBg = (TH1D*)hRareBg->GetHistogram(); 
  hhRareBg->SetAxisRange(4.9, 5.9, "X"); 
  setTitles(hhRareBg, "m_{#mu #mu} [GeV]", Form("Candidates/%4.3f GeV", hhRareBg->GetBinWidth(1)), 0.06, 0.9, 1.5);
  legg->SetHeader("CMS simulation"); 
  legg->Draw(); 
  hhRareBg->Draw("same");
  stamp(0.18, fStampString, 0.67, fStampCms); 
  double size = tl->GetTextSize();
  tl->SetTextSize(0.07); 
  tl->DrawLatex(0.25, 0.8, (chan == 0?"Barrel":"Endcap"));   
  tl->SetTextSize(size); 

  string pdfname;
  c0->SaveAs(Form("%s/%s-%s_rare_chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), (fDoUseBDT?"bdt":"cnc"), chan));

  hslRareBg->SetMaximum(1.3*hslRareBg->GetMaximum()); 
  hslRareBg->Draw();
  TH1D *hhslRareBg = (TH1D*)hslRareBg->GetHistogram(); 
  hhslRareBg->SetAxisRange(4.9, 5.9, "X"); 
  setTitles(hhslRareBg, "m_{#mu #mu} [GeV]", Form("Candidates/%4.3f GeV", hhslRareBg->GetBinWidth(1)), 0.06, 0.9, 1.5);
  leggsl->SetHeader("CMS simulation"); 
  leggsl->Draw(); 
  hhslRareBg->Draw("same");
  stamp(0.18, fStampString, 0.67, fStampCms); 
  size = tl->GetTextSize();
  tl->SetTextSize(0.07); 
  tl->DrawLatex(0.25, 0.8, (chan == 0?"Barrel":"Endcap"));   
  tl->SetTextSize(size); 

  c0->SaveAs(Form("%s/%s-%s_slrare_chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), (fDoUseBDT?"bdt":"cnc"), chan));

  TDirectory *pD = gDirectory;
  TFile *f = TFile::Open("hist.root", "UPDATE"); 
  TH1D *bCopy = (TH1D*)bslRare->Clone(Form("bslRare_chan%d", chan)); 
  bCopy->SetDirectory(f); 
  bCopy->Write();
  bCopy = (TH1D*)bRare->Clone(Form("bRare_chan%d", chan)); 
  bCopy->SetDirectory(f); 
  bCopy->Write();
  f->Close();
  gDirectory = pD; 

  fNumbersBs[chan]->bsRare = rareBs; 
  fNumbersBs[chan]->bsRareE= rareBsE; 
  fNumbersBs[chan]->bdRare = rareBd; 
  fNumbersBs[chan]->bdRareE= rareBdE; 

  double relErr(0.05); 
  
  fNumbersBs[chan]->offLoRare = bRare->Integral(bRare->FindBin(4.9001), bRare->FindBin(5.1999));
  fNumbersBs[chan]->offLoRareE= (vRareE[0]/vRare[0])*fNumbersBs[chan]->offLoRare; 
  fNumbersBs[chan]->offHiRare = bRare->Integral(bRare->FindBin(5.4501), bRare->FindBin(5.8999));
  fNumbersBs[chan]->offHiRareE= (vRareE[3]/vRare[3])*fNumbersBs[chan]->offHiRare; 

  // -- new numbers
  fNumbersBs[chan]->fBgPeakLo   = bRare->Integral(bRare->FindBin(4.9001), bRare->FindBin(5.1999));
  fNumbersBs[chan]->fBgPeakLoE1 = relErr * fNumbersBs[chan]->fBgPeakLo; 
  fNumbersBs[chan]->fBgPeakLoE2 = (vRareE[0]/vRare[0]) * fNumbersBs[chan]->fBgPeakLo; 
  fNumbersBs[chan]->fBgPeakBd   = bRare->Integral(bRare->FindBin(fCuts[chan]->mBdLo+fEpsilon), bRare->FindBin(fCuts[chan]->mBdHi-fEpsilon));
  fNumbersBs[chan]->fBgPeakBdE1 = relErr * fNumbersBs[chan]->fBgPeakBd; 
  fNumbersBs[chan]->fBgPeakBdE2 = (vRareE[1]/vRare[1]) * fNumbersBs[chan]->fBgPeakBd; 
  fNumbersBs[chan]->fBgPeakBs   = bRare->Integral(bRare->FindBin(fCuts[chan]->mBsLo+fEpsilon), bRare->FindBin(fCuts[chan]->mBsHi-fEpsilon));
  fNumbersBs[chan]->fBgPeakBsE1 = relErr * fNumbersBs[chan]->fBgPeakBs; 
  fNumbersBs[chan]->fBgPeakBsE2 = (vRareE[2]/vRare[2]) * fNumbersBs[chan]->fBgPeakBs;
  fNumbersBs[chan]->fBgPeakHi   = bRare->Integral(bRare->FindBin(5.4501), bRare->FindBin(5.8999));
  fNumbersBs[chan]->fBgPeakHiE1 = relErr * fNumbersBs[chan]->fBgPeakHi; 
  fNumbersBs[chan]->fBgPeakHiE2 = (vRareE[3]/vRare[3]) * fNumbersBs[chan]->fBgPeakHi;

  fNumbersBs[chan]->fBgRslLo   = bslRare->Integral(bslRare->FindBin(4.9001), bslRare->FindBin(5.1999));
  fNumbersBs[chan]->fBgRslLoE1 = relErr * fNumbersBs[chan]->fBgRslLo; 
  fNumbersBs[chan]->fBgRslLoE2 = (vRareSlE[0]/vRareSl[0]) * fNumbersBs[chan]->fBgRslLo; 
  fNumbersBs[chan]->fBgRslBd   = bslRare->Integral(bslRare->FindBin(fCuts[chan]->mBdLo+fEpsilon), bslRare->FindBin(fCuts[chan]->mBdHi-fEpsilon));
  fNumbersBs[chan]->fBgRslBdE1 = relErr * fNumbersBs[chan]->fBgRslBd; 
  fNumbersBs[chan]->fBgRslBdE2 = (vRareSlE[1]/vRareSl[1]) * fNumbersBs[chan]->fBgRslBd; 
  fNumbersBs[chan]->fBgRslBs   = bslRare->Integral(bslRare->FindBin(fCuts[chan]->mBsLo+fEpsilon), bslRare->FindBin(fCuts[chan]->mBsHi-fEpsilon));
  fNumbersBs[chan]->fBgRslBsE1 = relErr * fNumbersBs[chan]->fBgRslBs; 
  fNumbersBs[chan]->fBgRslBsE2 = (vRareSlE[2]/vRareSl[2]) * fNumbersBs[chan]->fBgRslBs;
  fNumbersBs[chan]->fBgRslHi   = bslRare->Integral(bslRare->FindBin(5.4501), bslRare->FindBin(5.8999));
  fNumbersBs[chan]->fBgRslHiE1 = relErr * fNumbersBs[chan]->fBgRslHi; 
  fNumbersBs[chan]->fBgRslHiE2 = (vRareSlE[3]/vRareSl[3]) * fNumbersBs[chan]->fBgRslHi;

  fNumbersBs[chan]->fBgRareLo   = fNumbersBs[chan]->fBgRslLo + fNumbersBs[chan]->fBgPeakLo;
  fNumbersBs[chan]->fBgRareLoE1 = quadraticSum(2, fNumbersBs[chan]->fBgRslLoE1, fNumbersBs[chan]->fBgPeakLoE1);
  fNumbersBs[chan]->fBgRareLoE2 = quadraticSum(2, fNumbersBs[chan]->fBgRslLoE2, fNumbersBs[chan]->fBgPeakLoE2);
  fNumbersBs[chan]->fBgRareBs   = fNumbersBs[chan]->fBgRslBs + fNumbersBs[chan]->fBgPeakBs;
  fNumbersBs[chan]->fBgRareBsE1 = quadraticSum(2, fNumbersBs[chan]->fBgRslBsE1, fNumbersBs[chan]->fBgPeakBsE1);
  fNumbersBs[chan]->fBgRareBsE2 = quadraticSum(2, fNumbersBs[chan]->fBgRslBsE2, fNumbersBs[chan]->fBgPeakBsE2);
  fNumbersBs[chan]->fBgRareBd   = fNumbersBs[chan]->fBgRslBd + fNumbersBs[chan]->fBgPeakBd;
  fNumbersBs[chan]->fBgRareBdE1 = quadraticSum(2, fNumbersBs[chan]->fBgRslBdE1, fNumbersBs[chan]->fBgPeakBdE1);
  fNumbersBs[chan]->fBgRareBdE2 = quadraticSum(2, fNumbersBs[chan]->fBgRslBdE2, fNumbersBs[chan]->fBgPeakBdE2);
  fNumbersBs[chan]->fBgRareHi   = fNumbersBs[chan]->fBgRslHi + fNumbersBs[chan]->fBgPeakHi;
  fNumbersBs[chan]->fBgRareHiE1 = quadraticSum(2, fNumbersBs[chan]->fBgRslHiE1, fNumbersBs[chan]->fBgPeakHiE1);
  fNumbersBs[chan]->fBgRareHiE2 = quadraticSum(2, fNumbersBs[chan]->fBgRslHiE2, fNumbersBs[chan]->fBgPeakHiE2);


// -- This comparison agrees perfectly!
  cout << "COMPARISON: " 
       << fNumbersBs[chan]->fBgPeakLo 
       << " vs hh sum: " << (vRarePiPi[0]+vRareKK[0]+vRareKPi[0]+vRareKP[0]+vRarePiP[0])
       << " vs hh: " << vRareHh[0]
       << " +/- " << vRareHhE[0]
       << endl;

  cout << "COMPARISON: " 
       << fNumbersBs[chan]->fBgPeakBd 
       << " vs hh sum: " << (vRarePiPi[1]+vRareKK[1]+vRareKPi[1]+vRareKP[1]+vRarePiP[1])
       << " vs hh: " << vRareHh[1]
       << " +/- " << vRareHhE[1]
       << endl;

  cout << "COMPARISON: " 
       << fNumbersBs[chan]->fBgPeakBs
       << " vs hh sum: " << (vRarePiPi[2]+vRareKK[2]+vRareKPi[2]+vRareKP[2]+vRarePiP[2])
       << " vs hh: " << vRareHh[2]
       << " +/- " << vRareHhE[2]
       << endl;

  cout << "COMPARISON: " 
       << fNumbersBs[chan]->fBgPeakHi 
       << " vs hh sum: " << (vRarePiPi[3]+vRareKK[3]+vRareKPi[3]+vRareKP[3]+vRarePiP[3])
       << " vs hh: " << vRareHh[3]
       << " +/- " << vRareHhE[3]
       << endl;


  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << "% --- rare Background summary numbers" << endl;  
  fTEX <<  Form("\\vdef{%s:bsRare%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, rareBs) << endl;
  fTEX <<  Form("\\vdef{%s:bsRare%dE}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, rareBsE) << endl;
  fTEX <<  Form("\\vdef{%s:bdRare%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, rareBd) << endl;
  fTEX <<  Form("\\vdef{%s:bdRare%dE}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, rareBdE) << endl;

  fSuffix = cache; 

}



// ----------------------------------------------------------------------
void plotResults::calculateSgNumbers(int chan) {

  // -- efficiency and acceptance
  fSetup = "SgMc"; 
  numbersFromHist(chan, 0, fNumbersBs[chan]); 
  numbersFromHist(chan, 1, fNumbersBd[chan]); 

  fSetup = "SgData"; 
  string name(""); 
  
  string modifier = (fDoUseBDT?"bdt":"cnc"); 
  name = Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), 5, chan);
  TH1D *h1 = (TH1D*)fHistFile->Get(name.c_str());
  TH1D *h1u = (TH1D*)h1->Clone(Form("h1u%d", chan)); 

  name = Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), 5, 1-chan);
  TH1D *ho = (TH1D*)fHistFile->Get(name.c_str());

  name = Form("hMassWithAllCutsManyBins_%s_%d_chan%d", modifier.c_str(), 5, chan);
  TH1D *h = (TH1D*)fHistFile->Get(name.c_str());

  // -- blinded version
  setHist(h1, kBlack, 20, 1.); 
  setTitles(h1, "m_{#mu#mu} [GeV]", Form("Candidates / %3.3f GeV", h1->GetBinWidth(1)), 0.05, 1.2, 1.3); 
  h1->SetMinimum(0.01); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h1->SetAxisRange(4.9, 5.9-fEpsilon, "X"); 
  h1->SetMaximum(h1->GetBinContent(h1->GetMaximumBin())+2); 

  double ymax0 = (h1->GetBinContent(h1->GetMaximumBin()) > ho->GetBinContent(ho->GetMaximumBin()) ? 
		  h1->GetBinContent(h1->GetMaximumBin())+3:ho->GetBinContent(ho->GetMaximumBin())+3);
  h1->SetMaximum(ymax0); 
  //  h1->Draw();
  //       drawArrow(0.5, 2); 
  //       drawArrow(0.5, 1); 
  bgBlind(h1, 5, 5.45, 5.9);
  drawBox(2, 0.5); 
  drawBox(1, 0.5); 

  TH1D *dummy1 = new TH1D("dummy1", "", 10, 0., 10.); setFilledHist(dummy1, kBlue, kBlue, 3356); 
  TH1D *dummy2 = new TH1D("dummy2", "", 10, 0., 10.); setFilledHist(dummy2, kRed, kRed, 3365); 
      
  newLegend(0.4, 0.65, 0.8, 0.8); 
  legg->SetTextSize(0.04);  
  
  legg->AddEntry(dummy1, "B_{s}^{0} signal window", "f"); 
  legg->AddEntry(dummy2, "B^{0} signal window", "f"); 
  if (0 == chan) {
    legg->SetHeader("Barrel");   
  } else {
    legg->SetHeader("Endcap");   
  }      
  
  legg->Draw();
  
  
  stamp(0.18, fStampString, 0.67, fStampCms); 
  if (fDoPrint)  {
    if (fDoUseBDT) c0->SaveAs(Form("%s/%s-bdtsig-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan));
    else c0->SaveAs(Form("%s/%s-cncsig-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan));
  }



  // -- unblinded version
  h1u->SetAxisRange(4.9, 5.9-fEpsilon, "X"); 
  h1u->Draw();
  double yleg = ymax0 - 0.3;
  h1u->SetMaximum(ymax0); 
  
  drawArrow(0.6, 1, yleg); 
  drawArrow(0.4, 2, yleg-0.5); 
  //   tl->DrawLatex(0.4, yleg, Form("obs: %3.0f", bs));
  //   tl->DrawLatex(0.4, yleg-0.35, Form("obs: %3.0f", b0));
  //   tl->DrawLatex(0.2, 0.92, (0 == i?"Barrel":"Endcap"));

  stamp(0.18, fStampString, 0.67, fStampCms); 
  if (fDoPrint)  {
    if (fDoUseBDT) c0->SaveAs(Form("%s/%s-bdtunblinded-sig-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan));
    else c0->SaveAs(Form("%s/%s-cncunblinded-sig-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan));
  }
  

  // -- unblinded version
  h->SetAxisRange(4.9, 5.9, "X"); 
  h->Draw();
  
  drawArrow(0.6, 1, yleg); 
  drawArrow(0.4, 2, yleg-0.5); 

  stamp(0.18, fStampString, 0.67, fStampCms); 
  if (fDoPrint)  {
    if (fDoUseBDT) c0->SaveAs(Form("%s/%s-bdtunblinded-manybins-sig-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan));
    else c0->SaveAs(Form("%s/%s-cncunblinded-manybins-sig-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan));
  }
  

  double eps(0.00001);
  fNumbersBs[chan]->bsObs = h->Integral(h->FindBin(fCuts[chan]->mBsLo+eps), h->FindBin(fCuts[chan]->mBsHi-eps)); 
  fNumbersBs[chan]->bdObs = h->Integral(h->FindBin(fCuts[chan]->mBdLo+eps), h->FindBin(fCuts[chan]->mBdHi-eps)); 

  fNumbersBs[chan]->offHi = fBgHistHi;
  fNumbersBs[chan]->offLo = fBgHistLo; 
  fNumbersBs[chan]->bgObs = fBgHistHi + fBgHistLo;

  // -- non-peaking background
  fNumbersBs[chan]->bgBsExp    = fBsBgExp;	
  fNumbersBs[chan]->bgBsExpE   = fBsBgExpE;
  fNumbersBs[chan]->bgBdExp    = fBdBgExp; 
  fNumbersBs[chan]->bgBdExpE   = fBdBgExpE;
  
  double relCombError = (fBgHistHi>0?1./TMath::Sqrt(fBgHistHi): 1.);
  cout << "^^^^^^^^^^^ relative statistical error on combinatorial bg: " << relCombError << " from " << fBgHistHi << endl;

  // -- new numbers: combined combinatorial and scaled rare sl bg  (E2 means stat + syst here!)

  // -- scaled sl rare bg
  double err1(0.), err2(0.); 
  err1 = (fNumbersBs[chan]->fBgRslLo>0?fNumbersBs[chan]->fBgRslLoE1/fNumbersBs[chan]->fBgRslLo:0.); 
  err2 = (fNumbersBs[chan]->fBgRslLo>0?fNumbersBs[chan]->fBgRslLoE2/fNumbersBs[chan]->fBgRslLo:0.); 
  fNumbersBs[chan]->fBgRslsLo   = fLoSlBgExp;
  fNumbersBs[chan]->fBgRslsLoE1 = err1*fNumbersBs[chan]->fBgRslsLo;
  fNumbersBs[chan]->fBgRslsLoE2 = err2*fNumbersBs[chan]->fBgRslsLo;

  err1 = (fNumbersBs[chan]->fBgRslBd>0?fNumbersBs[chan]->fBgRslBdE1/fNumbersBs[chan]->fBgRslBd:0.); 
  err2 = (fNumbersBs[chan]->fBgRslBd>0?fNumbersBs[chan]->fBgRslBdE2/fNumbersBs[chan]->fBgRslBd:0.); 
  fNumbersBs[chan]->fBgRslsBd   = fBdSlBgExp;
  fNumbersBs[chan]->fBgRslsBdE1 = err1*fNumbersBs[chan]->fBgRslsBd;
  fNumbersBs[chan]->fBgRslsBdE2 = err2*fNumbersBs[chan]->fBgRslsBd;

  err1 = (fNumbersBs[chan]->fBgRslBs>0?fNumbersBs[chan]->fBgRslBsE1/fNumbersBs[chan]->fBgRslBs:0.); 
  err2 = (fNumbersBs[chan]->fBgRslBs>0?fNumbersBs[chan]->fBgRslBsE2/fNumbersBs[chan]->fBgRslBs:0.); 
  fNumbersBs[chan]->fBgRslsBs   = fBsSlBgExp;
  fNumbersBs[chan]->fBgRslsBsE1 = err1*fNumbersBs[chan]->fBgRslsBs;
  fNumbersBs[chan]->fBgRslsBsE2 = err2*fNumbersBs[chan]->fBgRslsBs;

  err1 = (fNumbersBs[chan]->fBgRslHi>0?fNumbersBs[chan]->fBgRslHiE1/fNumbersBs[chan]->fBgRslHi:0.); 
  err2 = (fNumbersBs[chan]->fBgRslHi>0?fNumbersBs[chan]->fBgRslHiE2/fNumbersBs[chan]->fBgRslHi:0.); 
  fNumbersBs[chan]->fBgRslsHi   = fHiSlBgExp;
  fNumbersBs[chan]->fBgRslsHiE1 = err1*fNumbersBs[chan]->fBgRslsHi;
  fNumbersBs[chan]->fBgRslsHiE2 = err2*fNumbersBs[chan]->fBgRslsHi;

  // -- combinatorial bg: E2 means stat + syst here!
  double combSystRel = 0.1; 
  fNumbersBs[chan]->fBgCombLo   = fLoCoBgExp;
  fNumbersBs[chan]->fBgCombLoE1 = relCombError*fLoCoBgExp;
  fNumbersBs[chan]->fBgCombLoE2 = quadraticSum(2, relCombError*fLoCoBgExp, (fLoCoBgExp>0? combSystRel*fLoCoBgExp: combSystRel)); 

  fNumbersBs[chan]->fBgCombBd   = fBdCoBgExp;
  fNumbersBs[chan]->fBgCombBdE1 = relCombError*fBdCoBgExp;
  fNumbersBs[chan]->fBgCombBdE2 = quadraticSum(2, relCombError*fBdCoBgExp, (fBdCoBgExp>0? combSystRel*fBdCoBgExp: combSystRel));

  fNumbersBs[chan]->fBgCombBs   = fBsCoBgExp;
  fNumbersBs[chan]->fBgCombBsE1 = relCombError*fBsCoBgExp;
  fNumbersBs[chan]->fBgCombBsE2 = quadraticSum(2, relCombError*fBsCoBgExp, (fBsCoBgExp>0? combSystRel*fBsCoBgExp: combSystRel));

  fNumbersBs[chan]->fBgCombHi   = fHiCoBgExp;
  fNumbersBs[chan]->fBgCombHiE1 = relCombError*fHiCoBgExp;
  fNumbersBs[chan]->fBgCombHiE2 = quadraticSum(2, relCombError*fHiCoBgExp, (fHiCoBgExp>0? combSystRel*fHiCoBgExp: combSystRel));

  // -- non-peaking sum
  fNumbersBs[chan]->fBgNonpLo   = fLoBgExp;
  fNumbersBs[chan]->fBgNonpLoE1 = quadraticSum(2, fNumbersBs[chan]->fBgCombLoE1, fNumbersBs[chan]->fBgRslsLoE1);
  fNumbersBs[chan]->fBgNonpLoE2 = quadraticSum(2, fNumbersBs[chan]->fBgCombLoE2, fNumbersBs[chan]->fBgRslsLoE2);

  fNumbersBs[chan]->fBgNonpBd   = fBdBgExp;
  fNumbersBs[chan]->fBgNonpBdE1 = quadraticSum(2, fNumbersBs[chan]->fBgCombBdE1, fNumbersBs[chan]->fBgRslsBdE1);
  fNumbersBs[chan]->fBgNonpBdE2 = quadraticSum(2, fNumbersBs[chan]->fBgCombBdE2, fNumbersBs[chan]->fBgRslsBdE2);

  fNumbersBs[chan]->fBgNonpBs   = fBsBgExp;
  fNumbersBs[chan]->fBgNonpBsE1 = quadraticSum(2, fNumbersBs[chan]->fBgCombBsE1, fNumbersBs[chan]->fBgRslsBsE1);
  fNumbersBs[chan]->fBgNonpBsE2 = quadraticSum(2, fNumbersBs[chan]->fBgCombBdE2, fNumbersBs[chan]->fBgRslsBsE2);

  fNumbersBs[chan]->fBgNonpHi   = fHiBgExp;
  fNumbersBs[chan]->fBgNonpHiE1 = quadraticSum(2, fNumbersBs[chan]->fBgCombHiE1, fNumbersBs[chan]->fBgRslsHiE1);
  fNumbersBs[chan]->fBgNonpHiE2 = quadraticSum(2, fNumbersBs[chan]->fBgCombHiE2, fNumbersBs[chan]->fBgRslsHiE2);

  fNumbersBs[chan]->tauBs  = fBsBgExp/(fLoBgExp + fHiBgExp); 
  fNumbersBs[chan]->tauBsE = 0.04*fNumbersBs[chan]->tauBs; // it's only the SHAPE that is relevant here

  fNumbersBs[chan]->tauBd = fBdBgExp/(fLoBgExp + fHiBgExp); 
  fNumbersBs[chan]->tauBdE = 0.04*fNumbersBs[chan]->tauBd; // it's only the SHAPE that is relevant here


  // -- all backgrounds
  fNumbersBs[chan]->fBgTotLo   = fNumbersBs[chan]->fBgNonpLo + fNumbersBs[chan]->fBgPeakLo;
  fNumbersBs[chan]->fBgTotLoE2 = quadraticSum(2, fNumbersBs[chan]->fBgNonpLoE2, fNumbersBs[chan]->fBgPeakLoE2); 
  fNumbersBs[chan]->fBgTotBd   = fNumbersBs[chan]->fBgNonpBd + fNumbersBs[chan]->fBgPeakBd;
  fNumbersBs[chan]->fBgTotBdE2 = quadraticSum(2, fNumbersBs[chan]->fBgNonpBdE2,  fNumbersBs[chan]->fBgPeakBdE2); 
  fNumbersBs[chan]->fBgTotBs   = fNumbersBs[chan]->fBgNonpBs + fNumbersBs[chan]->fBgPeakBs;
  fNumbersBs[chan]->fBgTotBsE2 = quadraticSum(2, fNumbersBs[chan]->fBgNonpBsE2, fNumbersBs[chan]->fBgPeakBsE2);
  fNumbersBs[chan]->fBgTotHi   = fNumbersBs[chan]->fBgNonpHi + fNumbersBs[chan]->fBgPeakHi;
  fNumbersBs[chan]->fBgTotHiE2 = quadraticSum(2, fNumbersBs[chan]->fBgNonpHiE2, fNumbersBs[chan]->fBgPeakHiE2);


  // -- scale SM expectation to normalization and compute total signal boxes contents
  for (unsigned int i = 0; i < fNchan; ++i) {
    
    fNumbersBs[i]->mBdHi = fCuts[i]->mBdHi;
    fNumbersBs[i]->mBdLo = fCuts[i]->mBdLo;
    fNumbersBs[i]->mBsHi = fCuts[i]->mBsHi;
    fNumbersBs[i]->mBsLo = fCuts[i]->mBsLo;

    // -- Bs -> mu mu 
    name  = Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), 0, chan);
    h1    = (TH1D*)fHistFile->Get(name.c_str()); 
    scaledYield(fNumbersBs[i], fNumbersNo[i], "SgMc", fsfu);

    // -- B0 -> mu mu 
    name  = Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), 1, chan);
    h1    = (TH1D*)fHistFile->Get(name.c_str()); 
    scaledYield(fNumbersBd[i], fNumbersNo[i], "BdMc", 1.);

    // -- contents in boxes 
    fNumbersBs[i]->bsExpObs = fNumbersBs[i]->bgBsExp + fNumbersBs[i]->bsRare + fNumbersBs[i]->bsNoScaled;
    cout << "old approach: " <<     fNumbersBs[i]->bsExpObs; 

    fNumbersBs[i]->bsExpObs = fNumbersBs[i]->fBgNonpBs + fNumbersBs[i]->fBgPeakBs + fNumbersBs[i]->fSgBs + fNumbersBs[i]->fSgBd;
    cout << " new approach: " <<     fNumbersBs[i]->bsExpObs << endl;
  
    fNumbersBs[i]->bsExpObsE = TMath::Sqrt(fNumbersBs[i]->bgBsExpE*fNumbersBs[i]->bgBsExpE //??
					   + fNumbersBs[i]->bsRareE*fNumbersBs[i]->bsRareE
					   + fNumbersBs[i]->bsNoScaledE*fNumbersBs[i]->bsNoScaledE);
    
    fNumbersBs[i]->bdExpObs = fNumbersBs[i]->bgBdExp + fNumbersBs[i]->bdRare + fNumbersBd[i]->bdNoScaled;
    cout << "old approach: " <<     fNumbersBs[i]->bdExpObs; 
    
    fNumbersBs[i]->bdExpObs = fNumbersBs[i]->fBgNonpBd + fNumbersBs[i]->fBgPeakBd + fNumbersBs[i]->fBdBd + fNumbersBs[i]->fSgBd;
    cout << " new approach: " <<     fNumbersBs[i]->bdExpObs << endl;


    fNumbersBs[i]->bdExpObsE = TMath::Sqrt(fNumbersBs[i]->bgBdExpE*fNumbersBs[i]->bgBdExpE //??
					   + fNumbersBs[i]->bdRareE*fNumbersBs[i]->bdRareE
					   + fNumbersBd[i]->bdNoScaledE*fNumbersBd[i]->bdNoScaledE);
  }


  // -- new numbers: Bd yields
  fNumbersBs[chan]->fBdLo   = fNumbersBd[chan]->fBdLo;
  fNumbersBs[chan]->fBdLoE1 = fNumbersBd[chan]->fBdLoE1;
  fNumbersBs[chan]->fBdLoE2 = fNumbersBd[chan]->fBdLoE2;

  fNumbersBs[chan]->fBdBd   = fNumbersBd[chan]->fBdBd;
  fNumbersBs[chan]->fBdBdE1 = fNumbersBd[chan]->fBdBdE1;
  fNumbersBs[chan]->fBdBdE2 = fNumbersBd[chan]->fBdBdE1;

  fNumbersBs[chan]->fBdBs   = fNumbersBd[chan]->fBdBs;
  fNumbersBs[chan]->fBdBsE1 = fNumbersBd[chan]->fBdBsE1;
  fNumbersBs[chan]->fBdBsE2 = fNumbersBd[chan]->fBdBsE2;

  fNumbersBs[chan]->fBdHi   = fNumbersBd[chan]->fBdHi;
  fNumbersBs[chan]->fBdHiE1 = fNumbersBd[chan]->fBdHiE1;
  fNumbersBs[chan]->fBdHiE2 = fNumbersBd[chan]->fBdHiE2;

  // -- all backgrounds + signal + cross-feed
  fNumbersBs[chan]->fSgAndBgLo  = fNumbersBs[chan]->fBgNonpLo + fNumbersBs[chan]->fBgPeakLo + fNumbersBs[chan]->fSgLo + fNumbersBs[chan]->fBdLo;
  fNumbersBs[chan]->fSgAndBgLoE2= quadraticSum(4, fNumbersBs[chan]->fBgNonpLoE2, fNumbersBs[chan]->fBgPeakLoE2, fNumbersBs[chan]->fSgLoE2, fNumbersBs[chan]->fBdLoE2);
  fNumbersBs[chan]->fSgAndBgBd  = fNumbersBs[chan]->fBgNonpBd + fNumbersBs[chan]->fBgPeakBd + fNumbersBs[chan]->fSgBd + fNumbersBs[chan]->fBdBd;
  fNumbersBs[chan]->fSgAndBgBdE2= quadraticSum(4, fNumbersBs[chan]->fBgNonpBdE2, fNumbersBs[chan]->fBgPeakBdE2, fNumbersBs[chan]->fSgBdE2, fNumbersBs[chan]->fBdBdE2);
  fNumbersBs[chan]->fSgAndBgBs  = fNumbersBs[chan]->fBgNonpBs + fNumbersBs[chan]->fBgPeakBs + fNumbersBs[chan]->fSgBs + fNumbersBs[chan]->fBdBs;
  fNumbersBs[chan]->fSgAndBgBsE2= quadraticSum(4, fNumbersBs[chan]->fBgNonpBsE2, fNumbersBs[chan]->fBgPeakBsE2, fNumbersBs[chan]->fSgBsE2, fNumbersBs[chan]->fBdBsE2);
  fNumbersBs[chan]->fSgAndBgHi  = fNumbersBs[chan]->fBgNonpHi + fNumbersBs[chan]->fBgPeakHi + fNumbersBs[chan]->fSgHi + fNumbersBs[chan]->fBdHi;
  fNumbersBs[chan]->fSgAndBgHiE2= quadraticSum(4, fNumbersBs[chan]->fBgNonpHiE2, fNumbersBs[chan]->fBgPeakHiE2, fNumbersBs[chan]->fSgHiE2, fNumbersBs[chan]->fBdHiE2);

  computeErrors(fNumbersBs); 
  computeErrors(fNumbersBd); 
}


// ----------------------------------------------------------------------
void plotResults::fillAndSaveHistograms(int nevents) {

  // -- dump histograms
  string hfname  = fDirectory + "/anaBmm.plotResults." + fSuffix + ".root";
  cout << "fHistFile: " << hfname;
  fHistFile = TFile::Open(hfname.c_str(), "UPDATE");
  cout << " opened, running on " << nevents << " entries" << endl;

  TTree *t(0);

  fSaveSmallTree = true; 

  // -- DBX
  if (0) {
    resetHistograms();
    fSetup = "SgMcAcc";
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    otherNumbers(fSetup); 
    saveHists(fSetup);
    return;
  } // FIXME

  if (1) {
    // -- rare backgrounds
    resetHistograms();
    rareBgHists("nada", nevents); 
  }
    
  if (1) {
    // -- normalization
    resetHistograms();
    fSetup = "NoData"; 
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    saveHists(fSetup);

    resetHistograms();
    fSetup = "NoMc";
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    otherNumbers(fSetup); 
    saveHists(fSetup);
  }

  if (1) {
    // -- dimuons
    resetHistograms();
    fSetup = "SgData"; 
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    saveHists(fSetup);

    resetHistograms();
    fSetup = "SgMc";
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    otherNumbers(fSetup); 
    saveHists(fSetup);

    resetHistograms();
    fSetup = "BdMc";
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    otherNumbers(fSetup); 
    saveHists(fSetup);

    resetHistograms();
    fSetup = "SgMcAcc";
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    saveHists(fSetup);
  }


  if (1) {
    // -- control sample
    resetHistograms();
    fSetup = "CsData"; 
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    saveHists(fSetup);
    
    resetHistograms();
    fSetup = "CsMc";
    t = getTree(fSetup); 
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);
    otherNumbers(fSetup); 
    saveHists(fSetup);
  }

  // -- close the file
  fHistFile->cd();
  fHistFile->Close(); 
  
  fSaveSmallTree = false; 

}


// ----------------------------------------------------------------------
void plotResults::saveLargerTree(string mode) {

  TTree *t(0); 
  
  fSaveLargerTree = true; 
  fSaveSmallTree = true; 
  
  t = getTree(mode); 
  setupTree(t, mode); 
  loopOverTree(t, mode, 1);
  
  
}



// ----------------------------------------------------------------------
void plotResults::otherNumbers(string smode) {
  int ibin(0); 
  string accname = smode + "Acc"; 
  string directory("candAnaMuMu"); 
  numbers* aa[2]; 
  if (string::npos != fSetup.find("No"))  {
    directory = "candAnaBu2JpsiK";
    aa[0] = fNumbersNo[0]; 
    aa[1] = fNumbersNo[1]; 
  }

  if (string::npos != fSetup.find("Cs"))  {
    directory = "candAnaBs2JpsiPhi";
    aa[0] = fNumbersCs[0]; 
    aa[1] = fNumbersCs[1]; 
  }

  if (string::npos != fSetup.find("Sg"))  {
    directory = "candAnaMuMu";
    aa[0] = fNumbersBs[0]; 
    aa[1] = fNumbersBs[1]; 
  }

  if (string::npos != fSetup.find("Bd"))  {
    directory = "candAnaMuMu";
    aa[0] = fNumbersBd[0]; 
    aa[1] = fNumbersBd[1]; 
  }


  cout << "accname: " << accname << " directory: " << directory 
       << " numbers: " << aa[0]->name << " index = " << aa[0]->index 
       << endl;
  
  for (unsigned int i = 0; i < fNchan; ++i) {
    fChan = i; 
    accEffFromEffTree(accname, directory, *aa[i], *fCuts[i], -1);
    fF[fSetup]->cd(directory.c_str());
    double effFilter = fFilterEff[fSetup];
    if (effFilter < 1e-6) effFilter = 1.0;
    double genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
    ibin = 1; fhGenAndAccNumbers[i]->SetBinContent(ibin, effFilter); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "effFilter");
    ibin = 2; fhGenAndAccNumbers[i]->SetBinContent(ibin, genFileYield); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "genFileYield");
    ibin = 10;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->effPtReco); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "effPtReco");
    ibin = 11;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->effPtRecoE); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "effPtRecoE");
    ibin = 12;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->genAccFileYield); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "genAccFileYield");
    ibin = 13;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->genAccYield); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "genAccYield");
    ibin = 14;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->recoYield); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "recoYield");
    ibin = 15;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->muidYield); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "muidYield");
    ibin = 16;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->trigYield); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "trigYield");
    ibin = 17;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->candYield); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "candYield");
    ibin = 18;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->chanYield); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "chanYield");

    ibin = 20;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->acc); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "acc");
    ibin = 21;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->accE); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "accE");

    ibin = 30;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->effCand); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "effCand");
    ibin = 31;fhGenAndAccNumbers[i]->SetBinContent(ibin, aa[i]->effCandE); 
              fhGenAndAccNumbers[i]->GetXaxis()->SetBinLabel(ibin, "effCandE");
  }
}



// ----------------------------------------------------------------------
void plotResults::saveHists(string smode) {

  fHistFile->cd();

  int mode(0); 
  if ("SgMc" == smode)   mode = 0; 
  if ("BdMc" == smode)   mode = 1; 
  if ("SgMcAcc" == smode)mode = 2; 
  if ("SgData" == smode) mode = 5; 

  if ("NoMc" == smode)   mode = 10; 
  if ("NoData" == smode) mode = 15; 

  if ("CsMc" == smode)   mode = 20; 
  if ("CsData" == smode) mode = 25; 
 
  TH1D *h1(0); 
  TH2D *h2(0); 

  for (unsigned int i = 0; i < fNchan; ++i) {
    string modifier = (fDoUseBDT?"bdt":"cnc"); 
    h1 = (TH1D*)(fhGenAndAccNumbers[i]->Clone(Form("hGenAndAccNumbers_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h2 = (TH2D*)(fhAccAll[i]->Clone(Form("hAccAll_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h2->SetDirectory(fHistFile); 
    h2->Write();
    h2 = (TH2D*)(fhAccPass[i]->Clone(Form("hAccPass_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h2->SetDirectory(fHistFile); 
    h2->Write();

    h1 = (TH1D*)(fhAccEtaAll[i]->Clone(Form("hAccEtaAll_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();
    h1 = (TH1D*)(fhAccEtaPass[i]->Clone(Form("hAccEtaPass_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();

    h1 = (TH1D*)(fhAccPtAll[i]->Clone(Form("hAccPtAll_%s_%d_chan%d", modifier.c_str(), mode, i))); 
    h1->SetDirectory(fHistFile); 
    h1->Write();
    h1 = (TH1D*)(fhAccPtPass[i]->Clone(Form("hAccPtPass_%s_%d_chan%d", modifier.c_str(), mode, i))); 
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
void plotResults::loopFunction(int function, int mode) {
  if (1 == function) loopFunction1(mode);
  if (2 == function) loopFunction2(mode);
}


// ----------------------------------------------------------------------
void plotResults::loopFunction1(int mode) {

  // 173692/2704432277
  // 172992/459341492
  //  const Long64_t printRun(172992), printEvt(459341492); 
  //   const Long64_t printRun(170854), printEvt(280769371); 
  //   if (fb.run == printRun && fb.evt == printEvt) {
  //     cout << "BDT = " << fBDT << ": "; 
  //     printRedTreeEvt(fb); 
  //   }

  if (fChan < 0) return;
  
  double mass = fb.m; 

  bool bp2jpsikp(false), bs2jpsiphi(false); 
  if (10 == mode) bp2jpsikp = true; 
  if (15 == mode) bp2jpsikp = true; 

  if (20 == mode) bs2jpsiphi = true; 
  if (25 == mode) bs2jpsiphi = true; 
  
  fhMassAbsNoCuts[fChan]->Fill(mass);
  fhAccAll[fChan]->Fill(TMath::Abs(fb.eta), fb.pt); 
  fhAccPtAll[fChan]->Fill(fb.pt); 
  fhAccEtaAll[fChan]->Fill(TMath::Abs(fb.eta)); 

  if (!fGoodAcceptance) return;
  fhAccPass[fChan]->Fill(TMath::Abs(fb.eta), fb.pt); 
  fhAccPtPass[fChan]->Fill(fb.pt); 
  fhAccEtaPass[fChan]->Fill(TMath::Abs(fb.eta)); 

  //   if (fb.run == printRun && fb.evt == printEvt) cout << "after acceptance" << endl;

  // -- this is the base, after the raw acceptance cuts
  fhMassNoCuts[fChan]->Fill(mass);
  fhMassNoCutsManyBins[fChan]->Fill(mass); 

  if (fDoUseBDT) {
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "plotResults fGoodQ: " << fGoodQ << " from " 
    // 						       << fb.m1q << " " << fb.m2q << endl;
    if (!fGoodQ) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after goodq" << endl;
    if (!fGoodPvAveW8) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after pvw8" << endl;
    if (!fGoodTracks) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after gtqual" << endl;
    if (!fGoodTracksPt) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after gtpt" << endl;
    if (!fGoodTracksEta) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after gteta" << endl;
    if (!fGoodBdtPt) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after gmpt" << endl;
    if (!fGoodMuonsEta) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after gmeta" << endl;
    if (!fGoodJpsiCuts) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after jpsi" << endl;
    if (fBDT < fCuts[fChan]->bdt) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after maxbdt" << endl;
    if (fBDT > fCuts[fChan]->bdtMax) return;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "after minbdt" << endl;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "chan = " << fChan << endl;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "muonid = " << fGoodMuonsID << endl;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "hlt = " << fGoodHLT << endl;
    //     if (fb.run == printRun && fb.evt == printEvt) cout << "mass = " << mass << endl;
  } else {
    if (!fGoodQ) return;
    if (!fGoodMuonsPt) return;
    if (!fGoodMuonsEta) return;
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

  //   if (fBDT > 0.2 
  //       && fChan == 0 
  //       && fGoodMuonsID 
  //       && fGoodHLT 
  //       && 4.9 < mass && mass < 5.9) 
  //     printRedTreeEvt(fb); 
    
  fhMassWithAnaCuts[fChan]->Fill(mass); 
  fhMassWithAnaCutsManyBins[fChan]->Fill(mass); 

  // -- weighted with fake rate
  fhW8MassWithAllCuts[fChan]->Fill(mass, fW8MisId);
  if (5 == mode && !(5.2 < mass && mass < 5.45)) {
    fhW8MassWithAllCutsBlind[fChan]->Fill(mass, fW8MisId);
  }
  fhW8MassWithAllCutsManyBins[fChan]->Fill(mass, fW8MisId);

  // -- MUON ID
  if (false == fGoodMuonsID) return;
  fhMassWithMuonCuts[fChan]->Fill(mass); 
  fhMassWithMuonCutsManyBins[fChan]->Fill(mass); 

  // -- muon ID: Data PidTables AFTER THE REQUIREMENT!
  if (fW8DmuID > 0.) fhMuId[fChan]->Fill(fW8DmuID, 1./fW8DmuID); 
  // -- muon ID: MC PidTables
  if (fW8MmuID > 0.) fhMuIdMC[fChan]->Fill(fW8MmuID, 1./fW8MmuID); 

  //   PidTable *pT, *pT1, *pT2; 
  //   double w1(0.), w2(0.), am1eta(0.), am2eta(0.); 
  //   if (1) {
  //     if (fIsCowboy) {
  //       pT  = fptCbMMC; 
  //       pT1 = fptCbT1MC;
  //       pT2 = fptCbT2MC;
  //     } else {
  //       pT  = fptSgMMC; 
  //       pT1 = fptSgT1MC;
  //       pT2 = fptSgT2MC;
  //     }
  //     am1eta = TMath::Abs(fb.m1eta);
  //     am2eta = TMath::Abs(fb.m2eta);
  //     w1       = pT1->effD(fb.m1pt, am1eta, fb.m1phi) * pT2->effD(fb.m1pt, am1eta, fb.m1phi);
  //     w2       = pT1->effD(fb.m2pt, am2eta, fb.m2phi) * pT2->effD(fb.m2pt, am2eta, fb.m2phi);
  //   }
  
  //   // -- TRIGGER
  //   if (fChan == 1) cout << "fGoodHLT: " << fGoodHLT 
  // 		       << " fw8Mtrig = " << fW8Mtrig 
  // 		       << " w8 = " << w1*w2 
  // 		       << " pt:  "  << fb.m1pt << "/" << fb.m2pt  
  // 		       << " eta:  " << fb.m1eta << "/" << fb.m2eta  
  // 		       << " L1L2: " << pT1->effD(fb.m1pt, am1eta, fb.m1phi) << "/" <<  pT1->effD(fb.m2pt, am2eta, fb.m2phi)
  // 		       << " L3: "   << pT2->effD(fb.m1pt, am1eta, fb.m1phi) << "/" <<  pT2->effD(fb.m2pt, am2eta, fb.m2phi)
  // 		       << endl;

  if (false == fGoodHLT) return;
  fhMassWithTriggerCuts[fChan]->Fill(mass); 
  fhMassWithTriggerCutsManyBins[fChan]->Fill(mass); 

  // -- muon trigger: Data PidTables AFTER THE REQUIREMENT!
  if (fW8Dtrig > 0.) fhMuTr[fChan]->Fill(fW8Dtrig, 1./fW8Dtrig); 
  // -- muon trigger: MC PidTables
  if (fW8Mtrig > 0.) fhMuTrMC[fChan]->Fill(fW8Mtrig, 1./fW8Mtrig); 
  
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
void plotResults::resetHistograms() {

  for (unsigned int i = 0; i < fNchan; ++i) {
    fhGenAndAccNumbers[i]->Reset();

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

    fhW8MassWithAllCuts[i]->Reset();
    fhW8MassWithAllCutsBlind[i]->Reset();
    fhW8MassWithAllCutsManyBins[i]->Reset();

    fhMassWithMassCutsManyBins[i]->Reset();
    fhMassWithMassCuts[i]->Reset();

    fhNorm[i]->Reset();
    fhNormC[i]->Reset();
    fhDstarPi[i]->Reset();

    fhMuId[i]->Reset();
    fhMuIdMC[i]->Reset();
    fhMuTr[i]->Reset();
    fhMuTrMC[i]->Reset();

    fhBdtMass[i]->Reset();

    fhAccAll[i]->Reset(); 
    fhAccPass[i]->Reset();  
    fhAccPtAll[i]->Reset();  
    fhAccPtPass[i]->Reset(); 
    fhAccEtaAll[i]->Reset(); 
    fhAccEtaPass[i]->Reset();
   


  }

}

// ----------------------------------------------------------------------
void plotResults::loopFunction2(int mode) {
  cout << "loopfunction2: " << mode << endl;
}


// ----------------------------------------------------------------------
void plotResults::fitHists(int chan) {

  string hfname  = fDirectory + "/anaBmm.plotResults." + fSuffix + ".root";
  cout << "open fHistFile: " << hfname << " readonly" << endl;
  fHistFile = new TFile(hfname.c_str());

  string name; 
  TH1D *h(0); 

//   // -- normalization sample 
//   name = Form("hNorm_bdt_15_chan%d", chan);
//   h =  (TH1D*)fHistFile->Get(name.c_str()); 
//   normYield2(h, chan, 5.0); 

//   // -- control sample 
//   name = Form("hNorm_bdt_25_chan%d", chan);
//   h = (TH1D*)fHistFile->Get(name.c_str()); 
//   csYield(h, chan, 5.1); 

  // -- dimuon background 
  name = Form("hMassWithAllCuts_bdt_5_chan%d", chan);
  h = (TH1D*)fHistFile->Get(name.c_str()); 

  cout << "fSgLo: " << fSgLo << " fSgHi: " << fSgHi << endl;

  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 
  gStyle->SetOptTitle(0); 

  double xt(0.5), yt(0.8), ydec(0.05); 

  zone(1);
  bgBlind(h, 0, fBgLo, fBgHi); 
  cout << "==> mode 0: expected bg = " << fBsBgExp << "+/-" << fBsBgExpE 
       << " Lo: " << fBgHistLo << "/" << fLoBgExp
       << " Hi: " << fBgHistHi << "/" << fHiBgExp 
       << endl;
  tl->SetTextSize(0.03);
  tl->DrawLatex(xt, yt,        Form("expected: %4.2f+/-%4.2f", fBsBgExp, fBsBgExpE)); 
  tl->DrawLatex(xt, yt-ydec,   Form("Lo (obs/exp): %3.0f/%3.1f", fBgHistLo, fLoBgExp)); 
  tl->DrawLatex(xt, yt-2.*ydec,Form("Hi (obs/exp): %3.0f/%3.1f", fBgHistHi, fHiBgExp)); 
  c0->SaveAs(Form("%s/%s-bgestimate-mode0-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan)); 

  zone(1);
  bgBlind(h, 1, fBgLo, fBgHi); 
  cout << "==> mode 1: expected bg = " << fBsBgExp << "+/-" << fBsBgExpE 
       << " Lo: " << fBgHistLo << "/" << fLoBgExp  
       << " Hi: " << fBgHistHi << "/" << fHiBgExp
       << endl;
  tl->DrawLatex(xt, yt,        Form("expected: %4.2f+/-%4.2f", fBsBgExp, fBsBgExpE)); 
  tl->DrawLatex(xt, yt-ydec,   Form("Lo (obs/exp): %3.0f/%3.1f", fBgHistLo, fLoBgExp)); 
  tl->DrawLatex(xt, yt-2.*ydec,Form("Hi (obs/exp): %3.0f/%3.1f", fBgHistHi, fHiBgExp)); 
  c0->SaveAs(Form("%s/%s-bgestimate-mode1-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan)); 

  zone(1);
  bgBlind(h, 2, fBgLo, fBgHi); 
  cout << "==> mode 2: expected bg = " << fBsBgExp << "+/-" << fBsBgExpE 
       << " Lo: " << fBgHistLo << "/" << fLoBgExp
       << " Hi: " << fBgHistHi << "/" << fHiBgExp  
       << endl;
  tl->DrawLatex(xt, yt,        Form("expected: %4.2f+/-%4.2f", fBsBgExp, fBsBgExpE)); 
  tl->DrawLatex(xt, yt-ydec,   Form("Lo (obs/exp): %3.0f/%3.1f", fBgHistLo, fLoBgExp)); 
  tl->DrawLatex(xt, yt-2.*ydec,Form("Hi (obs/exp): %3.0f/%3.1f", fBgHistHi, fHiBgExp)); 
  c0->SaveAs(Form("%s/%s-bgestimate-mode2-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan)); 


  zone(1);
  bgBlind(h, 3, fBgLo, fBgHi); 
  cout << "==> mode 3: expected bg = " << fBsBgExp << "+/-" << fBsBgExpE 
       << " Lo: " << fBgHistLo << "/" << fLoBgExp
       << " Hi: " << fBgHistHi << "/" << fHiBgExp  
       << endl;
  tl->DrawLatex(xt, yt,        Form("expected: %4.2f+/-%4.2f", fBsBgExp, fBsBgExpE)); 
  tl->DrawLatex(xt, yt-ydec,   Form("Lo (obs/exp): %3.0f/%3.1f", fBgHistLo, fLoBgExp)); 
  tl->DrawLatex(xt, yt-2.*ydec,Form("Hi (obs/exp): %3.0f/%3.1f", fBgHistHi, fHiBgExp)); 
  c0->SaveAs(Form("%s/%s-bgestimate-mode3-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan)); 

  zone(1);
  bgBlind(h, 4, 5.4, 5.9); 
  cout << "==> mode 4: expected bg = " << fBsBgExp << "+/-" << fBsBgExpE 
       << " Lo: " << fBgHistLo << "/" << fLoBgExp
       << " Hi: " << fBgHistHi << "/" << fHiBgExp  
       << endl;
  tl->DrawLatex(xt, yt,        Form("expected: %4.2f+/-%4.2f", fBsBgExp, fBsBgExpE)); 
  tl->DrawLatex(xt, yt-ydec,   Form("Lo (obs/exp): %3.0f/%3.1f", fBgHistLo, fLoBgExp)); 
  tl->DrawLatex(xt, yt-2.*ydec,Form("Hi (obs/exp): %3.0f/%3.1f", fBgHistHi, fHiBgExp)); 
  c0->SaveAs(Form("%s/%s-bgestimate-mode4-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan)); 

  zone(1);
  bgBlind(h, 5, 5.4, 5.9); 
  cout << "==> mode 5: expected bg = " << fBsBgExp << "+/-" << fBsBgExpE 
       << " Lo: " << fBgHistLo << "/" << fLoBgExp  
       << " Hi: " << fBgHistHi << "/" << fHiBgExp
       << endl;
  tl->DrawLatex(xt, yt,        Form("expected: %4.2f+/-%4.2f", fBsBgExp, fBsBgExpE)); 
  tl->DrawLatex(xt, yt-ydec,   Form("Lo (obs/exp): %3.0f/%3.1f", fBgHistLo, fLoBgExp)); 
  tl->DrawLatex(xt, yt-2.*ydec,Form("Hi (obs/exp): %3.0f/%3.1f", fBgHistHi, fHiBgExp)); 
  c0->SaveAs(Form("%s/%s-bgestimate-mode5-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), chan)); 

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
      yield = scaledYield(fNumbersBs[i], fNumbersNo[i], "SgMc", fsfu);
      h = (TH1D*)(fhMassWithAllCuts[i]->Clone(Form("Bs_%s_chan%d", cache.c_str(), i)));  
      h->Scale(yield/tot);
    } 
    else if (1 == mode) {
      tot  = fhMassWithAllCuts[i]->Integral(0, fhMassWithAllCuts[i]->GetNbinsX()+1); 
      yield = scaledYield(fNumbersBd[i], fNumbersNo[i], "BdMc", 1.);
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
  for (unsigned int i = 0; i < lines.size(); ++i) {
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

    syst = sysMuid[i]*a[i]->effRhoMuidTNPMC;
    a[i]->effRhoMuidTNPMCTE = TMath::Sqrt(a[i]->effRhoMuidTNPMCE*a[i]->effRhoMuidTNPMCE + syst*syst); 

    syst = sysMuid[i]*a[i]->effRhoMuidTNP;
    a[i]->effRhoMuidTNPMCTE = TMath::Sqrt(a[i]->effRhoMuidTNPE*a[i]->effRhoMuidTNPE + syst*syst); 

    // -- trigger
    syst = sysTrig[i]*a[i]->effTrigMC;
    a[i]->effTrigTE = TMath::Sqrt(a[i]->effTrigMCE*a[i]->effTrigMCE + syst*syst); 

    syst = sysTrig[i]*a[i]->effTrigMC;
    a[i]->effTrigMCTE = TMath::Sqrt(a[i]->effTrigMCE*a[i]->effTrigMCE + syst*syst); 

    syst = sysTrig[i]*a[i]->effTrigTNP;
    a[i]->effTrigTNPTE = TMath::Sqrt(a[i]->effTrigTNPE*a[i]->effTrigTNPE + syst*syst); 

    syst = sysTrig[i]*a[i]->effTrigTNPMC;
    a[i]->effTrigTNPMCTE = TMath::Sqrt(a[i]->effTrigTNPMCE*a[i]->effTrigTNPMCE + syst*syst); 

    syst = sysTrig[i]*a[i]->effRhoTrigTNPMC;
    a[i]->effRhoTrigTNPMCTE = TMath::Sqrt(a[i]->effRhoTrigTNPMCE*a[i]->effRhoTrigTNPMCE + syst*syst); 

    syst = sysTrig[i]*a[i]->effRhoTrigTNP;
    a[i]->effRhoTrigTNPTE = TMath::Sqrt(a[i]->effRhoTrigTNPE*a[i]->effRhoTrigTNPE + syst*syst); 

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
  string cache = fSuffix; 
  fSuffix = (fDoUseBDT?"bdt":"cnc") + fSuffix;

  if (1) {

    //     printNumbers(*fNumbersNo[0], cout); 
    //     printNumbers(*fNumbersNo[1], cout); 

    for (unsigned int i = 0; i < fNchan; ++i) {
      OUT << "# -- NORMALIZATION " << i << endl;
      fTEX << "% -- NORMALIZATION " << i << endl;

      OUT << "#EFF_TOT_BPLUS\t" << i << "\t" << fNumbersNo[i]->effTot << endl;
      fTEX << formatTex(fNumbersNo[i]->effTot, Form("%s:N-EFF-TOT-BPLUS%i:val", fSuffix.c_str(), i), 5) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTotE, Form("%s:N-EFF-TOT-BPLUS%i:err", fSuffix.c_str(), i), 6) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTotTE, Form("%s:N-EFF-TOT-BPLUS%i:tot", fSuffix.c_str(), i), 5) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTot, fNumbersNo[i]->effTotTE, 
			    Form("%s:N-EFF-TOT-BPLUS%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;

      OUT << "#EFF_PRODMC_BPLUS\t" << i << "\t" << fNumbersNo[i]->effProdMC << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdMC, Form("%s:N-EFF-PRODMC-BPLUS%i:val", fSuffix.c_str(), i), 5) << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdMCE, Form("%s:N-EFF-PRODMC-BPLUS%i:err", fSuffix.c_str(), i), 6) << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdMCE, Form("%s:N-EFF-PRODMC-BPLUS%i:tot", fSuffix.c_str(), i), 5) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effProdMC, fNumbersNo[i]->effProdMCE, 
			    Form("%s:N-EFF-PRODMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;

      OUT << "#EFF_PRODTNP_BPLUS\t" << i << "\t" << fNumbersNo[i]->effProdTNP << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdTNP, Form("%s:N-EFF-PRODTNP-BPLUS%i:val", fSuffix.c_str(), i), 5) << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdTNPE, Form("%s:N-EFF-PRODTNP-BPLUS%i:err", fSuffix.c_str(), i), 6) << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdTNPE, Form("%s:N-EFF-PRODTNP-BPLUS%i:tot", fSuffix.c_str(), i), 5) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effProdTNP, fNumbersNo[i]->effProdTNPE, 
			    Form("%s:N-EFF-PRODTNP-BPLUS%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;

      OUT << "#EFF_PRODTNPMC_BPLUS\t" << i << "\t" << fNumbersNo[i]->effProdTNPMC << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdTNPMC, Form("%s:N-EFF-PRODTNPMC-BPLUS%i:val", fSuffix.c_str(), i), 5) << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdTNPMCE, Form("%s:N-EFF-PRODTNPMC-BPLUS%i:err", fSuffix.c_str(), i), 6) << endl;
      fTEX << formatTex(fNumbersNo[i]->effProdTNPMCE, Form("%s:N-EFF-PRODTNPMC-BPLUS%i:tot", fSuffix.c_str(), i), 5) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effProdMC, fNumbersNo[i]->effProdTNPMCE, 
			    Form("%s:N-EFF-PRODTNPMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;


      OUT << "ACC_BPLUS\t" << i << "\t" << fNumbersNo[i]->acc <<"\t" << fNumbersNo[i]->accTE  << endl;
      fTEX << formatTex(fNumbersNo[i]->acc, Form("%s:N-ACC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->accE, Form("%s:N-ACC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->accTE, Form("%s:N-ACC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->acc, fNumbersNo[i]->accTE, 
			    Form("%s:N-ACC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      cout << "1: " << fNumbersNo[i]->effMuidTNP << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNP, Form("%s:N-EFF-MU-PID-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPE, Form("%s:N-EFF-MU-PID-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPTE, Form("%s:N-EFF-MU-PID-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effMuidTNP, fNumbersNo[i]->effMuidTNPTE, 
			    Form("%s:N-EFF-MU-PID-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      cout << "2: " << fNumbersNo[i]->effRhoMuidTNP << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoMuidTNP, Form("%s:N-EFFRHO-MU-PID-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoMuidTNPE, Form("%s:N-EFFRHO-MU-PID-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoMuidTNPTE, Form("%s:N-EFFRHO-MU-PID-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effRhoMuidTNP, fNumbersNo[i]->effRhoMuidTNPTE, 
			    Form("%s:N-EFFRHO-MU-PID-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      cout << "3: " << fNumbersNo[i]->effMuidMC << endl;
      OUT << "EFF_MU_BPLUS\t" << i << "\t"  << fNumbersNo[i]->effMuidMC << "\t"  << "0." << endl;

      cout << "4: " << fNumbersNo[i]->effMuidTNPMC << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPMCE, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effMuidTNPMCTE, Form("%s:N-EFF-MU-PIDMC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effMuidTNPMC, fNumbersNo[i]->effMuidTNPMCTE, 
			    Form("%s:N-EFF-MU-PIDMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      cout << "5: " << fNumbersNo[i]->effRhoMuidTNPMC << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoMuidTNPMC, Form("%s:N-EFFRHO-MU-PIDMC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoMuidTNPMCE, Form("%s:N-EFFRHO-MU-PIDMC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoMuidTNPMCTE, Form("%s:N-EFFRHO-MU-PIDMC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effRhoMuidTNPMC, fNumbersNo[i]->effRhoMuidTNPMCTE, 
			    Form("%s:N-EFFRHO-MU-PIDMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

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

      fTEX << formatTex(fNumbersNo[i]->effRhoTrigTNP, Form("%s:N-EFFRHO-TRIG-PID-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoTrigTNPE, Form("%s:N-EFFRHO-TRIG-PID-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoTrigTNPTE, Form("%s:N-EFFRHO-TRIG-PID-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effRhoTrigTNP, fNumbersNo[i]->effRhoTrigTNPTE, 
			    Form("%s:N-EFFRHO-TRIG-PID-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      fTEX << formatTex(fNumbersNo[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigTNPMCE, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effTrigTNPMCTE, Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effTrigTNPMC, fNumbersNo[i]->effTrigTNPMCTE, 
			    Form("%s:N-EFF-TRIG-PIDMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

      fTEX << formatTex(fNumbersNo[i]->effRhoTrigTNPMC, Form("%s:N-EFFRHO-TRIG-PIDMC-BPLUS%i:val", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoTrigTNPMCE, Form("%s:N-EFFRHO-TRIG-PIDMC-BPLUS%i:err", fSuffix.c_str(), i), 3) << endl;
      fTEX << formatTex(fNumbersNo[i]->effRhoTrigTNPMCTE, Form("%s:N-EFFRHO-TRIG-PIDMC-BPLUS%i:tot", fSuffix.c_str(), i), 3) << endl;
      fTEX << scientificTex(fNumbersNo[i]->effRhoTrigTNPMC, fNumbersNo[i]->effRhoTrigTNPMCTE, 
			    Form("%s:N-EFFRHO-TRIG-PIDMC-BPLUS%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

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
  double scaledSig, scaledSigE(0.); 
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

    OUT << "#EFF_PRODMC_BSMM\t" << i << "\t" << fNumbersBs[i]->effProdMC << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdMC, Form("%s:N-EFF-PRODMC-BSMM%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdMCE, Form("%s:N-EFF-PRODMC-BSMM%i:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdMCE, Form("%s:N-EFF-PRODMC-BSMM%i:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effProdMC, fNumbersBs[i]->effProdMCE, 
			  Form("%s:N-EFF-PRODMC-BSMM%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;

    // -- efficiency ratio
    fTEX << formatTex(fNumbersNo[i]->effTot/fNumbersBs[i]->effTot, 
		      Form("%s:N-EFFRATIO-TOT-BSMM%i:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(dRatio(fNumbersNo[i]->effTot, fNumbersNo[i]->effTotE, fNumbersBs[i]->effTot, fNumbersBs[i]->effTotE), 
		      Form("%s:N-EFFRATIO-TOT-BSMM%i:err", fSuffix.c_str(), i), 3) << endl;

    fTEX << formatTex(fNumbersNo[i]->effProdMC/fNumbersBs[i]->effProdMC, 
		      Form("%s:N-EFFRATIO-PRODMC-BSMM%i:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(dRatio(fNumbersNo[i]->effProdMC, fNumbersNo[i]->effProdMCE, fNumbersBs[i]->effProdMC, fNumbersBs[i]->effProdMCE), 
		      Form("%s:N-EFFRATIO-PRODMC-BSMM%i:err", fSuffix.c_str(), i), 3) << endl;

    fTEX << formatTex(fNumbersNo[i]->effProdTNP/fNumbersBs[i]->effProdTNP, 
		      Form("%s:N-EFFRATIO-PRODTNP-BSMM%i:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(dRatio(fNumbersNo[i]->effProdTNP, fNumbersNo[i]->effProdTNPE, fNumbersBs[i]->effProdTNP, fNumbersBs[i]->effProdTNPE),
		      Form("%s:N-EFFRATIO-PRODTNP-BSMM%i:err", fSuffix.c_str(), i), 3) << endl;

    fTEX << formatTex(fNumbersNo[i]->effProdTNPMC/fNumbersBs[i]->effProdTNPMC, 
		      Form("%s:N-EFFRATIO-PRODTNPMC-BSMM%i:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(dRatio(fNumbersNo[i]->effProdTNPMC, fNumbersNo[i]->effProdTNPMCE, 
			     fNumbersBs[i]->effProdTNPMC, fNumbersBs[i]->effProdTNPMCE), 
		      Form("%s:N-EFFRATIO-PRODTNPMC-BSMM%i:err", fSuffix.c_str(), i), 3) << endl;

    
    OUT << "#EFF_PRODTNP_BSMM\t" << i << "\t" << fNumbersBs[i]->effProdTNP << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdTNP, Form("%s:N-EFF-PRODTNP-BSMM%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdTNPE, Form("%s:N-EFF-PRODTNP-BSMM%i:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdTNPE, Form("%s:N-EFF-PRODTNP-BSMM%i:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effProdTNP, fNumbersBs[i]->effProdTNPE, 
			  Form("%s:N-EFF-PRODTNP-BSMM%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;
    
    OUT << "#EFF_PRODTNPMC_BSMM\t" << i << "\t" << fNumbersBs[i]->effProdTNPMC << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdTNPMC, Form("%s:N-EFF-PRODTNPMC-BSMM%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdTNPMCE, Form("%s:N-EFF-PRODTNPMC-BSMM%i:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBs[i]->effProdTNPMCE, Form("%s:N-EFF-PRODTNPMC-BSMM%i:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effProdMC, fNumbersBs[i]->effProdTNPMCE, 
			  Form("%s:N-EFF-PRODTNPMC-BSMM%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;


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

    fTEX << formatTex(fNumbersBs[i]->effRhoMuidTNP, Form("%s:N-EFFRHO-MU-PID-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effRhoMuidTNPE, Form("%s:N-EFFRHO-MU-PID-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effRhoMuidTNPTE, Form("%s:N-EFFRHO-MU-PID-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effRhoMuidTNP, fNumbersBs[i]->effRhoMuidTNPTE, 
			  Form("%s:N-EFFRHO-MU-PID-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBs[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidTNPMCE, Form("%s:N-EFF-MU-PIDMC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effMuidTNPMCTE, Form("%s:N-EFF-MU-PIDMC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effMuidTNPMC, fNumbersBs[i]->effMuidTNPMCTE, 
			  Form("%s:N-EFF-MU-PIDMC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBs[i]->effRhoMuidTNPMC, Form("%s:N-EFFRHO-MU-PIDMC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effRhoMuidTNPMCE, Form("%s:N-EFFRHO-MU-PIDMC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effRhoMuidTNPMCTE, Form("%s:N-EFFRHO-MU-PIDMC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effRhoMuidTNPMC, fNumbersBs[i]->effRhoMuidTNPMCTE, 
			  Form("%s:N-EFFRHO-MU-PIDMC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

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

    fTEX << formatTex(fNumbersBs[i]->effRhoTrigTNP, Form("%s:N-EFFRHO-TRIG-PID-BSMM%i:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effRhoTrigTNPE, Form("%s:N-EFFRHO-TRIG-PID-BSMM%i:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effRhoTrigTNPTE, Form("%s:N-EFFRHO-TRIG-PID-BSMM%i:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effRhoTrigTNP, fNumbersBs[i]->effRhoTrigTNPTE, 
			  Form("%s:N-EFFRHO-TRIG-PID-BSMM%i:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBs[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigTNPMCE, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effTrigTNPMCTE, Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effTrigTNPMC, fNumbersBs[i]->effTrigTNPMCTE, 
			  Form("%s:N-EFF-TRIG-PIDMC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBs[i]->effRhoTrigTNPMC, Form("%s:N-EFFRHO-TRIG-PIDMC-BSMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effRhoTrigTNPMCE, Form("%s:N-EFFRHO-TRIG-PIDMC-BSMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBs[i]->effRhoTrigTNPMCTE, Form("%s:N-EFFRHO-TRIG-PIDMC-BSMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBs[i]->effRhoTrigTNPMC, fNumbersBs[i]->effRhoTrigTNPMCTE, 
			  Form("%s:N-EFFRHO-TRIG-PIDMC-BSMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

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
    fTEX << formatTex(fNumbersBd[i]->effTotTE, Form("%s:N-EFF-TOT-BDMM%d:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effTot, fNumbersBd[i]->effTotTE, 
			  Form("%s:N-EFF-TOT-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    OUT << "#EFF_PRODMC_BDMM\t" << i << "\t" << fNumbersBd[i]->effProdMC << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdMC, Form("%s:N-EFF-PRODMC-BDMM%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdMCE, Form("%s:N-EFF-PRODMC-BDMM%i:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdMCE, Form("%s:N-EFF-PRODMC-BDMM%i:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effProdMC, fNumbersBd[i]->effProdMCE, 
			  Form("%s:N-EFF-PRODMC-BDMM%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;
    
    OUT << "#EFF_PRODTNP_BDMM\t" << i << "\t" << fNumbersBd[i]->effProdTNP << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdTNP, Form("%s:N-EFF-PRODTNP-BDMM%i:val", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdTNPE, Form("%s:N-EFF-PRODTNP-BDMM%i:err", fSuffix.c_str(), i), 4) << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdTNPE, Form("%s:N-EFF-PRODTNP-BDMM%i:tot", fSuffix.c_str(), i), 4) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effProdTNP, fNumbersBd[i]->effProdTNPE, 
			  Form("%s:N-EFF-PRODTNP-BDMM%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;
    
    OUT << "#EFF_PRODTNPMC_BDMM\t" << i << "\t" << fNumbersBd[i]->effProdTNPMC << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdTNPMC, Form("%s:N-EFF-PRODTNPMC-BDMM%i:val", fSuffix.c_str(), i), 5) << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdTNPMCE, Form("%s:N-EFF-PRODTNPMC-BDMM%i:err", fSuffix.c_str(), i), 6) << endl;
    fTEX << formatTex(fNumbersBd[i]->effProdTNPMCE, Form("%s:N-EFF-PRODTNPMC-BDMM%i:tot", fSuffix.c_str(), i), 5) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effProdMC, fNumbersBd[i]->effProdTNPMCE, 
			  Form("%s:N-EFF-PRODTNPMC-BDMM%i:all", fSuffix.c_str(), i), 1e-3, 2) << endl;

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

    fTEX << formatTex(fNumbersBd[i]->effRhoMuidTNP, Form("%s:N-EFFRHO-MU-PID-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effRhoMuidTNPE, Form("%s:N-EFFRHO-MU-PID-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effRhoMuidTNPTE, Form("%s:N-EFFRHO-MU-PID-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effRhoMuidTNP, fNumbersBd[i]->effRhoMuidTNPTE, 
			  Form("%s:N-EFFRHO-MU-PID-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBd[i]->effMuidTNPMC, Form("%s:N-EFF-MU-PIDMC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidTNPMCE, Form("%s:N-EFF-MU-PIDMC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effMuidTNPMCTE, Form("%s:N-EFF-MU-PIDMC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effMuidTNPMC, fNumbersBd[i]->effMuidTNPMCTE, 
			  Form("%s:N-EFF-MU-PIDMC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBd[i]->effRhoMuidTNPMC, Form("%s:N-EFFRHO-MU-PIDMC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effRhoMuidTNPMCE, Form("%s:N-EFFRHO-MU-PIDMC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effRhoMuidTNPMCTE, Form("%s:N-EFFRHO-MU-PIDMC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effRhoMuidTNPMC, fNumbersBd[i]->effRhoMuidTNPMCTE, 
			  Form("%s:N-EFFRHO-MU-PIDMC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

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

    fTEX << formatTex(fNumbersBd[i]->effRhoTrigTNP, Form("%s:N-EFFRHO-TRIG-PID-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effRhoTrigTNPE, Form("%s:N-EFFRHO-TRIG-PID-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effRhoTrigTNPTE, Form("%s:N-EFFRHO-TRIG-PID-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effRhoTrigTNP, fNumbersBd[i]->effRhoTrigTNPTE, 
			  Form("%s:N-EFFRHO-TRIG-PID-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBd[i]->effTrigTNPMC, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigTNPMCE, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effTrigTNPMCTE, Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effTrigTNPMC, fNumbersBd[i]->effTrigTNPMCTE, 
			  Form("%s:N-EFF-TRIG-PIDMC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

    fTEX << formatTex(fNumbersBd[i]->effRhoTrigTNPMC, Form("%s:N-EFFRHO-TRIG-PIDMC-BDMM%d:val", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effRhoTrigTNPMCE, Form("%s:N-EFFRHO-TRIG-PIDMC-BDMM%d:err", fSuffix.c_str(), i), 3) << endl;
    fTEX << formatTex(fNumbersBd[i]->effRhoTrigTNPMCTE, Form("%s:N-EFFRHO-TRIG-PIDMC-BDMM%d:tot", fSuffix.c_str(), i), 3) << endl;
    fTEX << scientificTex(fNumbersBd[i]->effRhoTrigTNPMC, fNumbersBd[i]->effRhoTrigTNPMCTE, 
			  Form("%s:N-EFFRHO-TRIG-PIDMC-BDMM%d:all", fSuffix.c_str(), i), 1e-2, 2) << endl;

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

    cout << "=====================> fBgPeakBs[" << i << "] = " << fNumbersBs[i]->fBgPeakBs << endl;
    cout << "=====================> fBgPeakBs[" << i << "] = " 
	 << formatTex(fNumbersBs[i]->fBgPeakBs, Form("%s:N-PEAK-BKG-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    OUT << "PEAK_BKG_BS\t" << i << "\t" << fNumbersBs[i]->fBgPeakBs << "\t" << fNumbersBs[i]->fBgPeakBsE2 << endl;
    fTEX << formatTex(fNumbersBs[i]->fBgPeakBs, Form("%s:N-PEAK-BKG-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->fBgPeakBsE2,Form("%s:N-PEAK-BKG-BS%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "PEAK_BKG_BD\t" << i << "\t"	<< fNumbersBs[i]->fBgPeakBd << "\t" << fNumbersBs[i]->fBgPeakBdE2 << endl;
    fTEX << formatTex(fNumbersBs[i]->fBgPeakBd, Form("%s:N-PEAK-BKG-BD%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->fBgPeakBdE2,Form("%s:N-PEAK-BKG-BD%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "TAU_BS\t" << i << "\t" << fNumbersBs[i]->tauBs << "\t" << fNumbersBs[i]->tauBsE << endl;
    fTEX << formatTex(fNumbersBs[i]->tauBs, Form("%s:N-TAU-BS%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->tauBsE, Form("%s:N-TAU-BS%d:err", fSuffix.c_str(), i), 2) << endl;

    OUT << "TAU_BD\t" << i << "\t" << fNumbersBs[i]->tauBd << "\t" << fNumbersBs[i]->tauBdE << endl;
    fTEX << formatTex(fNumbersBs[i]->tauBd, Form("%s:N-TAU-BD%d:val", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fNumbersBs[i]->tauBdE, Form("%s:N-TAU-BD%d:err", fSuffix.c_str(), i), 2) << endl;


    fTEX << formatTex(fNumbersBs[i]->offHi, Form("%s:N-OBS-OFFHI%d:val", fSuffix.c_str(), i), 0) << endl;
    fTEX << formatTex(fNumbersBs[i]->offLo, Form("%s:N-OBS-OFFLO%d:val", fSuffix.c_str(), i), 0) << endl;

    double sb = fNumbersBs[i]->bsNoScaled/(fNumbersBs[i]->bsRare + fNumbersBs[i]->bgBsExp);
    fTEX << formatTex(sb, Form("%s:N-EXP-SoverB%d:val", fSuffix.c_str(), i), 2) << endl;
    double ssb = fNumbersBs[i]->bsNoScaled/TMath::Sqrt(fNumbersBs[i]->bsNoScaled + fNumbersBs[i]->bsRare + fNumbersBs[i]->bgBsExp);
    fTEX << formatTex(ssb, Form("%s:N-EXP-SoverSplusB%d:val", fSuffix.c_str(), i), 2) << endl;
  }


  for (unsigned int chan = 0; chan < fNchan; ++chan) {

    fTEX << "% ----------------------------------------------------------------------" << endl;
    fTEX << "% --- signal summary numbers (new numbers)" << endl;  

    fTEX <<  Form("\\vdef{%s:SgBd%d:val}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgBd) << endl;
    fTEX <<  Form("\\vdef{%s:SgBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:SgBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgBdE2) << endl;

    fTEX <<  Form("\\vdef{%s:SgBs%d:val}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgBs) << endl;
    fTEX <<  Form("\\vdef{%s:SgBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:SgBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgBsE2) << endl;

    fTEX <<  Form("\\vdef{%s:SgLo%d:val}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgLo) << endl;
    fTEX <<  Form("\\vdef{%s:SgLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:SgLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgLoE2) << endl;

    fTEX <<  Form("\\vdef{%s:SgHi%d:val}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgHi) << endl;
    fTEX <<  Form("\\vdef{%s:SgHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:SgHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgHiE2) << endl;

    fTEX <<  Form("\\vdef{%s:BdBd%d:val}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdBd) << endl;
    fTEX <<  Form("\\vdef{%s:BdBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:BdBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdBdE2) << endl;

    fTEX <<  Form("\\vdef{%s:BdBs%d:val}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdBs) << endl;
    fTEX <<  Form("\\vdef{%s:BdBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:BdBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdBsE2) << endl;

    fTEX <<  Form("\\vdef{%s:BdLo%d:val}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdLo) << endl;
    fTEX <<  Form("\\vdef{%s:BdLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:BdLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdLoE2) << endl;

    fTEX <<  Form("\\vdef{%s:BdHi%d:val}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdHi) << endl;
    fTEX <<  Form("\\vdef{%s:BdHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:BdHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBd[chan]->fBdHiE2) << endl;


    fTEX << "% ----------------------------------------------------------------------" << endl;
    fTEX << "% --- background summary numbers (new numbers)" << endl;  

    fTEX <<  Form("\\vdef{%s:BgPeakLo%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakLo) << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakLoE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgPeakBd%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakBd) << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakBdE2) << endl;

    cout << "=====================> fBgPeakBs[" << chan << "] = " << fNumbersBs[chan]->fBgPeakBs << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakBs) << endl;
    cout << "=====================> fBgPeakBs[" << chan << "] = " 
	 << Form("\\vdef{%s:BgPeakBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakBs) << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakBsE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgPeakHi%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakHi) << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgPeakHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgPeakHiE2) << endl;


    fTEX <<  Form("\\vdef{%s:BgRslLo%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslLo) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslLoE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgRslBd%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslBd) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslBdE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgRslBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslBs) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslBsE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgRslHi%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslHi) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslHiE2) << endl;


    fTEX <<  Form("\\vdef{%s:BgRareLo%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareLo) << endl;
    fTEX <<  Form("\\vdef{%s:BgRareLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRareLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareLoE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgRareBd%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareBd) << endl;
    fTEX <<  Form("\\vdef{%s:BgRareBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRareBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareBdE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgRareBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareBs) << endl;
    fTEX <<  Form("\\vdef{%s:BgRareBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRareBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareBsE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgRareHi%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareHi) << endl;
    fTEX <<  Form("\\vdef{%s:BgRareHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRareHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRareHiE2) << endl;


    fTEX <<  Form("\\vdef{%s:BgRslsLo%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsLo) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslsLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslsLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsLoE2) << endl;
				   
    fTEX <<  Form("\\vdef{%s:BgRslsBd%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsBd) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslsBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslsBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsBdE2) << endl;
				   
    fTEX <<  Form("\\vdef{%s:BgRslsBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsBs) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslsBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslsBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsBsE2) << endl;
				   
    fTEX <<  Form("\\vdef{%s:BgRslsHi%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsHi) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslsHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgRslsHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgRslsHiE2) << endl;


    fTEX <<  Form("\\vdef{%s:BgCombLo%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombLo) << endl;
    fTEX <<  Form("\\vdef{%s:BgCombLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgCombLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombLoE2) << endl;
				   
    fTEX <<  Form("\\vdef{%s:BgCombBd%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombBd) << endl;
    fTEX <<  Form("\\vdef{%s:BgCombBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgCombBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombBdE2) << endl;
				   
    fTEX <<  Form("\\vdef{%s:BgCombBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombBs) << endl;
    fTEX <<  Form("\\vdef{%s:BgCombBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgCombBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombBsE2) << endl;
				   
    fTEX <<  Form("\\vdef{%s:BgCombHi%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombHi) << endl;
    fTEX <<  Form("\\vdef{%s:BgCombHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgCombHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgCombHiE2) << endl;


    fTEX <<  Form("\\vdef{%s:BgNonpLo%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpLo) << endl;
    fTEX <<  Form("\\vdef{%s:BgNonpLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgNonpLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpLoE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgNonpBd%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpBd) << endl;
    fTEX <<  Form("\\vdef{%s:BgNonpBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgNonpBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpBdE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgNonpBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpBs) << endl;
    fTEX <<  Form("\\vdef{%s:BgNonpBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgNonpBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpBsE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgNonpHi%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpHi) << endl;
    fTEX <<  Form("\\vdef{%s:BgNonpHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgNonpHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgNonpHiE2) << endl;


    fTEX <<  Form("\\vdef{%s:BgTotLo%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotLo) << endl;
    fTEX <<  Form("\\vdef{%s:BgTotLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgTotLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotLoE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgTotBd%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotBd) << endl;
    fTEX <<  Form("\\vdef{%s:BgTotBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgTotBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotBdE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgTotBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotBs) << endl;
    fTEX <<  Form("\\vdef{%s:BgTotBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgTotBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotBsE2) << endl;

    fTEX <<  Form("\\vdef{%s:BgTotHi%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotHi) << endl;
    fTEX <<  Form("\\vdef{%s:BgTotHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:BgTotHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fBgTotHiE2) << endl;


    fTEX <<  Form("\\vdef{%s:SgAndBgLo%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgLo) << endl;
    fTEX <<  Form("\\vdef{%s:SgAndBgLo%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgLoE1) << endl;
    fTEX <<  Form("\\vdef{%s:SgAndBgLo%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgLoE2) << endl;

    fTEX <<  Form("\\vdef{%s:SgAndBgBd%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgBd) << endl;
    fTEX <<  Form("\\vdef{%s:SgAndBgBd%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgBdE1) << endl;
    fTEX <<  Form("\\vdef{%s:SgAndBgBd%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgBdE2) << endl;

    fTEX <<  Form("\\vdef{%s:SgAndBgBs%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgBs) << endl;
    fTEX <<  Form("\\vdef{%s:SgAndBgBs%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgBsE1) << endl;
    fTEX <<  Form("\\vdef{%s:SgAndBgBs%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgBsE2) << endl;

    fTEX <<  Form("\\vdef{%s:SgAndBgHi%d:val}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgHi) << endl;
    fTEX <<  Form("\\vdef{%s:SgAndBgHi%d:e1}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgHiE1) << endl;
    fTEX <<  Form("\\vdef{%s:SgAndBgHi%d:e2}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), chan, fNumbersBs[chan]->fSgAndBgHiE2) << endl;


    fTEX << "% ----------------------------------------------------------------------" << endl;
    fTEX << "% --- scaled/expected background summary numbers (new numbers)" << endl;  



  }


  OUT.close();

  fSuffix = cache; 

}


// ----------------------------------------------------------------------
void plotResults::printCsBFNumbers() {

  string cache = fSuffix; 
  fSuffix = (fDoUseBDT? "bdt":"cnc") + fSuffix; 
  
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
double  plotResults::scaledYield(numbers *a, numbers *no, string chan, double lfsfu) {

  bool verbose(true); 

  double chanbf = fBF[chan]; 

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
  double yield  = (chanbf/6.0e-5) * (lfsfu) * (a->effTot/no->effTot) * no->fitYield; 

  if (string::npos != a->name.find("signal Bs2MuMu")) {
    cout << "SCALING SIGNAL BS2MUMU" << endl;
    a->bsNoScaled  = a->pss * yield;
    a->bsNoScaledE = relError*a->bsNoScaled; 

    // -- new numbers
    a->fSgLo   = a->pls * yield;
    a->fSgLoE1 = TMath::Sqrt(a->fSgLo); 
    a->fSgLoE2 = relError*a->fSgLo; 

    a->fSgBd   = a->pds * yield;
    a->fSgBdE1 = TMath::Sqrt(a->fSgBd); 
    a->fSgBdE2 = relError*a->fSgBd; 

    a->fSgBs   = a->pss * yield;
    a->fSgBsE1 = TMath::Sqrt(a->fSgBs); 
    a->fSgBsE2 = relError*a->fSgBs; 

    a->fSgHi   = a->phs * yield;
    a->fSgHiE1 = TMath::Sqrt(a->fSgHi); 
    a->fSgHiE2 = relError*a->fSgHi; 
  } else if (string::npos != a->name.find("signal Bd2MuMu")) {
    cout << "SCALING SIGNAL BD2MUMU for numbers: " << a->name << " with ";
    a->bdNoScaled  = a->pdd * yield;
    a->bdNoScaledE = relError*a->bdNoScaled; 

    // -- new numbers
    a->fBdLo   = a->pld * yield;
    a->fBdLoE1 = TMath::Sqrt(a->fBdLo); 
    a->fBdLoE2 = relError*a->fBdLo; 

    a->fBdBd   = a->pdd * yield;
    a->fBdBdE1 = TMath::Sqrt(a->fBdBd); 
    a->fBdBdE2 = relError*a->fBdBd; 

    a->fBdBs   = a->psd * yield;
    a->fBdBsE1 = TMath::Sqrt(a->fBdBs); 
    a->fBdBsE2 = relError*a->fBdBs; 

    a->fBdHi   = a->phd * yield;
    a->fBdHiE1 = TMath::Sqrt(a->fBdHi); 
    a->fBdHiE2 = relError*a->fBdHi; 

    cout << "lo/Bd/Bs/hi: " << a->fBdLo << "/" << a->fBdBd << "/" << a->fBdBs << "/" << a->fBdHi << endl;
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
void plotResults::numbersFromHist(int chan, int mode, numbers *aa) {
  // -- efficiency and acceptance
  string modifier = (fDoUseBDT?"bdt":"cnc");
  TH1D *hAcceptance              = (TH1D*)fHistFile->Get(Form("hGenAndAccNumbers_%s_%d_chan%d", modifier.c_str(), mode, chan)); 

  TH1D *hMassAbsNoCuts           = (TH1D*)fHistFile->Get(Form("hMassAbsNoCuts_%s_%d_chan%d", modifier.c_str(), mode, chan));
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

  double effFilter    = getValueByLabel(hAcceptance, "effFilter");
  if (effFilter < 1e-6) {
    cout << "resetting effFilter to 1" << endl;
    effFilter = 1.0;
  }
  double genFileYield = getValueByLabel(hAcceptance, "genFileYield");

  aa->effPtReco       = getValueByLabel(hAcceptance, "effPtReco");
  aa->effPtRecoE      = getValueByLabel(hAcceptance, "effPtRecoE");

  aa->genAccFileYield = getValueByLabel(hAcceptance, "genAccFileYield");
  aa->genAccYield     = getValueByLabel(hAcceptance, "genAccYield");
  aa->recoYield       = getValueByLabel(hAcceptance, "recoYield");
  aa->muidYield       = getValueByLabel(hAcceptance, "muidYield");
  aa->trigYield       = getValueByLabel(hAcceptance, "trigYield");
  aa->candYield       = getValueByLabel(hAcceptance, "candYield");

  aa->acc             = getValueByLabel(hAcceptance, "acc");
  aa->accE            = getValueByLabel(hAcceptance, "accE");

  aa->effCand         = getValueByLabel(hAcceptance, "effCand");
  aa->effCandE        = getValueByLabel(hAcceptance, "effCandE");


  double a = hMassNoCuts->Integral(0, hMassNoCuts->GetNbinsX()+1); 
  double b = hMassWithAnaCuts->Integral(0, hMassWithAnaCuts->GetNbinsX()+1);
  double c = hMassWithMuonCuts->Integral(0, hMassWithMuonCuts->GetNbinsX()+1);
  double d = hMassWithTriggerCuts->Integral(0, hMassWithTriggerCuts->GetNbinsX()+1);
  double e = hMassWithAllCuts->Integral(0, hMassWithAllCuts->GetNbinsX()+1);
  double f = hMassWithMassCuts->Integral(0, hMassWithMassCuts->GetNbinsX()+1);
  aa->absNoCutsYield   = hMassAbsNoCuts->Integral(0, hMassAbsNoCuts->GetNbinsX()+1); 
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
  aa->effAna           = b/a;
  aa->effAnaE          = dEff(static_cast<int>(b), static_cast<int>(a)); 
  aa->effMuidMC        = c/b;
  cout << "*************************** name: " << aa->name << " effMuidMC: " << aa->effMuidMC << endl;
  aa->effMuidMCE       = dEff(static_cast<int>(c), static_cast<int>(b));
  aa->effMuidTNP       = hMuId->GetMean();
  aa->effMuidTNPE      = hMuId->GetMeanError();
  aa->effMuidTNPMC     = hMuIdMC->GetMean();
  aa->effMuidTNPMCE    = hMuIdMC->GetMeanError();
  // -- rho factor and application
  aa->rhoMuidTNPMC     = aa->effMuidMC/aa->effMuidTNPMC;
  aa->rhoMuidTNPMCE    = dRatio(aa->effMuidMC, aa->effMuidMCE, aa->effMuidTNPMC, aa->effMuidTNPMCE); 
  aa->effRhoMuidTNPMC  = aa->effMuidTNPMC*aa->rhoMuidTNPMC;
  aa->effRhoMuidTNPMCE = aa->effMuidTNPMCE*aa->rhoMuidTNPMC; // FIXME?
  aa->effRhoMuidTNP    = aa->effMuidTNP*aa->rhoMuidTNPMC;
  aa->effRhoMuidTNPE   = aa->effMuidTNPE*aa->rhoMuidTNPMC; // FIXME?


  aa->effTrigMC        = d/c;
  aa->effTrigMCE       = dEff(static_cast<int>(d), static_cast<int>(c));
  aa->effTrigTNP       = hMuTr->GetMean();
  aa->effTrigTNPE      = hMuTr->GetMeanError();
  aa->effTrigTNPMC     = hMuTrMC->GetMean();
  aa->effTrigTNPMCE    = hMuTrMC->GetMeanError();
  aa->genFileYield     = genFileYield;
  aa->effGenFilter     = effFilter; 
  aa->genYield         = aa->genFileYield/aa->effGenFilter;
  aa->effTot           = e/aa->genYield;
  aa->effTotE          = dEff(static_cast<int>(e), static_cast<int>(aa->genYield));

  // -- rho factor and application
  aa->rhoTrigTNPMC     = aa->effTrigMC/aa->effTrigTNPMC;
  aa->rhoTrigTNPMCE    = dRatio(aa->effTrigMC, aa->effTrigMCE, aa->effTrigTNPMC, aa->effTrigTNPMCE); 
  aa->effRhoTrigTNPMC  = aa->effTrigTNPMC*aa->rhoTrigTNPMC;
  aa->effRhoTrigTNPMCE = aa->effTrigTNPMCE*aa->rhoTrigTNPMC; // FIXME?
  aa->effRhoTrigTNP    = aa->effTrigTNP*aa->rhoTrigTNPMC;
  aa->effRhoTrigTNPE   = aa->effTrigTNPE*aa->rhoTrigTNPMC; // FIXME?

  // -- prod efficiencies
  aa->effProdMC        = aa->acc * aa->effCand * aa->effAna * aa->effMuidMC * aa->effTrigMC;
  aa->effProdMCE       = aa->effProdMC * 
    quadraticSum(4, aa->accE/aa->acc, aa->effAnaE/aa->effAna, aa->effMuidMCE/aa->effMuidMC, aa->effTrigMCE/aa->effTrigMC);
  aa->effProdTNP       = aa->acc * aa->effCand * aa->effAna * aa->effMuidTNP * aa->effTrigTNP;
  aa->effProdTNPE      = aa->effProdTNP * 
    quadraticSum(4, aa->accE/aa->acc, aa->effAna/aa->effAnaE, aa->effMuidMCE/aa->effMuidMC, aa->effTrigTNPE/aa->effTrigTNP);
  aa->effProdTNPMC     = aa->acc * aa->effCand * aa->effAna * aa->effMuidTNPMC * aa->rhoMuidTNPMC * aa->effTrigTNPMC  * aa->rhoTrigTNPMC;
  aa->effProdTNPMCE    = aa->effProdTNPMC *
    quadraticSum(4, aa->accE/aa->acc, aa->effAnaE/aa->effAna, aa->effMuidMCE/aa->effMuidMC, aa->effTrigTNPMCE/aa->effTrigTNPMC);


  aa->combGenYield     = e/(aa->effProdMC);
  aa->prodGenYield     = e/(aa->effTot); 


  if (0 == mode) {
    double tot   = hMassWithAllCutsManyBins->Integral(0, hMassWithAllCutsManyBins->GetNbinsX()+1);
    double lo    = hMassWithAllCutsManyBins->Integral(hMassWithAllCutsManyBins->FindBin(4.9), 
						      hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBdLo));
    double hi    = hMassWithAllCutsManyBins->Integral(hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBsHi), 
						      hMassWithAllCutsManyBins->FindBin(5.9));
    double bd    = hMassWithAllCutsManyBins->Integral(hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBdLo), 
						      hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBdHi));
    double bs    = hMassWithAllCutsManyBins->Integral(hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBsLo), 
						      hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBsHi));
    
    aa->pss = bs/tot;
    aa->pssE = dEff(static_cast<int>(bs), static_cast<int>(tot)); 
    aa->pds = bd/tot;
    aa->pdsE = dEff(static_cast<int>(bd), static_cast<int>(tot)); 
    aa->pls = lo/tot;
    aa->plsE = dEff(static_cast<int>(lo), static_cast<int>(tot)); 
    aa->phs = hi/tot;
    aa->phsE = dEff(static_cast<int>(hi), static_cast<int>(tot)); 
  } 
    
  if (1 == mode) {
    double tot   = hMassWithAllCutsManyBins->Integral(0, hMassWithAllCutsManyBins->GetNbinsX()+1);
    double lo    = hMassWithAllCutsManyBins->Integral(hMassWithAllCutsManyBins->FindBin(4.9), 
						      hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBdLo));
    double hi    = hMassWithAllCutsManyBins->Integral(hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBsHi), 
						      hMassWithAllCutsManyBins->FindBin(5.9));
    double bd    = hMassWithAllCutsManyBins->Integral(hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBdLo), 
						      hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBdHi));
    double bs    = hMassWithAllCutsManyBins->Integral(hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBsLo), 
						      hMassWithAllCutsManyBins->FindBin(fCuts[chan]->mBsHi));
    
    aa->pdd  = bd/tot;
    aa->pddE = dEff(static_cast<int>(bd), static_cast<int>(tot)); 
    aa->psd  = bs/tot;
    aa->psdE = dEff(static_cast<int>(bs), static_cast<int>(tot)); 
    aa->pld  = lo/tot;
    aa->pldE = dEff(static_cast<int>(lo), static_cast<int>(tot)); 
    aa->phd  = hi/tot;
    aa->phdE = dEff(static_cast<int>(hi), static_cast<int>(tot)); 
  } 

  printNumbers(*aa, cout); 
    
}


// ----------------------------------------------------------------------
void plotResults::rareBgHists(std::string mode, int nevents) {

  fSetup = "BgRareMc"; 

  TH1D *h1(0); 
  TTree *t(0); 
  string smode("cnc"); 
  if (fDoUseBDT) smode = "bdt"; 

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

    t = (TTree*)(fF[imap->first]->Get("candAnaMuMu/events"));
    cout << "======> " << imap->first  << " with the tree at " << t << endl;

    fRareName = imap->first; 
    resetHistograms();
    setupTree(t, fSetup); 
    loopOverTree(t, fSetup, 1, nevents);

    for (unsigned int ichan = 0; ichan < fNchan; ++ichan) {
      h1 = (TH1D*)(fhMassWithAllCutsManyBins[ichan]->Clone(Form("hMassWithAllCutsManyBins_%s_%s_chan%d", 
								smode.c_str(), fRareName.c_str(), ichan)));  
      h1->SetTitle(Form("hMassWithAllCutsManyBins_%s_%s_chan%d", smode.c_str(), fRareName.c_str(), ichan));
      h1->SetDirectory(fHistFile); 
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCuts[ichan]->Clone(Form("hMassWithAllCuts_%s_%s_chan%d", smode.c_str(), fRareName.c_str(), ichan)));  
      h1->SetTitle(Form("hMassWithAllCuts_%s_%s_chan%d", smode.c_str(), fRareName.c_str(), ichan));
      h1->SetDirectory(fHistFile); 
      h1->Write();

      h1 = (TH1D*)(fhW8MassWithAllCutsManyBins[ichan]->Clone(Form("hW8MassWithAllCutsManyBins_%s_%s_chan%d", 
								smode.c_str(), fRareName.c_str(), ichan)));  
      h1->SetTitle(Form("hW8MassWithAllCutsManyBins_%s_%s_chan%d", smode.c_str(), fRareName.c_str(), ichan));
      h1->SetDirectory(fHistFile); 
      h1->Write();

      h1 = (TH1D*)(fhW8MassWithAllCuts[ichan]->Clone(Form("hW8MassWithAllCuts_%s_%s_chan%d", smode.c_str(), fRareName.c_str(), ichan)));  
      h1->SetTitle(Form("hW8MassWithAllCuts_%s_%s_chan%d", smode.c_str(), fRareName.c_str(), ichan));
      h1->SetDirectory(fHistFile); 
      h1->Write();
    }
  }

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

  for (unsigned int i = 0; i < processes.size(); ++i) {
    int mode = processes[i];

    if (0 == mode) fF["SgMcAcc"]->cd();
    else if (10 == mode) fF["NoMcAcc"]->cd();
    else break;

    for (int chan = 0; chan < 2; ++chan) {

      fTEX << "% -- acceptancePerProcess: " << mode << " " << chan << endl;
      
      string cuts; 
      
      if (0 == mode) {
	if (0 == chan) {
	  cuts = "g1pt>3.5&&g2pt>3.5&&abs(g1eta)<1.4&&abs(g2eta)<1.4";
	  cuts += "&&m1pt>3.5&&m2pt>3.5&&abs(m1eta)<1.4&&abs(m2eta)<1.4&&m1gt&&m2gt";
	} else {
	  cuts = "g1pt>3.5&&g2pt>3.5&&(abs(g1eta)>1.4||abs(g2eta)>1.4)&&abs(g1eta)<2.5&&abs(g2eta)<2.5";
	  cuts += "&&m1pt>3.5&&m2pt>3.5&&(abs(m1eta)>1.4||abs(m2eta)>1.4)&&abs(m1eta)<2.4&&abs(m2eta)<2.4&&m1gt&&m2gt";
	}
      } else if (10 == mode) {
	if (0 == chan) {
	  cuts = "g1pt>3.5&&g2pt>3.5&&g3pt>0.4&&abs(g1eta)<1.4&&abs(g2eta)<1.4&&abs(g3eta)<2.5";
	  cuts += "&&m1pt>3.5&&m2pt>3.5&&k1pt>0.5&&abs(m1eta)<1.4&&abs(m2eta)<1.4&&abs(k1eta)<2.4&&m1gt&&m2gt&&k1gt";
	} else {
	  cuts = "g1pt>3.5&&g2pt>3.5&&g3pt>0.4&&(abs(g1eta)>1.4||abs(g2eta)>1.4)&&abs(g1eta)<2.5&&abs(g2eta)<2.5&&abs(g3eta)<2.5";
	  cuts += "&&m1pt>3.5&&m2pt>3.5&&k1pt>0.5&&(abs(m1eta)>1.4||abs(m2eta)>1.4)&&abs(m1eta)<2.4&&abs(m2eta)<2.4&&abs(k1eta)<2.4&&m1gt&&m2gt&&k1gt";
	}
      }
      
      TTree *t(0); 
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

  fNumbersNo[0]->effTot = 0.00166;
  fNumbersNo[1]->effTot = 0.00037;
  fNumbersNo[0]->fitYield = 129806;
  fNumbersNo[1]->fitYield = 25578;

}
