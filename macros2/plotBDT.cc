#include "plotBDT.hh"

#include <TMath.h>
#include <TGaxis.h>

using namespace std; 

ClassImp(plotBDT)

// ----------------------------------------------------------------------
plotBDT::plotBDT(const char *files, const char *dir, const char *cuts, int mode) : plotClass(files, dir, cuts, mode) { 
  fNumbersFileName = fDirectory + "/anaBmm.plotBDT." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  cout << "==> plotBDT files: " << files << " dir: " << dir << " cuts: " << cuts << endl;

  string hfname  = fDirectory + "/anaBmm.plotBDT." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  //  if (fHistFile) fHistFile->Close();
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

  fIsMC = false;

  fDoUseBDT = true; 
  fSaveSmallTree = false; 

  int NBINS = static_cast<int>((fMassHi - fMassLo)/0.025);
  TH1D *h(0); 
  TProfile *p(0); 
  for (unsigned int i = 0; i < fNchan; ++i) {
    h = new TH1D(Form("hMassNoCuts%d", i), Form("hMassNoCuts%d", i), NBINS, fMassLo, fMassHi);
    setTitles(h, "m [GeV]", "mean(b)", 0.05, 1.2, 1.5);
    fhMassNoCuts.push_back(h); 

    h = new TH1D(Form("hMass%d", i), Form("hMass%d", i), NBINS, fMassLo, fMassHi);
    setTitles(h, "m [GeV]", "mean(b)", 0.05, 1.2, 1.5);
    fhMass.push_back(h); 
    
    // -- BDT output distributions for all, lo, hi, and signal mass ranges
    h = new TH1D(Form("hBDT%d", i), Form("hBDT%d", i), 200, -1., 1.);
    fhBDT.push_back(h); 

    h = new TH1D(Form("hLoBDT%d", i), Form("hLoBDT%d", i), 200, -1., 1.);
    fhLoBDT.push_back(h); 

    h = new TH1D(Form("hHiBDT%d", i), Form("hHiBDT%d", i), 200, -1., 1.);
    fhHiBDT.push_back(h); 

    h = new TH1D(Form("hInBDT%d", i), Form("hInBDT%d", i), 200, -1., 1.);
    fhInBDT.push_back(h); 


    p = new TProfile(Form("pMassBDT%d", i), Form("pMassBDT%d", i), NBINS/2, fMassLo, fMassHi);
    fpMassBDT.push_back(p); 
    p = new TProfile(Form("pMassAcBDT%d", i), Form("pMassAcBDT%d", i), NBINS/2, fMassLo, fMassHi);
    fpMassAcBDT.push_back(p); 
    p = new TProfile(Form("pMassAdBDT%d", i), Form("pMassAdBDT%d", i), NBINS/2, fMassLo, fMassHi);
    fpMassAdBDT.push_back(p); 

    p = new TProfile(Form("pNpvBDT%d", i), Form("pNpvBDT%d", i), 30, 0., 30.);
    fpNpvBDT.push_back(p); 
    p = new TProfile(Form("pNpvAcBDT%d", i), Form("pNpvAcBDT%d", i), 30, 0., 30.);
    fpNpvAcBDT.push_back(p);    
    p = new TProfile(Form("pNpvAdBDT%d", i), Form("pNpvAdBDT%d", i), 30, 0., 30.);
    fpNpvAdBDT.push_back(p);    
  }


  fmeanNpvBDTchan0   = new TH1D("fmeanNpvBDTchan0", "fmeanNpvBDTchan0", 50, 0., 100.);     
  setTitles(fmeanNpvBDTchan0, "N_{PV}", "mean(b)", 0.05, 1.2, 1.5);
  fmeanNpvBDTchan1   = new TH1D("fmeanNpvBDTchan1", "fmeanNpvBDTchan0", 50, 0., 100.);     
  setTitles(fmeanNpvBDTchan1, "N_{PV}", "mean(b)", 0.05, 1.2, 1.5);
  fmeanNpvAcBDTchan0 = new TH1D("fmeanNpvAcBDTchan0", "fmeanNpvAcBDTchan0", 50, 0., 100.); 
  setTitles(fmeanNpvAcBDTchan0, "N_{PV}", "mean(b)", 0.05, 1.2, 1.5);
  fmeanNpvAcBDTchan1 = new TH1D("fmeanNpvAcBDTchan1", "fmeanNpvAcBDTchan0", 50, 0., 100.); 
  setTitles(fmeanNpvAcBDTchan1, "N_{PV}", "mean(b)", 0.05, 1.2, 1.5);
  fmeanNpvAdBDTchan0 = new TH1D("fmeanNpvAdBDTchan0", "fmeanNpvAdBDTchan0", 50, 0., 100.); 
  setTitles(fmeanNpvAdBDTchan0, "N_{PV}", "mean(b)", 0.05, 1.2, 1.5);
  fmeanNpvAdBDTchan1 = new TH1D("fmeanNpvAdBDTchan1", "fmeanNpvAdBDTchan0", 50, 0., 100.); 
  setTitles(fmeanNpvAdBDTchan1, "N_{PV}", "mean(b)", 0.05, 1.2, 1.5);

  frmsNpvBDTchan0   = new TH1D("frmsNpvBDTchan0", "frmsNpvBDTchan0", 50, 0., 100.);        
  setTitles(frmsNpvBDTchan0, "N_{PV}", "rms(b)", 0.05, 1.2, 1.5);  
  frmsNpvBDTchan1   = new TH1D("frmsNpvBDTchan1", "frmsNpvBDTchan0", 50, 0., 100.); 	   
  setTitles(frmsNpvBDTchan1, "N_{PV}", "rms(b)", 0.05, 1.2, 1.5);  
  frmsNpvAcBDTchan0 = new TH1D("frmsNpvAcBDTchan0", "frmsNpvAcBDTchan0", 50, 0., 100.);    
  setTitles(frmsNpvAcBDTchan0, "N_{PV}", "rms(b)", 0.05, 1.2, 1.5);
  frmsNpvAcBDTchan1 = new TH1D("frmsNpvAcBDTchan1", "frmsNpvAcBDTchan0", 50, 0., 100.);    
  setTitles(frmsNpvAcBDTchan1, "N_{PV}", "rms(b)", 0.05, 1.2, 1.5);
  frmsNpvAdBDTchan0 = new TH1D("frmsNpvAdBDTchan0", "frmsNpvAdBDTchan0", 50, 0., 100.);    
  setTitles(frmsNpvAdBDTchan0, "N_{PV}", "rms(b)", 0.05, 1.2, 1.5);
  frmsNpvAdBDTchan1 = new TH1D("frmsNpvAdBDTchan1", "frmsNpvAdBDTchan0", 50, 0., 100.);    
  setTitles(frmsNpvAdBDTchan1, "N_{PV}", "rms(b)", 0.05, 1.2, 1.5);

  for (unsigned int pv = 0; pv < 50; ++pv) {
    h = new TH1D(Form("hNpvBDTChan0_%d", pv), Form("hNpvBDTChan0_%d", pv), 200, -1., 1.);
    fhNpvBDTchan0.push_back(h);
    h = new TH1D(Form("hNpvBDTChan1_%d", pv), Form("hNpvBDTChan1_%d", pv), 200, -1., 1.);
    fhNpvBDTchan1.push_back(h);

    h = new TH1D(Form("hNpvAcBDTChan0_%d", pv), Form("hNpvAcBDTChan0_%d", pv), 200, -1., 1.);
    fhNpvAcBDTchan0.push_back(h);
    h = new TH1D(Form("hNpvAcBDTChan1_%d", pv), Form("hNpvAcBDTChan1_%d", pv), 200, -1., 1.);
    fhNpvAcBDTchan1.push_back(h);

    h = new TH1D(Form("hNpvAdBDTChan0_%d", pv), Form("hNpvAdBDTChan0_%d", pv), 200, -1., 1.);
    fhNpvAdBDTchan0.push_back(h);
    h = new TH1D(Form("hNpvAdBDTChan1_%d", pv), Form("hNpvAdBDTChan1_%d", pv), 200, -1., 1.);
    fhNpvAdBDTchan1.push_back(h);
  }


  double xmin(-0.2), xmax(0.6); 
  int nbins; 
  string htitle;
  for (unsigned int i = 0; i < 3; ++i) {
    // -- hackedMC1
    nbins = static_cast<int>((xmax - xmin)/0.025);

    h = new TH1D(Form("hMcBDT%d", i), Form("hMcBDT%d", i), nbins, xmin, xmax);
    if (0 == i) setHist(h, kBlack, 24); 
    if (1 == i) setHist(h, kRed, 25); 
    if (2 == i) setHist(h, kBlue, 26); 
    fhMcBDT.push_back(h); 
    h = new TH1D(Form("hMcMass%d", i), Form("hMcMass%d", i), 100, 4.5, 6.0);
    if (0 == i) setHist(h, kBlack, 24); 
    if (1 == i) setHist(h, kRed, 25); 
    if (2 == i) setHist(h, kBlue, 26); 
    fhMcMass.push_back(h); 
    h = new TH1D(Form("hMcRatio%d", i), Form("hMcRatio%d", i), nbins, xmin, xmax);
    h->Sumw2();
    if (0 == i) setHist(h, kBlack, 20); 
    if (1 == i) setHist(h, kRed, 20); 
    if (2 == i) setHist(h, kBlue, 20); 
    fhMcRatio.push_back(h); 

    // -- hackedMC2

    if (0 == i) htitle = "B_{s} #rightarrow #mu #mu";
    if (1 == i) htitle = "B_{x} #rightarrow #mu #mu";
    if (2 == i) htitle = "B_{y} #rightarrow #mu #mu";
    int color(kBlack); 
    if (1 == i) color = kRed; 
    if (2 == i) color = kBlue; 
    nbins = 200; 
    h = new TH1D(Form("hMcMass0_%d", i), htitle.c_str(), nbins, 4.5, 6.0); setHist(h, color, 20); 
    fhMcMass0.push_back(h); 
    h = new TH1D(Form("hMcMass1_%d", i), htitle.c_str(), nbins, 4.5, 6.0); setHist(h, color, 20); 
    fhMcMass1.push_back(h); 
    h = new TH1D(Form("hMcMass2_%d", i), htitle.c_str(), nbins, 4.5, 6.0); setHist(h, color, 24); 
    fhMcMass2.push_back(h); 

    h = new TH1D(Form("hMcRatio1_%d", i), "ratio 1/0", nbins, 4.5, 6.0); h->Sumw2(); setHist(h, color, 20); 
    fhMcRatio1.push_back(h); 
    h = new TH1D(Form("hMcRatio2_%d", i), "ratio 2/0", nbins, 4.5, 6.0); h->Sumw2(); setHist(h, color, 24); 
    fhMcRatio2.push_back(h); 

  }


  nbins = static_cast<int>((xmax - xmin)/0.025);
  h = new TH1D("hMcBdt0", "Bs2MuMu", nbins, xmin, xmax);     h->Sumw2();  setHist(h, kBlack);               fhMcBDT5.push_back(h); 
  h = new TH1D("hMcBdt1", "Bx2MuMu", nbins, xmin, xmax);     h->Sumw2();  setHist(h, kBlue);                fhMcBDT5.push_back(h); 
  h = new TH1D("hMcBdt2", "By2MuMu", nbins, xmin, xmax);     h->Sumw2();  setHist(h, kRed);                 fhMcBDT5.push_back(h); 
  h = new TH1D("hMcBdt3", "Bs2JpsiPhi", nbins, xmin, xmax);  h->Sumw2();  setFilledHist(h, 35, 35, 3365);   fhMcBDT5.push_back(h); 
  h = new TH1D("hMcBdt4", "Bu2JpsiK", nbins, xmin, xmax);    h->Sumw2();  setFilledHist(h, 45, 45, 3356);   fhMcBDT5.push_back(h);  

}


// ----------------------------------------------------------------------
plotBDT::~plotBDT() { 
 
  if (fHistFile) {
    fHistFile->Write();
    fHistFile->Close();
  }
}

// ----------------------------------------------------------------------
void plotBDT::makeAll(int channels) { 

  if (channels & 1) {
    //    npvSpecial("SgData");
    //    npvSpecial("NoData");


    cout << "--> tmvaControlPlots()" << endl;
    tmvaControlPlots();
    cout << "--> plotSSB()" << endl;
    plotSSB();
    cout << "--> overlayBdtOutput()" << endl;
    overlayBdtOutput(); 
    cout << "--> validateAllDistributions()" << endl;
    validateAllDistributions();
    cout << "--> allCorrelationPlots(...)" << endl;
    allCorrelationPlots(0., fCuts[0]->xmlFile.c_str());
    allCorrelationPlots(0., fCuts[1]->xmlFile.c_str());
    cout << "--> bdtScan()" << endl;
    bdtScan();
    cout << "--> plotEffVsBg(...)" << endl;
    for (int i = 0; i < fNchan; ++i) {
      fBdtString = Form("%s", fCuts[i]->xmlFile.c_str());
      plotEffVsBg(i);   
    }

    cout << "--> bdtDependencies(\"SgData\")" << endl;
    bdtDependencies("SgData");

    cout << "--> xmlParsing()" << endl;
    xmlParsing(); 

    cout << "--> hackedMC(...)" << endl;
    if (0) {
      hackedMC(0); 
      hackedMC(1); 
    }
  }

  if (channels & 2) {
    hackedMC2(0.0, 0.2); 
    hackedMC2(0.0, 0.3); 
    hackedMC1("abs(m1eta)<1.4&&abs(m2eta)<1.4", -0.3, 0.5, "chan0");
    hackedMC1("!(abs(m1eta)<1.4&&abs(m2eta)<1.4)", -0.3, 0.5, "chan1");
    hackedMC3(0); 
    hackedMC3(1); 
  }

  if (channels & 8) {
    illustrateAntiMuonSample();
  }


}



// ----------------------------------------------------------------------
void plotBDT::loopFunction(int function, int mode) {
  if (1 == function) loopFunction1(mode);
  if (2 == function) loopFunction2(mode);
}


// ----------------------------------------------------------------------
void plotBDT::loopFunction1(int mode) {

  if (fChan < 0) return;

  // -- gen-level acceptance cuts
  if (fIsMC) {
    if (TMath::Abs(fb.g1eta) > 2.5) return;
    if (TMath::Abs(fb.g2eta) > 2.5) return;
    if (fb.g1pt < 1.0) return;
    if (fb.g2pt < 1.0) return;
  }

  // -- immutable cuts: require basic muon and trackQual cuts
  if (false == fb.gtqual)  return;
  if (TMath::Abs(fb.m1eta) > 2.4) return;
  if (TMath::Abs(fb.m2eta) > 2.4) return;
  
  if (fb.m1pt < 1.0) return;
  if (fb.m2pt < 1.0) return;
  
  // -- no cuts for normalization
  fhMassNoCuts[fChan]->Fill(fb.m); 

  if (fBDT > fCuts[fChan]->bdt) {
    fhMass[fChan]->Fill(fb.m); 
  }

  // -- the good events must also pass GMUID and HLT!
  //  if (fb.gmuid && fb.hlt && fb.gtqual) {
  if (!fGoodAcceptance) return;
  if (!fGoodTracks) return;
  if (!fGoodTracksPt) return;
  if (!fGoodTracksEta) return;
  if (!fGoodMuonsPt) return;
  if (!fGoodMuonsEta) return;
  if (!fGoodJpsiCuts) return;
  if (!fGoodMuonsID) return;
  if (!fGoodHLT) return;

  fhBDT[fChan]->Fill(fBDT);
  if (fb.m > 5.45) fhHiBDT[fChan]->Fill(fBDT);
  if (fb.m < 5.20) fhLoBDT[fChan]->Fill(fBDT);
  if (fb.m > 5.20 && fb.m < 5.45) fhInBDT[fChan]->Fill(fBDT); 
  
  if (fBDT > -1.)  fpMassBDT[fChan]->Fill(fb.m, fBDT);
  if (fBDT > 0.00) fpMassAcBDT[fChan]->Fill(fb.m, fBDT);
  if (fBDT > 0.10) fpMassAdBDT[fChan]->Fill(fb.m, fBDT);
  
  if (fBDT > -1.)  fpNpvBDT[fChan]->Fill(fb.pvn, fBDT);
  if (fBDT > 0.00) fpNpvAcBDT[fChan]->Fill(fb.pvn, fBDT);
  if (fBDT > 0.10) fpNpvAdBDT[fChan]->Fill(fb.pvn, fBDT);

  int pvbin = fb.pvn/2; 
  if (pvbin > 49) pvbin = 49; 
  if (0 == fChan) {
    if (fBDT > -1.)  fhNpvBDTchan0[pvbin]->Fill(fBDT); 
    if (fBDT > 0.00) fhNpvAcBDTchan0[pvbin]->Fill(fBDT); 
    if (fBDT > 0.10) fhNpvAdBDTchan0[pvbin]->Fill(fBDT); 
  } 
  if (1 == fChan) {
    if (fBDT > -1.)  fhNpvBDTchan1[pvbin]->Fill(fBDT); 
    if (fBDT > 0.00) fhNpvAcBDTchan1[pvbin]->Fill(fBDT); 
    if (fBDT > 0.10) fhNpvAdBDTchan1[pvbin]->Fill(fBDT); 
  } 

}


// ----------------------------------------------------------------------
void  plotBDT::resetHistograms() {
  for (unsigned int i = 0; i < fNchan; ++i) {
    fhBDT[i]->Reset();
    fhLoBDT[i]->Reset();
    fhInBDT[i]->Reset();
    fhHiBDT[i]->Reset();

    fhMass[i]->Reset();
    fhMassNoCuts[i]->Reset();

    fpMassBDT[i]->Reset();
    fpMassAcBDT[i]->Reset();
    fpMassAdBDT[i]->Reset();
    fpNpvBDT[i]->Reset();
    fpNpvAcBDT[i]->Reset();
    fpNpvAdBDT[i]->Reset();

  }    

  fmeanNpvBDTchan0->Reset();  
  fmeanNpvBDTchan1->Reset();  
  fmeanNpvAcBDTchan0->Reset();
  fmeanNpvAcBDTchan1->Reset();
  fmeanNpvAdBDTchan0->Reset();
  fmeanNpvAcBDTchan1->Reset();

  frmsNpvBDTchan0->Reset();  
  frmsNpvBDTchan1->Reset();  
  frmsNpvAcBDTchan0->Reset();
  frmsNpvAcBDTchan1->Reset();
  frmsNpvAdBDTchan0->Reset();
  frmsNpvAcBDTchan1->Reset();

  for (unsigned int i = 0; i < 50; ++i) {
    fhNpvBDTchan0[i]->Reset();
    fhNpvBDTchan1[i]->Reset();

    fhNpvAcBDTchan0[i]->Reset();
    fhNpvAcBDTchan1[i]->Reset();

    fhNpvAdBDTchan0[i]->Reset();
    fhNpvAdBDTchan1[i]->Reset();
  }
    
  for (unsigned int i = 0; i < 3; ++i) {
    fhMcMass[i]->Reset();
    fhMcBDT[i]->Reset();
    fhMcRatio[i]->Reset();

    fhMcMass0[i]->Reset();
    fhMcMass1[i]->Reset();
    fhMcMass2[i]->Reset();
    fhMcRatio1[i]->Reset();
    fhMcRatio2[i]->Reset();
  }

  for (unsigned int i = 0; i < 3; ++i) {
    fhMcBDT5[i]->Reset();
  }

}

// ----------------------------------------------------------------------
void plotBDT::testLoop(string mode) { 

}



// ----------------------------------------------------------------------
void  plotBDT::overlayBdtOutput() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  vector<string> type; 
  type.push_back("sg"); 
  type.push_back("bg"); 

  TH1D *hte0(0), *htr0(0), *hte1(0), *htr1(0), *hte2(0), *htr2(0); 
  TH1D *hap0(0), *hap1(0), *hap2(0);
  for (unsigned int j = 0; j < type.size(); ++j) {
    for (unsigned int i = 0; i < fNchan; ++i) {
      string rootfile = "weights/" + fCuts[i]->xmlFile + "-combined.root";
      fRootFile = TFile::Open(rootfile.c_str());
      cout << "fRootFile: " << rootfile << endl;
      
      hte0 = (TH1D*)fRootFile->Get(Form("te0%sBDT", type[j].c_str())); 
      hte1 = (TH1D*)fRootFile->Get(Form("te1%sBDT", type[j].c_str())); 
      hte2 = (TH1D*)fRootFile->Get(Form("te2%sBDT", type[j].c_str())); 

      htr0 = (TH1D*)fRootFile->Get(Form("tr0%sBDT", type[j].c_str())); 
      htr1 = (TH1D*)fRootFile->Get(Form("tr1%sBDT", type[j].c_str())); 
      htr2 = (TH1D*)fRootFile->Get(Form("tr2%sBDT", type[j].c_str())); 

      hap0 = (TH1D*)fRootFile->Get(Form("ap0%sBDT", type[j].c_str())); 
      hap1 = (TH1D*)fRootFile->Get(Form("ap1%sBDT", type[j].c_str())); 
      hap2 = (TH1D*)fRootFile->Get(Form("ap2%sBDT", type[j].c_str())); 

      double hmax = htr0->GetMaximum(); 
      if (htr1->GetMaximum() > hmax) hmax = htr1->GetMaximum();
      if (htr2->GetMaximum() > hmax) hmax = htr2->GetMaximum();
      htr0->SetMaximum(1.1*hmax);
    
      SetSignalAndBackgroundStyle(htr0, htr1, htr2);            
      SetSignalAndBackgroundStyle(hte0, hte1, hte2);            
      setTitles(htr0, "b", ""); 

      c0->cd();
      shrinkPad(0.15, 0.15); 
      setTitles(htr0, "b", "candidates", 0.05, 1.2, 1.5); 
      htr0->Draw("p");
      hte0->Draw("same");
    
      htr1->Draw("psame");
      hte1->Draw("same");
    
      htr2->Draw("psame");
      hte2->Draw("same");
    
      double x = (j==0?0.25:0.50); 
      newLegend(x, 0.6, x+0.3, 0.85, (i==0?(j==0?"Barrel signal":"Barrel background"):(j==0)?"Endcap signal":"Endcap background")); 
      legg->AddEntry(htr0, "BDT 0 (train)", "p"); 
      legg->AddEntry(hte0, "BDT 0 (test) ", "l"); 
      legg->AddEntry(htr1, "BDT 1 (train)", "p"); 
      legg->AddEntry(htr1, "BDT 1 (test) ", "l"); 
      legg->AddEntry(htr2, "BDT 2 (train)", "p"); 
      legg->AddEntry(htr2, "BDT 2 (test) ", "l"); 
      legg->Draw();

      c0->SaveAs(Form("%s/%s-b-overlays-%s.pdf", fDirectory.c_str(), fCuts[i]->xmlFile.c_str(), type[j].c_str())); 

      SetSignalAndBackgroundStyle(hap0, hap1, hap2);            
      setTitles(hap0, "b", "candidates", 0.05, 1.2, 1.5); 
      hap0->Draw("");
      hap1->Draw("same");
      hap2->Draw("same");
      newLegend(x, 0.6, x+0.3, 0.85, (i==0?(j==0?"Barrel signal":"Barrel background"):(j==0)?"Endcap signal":"Endcap background")); 
      legg->AddEntry(hap0, "BDT 0", "p"); 
      legg->AddEntry(hap1, "BDT 1", "p"); 
      legg->AddEntry(hap2, "BDT 2", "p"); 
      legg->Draw();
      c0->SaveAs(Form("%s/%s-b-all-overlays-%s.pdf", fDirectory.c_str(), fCuts[i]->xmlFile.c_str(), type[j].c_str())); 

    }
  }

}




// ----------------------------------------------------------------------
void plotBDT::overlap() {


}


// ----------------------------------------------------------------------
void plotBDT::bdtDependencies(string mode) {

  //  string mode = "SgData";
  cout << "bdtDependencies: running over " << mode << endl;
  TTree *t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  resetHistograms();
  setupTree(t, mode);
  loopOverTree(t, mode, 1); 

  // -- rescale errors for low-stats entries
  for (unsigned int i = 0; i < fNchan; ++i) {
    for (int j = 1; j < fpMassAcBDT[i]->GetNbinsX(); ++j) {
      if (fpMassAcBDT[i]->GetBinEntries(j) < 5) fpMassAcBDT[i]->SetBinEntries(j, 0); 
    }

    for (int j = 1; j < fpMassBDT[i]->GetNbinsX(); ++j) {
      if (fpMassBDT[i]->GetBinEntries(j) < 5) fpMassBDT[i]->SetBinEntries(j, 0); 
    }

    for (int j = 1; j < fpNpvBDT[i]->GetNbinsX(); ++j) {
      if (fpNpvBDT[i]->GetBinEntries(j) < 5) fpNpvBDT[i]->SetBinEntries(j, 0); 
    }

    for (int j = 1; j < fpNpvAcBDT[i]->GetNbinsX(); ++j) {
      if (fpNpvAcBDT[i]->GetBinEntries(j) < 5) fpNpvAcBDT[i]->SetBinEntries(j, 0); 
    }

    for (int j = 1; j < fpMassAcBDT[i]->GetNbinsX(); ++j) {
      if (fpMassAcBDT[i]->GetBinEntries(j) < 5) fpMassAcBDT[i]->SetBinEntries(j, 0); 
    }

  }

  // -- produce summary plots
  int nmin(10); 
  for (int j = 0; j < 50; ++j) {
    if (fhNpvBDTchan0[j]->GetEntries() > nmin) {
      fmeanNpvBDTchan0->SetBinContent(j, fhNpvBDTchan0[j]->GetMean()); fmeanNpvBDTchan0->SetBinError(j, fhNpvBDTchan0[j]->GetMeanError());
    }
    if (fhNpvBDTchan1[j]->GetEntries() > nmin) {
      fmeanNpvBDTchan1->SetBinContent(j, fhNpvBDTchan1[j]->GetMean()); fmeanNpvBDTchan1->SetBinError(j, fhNpvBDTchan1[j]->GetMeanError());
    }

    if (fhNpvAcBDTchan0[j]->GetEntries() > nmin) {
      fmeanNpvAcBDTchan0->SetBinContent(j, fhNpvAcBDTchan0[j]->GetMean()); fmeanNpvAcBDTchan0->SetBinError(j, fhNpvAcBDTchan0[j]->GetMeanError());
    }
    if (fhNpvAcBDTchan1[j]->GetEntries() > nmin) {
      fmeanNpvAcBDTchan1->SetBinContent(j, fhNpvAcBDTchan1[j]->GetMean()); fmeanNpvAcBDTchan1->SetBinError(j, fhNpvAcBDTchan1[j]->GetMeanError());
    }

    if (fhNpvAdBDTchan0[j]->GetEntries() > nmin) {
      fmeanNpvAdBDTchan0->SetBinContent(j, fhNpvAdBDTchan0[j]->GetMean()); fmeanNpvAdBDTchan0->SetBinError(j, fhNpvAdBDTchan0[j]->GetMeanError());
    }
    if (fhNpvAdBDTchan1[j]->GetEntries() > nmin) {
      fmeanNpvAdBDTchan1->SetBinContent(j, fhNpvAdBDTchan1[j]->GetMean()); fmeanNpvAdBDTchan1->SetBinError(j, fhNpvAdBDTchan1[j]->GetMeanError());
    }

    if (fhNpvBDTchan0[j]->GetEntries() > nmin) {
      frmsNpvBDTchan0->SetBinContent(j, fhNpvBDTchan0[j]->GetRMS()); frmsNpvBDTchan0->SetBinError(j, fhNpvBDTchan0[j]->GetRMSError());
    }
    if (fhNpvBDTchan1[j]->GetEntries() > nmin) {
      frmsNpvBDTchan1->SetBinContent(j, fhNpvBDTchan1[j]->GetRMS()); frmsNpvBDTchan1->SetBinError(j, fhNpvBDTchan1[j]->GetRMSError());
    }

    if (fhNpvAcBDTchan0[j]->GetEntries() > nmin) {
      frmsNpvAcBDTchan0->SetBinContent(j, fhNpvAcBDTchan0[j]->GetRMS()); frmsNpvAcBDTchan0->SetBinError(j, fhNpvAcBDTchan0[j]->GetRMSError());
    }

    if (fhNpvAcBDTchan1[j]->GetEntries() > nmin) {
      frmsNpvAcBDTchan1->SetBinContent(j, fhNpvAcBDTchan1[j]->GetRMS()); frmsNpvAcBDTchan1->SetBinError(j, fhNpvAcBDTchan1[j]->GetRMSError());
    }

    if (fhNpvAdBDTchan0[j]->GetEntries() > nmin) {
      frmsNpvAdBDTchan0->SetBinContent(j, fhNpvAdBDTchan0[j]->GetRMS()); frmsNpvAdBDTchan0->SetBinError(j, fhNpvAdBDTchan0[j]->GetRMSError());
    }

    if (fhNpvAdBDTchan1[j]->GetEntries() > nmin) {
      frmsNpvAdBDTchan1->SetBinContent(j, fhNpvAdBDTchan1[j]->GetRMS()); frmsNpvAdBDTchan1->SetBinError(j, fhNpvAdBDTchan1[j]->GetRMSError());
    }
  }


  double npvmax(30.); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 
  gStyle->SetOptTitle(0); 
  tl->SetTextSize(0.03);
  for (int j = 0; j < 2; ++j) {
    zone(1);
    shrinkPad(0.15, 0.15); 
    // -- mean 
    fmeanNpvBDTchan0->SetAxisRange(0., npvmax, "X");
    fmeanNpvBDTchan0->Fit(Form("pol%d", j));  
    if (fmeanNpvBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  fmeanNpvBDTchan0->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  fmeanNpvBDTchan0->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (fmeanNpvBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  fmeanNpvBDTchan0->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  fmeanNpvBDTchan0->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvmean-nocuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[0]->xmlFile.c_str()));

    fmeanNpvAcBDTchan0->SetAxisRange(0., npvmax, "X");
    fmeanNpvAcBDTchan0->Fit(Form("pol%d", j));  
    if (fmeanNpvAcBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  fmeanNpvAcBDTchan0->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  fmeanNpvAcBDTchan0->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (fmeanNpvAcBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  fmeanNpvAcBDTchan0->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  fmeanNpvAcBDTchan0->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvmean-accuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[0]->xmlFile.c_str()));

    fmeanNpvAdBDTchan0->SetAxisRange(0., npvmax, "X");
    fmeanNpvAdBDTchan0->Fit(Form("pol%d", j));  
    if (fmeanNpvAdBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  fmeanNpvAdBDTchan0->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  fmeanNpvAdBDTchan0->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (fmeanNpvAdBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  fmeanNpvAdBDTchan0->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  fmeanNpvAdBDTchan0->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvmean-adcuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[0]->xmlFile.c_str()));

  
    fmeanNpvBDTchan1->SetAxisRange(0., npvmax, "X");
    fmeanNpvBDTchan1->Fit(Form("pol%d", j));  
    if (fmeanNpvBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  fmeanNpvBDTchan1->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  fmeanNpvBDTchan1->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (fmeanNpvBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  fmeanNpvBDTchan1->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  fmeanNpvBDTchan1->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvmean-nocuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[1]->xmlFile.c_str()));

    fmeanNpvAcBDTchan1->SetAxisRange(0., npvmax, "X");
    fmeanNpvAcBDTchan1->Fit(Form("pol%d", j));  
    if (fmeanNpvAcBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  fmeanNpvAcBDTchan1->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  fmeanNpvAcBDTchan1->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (fmeanNpvAcBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  fmeanNpvAcBDTchan1->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  fmeanNpvAcBDTchan1->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvmean-accuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[1]->xmlFile.c_str()));

    fmeanNpvAdBDTchan1->SetAxisRange(0., npvmax, "X");
    fmeanNpvAdBDTchan1->Fit(Form("pol%d", j));  
    if (fmeanNpvAdBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  fmeanNpvAdBDTchan1->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  fmeanNpvAdBDTchan1->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (fmeanNpvAdBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  fmeanNpvAdBDTchan1->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  fmeanNpvAdBDTchan1->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvmean-adcuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[1]->xmlFile.c_str()));


    // -- RMS
    frmsNpvBDTchan0->SetAxisRange(0., npvmax, "X");
    frmsNpvBDTchan0->Fit(Form("pol%d", j));  
    if (frmsNpvBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  frmsNpvBDTchan0->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  frmsNpvBDTchan0->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (frmsNpvBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  frmsNpvBDTchan0->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  frmsNpvBDTchan0->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvrms-nocuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[0]->xmlFile.c_str()));

    frmsNpvAcBDTchan0->SetAxisRange(0., npvmax, "X");
    frmsNpvAcBDTchan0->Fit(Form("pol%d", j));  
    if (frmsNpvAcBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  frmsNpvAcBDTchan0->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  frmsNpvAcBDTchan0->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (frmsNpvAcBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  frmsNpvAcBDTchan0->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  frmsNpvAcBDTchan0->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvrms-accuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[0]->xmlFile.c_str()));

    frmsNpvAdBDTchan0->SetAxisRange(0., npvmax, "X");
    frmsNpvAdBDTchan0->Fit(Form("pol%d", j));  
    if (frmsNpvAdBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  frmsNpvAdBDTchan0->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  frmsNpvAdBDTchan0->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (frmsNpvAdBDTchan0->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  frmsNpvAdBDTchan0->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  frmsNpvAdBDTchan0->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvrms-adcuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[0]->xmlFile.c_str()));

    
    frmsNpvBDTchan1->SetAxisRange(0., npvmax, "X");
    frmsNpvBDTchan1->Fit(Form("pol%d", j));  
    if (frmsNpvBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  frmsNpvBDTchan1->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  frmsNpvBDTchan1->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (frmsNpvBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  frmsNpvBDTchan1->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  frmsNpvBDTchan1->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvrms-nocuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[1]->xmlFile.c_str()));

    frmsNpvAcBDTchan1->SetAxisRange(0., npvmax, "X");
    frmsNpvAcBDTchan1->Fit(Form("pol%d", j));  
    if (frmsNpvAcBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  frmsNpvAcBDTchan1->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  frmsNpvAcBDTchan1->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (frmsNpvAcBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  frmsNpvAcBDTchan1->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  frmsNpvAcBDTchan1->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvrms-accuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[1]->xmlFile.c_str()));

    frmsNpvAdBDTchan1->SetAxisRange(0., npvmax, "X");
    frmsNpvAdBDTchan1->Fit(Form("pol%d", j));  
    if (frmsNpvAdBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.2, 0.92, Form("p%d = %5.4f #pm %5.4f", j,
				  frmsNpvAdBDTchan1->GetFunction(Form("pol%d", j))->GetParameter(j), 
				  frmsNpvAdBDTchan1->GetFunction(Form("pol%d", j))->GetParError(j)));
    if (frmsNpvAdBDTchan1->GetFunction(Form("pol%d", j))) tl->DrawLatex(0.6, 0.92, Form("#chi^{2}/dof = %3.1f/%3d", 
				  frmsNpvAdBDTchan1->GetFunction(Form("pol%d", j))->GetChisquare(), 
				  frmsNpvAdBDTchan1->GetFunction(Form("pol%d", j))->GetNDF()));
    c0->SaveAs(Form("%s/dep-bdt-%s-npvrms-adcuts-pol%d-%s.pdf", fDirectory.c_str(), mode.c_str(), j, fCuts[1]->xmlFile.c_str()));
  }
  
  for (unsigned int i = 0; i < fNchan; ++i) {
    fpMassBDT[i]->SetAxisRange(4.9, 5.9, "X");
    fpMassBDT[i]->Fit("pol0");  
    c0->SaveAs(Form("%s/dep-bdt-%s-%s-mass-nocuts%d.pdf", fDirectory.c_str(), fCuts[i]->xmlFile.c_str(), mode.c_str(), i));

    fpMassAcBDT[i]->SetAxisRange(4.9, 5.9, "X");
    fpMassAcBDT[i]->Fit("pol0");  
    c0->SaveAs(Form("%s/dep-bdt-%s-%s-mass-aftercuts%d.pdf", fDirectory.c_str(), fCuts[1]->xmlFile.c_str(), mode.c_str(), i));

  }

}


// ----------------------------------------------------------------------
void plotBDT::bdtScan() { 

  TH1D *h(0);
  string mode = "SgMc";
  fIsMC = true; 
  cout << "bdtScan: running over " << mode << endl;
  TTree *t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  resetHistograms();
  setupTree(t, mode);
  loopOverTree(t, mode, 1); 
  //  loopOverTree(t, mode, 1, 50000); 
  fhSgBDT.push_back((TH1D*)fhBDT[0]->Clone("hSgBDT0"));
  fhSgBDT.push_back((TH1D*)fhBDT[1]->Clone("hSgBDT1"));

  // -- efficiency histogram (includes muid and HLT!)
  double pass(0.), all(0.), eff(0.);
  for (unsigned int i = 0; i < fNchan; ++i) {
    cout << "channel " << i << endl;
    h = new TH1D(Form("hSgEff%d", i), Form("hSgEff%d", i), 200, -1., 1.); 
    h->SetDirectory(fHistFile); 

    all  = fhMassNoCuts[i]->GetSumOfWeights();
    for (unsigned int j = 1; j < 201; ++j) {
      pass = fhSgBDT[i]->Integral(j, 200);
      eff = pass/all;
      h->SetBinContent(j, eff); 
    }

  }

  mode = "SgData";
  fIsMC = false;
  cout << "bdtScan: running over " << mode << endl;
  t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  resetHistograms();
  setupTree(t, mode);
  loopOverTree(t, mode, 1); 
  fhBgBDT.push_back((TH1D*)fhBDT[0]->Clone("hBgBDT0"));
  fhBgBDT.push_back((TH1D*)fhBDT[1]->Clone("hBgBDT1"));
  
  // -- background count histogram
  for (unsigned int i = 0; i < fNchan; ++i) {
    cout << "channel " << i << endl;
    h = new TH1D(Form("hBgEvts%d", i), Form("hBgEvts%d", i), 200, -1., 1.); 
    h->SetDirectory(fHistFile); 

    for (int j = 1; j < 201; ++j) {
      pass = fhBgBDT[i]->Integral(j, 200);
      h->SetBinContent(j, pass); 
    }

    // -- all mass range
    h = new TH1D(Form("hBgEff%d", i), Form("hBgEff%d", i), 200, -1., 1.); 
    h->SetDirectory(fHistFile); 
    
    all  = fhMassNoCuts[i]->GetSumOfWeights();
    for (int j = 1; j < 201; ++j) {
      pass = fhBgBDT[i]->Integral(j, 200);
      eff = pass/all;
      h->SetBinContent(j, eff); 
    }

    // -- low mass range
    h = new TH1D(Form("hBgLoEff%d", i), Form("hBgLoEff%d", i), 200, -1., 1.); 
    h->SetDirectory(fHistFile); 
    
    all  = fhMassNoCuts[i]->Integral(fhMassNoCuts[i]->FindBin(4.9), fhMassNoCuts[i]->FindBin(5.2));
    for (int j = 1; j < 201; ++j) {
      pass = fhLoBDT[i]->Integral(j, 200);
      eff = pass/all;
      h->SetBinContent(j, eff); 
    }

    // -- high mass range
    h = new TH1D(Form("hBgHiEff%d", i), Form("hBgHiEff%d", i), 200, -1., 1.); 
    h->SetDirectory(fHistFile); 
    
    all  = fhMassNoCuts[i]->Integral(fhMassNoCuts[i]->FindBin(5.45), fhMassNoCuts[i]->FindBin(5.9));
    for (int j = 1; j < 201; ++j) {
      pass = fhHiBDT[i]->Integral(j, 200);
      eff = pass/all;
      h->SetBinContent(j, eff); 
    }
    
  }
  
  fHistFile->cd();
  
  ((TH1D*)gDirectory->Get("hBgEvts0"))->Write();
  ((TH1D*)gDirectory->Get("hBgEvts1"))->Write();
  ((TH1D*)gDirectory->Get("hSgEff0"))->Write();
  ((TH1D*)gDirectory->Get("hSgEff1"))->Write();

  ((TH1D*)gDirectory->Get("hBgLoEff0"))->Write();
  ((TH1D*)gDirectory->Get("hBgLoEff1"))->Write();

  ((TH1D*)gDirectory->Get("hBgHiEff0"))->Write();
  ((TH1D*)gDirectory->Get("hBgHiEff1"))->Write();

  // -- overlay low and high mass efficiencies
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  zone();
  TH1D *H1 = (TH1D*)gDirectory->Get("hBgLoEff0");
  setTitles(H1, "b > ", "relative efficiency");
  cout << H1->GetTitle() << " " << 1./H1->GetBinContent(1) << endl;
  TH1D *H2 = (TH1D*)gDirectory->Get("hBgHiEff0");
  cout << H2->GetTitle() << " " << 1./H2->GetBinContent(1) << endl;
  H1->Scale(1./H1->GetBinContent(1));
  H2->Scale(1./H2->GetBinContent(1));
  setHist(H1, kBlue); 
  setHist(H2, kRed); H2->SetLineStyle(kDashed);
  gPad->SetLogy(0); 
  H1->Draw();
  H2->Draw("same");
  newLegend(0.55, 0.7, 0.85, 0.85, "Barrel"); 
  legg->AddEntry(H1, "low sideband", "l"); 
  legg->AddEntry(H2, "high sideband", "l"); 
  legg->Draw();

  c0->SaveAs(Form("%s/%s-lin-Bg-HiLo-Efficiency-chan0.pdf", fDirectory.c_str(), fCuts[0]->xmlFile.c_str())); 
  gPad->SetLogy(1); 
  H1->Draw();
  H2->Draw("same");
  c0->SaveAs(Form("%s/%s-log-Bg-HiLo-Efficiency-chan0.pdf", fDirectory.c_str(), fCuts[1]->xmlFile.c_str())); 

  zone();
  H1 = (TH1D*)gDirectory->Get("hBgLoEff1"); 
  setTitles(H1, "b > ", "relative efficiency");
  cout << H1->GetTitle() << " " << 1./H1->GetBinContent(1) << endl;
  H2 = (TH1D*)gDirectory->Get("hBgHiEff1");
  cout << H2->GetTitle() << " " << 1./H2->GetBinContent(1) << endl;
  H1->Scale(1./H1->GetBinContent(1));
  H2->Scale(1./H2->GetBinContent(1));
  setHist(H1, kBlue);
  setHist(H2, kRed); H2->SetLineStyle(kDashed);
  gPad->SetLogy(0); 
  H1->Draw();
  H2->Draw("same");
  newLegend(0.55, 0.7, 0.85, 0.85, "Endcap"); 
  legg->AddEntry(H1, "low sideband", "l"); 
  legg->AddEntry(H2, "high sideband", "l"); 
  legg->Draw();

  c0->SaveAs(Form("%s/%s-lin-Bg-HiLo-Efficiency-chan1.pdf", fDirectory.c_str(), fCuts[0]->xmlFile.c_str())); 
  gPad->SetLogy(1);
  H1->Draw();
  H2->Draw("same");
  c0->SaveAs(Form("%s/%s-log-Bg-HiLo-Efficiency-chan1.pdf", fDirectory.c_str(), fCuts[1]->xmlFile.c_str())); 

}


// ----------------------------------------------------------------------
void plotBDT::tmvaControlPlots() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  string type[3] = {"Events0", "Events1", "Events2"};
  string XmlName;
  for (int j = 0; j < 3; ++j) {

    for (unsigned int i = 0; i < fNchan; ++i) {
    
      XmlName = "weights/" + fCuts[i]->xmlFile + Form("-%s_BDT.weights.xml", type[j].c_str()); 
      string rootfile = XmlName; 
      replaceAll(rootfile, "_BDT.weights.xml", ".root"); 
      fRootFile = TFile::Open(rootfile.c_str());
      cout << "fRootFile: " << rootfile << endl;

      fBdtString = fRootFile->GetName(); 
      fBdtString = fBdtString.substr(0, fBdtString.find(".root"));
      fBdtString = fBdtString.substr(fBdtString.find("/")+1);
      fBdtString = Form("%s-%s", fCuts[i]->xmlFile.c_str(), type[j].c_str()); 
      cout << "fBdtString: " << fBdtString << endl;
      
      dumpParameters();
      tmvaPlots(type[j]);
      fRootFile->Close();

      if (0 == j) {
	rootfile = "weights/" + fCuts[i]->xmlFile + "-combined.root";
	fRootFile = TFile::Open(rootfile.c_str());
	cout << "fRootFile: " << rootfile << endl;

	fBdtString = fRootFile->GetName(); 
	fBdtString = fBdtString.substr(0, fBdtString.find(".root"));
	fBdtString = fBdtString.substr(fBdtString.find("/")+1);
	fBdtString = Form("TMVA-%d", i); 
	cout << "fBdtString: " << fBdtString << endl;

	variableRanking();

	c0->cd();
	gPad->SetLeftMargin(0.20); 
	TH2F* frame(0); 
	frame = new TH2F("frame", "BDT output distributions", 100, 0., 0.65, 100, 0.99999, 1.000001);
	frame->GetXaxis()->SetTitle(" #epsilon_{S}");
	frame->GetYaxis()->SetTitle(" 1 - #epsilon_{B}");
	frame->Draw();  

	TGraph *g = (TGraph*)fRootFile->Get("groc"); 
	g->Draw("p");
	TGraph *g0 = (TGraph*)fRootFile->Get("groc0"); 
	g0->SetMarkerColor(kBlue);
	g0->SetLineColor(kBlue);
	g0->SetLineWidth(2);
	g0->Draw("l");
	TGraph *g1 = (TGraph*)fRootFile->Get("groc1"); 
	g1->SetMarkerColor(kRed);
	g1->SetLineColor(kRed);
	g1->SetLineWidth(2);
	g1->Draw("l");
	TGraph *g2 = (TGraph*)fRootFile->Get("groc2"); 
	g2->SetMarkerColor(kBlack);
	g2->SetLineColor(kBlack);
	g2->SetLineWidth(2);
	g2->Draw("l");

	TGraph *grocop = (TGraph*)fRootFile->Get("grocop"); 
	grocop->SetMarkerColor(kBlack);
	grocop->SetLineColor(kBlack);
	grocop->SetLineWidth(2);
	grocop->SetMarkerStyle(34); 
	grocop->SetMarkerSize(3.); 
	grocop->Draw("p");


	newLegend(0.25, 0.3, 0.50, 0.5); 
	legg->AddEntry(g, "combined", "p"); 
	legg->AddEntry(g0, "BDT0", "l"); 
	legg->AddEntry(g1, "BDT1", "l"); 
	legg->AddEntry(g2, "BDT2", "l"); 
	legg->Draw();

	float integral(-1.); 
	sscanf(g->GetTitle(), "integral = %f", &integral); 
	fTEX << formatTex(integral, Form("%s:%s:ROCintegral",  fSuffix.c_str(), fBdtString.c_str()), 3) << endl;
	  
	c0->SaveAs(Form("%s/%s-roc.pdf", fDirectory.c_str(), fCuts[i]->xmlFile.c_str())); 


	fRootFile->Close();
      }

    }
  }

}


// ----------------------------------------------------------------------
void plotBDT::plotSSB() {
  cout << "fNchan: " << fNchan << "   " << fCuts[1]->xmlFile << endl;
  for (unsigned int i = 0; i < fNchan; ++i) {
    string rootfile = "weights/" + fCuts[i]->xmlFile + "-combined.root";
    cout << "fRootFile: " << rootfile << endl;
    fRootFile = TFile::Open(rootfile.c_str());
    
    fBdtString = fRootFile->GetName(); 
    fBdtString = fBdtString.substr(0, fBdtString.find(".root"));
    fBdtString = fBdtString.substr(fBdtString.find("/")+1);
    fBdtString = Form("%s", fCuts[i]->xmlFile.c_str()); 
    cout << "fBdtString: " << fBdtString << endl;
    
    ssb();
    fRootFile->Close();
  }
}


// ----------------------------------------------------------------------
void plotBDT::dumpParameters() {
  TH1D *h = (TH1D*)fRootFile->Get("hSetup");
  for (int ibin = 1; ibin < h->GetNbinsX(); ++ibin) {
    if (!strcmp("NTrees", h->GetXaxis()->GetBinLabel(ibin))) { 
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s:NTrees",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
    }
    if (!strcmp("nEventsMin", h->GetXaxis()->GetBinLabel(ibin))) { 
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s:nEventsMin",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
    }
    if (!strcmp("MaxDepth", h->GetXaxis()->GetBinLabel(ibin))) { 
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s:MaxDepth",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
    }
    if (!strcmp("nCuts", h->GetXaxis()->GetBinLabel(ibin)))  {
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s:nCuts",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
    }
    if (!strcmp("AdaBoostBeta", h->GetXaxis()->GetBinLabel(ibin))) {
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s:AdaBoostBeta",  fSuffix.c_str(), fBdtString.c_str()), 2) << endl;
    }
    if (!strcmp("NNodesMax", h->GetXaxis()->GetBinLabel(ibin))) {
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s:NNodesMax",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
    }

    if (!strcmp("sgcnt", h->GetXaxis()->GetBinLabel(ibin))) {
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s:sgcnt",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
    }

    if (!strcmp("bgcnt", h->GetXaxis()->GetBinLabel(ibin))) {
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s:bgcnt",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
    }
  }


  h = getPreselectionNumbers();
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    if (!strcmp("", h->GetXaxis()->GetBinLabel(i))) continue;
    fTEX << formatTex(h->GetBinContent(i), Form("%s:%s:%s",  fSuffix.c_str(), fBdtString.c_str(), h->GetXaxis()->GetBinLabel(i)), 2) << endl;
  }
}


// ----------------------------------------------------------------------
void plotBDT::validateAllDistributions() {

  validateDistributions(0, "train", 0);
  validateDistributions(0, "train", 1);

  validateDistributions(1, "train", 0);
  validateDistributions(1, "train", 1);

  validateDistributions(0, "test", 0);
  validateDistributions(0, "test", 1);

  validateDistributions(1, "test", 0);
  validateDistributions(1, "test", 1);

}



// ----------------------------------------------------------------------
void plotBDT::validateDistributions(int channel, const char *type, int classID) {

  string fname0 = Form("weights/%s-Events0.root", fCuts[channel]->xmlFile.c_str()); 
  string fname1 = Form("weights/%s-Events1.root", fCuts[channel]->xmlFile.c_str());
  string fname2 = Form("weights/%s-Events2.root", fCuts[channel]->xmlFile.c_str());

  //  string sname  = Form("%s-%s-%s", fCuts[channel]->xmlFile.c_str(), type, (classID == 0?"sg":"bg")); 
  string sname  = Form("%s-%s-%s", fCuts[channel]->xmlFile.c_str(), type, (classID == 0?"sg":"bg")); 

  TFile *fEvt0 = TFile::Open(fname0.c_str()); 
  TFile *fEvt1 = TFile::Open(fname1.c_str()); 
  TFile *fEvt2 = TFile::Open(fname2.c_str()); 

  vector<string> vNames, vA;
  vector<double> vMin, vMax; 
  vector<int> vNbins; 
  vector<int> vLog; 
  vNames.push_back("m1pt"); vMin.push_back(0.); vMax.push_back(40.); vNbins.push_back(100); vLog.push_back(0); vA.push_back("p_{T,#mu1} [GeV]");
  vNames.push_back("m2pt"); vMin.push_back(0.); vMax.push_back(20.); vNbins.push_back(100); vLog.push_back(0); vA.push_back("p_{T,#mu2} [GeV]");
  vNames.push_back("m1eta"); vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); vLog.push_back(0); vA.push_back("#eta_{#mu1}");
  vNames.push_back("m2eta"); vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); vLog.push_back(0); vA.push_back("#eta_{#mu2}");
  vNames.push_back("pt");    vMin.push_back(0.); vMax.push_back(40.); vNbins.push_back(100); vLog.push_back(0); vA.push_back("p_{T} [GeV]");
  vNames.push_back("eta");   vMin.push_back(-2.5); vMax.push_back(2.5); vNbins.push_back(100); vLog.push_back(0); vA.push_back("#eta");
  vNames.push_back("fls3d"); vMin.push_back(0.); vMax.push_back(100.); vNbins.push_back(120); vLog.push_back(1); vA.push_back("l/#sigma(l)"); 
  vNames.push_back("alpha"); vMin.push_back(0.); vMax.push_back(0.3); vNbins.push_back(100); vLog.push_back(0);  vA.push_back("#alpha"); 
  vNames.push_back("maxdoca"); vMin.push_back(0.); vMax.push_back(0.03); vNbins.push_back(100); vLog.push_back(0); vA.push_back("d^{max} [cm]"); 
  vNames.push_back("pvip"); vMin.push_back(0.); vMax.push_back(0.05); vNbins.push_back(100); vLog.push_back(0); vA.push_back("#delta_{3D}"); 
  vNames.push_back("pvips"); vMin.push_back(0.); vMax.push_back(5); vNbins.push_back(100); vLog.push_back(0); vA.push_back("#delta_{3D}/#sigma(#delta_{3D})"); 
  vNames.push_back("iso"); vMin.push_back(0.6); vMax.push_back(1.01); vNbins.push_back(41); vLog.push_back(0); vA.push_back("isolation"); 
  vNames.push_back("closetrk"); vMin.push_back(0.); vMax.push_back(21); vNbins.push_back(21); vLog.push_back(0); vA.push_back("N_{trk}^{close}"); 
  vNames.push_back("docatrk"); vMin.push_back(0.); vMax.push_back(0.1); vNbins.push_back(100); vLog.push_back(0); vA.push_back("d_{ca}^{0} [cm]"); 
  vNames.push_back("chi2dof"); vMin.push_back(0.); vMax.push_back(5); vNbins.push_back(100); vLog.push_back(0); vA.push_back("#chi^{2}/dof"); 

  
  TTree *tEvt0(0), *tEvt1(0), *tEvt2(0); 
  if (!strcmp("test", type)) {
    tEvt0 = (TTree*)fEvt0->Get("TestTree");
    tEvt1 = (TTree*)fEvt1->Get("TestTree");
    tEvt2 = (TTree*)fEvt2->Get("TestTree");
  } else {
    tEvt0 = (TTree*)fEvt0->Get("TrainTree");
    tEvt1 = (TTree*)fEvt1->Get("TrainTree");
    tEvt2 = (TTree*)fEvt2->Get("TrainTree");
  }

  TH1D *h0, *h1, *h2; 
  string h0Name, h1Name, h2Name; 

  // -- KS probabilities
  TH1D *sgKS = new TH1D("sgKS", "", 20, 0., 1.); setFilledHist(sgKS, kBlack, kYellow, 1000); 
  TH1D *bgKS = new TH1D("bgKS", "", 20, 0., 1.); setFilledHist(bgKS, kBlack, kYellow, 1000); 
  
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  for (unsigned int i = 0; i < vNames.size(); ++i) {
    
    if (fReaderVariables.end() == find(fReaderVariables.begin(), fReaderVariables.end(), vNames[i]))  {
      //      cout << " =====> " << vNames[i] << " not found in fReaderVariables" << endl;
      continue;
    }

    h0Name = Form("h0_%s", vNames[i].c_str());
    h0  = new TH1D(h0Name.c_str(), vNames[i].c_str(), vNbins[i], vMin[i], vMax[i]);

    h1Name = Form("h1_%s", vNames[i].c_str());
    h1  = new TH1D(h1Name.c_str(), vNames[i].c_str(), vNbins[i], vMin[i], vMax[i]);

    h2Name = Form("h2_%s", vNames[i].c_str());
    h2 = new TH1D(h2Name.c_str(), h2Name.c_str(), vNbins[i], vMin[i], vMax[i]);

    setTitles(h0, vA[i].c_str(), "#candidates/bin"); 
    SetSignalAndBackgroundStyle(h0, h1, h2);            

    tEvt0->Draw(Form("%s>>%s", vNames[i].c_str(), h0Name.c_str()), Form("classID==%d", classID)); 
    tEvt1->Draw(Form("%s>>%s", vNames[i].c_str(), h1Name.c_str()), Form("classID==%d", classID)); 
    tEvt2->Draw(Form("%s>>%s", vNames[i].c_str(), h2Name.c_str()), Form("classID==%d", classID)); 

    double hmax = h0->GetMaximum(); 
    if (h1->GetMaximum() > hmax) hmax = h1->GetMaximum();
    if (h2->GetMaximum() > hmax) hmax = h2->GetMaximum();
    h0->SetMaximum(1.3*hmax);
    
    h0->Draw();
    h1->Draw("same");
    h2->Draw("same");
    double ks01 = h0->KolmogorovTest(h1);
    double ks12 = h1->KolmogorovTest(h2);
    double ks20 = h2->KolmogorovTest(h0);
    tl->DrawLatex(0.15, 0.92, Form("P(KS)= %4.3f/%4.3f/%4.3f", ks01, ks12, ks20)); 

    if (0 == classID) {
      sgKS->Fill(ks01); 
      sgKS->Fill(ks12); 
      sgKS->Fill(ks20); 
    } else {
      bgKS->Fill(ks01); 
      bgKS->Fill(ks12); 
      bgKS->Fill(ks20); 
    }
    
    if (1 == vLog[i]) {
      gPad->SetLogy(1); 
    } else {
      gPad->SetLogy(0); 
    }

    c0->SaveAs(Form("%s/ks-%s_%s.pdf", 
		    fDirectory.c_str(), sname.c_str(), vNames[i].c_str())); 
  }


  if (0 == classID) {
    sgKS->Draw();
    tl->DrawLatex(0.15, 0.92, Form("mean = %4.3f", sgKS->GetMean())); 
    tl->DrawLatex(0.60, 0.92, Form("RMS = %4.3f",  sgKS->GetRMS())); 
    
    c0->SaveAs(Form("%s/ks-%s-sg-probs.pdf", fDirectory.c_str(), sname.c_str())); 
  }

  if (1 == classID) {
    bgKS->Draw();
    tl->DrawLatex(0.15, 0.92, Form("mean = %4.3f", bgKS->GetMean())); 
    tl->DrawLatex(0.60, 0.92, Form("RMS = %4.3f",  bgKS->GetRMS())); 
    c0->SaveAs(Form("%s/ks-%s-bg-probs.pdf", 
		    fDirectory.c_str(), sname.c_str())); 
  }

} 



// ----------------------------------------------------------------------
void plotBDT::tmvaPlots(string type) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  zone();
  shrinkPad(0.15, 0.15, 0.15); 
  TH2 *h2 = (TH2*)fRootFile->Get("CorrelationMatrixS");
  h2->Draw("colztext");
  c0->SaveAs(Form("%s/%s-CorrelationMatrixS.pdf", fDirectory.c_str(), fBdtString.c_str()));

  h2 = (TH2*)fRootFile->Get("CorrelationMatrixB");
  h2->Draw("colztext");
  c0->SaveAs(Form("%s/%s-CorrelationMatrixB.pdf", fDirectory.c_str(), fBdtString.c_str()));

  // -- from "mvas"
  // --------------
  Int_t width = 600;   // size of canvas
  
  // this defines how many canvases we need
  TCanvas *c = new TCanvas( Form("canvas%d", 1), "canvas1",  200, 20, width, static_cast<int>(width*0.78) ); 


  // search for the right histograms in full list of keys
  TIter next(fRootFile->GetListOfKeys());
  TKey *key(0);   
  while ((key = (TKey*)next())) {
    
    if (!TString(key->GetName()).BeginsWith("Method_")) continue;
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;
    
    TString methodName;
    GetMethodName(methodName,key);
    
    TDirectory* mDir = (TDirectory*)key->ReadObj();
    
    TIter keyIt(mDir->GetListOfKeys());
    TKey *titkey;
    while ((titkey = (TKey*)keyIt())) {
      
      if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;
      c->Clear();

      TDirectory *titDir = (TDirectory *)titkey->ReadObj();
      TString methodTitle;
      GetMethodTitle(methodTitle,titDir);
      
      cout << "--- Found directory for method: " << methodName << "::" << methodTitle << flush;
      TString hname = "MVA_" + methodTitle;
      TH1* sig = dynamic_cast<TH1*>(titDir->Get( hname + "_S" ));
      TH1* bgd = dynamic_cast<TH1*>(titDir->Get( hname + "_B" ));
      
      
      cout << " containing " << hname << "_S/_B" << endl;
      // chop off useless stuff
      sig->SetTitle( Form("TMVA overtraining check for classifier: %s", methodTitle.Data()) );
      
      // create new canvas
      TString ctitle = Form("TMVA comparison %s",methodTitle.Data()) ;
      
      // set the histogram style
      SetSignalAndBackgroundStyle( sig, bgd );
      
      // normalise both signal and background
      NormalizeHists( sig, bgd );
      
      // frame limits (choose judicuous x range)
      Float_t nrms = 10;
      cout << "--- Mean and RMS (S): " << sig->GetMean() << ", " << sig->GetRMS() << endl;
      cout << "--- Mean and RMS (B): " << bgd->GetMean() << ", " << bgd->GetRMS() << endl;
      Float_t xmin = TMath::Max( TMath::Min(sig->GetMean() - nrms*sig->GetRMS(), 
					    bgd->GetMean() - nrms*bgd->GetRMS() ),
				 sig->GetXaxis()->GetXmin() );
      Float_t xmax = TMath::Min( TMath::Max(sig->GetMean() + nrms*sig->GetRMS(), 
					    bgd->GetMean() + nrms*bgd->GetRMS() ),
				 sig->GetXaxis()->GetXmax() );
      Float_t ymin = 0;
      Float_t maxMult = 2.0;
      Float_t ymax = TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*maxMult;
      
      xmin = -1.;
      xmax = 1.;
      
      // build a frame
      Int_t nb = 500;
      TString hFrameName(TString("frame") + methodTitle);
      TObject *o = gROOT->FindObject(hFrameName);
      if(o) delete o;
      TH2F* frame = new TH2F( hFrameName, sig->GetTitle(), 
			      nb, xmin, xmax, nb, ymin, ymax );
      frame->GetXaxis()->SetTitle( methodTitle + " response"  );
      frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
      SetFrameStyle( frame );
      
      // eventually: draw the frame
      frame->Draw();  
      
      c->GetPad(0)->SetLeftMargin( 0.105 );
      frame->GetYaxis()->SetTitleOffset( 1.2 );
      
      // Draw legend               
      TLegend *legend= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin() - 0.12, 
				    c->GetLeftMargin() +  0.40, 1 - c->GetTopMargin() );
      legend->SetFillStyle( 1 );
      legend->AddEntry(sig,TString("Signal")     + " (test sample)", "F");
      legend->AddEntry(bgd,TString("Background") + " (test sample)", "F");
      legend->SetBorderSize(1);
      legend->SetMargin(0.2);
      legend->Draw("same");
      
      // overlay signal and background histograms
      sig->Draw("samehist");
      bgd->Draw("samehist");
      
      TH1* sigOv = 0;
      TH1* bgdOv = 0;
      
      TString ovname = hname += "_Train";
      sigOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_S" ));
      bgdOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_B" ));
      
      if (sigOv == 0 || bgdOv == 0) {
	cout << "+++ Problem in \"mvas.C\": overtraining check histograms do not exist" << endl;
      }
      else {
	cout << "--- Found comparison histograms for overtraining check" << endl;
	
	TLegend *legend2= new TLegend( 1 - c->GetRightMargin() - 0.42, 1 - c->GetTopMargin() - 0.12,
				       1 - c->GetRightMargin(), 1 - c->GetTopMargin() );
	legend2->SetFillStyle( 1 );
	legend2->SetBorderSize(1);
	legend2->AddEntry(sigOv,"Signal (training sample)","P");
	legend2->AddEntry(bgdOv,"Background (training sample)","P");
	legend2->SetMargin( 0.1 );
	legend2->Draw("same");
      }
      // normalise both signal and background
      NormalizeHists( sigOv, bgdOv );
      
      Int_t col = sig->GetLineColor();
      sigOv->SetMarkerColor( col );
      sigOv->SetMarkerSize( 0.7 );
      sigOv->SetMarkerStyle( 20 );
      sigOv->SetLineWidth( 1 );
      sigOv->SetLineColor( col );
      sigOv->Draw("e1same");
      
      col = bgd->GetLineColor();
      bgdOv->SetMarkerColor( col );
      bgdOv->SetMarkerSize( 0.7 );
      bgdOv->SetMarkerStyle( 20 );
      bgdOv->SetLineWidth( 1 );
      bgdOv->SetLineColor( col );
      bgdOv->Draw("e1same");
      
      ymax = TMath::Max(ymax, 
			TMath::Max( static_cast<Float_t>(sigOv->GetMaximum()), static_cast<Float_t>(bgdOv->GetMaximum()) )*maxMult);
      frame->GetYaxis()->SetLimits( 0, ymax );
      
      // for better visibility, plot thinner lines
      sig->SetLineWidth( 1 );
      bgd->SetLineWidth( 1 );
      
      // perform K-S test
      cout << "--- Perform Kolmogorov-Smirnov tests" << endl;
      Double_t kolS = sig->KolmogorovTest( sigOv );
      Double_t kolB = bgd->KolmogorovTest( bgdOv );
      cout << "--- Goodness of signal (background) consistency: " << kolS << " (" << kolB << ")" << endl;
      fTEX << formatTex(kolS, Form("%s:%s:KS:Sg",  fSuffix.c_str(), fBdtString.c_str()), 3) << endl;
      fTEX << formatTex(kolB, Form("%s:%s:KS:Bg",  fSuffix.c_str(), fBdtString.c_str()), 3) << endl;
      
      TString probatext = Form( "Kolmogorov-Smirnov test: signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
      TText* tt = new TText( 0.12, 0.91, probatext );
      tt->SetNDC(); tt->SetTextSize( 0.032 ); tt->AppendPad(); 
      
      
      
      
      // redraw axes
      frame->Draw("sameaxis");
      
      // text for overflows
      Int_t    nbin = sig->GetNbinsX();
      Double_t dxu  = sig->GetBinWidth(0);
      Double_t dxo  = sig->GetBinWidth(nbin+1);
      TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%", 
			     sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
			     sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );
      TText* t = new TText( 0.975, 0.115, uoflow );
      t->SetNDC();
      t->SetTextSize( 0.030 );
      t->SetTextAngle( 90 );
      t->AppendPad();    
      
      c->Update();
      
      // save canvas to file
      c->SaveAs(Form("%s/%s-overtrain0.pdf", fDirectory.c_str(), fBdtString.c_str())); 
      frame->GetYaxis()->SetLimits(1.e-5, 5.e1);
      c->SetLogy(1);
      c->SaveAs(Form("%s/%s-overtrain1.pdf", fDirectory.c_str(), fBdtString.c_str())); 
      
    }
    cout << "";
  }


  
  // -- from "variables"
  // -------------------

  c->SetLogy(0);
      
  TString title = "TMVA Input Variables";
  // obtain shorter histogram title 
  TString htitle = title; 
  htitle.ReplaceAll("variables ","variable");
  htitle.ReplaceAll("and target(s)","");
  htitle.ReplaceAll("(training sample)","");

  TString dirName = "Method_BDT/BDT";
  TDirectory* dir = (TDirectory*)fRootFile->Get(dirName);
  if (dir==0) {
    cout << "No information about " << title << " available in directory " << dirName << " of file " << fRootFile << endl;
    return;
  }
  dir->cd();


  // loop over all objects in directory
  TIter next1(dir->GetListOfKeys());
  while ((key = (TKey*)next1())) {
    // make sure, that we only look at histograms
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1 *sig = (TH1*)key->ReadObj();
    TString hname(sig->GetName());
    if (hname.Contains("__Background")) continue;
    TString htitle(sig->GetTitle());
    if (htitle.Contains("Input Variables")) continue;
    if (htitle.Contains(" ")) continue;

    if (htitle.Contains("MVA_BDT")) {
      //  cout << "--> SKIPPING " << htitle << endl;
      continue;
    }
    
    // find the corredponding backgrouns histo
    TString bgname = hname;
    bgname.ReplaceAll("__Signal","__Background");
    TH1 *bgd = (TH1*)dir->Get(bgname);
    if (bgd == NULL) {
      cout << "ERROR!!! couldn't find background histo for" << hname << endl;
      return;
    }

    // this is set but not stored during plot creation in MVA_Factory
    SetSignalAndBackgroundStyle(sig, bgd);            
    
    sig->SetTitle(TString( htitle ) + ": " + sig->GetTitle() );
    SetFrameStyle(sig, 1.2);

    // normalise both signal and background
    NormalizeHists(sig, bgd);
    
    // finally plot and overlay
    Float_t sc = 1.4;
    sig->SetMaximum( TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*sc );
    sig->Draw( "hist" );
    gPad->SetLeftMargin( 0.17 );
    
    sig->GetYaxis()->SetTitleOffset( 1.70 );
    bgd->Draw("histsame");
    TString ytit = TString("(1/N) ") + sig->GetYaxis()->GetTitle();
    sig->GetYaxis()->SetTitle( ytit ); // histograms are normalised
    
    // Draw legend
    TLegend *legend= new TLegend( gPad->GetLeftMargin(), 
				  1-gPad->GetTopMargin()-.15, 
				  gPad->GetLeftMargin()+.4, 
				  1-gPad->GetTopMargin() );
    legend->SetFillStyle(1);
    legend->AddEntry(sig,"Signal","F");
    legend->AddEntry(bgd,"Background","F");
    legend->SetBorderSize(1);
    legend->SetMargin( 0.3 );
    legend->Draw("same");
    
    // redraw axes
    sig->Draw("sameaxis");
    
    // text for overflows
    Int_t    nbin = sig->GetNbinsX();
    Double_t dxu  = sig->GetBinWidth(0);
    Double_t dxo  = sig->GetBinWidth(nbin+1);
    TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%", 
			   sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
			   sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );
    
    TText* t = new TText( 0.98, 0.14, uoflow );
    t->SetNDC();
    t->SetTextSize( 0.040 );
    t->SetTextAngle( 90 );
    t->AppendPad();    
    
    c->SaveAs(Form("%s/%s-%s.pdf", fDirectory.c_str(), fBdtString.c_str(), htitle.Data())); 

    delete legend;
  }
  
  

  
  // "BoostMonitor","BoostWeight","BoostWeightVsTree","ErrFractHist","NodesBeforePruning"
  dirName = "Method_BDT/BDT";
  dir = (TDirectory*)fRootFile->Get(dirName);
  if (dir==0) {
    cout << "No information about " << title << " available in directory " << dirName << " of file " << fRootFile << endl;
    return;
  }
  dir->cd();

  TCanvas *cc = new TCanvas("cc", "", 300, 200, 1000, 400);

  cc->cd();
  shrinkPad(0.20, 0.15, 0.1, 0.);
  TH1F *hf = (TH1F*)((TH1F*)dir->Get("BoostWeightVsTree"))->Clone("hf");
  setTitles(hf, hf->GetXaxis()->GetTitle(), hf->GetYaxis()->GetTitle(), 0.09, 1.02, 0.8, 0.09); 
  hf->Draw(); 
  cc->SaveAs(Form("%s/%s-BoostWeightVsTree.pdf", fDirectory.c_str(), fBdtString.c_str())); 

  cc->Clear();
  shrinkPad(0.20, 0.15, 0.1, 0.);
  delete hf; 
  hf = (TH1F*)((TH1F*)dir->Get("ErrFractHist"))->Clone("hf");
  setTitles(hf, hf->GetXaxis()->GetTitle(), hf->GetYaxis()->GetTitle(), 0.09, 1.02, 0.8, 0.09); 
  hf->Draw(); 
  cc->SaveAs(Form("%s/%s-ErrFractHist.pdf", fDirectory.c_str(), fBdtString.c_str())); 

  cc->Clear();
  shrinkPad(0.20, 0.15, 0.1, 0.);
  delete hf; 
  hf = (TH1F*)((TH1F*)dir->Get("NodesBeforePruning"))->Clone("hf");
  setTitles(hf, hf->GetXaxis()->GetTitle(), hf->GetYaxis()->GetTitle(), 0.09, 1.02, 0.8, 0.09); 
  hf->Draw();
  cc->SaveAs(Form("%s/%s-NodesBeforePruning.pdf", fDirectory.c_str(), fBdtString.c_str())); 
  

  
}



// ----------------------------------------------------------------------
string plotBDT::replaceLabelWithTex(string label) {
  if ("alpha" == label) return "\\ensuremath{\\alpha}";
  if ("closetrk" == label) return "\\closetrk";
  if ("docatrk" == label) return "\\docatrk";
  if ("fls3d" == label) return "\\fls";
  if ("iso" == label) return "I";
  if ("pvips" == label) return "\\ips";
  if ("pvip" == label) return "\\ip";
  if ("maxdoca" == label) return "\\dca";
  if ("pt" == label) return "\\ptb";
  if ("eta" == label) return "\\etab";
  if ("m1pt" == label) return "\\ptmuone";
  if ("m2pt" == label) return "\\ptmutwo";
  if ("m1eta" == label) return "\\etamuone";
  if ("m2eta" == label) return "\\etamutwo";
  if ("chi2dof" == label) return "\\chidof";
  return "unknown";
}

// ----------------------------------------------------------------------
void plotBDT::xmlResetHists() {
  fhBdtNodes->Reset(); 
  fhBdtNodesW8->Reset(); 

  fhBdtVariables->Reset(); 
  fhBdtVariablesW8->Reset(); 
  for (unsigned int i = 0; i < fhBdtVariableCuts.size(); ++i) {
    fhBdtVariableCuts[i]->Reset(); 
    fhBdtVariableCutsW8[i]->Reset();
  }

}


// ----------------------------------------------------------------------
void plotBDT::xmlParsing() {

  vector<string> etype; 
  etype.push_back("-Events0"); 
  etype.push_back("-Events1"); 
  etype.push_back("-Events2"); 

  string xmlfile = "weights/" + fCuts[0]->xmlFile + etype[0] + "_BDT.weights.xml";
  xmlParsingVariables(xmlfile);
  
  for (unsigned int ichan = 0; ichan < fNchan; ++ichan) {
    for (unsigned int i = 0; i < etype.size(); ++i) {
      xmlfile = "weights/" + fCuts[ichan]->xmlFile + etype[i] + "_BDT.weights.xml";
      fBdtString = fCuts[ichan]->xmlFile;
      cout << xmlfile << endl;

      xmlResetHists();
      xmlParsingReadTree(xmlfile);


    }    
  }


}

// ----------------------------------------------------------------------
void plotBDT::xmlParsingVariables(string weightfile) {

  fhBdtVariables   = new TH1D("bdtVariables", "", 15, 0., 15.); 
  fhBdtVariablesW8 = new TH1D("bdtVariablesW8", "", 15, 0., 15.); 

  fhBdtNodes       = new TH1D("bdtNodes", "", 10, 0., 10.); 
  fhBdtNodesW8     = new TH1D("bdtNodesW8", "", 10, 0., 10.); 

  // -- read in variables from weight file
  vector<string> allLines; 
  char  buffer[2000];
  ifstream is(weightfile.c_str()); 
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  int nvars(-1), ivar(-1); 
  string::size_type m0, m1, m2;
  string stype, varidx; 
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      cout << "  " << stype << " variables" << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m0 = allLines[j].find("VarIndex=\"");
	m1 = allLines[j].find("Expression=\"")+10; 
	m2 = allLines[j].find("\" Label=\"");
	varidx = allLines[j].substr(m0+10, m1-m0-22); 
	ivar = atoi(varidx.c_str()); 
	stype = allLines[j].substr(m1+2, m2-m1-2); 
	fVariableMap.insert(make_pair(ivar, stype)); 
	fhBdtVariables->GetXaxis()->SetBinLabel(ivar+1, stype.c_str()); 
	fhBdtVariablesW8->GetXaxis()->SetBinLabel(ivar+1, stype.c_str()); 
      }
      break;
    }
  }
  

  TH1D *h1(0); 
  for (map<int, string>::iterator imap = fVariableMap.begin(); imap != fVariableMap.end(); ++imap) {  
    cout << "index: " << imap->first << " name: " << imap->second << endl;
    if (imap->second == "m1pt")     h1 = new TH1D("m1pt", "m1pt", 100, 0., 40.); 
    if (imap->second == "m2pt")     h1 = new TH1D("m2pt", "m2pt", 100, 0., 20.); 
    if (imap->second == "m1eta")    h1 = new TH1D("m1eta", "m1eta", 100, -2.4, 2.4); 
    if (imap->second == "m2eta")    h1 = new TH1D("m2eta", "m2eta", 100, -2.4, 2.4); 
    if (imap->second == "pt")       h1 = new TH1D("pt", "pt", 100, 0., 50.); 
    if (imap->second == "eta")      h1 = new TH1D("eta", "eta", 100, -2.4, 2.4); 
    if (imap->second == "fls3d")    h1 = new TH1D("fls3d", "fls3d", 100, 0., 100.); 
    if (imap->second == "alpha")    h1 = new TH1D("alpha", "alpha", 100, 0., 0.3); 
    if (imap->second == "maxdoca")  h1 = new TH1D("maxdoca", "maxdoca", 100, 0., 0.05); 
    if (imap->second == "pvip")     h1 = new TH1D("pvip", "pvip", 100, 0., 0.05); 
    if (imap->second == "pvips")    h1 = new TH1D("pvips", "pvips", 100, 0., 5.); 
    if (imap->second == "iso")      h1 = new TH1D("iso", "iso", 101, 0., 1.01); 
    if (imap->second == "closetrk") h1 = new TH1D("closetrk", "closetrk", 21, 0., 21.); 
    if (imap->second == "docatrk")  h1 = new TH1D("docatrk", "docatrk", 100, 0., 0.2); 
    if (imap->second == "chi2/dof") h1 = new TH1D("chi2dof", "chi2dof", 100, 0., 10.); 
    h1->SetLineColor(kRed); h1->SetLineStyle(kDashed); 
    fhBdtVariableCuts.push_back(h1); 

    h1 = (TH1D*)h1->Clone(Form("w8_%s", h1->GetName())); 
    h1->SetLineColor(kBlue); h1->SetLineStyle(kSolid);
    fhBdtVariableCutsW8.push_back(h1); 
  }  

}



//  = node type: 1 signal node, -1 bkg leave, 0 intermediate Node 
//     <BinaryTree type="DecisionTree" boostWeight="1.6600676679081650e-01" itree="10">
//       <Node pos="s" depth="0" NCoef="0" IVar="14" Cut="3.80e+00" cType="1" res="-9.9+01" rms="0.00e+00" purity="4.68e-01" nType="0">
//         <Node pos="l" depth="1" NCoef="0" IVar="-1" Cut="0.00e+00" cType="1" res="-9.9e+01" rms="0.00e+00" purity="4.93e-01" nType="-1"/>
//         <Node pos="r" depth="1" NCoef="0" IVar="12" Cut="4.66e-02" cType="0" res="-9.9e+01" rms="0.00e+00" purity="2.32e-01" nType="0">
//           <Node pos="l" depth="2" NCoef="0" IVar="-1" Cut="0.00e+00" cType="1" res="-9.9e+01" rms="0.00e+00" purity="1.00e+00" nType="1"/>
//           <Node pos="r" depth="2" NCoef="0" IVar="-1" Cut="0.00e+00" cType="1" res="-9.9e+01" rms="0.00e+00" purity="2.70e-01" nType="-1"/>
//         </Node>
//       </Node>
//     </BinaryTree>
// ----------------------------------------------------------------------
void plotBDT::xmlParsingReadTree(string xmlfile) {
  
  vector<string> allLines; 
  char  buffer[2000];
  ifstream is(xmlfile.c_str()); 
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  float w8(0.), cut(0.); 
  int ivar(-1), nnodes(0); 
  string::size_type m0, m1;
  string stype, varidx; 
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("<BinaryTree type=\"DecisionTree\" boostWeight=")) {
      fhBdtNodes->Fill(nnodes); 
      fhBdtNodesW8->Fill(nnodes, w8); 
      nnodes = 0; 
      m0 = allLines[i].find("boostWeight=\"") + 13; 
      m1 = allLines[i].find("\" itree");
      stype = allLines[i].substr(m0, m1-m0); 
      w8 = atof(stype.c_str()); 
      //cout << "boost weight: " << w8 << " ->" << stype << "<-" 
      //   << " from m0: " << m0 << " m1: " << m1 
      //   << endl;
    }

    if (string::npos != allLines[i].find("nType=\"0\"")) {
      ++nnodes; 
      m0 = allLines[i].find("IVar=\"") + 6; 
      m1 = allLines[i].find("\" Cut");
      stype = allLines[i].substr(m0, m1-m0); 
      ivar = atoi(stype.c_str()); 
      //cout << "  ivar: ->" << stype << "<- " << ivar << endl;
      m0 = allLines[i].find("Cut=\"") + 5; 
      m1 = allLines[i].find("\" cType=");
      stype = allLines[i].substr(m0, m1-m0); 
      cut = atof(stype.c_str()); 
      //cout << "  cut: ->" << stype << "<- " << cut << endl;

      fhBdtVariables->Fill(ivar); 
      fhBdtVariablesW8->Fill(ivar, w8); 

      fhBdtVariableCuts[ivar]->Fill(cut); 
      fhBdtVariableCutsW8[ivar]->Fill(cut, w8); 
    }
  }

  zone(3,5); 
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  for (unsigned int i = 0; i < fhBdtVariableCuts.size(); ++i) {
    c0->cd(i+1); shrinkPad(0.18, 0.15);
    setTitles(fhBdtVariableCuts[i], Form("%s cut", fhBdtVariableCuts[i]->GetTitle()), "Number of decision trees", 0.07, 1.1, 1.1);
    fhBdtVariableCuts[i]->Draw("hist"); 
    fhBdtVariableCutsW8[i]->Scale(fhBdtVariableCuts[i]->GetSumOfWeights()/fhBdtVariableCutsW8[i]->GetSumOfWeights()); 
    fhBdtVariableCutsW8[i]->Draw("samehist"); 
  }
  c0->SaveAs(Form("%s/%s-bdtVariableCuts.pdf", fDirectory.c_str(), fBdtString.c_str())); 

  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 

  string pdfname = xmlfile; 
  rmSubString(pdfname, "weights/");
  rmSubString(pdfname, "_BDT.weights.xml");

  c0->Clear();
  shrinkPad(0.20, 0.15); 
  setTitles(fhBdtVariables, "BDT Variables", "Number of decision trees", 0.06, 1.1, 1.3); 
  setHist(fhBdtVariables, kRed, 20, 1, 3); fhBdtVariables->SetLineStyle(kDashed); 
  fhBdtVariables->Draw();

  fhBdtVariablesW8->Scale(fhBdtVariables->GetSumOfWeights()/fhBdtVariablesW8->GetSumOfWeights()); 
  setHist(fhBdtVariablesW8, kBlue, 20, 1, 3);  
  fhBdtVariablesW8->Draw("same");
  c0->Modified();  c0->Update();
  c0->SaveAs(Form("%s/%s-bdtVariables.pdf", fDirectory.c_str(), fBdtString.c_str())); 

  c0->Clear();
  shrinkPad(0.20, 0.15); 
  setTitles(fhBdtNodes, "BDT Nodes", "Number of decision trees", 0.06, 1.1, 1.3); 
  setHist(fhBdtNodes, kRed, 20, 1, 3); fhBdtNodes->SetLineStyle(kDashed); 
  fhBdtNodes->Draw();

  setHist(fhBdtNodesW8, kBlue, 20, 1, 3); 
  fhBdtNodesW8->Scale(fhBdtNodes->GetSumOfWeights()/fhBdtNodesW8->GetSumOfWeights()); 
  fhBdtNodesW8->Draw("same");
  c0->Modified();  c0->Update();
  c0->SaveAs(Form("%s/%s-bdtNodes.pdf", fDirectory.c_str(), fBdtString.c_str())); 
}





// ----------------------------------------------------------------------
void plotBDT::variableRanking() {
  string which = "IdTransformation";
  TH1D *h(0); 
  vector<string> splits; 
  splits.push_back("events0"); 
  splits.push_back("events1"); 
  splits.push_back("events2"); 

  for (unsigned int i = 0; i < splits.size(); ++i) {
    which = "IdTransformation";
    h = (TH1D*)fRootFile->Get(Form("rank_%s_%s", splits[i].c_str(), which.c_str()));
    if (0 == h) {
      cout << "==> histogram " << Form("rank_%s_%s", splits[i].c_str(), which.c_str()) << " NOT FOUND" << endl;
      break;
    }
    for (int ibin = 1; ibin < h->GetNbinsX(); ++ibin) {
      
      string label = h->GetXaxis()->GetBinLabel(ibin);
      if ("" == label) break;
      replaceAll(label, " ", ""); 
      string texlabel = replaceLabelWithTex(label);
      fTEX << Form("\\vdef{%s:%s_%s_%s:%d:name} {%s}", 
		   fSuffix.c_str(), fBdtString.c_str(), splits[i].c_str(), which.c_str(), ibin, texlabel.c_str()) 
	   << endl;
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s_%s_%s:%d:val", 
						     fSuffix.c_str(), fBdtString.c_str(), splits[i].c_str(), which.c_str(), ibin), 4)
	   << endl;
    }
    
    which = "BDT";
    h = (TH1D*)fRootFile->Get(Form("rank_%s_%s", splits[i].c_str(), which.c_str()));
    for (int ibin = 1; ibin < h->GetNbinsX(); ++ibin) {
      string label = h->GetXaxis()->GetBinLabel(ibin);
      if ("" == label) break;
      replaceAll(label, " ", ""); 
      if (0 == i) fReaderVariables.push_back(label); 
      string texlabel = replaceLabelWithTex(label);
      fTEX << Form("\\vdef{%s:%s_%s_%s:%d:name} {%s}", 
		   fSuffix.c_str(), fBdtString.c_str(), splits[i].c_str(), which.c_str(), ibin, texlabel.c_str()) 
	   << endl;
      fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s_%s_%s:%d:val", 
						     fSuffix.c_str(), fBdtString.c_str(), splits[i].c_str(), which.c_str(), ibin), 4)
	   << endl;
      //       fTEX << Form("\\vdef{%s:%s_%s:%d:name} {%s}", fSuffix.c_str(), fBdtString.c_str(), which.c_str(), ibin, texlabel.c_str()) << endl;
      //       fTEX << formatTex(h->GetBinContent(ibin), Form("%s:%s_%s:%d:val", fSuffix.c_str(), fBdtString.c_str(), which.c_str(), ibin), 4)
      // 	   << endl;
    }
  }

  cout << "==> Read BDT variables" << endl;
  for (unsigned int i = 0; i < fReaderVariables.size(); ++i) {
    cout << "  " << fReaderVariables[i] << endl;
  }

}


// ----------------------------------------------------------------------
void plotBDT::ssb() {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  TTree *t(0); 
  t = (TTree*)fRootFile->Get("bdtTree");
  if (0 == t) {
    cout << "%%%%%%% bdtTree not found in " << fRootFile->GetName() << endl;
    return;
  }

  int bdtBins(200); 
  double bdtMin(-1.0), bdtMax(1.0); 
  TH1D *sm = new TH1D("sm", "m(signal)", 100, 4.9, 5.9); sm->Sumw2();
  TH1D *dm = new TH1D("dm", "m(data)", 100, 4.9, 5.9); dm->Sumw2();
  TH1D *h  = new TH1D("s1", "S/Sqrt(S+B)", bdtBins, bdtMin, bdtMax); h->Sumw2();
  setTitles(h, "b >", "S/#sqrt{S+B}"); 

  double bdt, m, w8; 
  int classID;
  bool gmuid, hlt;
  t->SetBranchAddress("bdt", &bdt);
  t->SetBranchAddress("classID", &classID);
  t->SetBranchAddress("m", &m);
  t->SetBranchAddress("weight", &w8);
  t->SetBranchAddress("hlt", &hlt);
  t->SetBranchAddress("gmuid", &gmuid);

  // -- compute S and B
  double bdtCut, maxSSB(-1.), maxBDT(-1.), maxSSBsimple(-1.); 
  int nEvent(0); 
  //    for (int ibin = 80; ibin < 82; ++ibin) {
  for (int ibin = 0; ibin < bdtBins; ++ibin) {
    bdtCut = bdtMin + ibin*(bdtMax-bdtMin)/bdtBins;
    nEvent = t->GetEntries();
    cout << "=============>  bin " << ibin << " cutting at bdt > " << bdtCut << " with " << nEvent << " entries" << endl;
    sm->Reset();
    dm->Reset();
    for (Long64_t ievt=0; ievt<nEvent; ievt++) {
      t->GetEntry(ievt);
      if (false == hlt) continue;
      if (false == gmuid) continue;
      if (bdt < bdtCut) continue;
      
      if (1 == classID && 5.2<m&&m<5.45) continue;

      if (0 == classID) {
	if (fYear == 2012) {
	  sm->Fill(m, w8*0.6); // FIXME
	} else {
	  sm->Fill(m, w8); // FIXME
	}	  
      } else {
	dm->Fill(m, w8); 
      }
      
    }

    double s = sm->Integral(sm->FindBin(5.3), sm->FindBin(5.45));
    double pbg = 0.07*s;
    double d0 = dm->Integral(dm->FindBin(4.9), dm->FindBin(5.3));
    double d1 = dm->Integral(dm->FindBin(5.45), dm->GetNbinsX());
    double d = dm->Integral(0, dm->GetNbinsX());
    double bsimple = d*(5.45-5.30)/(5.9-4.9-0.25);

    bgBlind(dm, 3, 4.9, 5.9);
    double b = fBsBgExp;
    if (s+b >0) {
      double ssb = s/TMath::Sqrt(s+b+pbg);
      double ssbsimple = s/TMath::Sqrt(s+bsimple);
      if (ssbsimple > maxSSBsimple) {
	maxSSBsimple = ssbsimple; 
      }
      if (ssb > maxSSB) {
	maxSSB = ssb; 
	maxBDT = bdtCut;
      }
      h->SetBinContent(ibin, ssb); 
      //      h->SetBinError(ibin, TMath::Abs(ssb-ssbsimple)); 
    } else {
      h->SetBinContent(ibin, 0); 
    }
    //    cout << "S = " << s << " B = " << b << " => S/sqrt(S+B) = " << s/TMath::Sqrt(s+b) << endl;
    c0->Clear();
    dm->SetTitle(Form("b >%f", bdtCut));
    dm->Draw("e");
    sm->Draw("samehist");
    c0->Modified();
    c0->Update();
    if (s < 1e-6) break;
  }
  fTEX << formatTex(maxSSB, Form("%s:%s:maxSSB:val",  fSuffix.c_str(), fBdtString.c_str()), 2) << endl;
  //  fTEX << formatTex(maxSSBsimple, Form("%s:%s:maxSSBsimple:val",  fSuffix.c_str(), fBdtString.c_str()), 2) << endl;
  fTEX << formatTex(maxBDT, Form("%s:%s:maxSSB:bdt",  fSuffix.c_str(), fBdtString.c_str()), 2) << endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  

  h->SetMaximum(maxSSB+0.5); 
  h->Draw();
  //  tl->DrawLatex(0.25, 0.75, Form("S_{max} = %4.3f (%4.3f)", maxSSB, maxSSBsimple));
  tl->DrawLatex(0.25, 0.75, Form("S_{max} = %4.3f", maxSSB));
  tl->DrawLatex(0.25, 0.68, Form("B_{max} = %4.3f", maxBDT));
  c0->SaveAs(Form("%s/%s-ssb.pdf", fDirectory.c_str(), fBdtString.c_str())); 

  double xmax(0.), xmin(0.); 
  int nbins(0); 
  double maxVal = -1.;
  for (int i = 1; i < h->GetNbinsX(); ++i) 
    if (h->GetBinContent(i) > maxVal) maxVal = h->GetBinContent(i);

  for (int i = 1; i < h->GetNbinsX(); ++i) {
    if (h->GetBinContent(i) > 0.7*maxVal) {
      xmax = h->GetBinCenter(i); 
      nbins = i - h->GetMaximumBin(); 
    }
    h->SetBinError(i, 0.03*h->GetBinContent(i)); 
  }
  xmin = h->GetBinCenter(h->GetMaximumBin() - TMath::Abs(nbins)); 

  cout << "maxval: " << h->GetMaximum() << endl;
  cout << "maxbin: " << h->GetMaximumBin() << endl;
  cout << "xmax: " << xmax << endl;
  cout << "xmin: " << xmin << endl;
  cout << "nbins: " << nbins << endl;
  
  TF1 *f1 = fpFunc->pol2local(h, 0.05); 
  h->Fit(f1, "r", "", xmin, xmax); 
  double maxfitssbX = h->GetFunction("iF_pol2local")->GetParameter(2); 
  double maxfitssbY = h->GetFunction("iF_pol2local")->GetParameter(0); 

  fTEX << formatTex(maxfitssbX, Form("%s:%s:maxfitSSB:val",  fSuffix.c_str(), fBdtString.c_str()), 2) << endl;
  fTEX << formatTex(maxfitssbY, Form("%s:%s:maxfitSSB:bdt",  fSuffix.c_str(), fBdtString.c_str()), 2) << endl;

  c0->SaveAs(Form("%s/%s-fit-ssb.pdf", fDirectory.c_str(), fBdtString.c_str())); 


  if (fYear == 2012) {
    cout << "FIXME FIXME: patched MC weight with 0.6 to account for gen-level filters" << endl;
  }
}


// ----------------------------------------------------------------------
void plotBDT::plotEffVsBg(int offset) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  
  zone(1);
  shrinkPad(0.15, 0.15, 0.18); 
  double sf(1.1);
  
  //  TFile *f = TFile::Open("default/anaBmm.plotBDT.default-11.root");

  TH1D *h0 = (TH1D*)(gDirectory->Get(Form("hBgEvts%d", offset)))->Clone("h0"); 
  TH1D *h1 = (TH1D*)(gDirectory->Get(Form("hSgEff%d", offset)))->Clone("h1"); 

  h0->SetMaximum(100);
  h0->SetMinimum(0.5);
  double xmax = 1.0;
  if (h0->FindLastBinAbove(0.0) < 0.7) {
    cout << "reset last bin because " << h1->FindLastBinAbove(0.0) << endl;
    xmax = 0.8;
  }

  h0->SetAxisRange(-0.2, xmax, "X");
  setTitles(h0, "b >", "# background events", 0.05, 1.2, 1.5); 
  h0->Draw("e");
  c0->Update();

  double rightmax = sf*h1->GetMaximum();
  double scale = gPad->GetUymax()/rightmax;
  cout << rightmax << " " << gPad->GetUymax() << " " << scale << endl;
  h1->Scale(scale);
  h1->SetLineColor(kRed);
  h1->Draw("samehist");

  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetTitle("Efficiency");
  axis->SetTitleColor(kRed);
  axis->SetTitleOffset(1.8);
  axis->SetTitleFont(42);
  axis->SetLabelFont(42);
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetLabelColor(kRed);
  axis->Draw();

  
  c0->SaveAs(Form("%s/bdt-SgEff-BgEvts-%s.pdf", fDirectory.c_str(), fBdtString.c_str())); 

  // -- dump table of numbers
  h0 = (TH1D*)(gDirectory->Get(Form("hBgEvts%d", offset)))->Clone("H0"); 
  h1 = (TH1D*)(gDirectory->Get(Form("hSgEff%d", offset)))->Clone("H1"); 
  int icnt(0); 
  double eff(h1->GetBinContent(1)), step(0.002); 
  for (int i = 1; i < h0->GetNbinsX(); ++i) {
    if (h1->GetBinContent(i) > eff - icnt*step) continue;
    cout << " icnt = " << icnt << " i = " << i << " " << h0->GetBinContent(i) << " " << h1->GetBinContent(i) << endl;
    fTEX << formatTex(h0->GetBinLowEdge(i), Form("%s:bdtChan%d:Bdt:%d",  fSuffix.c_str(), offset, icnt), 2) << endl;
    fTEX << formatTex(h0->GetBinContent(i), Form("%s:bdtChan%d:BgCnt:%d",  fSuffix.c_str(), offset, icnt), 0) << endl;
    fTEX << formatTex(h1->GetBinContent(i), Form("%s:bdtChan%d:SgEff:%d",  fSuffix.c_str(), offset, icnt), 3) << endl;
    ++icnt;
    if (h1->GetBinContent(i) < 1e-5) break;
  }
}



// ----------------------------------------------------------------------
void plotBDT::hackedMC(int chan) {
  if (0 == chan) fSetup = "B";
  if (1 == chan) fSetup = "E";
  
  string mode; 
  TTree *t(0); 

  resetHistograms();
  
  // -- Bs2MuMu
  mode = "SgBs2MuMu";
  fMode = 0; 
  t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  setupTree(t, mode);
  loopOverTree(t, mode, 2); 

  // -- Bx2MuMu
  mode = "SgBx2MuMu";
  fMode = 1; 
  t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  setupTree(t, mode);
  loopOverTree(t, mode, 2); 

  // -- By2MuMu
  mode = "SgBy2MuMu";
  fMode = 2; 
  t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  setupTree(t, mode);
  loopOverTree(t, mode, 2); 

  // -- Bs2JpsiPhi
  mode = "CsMc";
  fMode = 3; 
  t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  setupTree(t, mode);
  loopOverTree(t, mode, 2); 

  // -- Bu2JpsiK
  mode = "NoMc";
  fMode = 4; 
  t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  setupTree(t, mode);
  loopOverTree(t, mode, 2); 

  // ----------------------------------------------------------------------
  // -- hacked MC1 analysis/display
  // ----------------------------------------------------------------------
  string fitstring("pol0"); 
  gStyle->SetOptTitle(0); 
  zone(2,2);
  fhMcBDT[1]->Scale(fhMcBDT[0]->GetSumOfWeights()/fhMcBDT[1]->GetSumOfWeights()); 
  fhMcBDT[2]->Scale(fhMcBDT[0]->GetSumOfWeights()/fhMcBDT[2]->GetSumOfWeights()); 

  fhMcMass[1]->Scale(fhMcMass[0]->GetSumOfWeights()/fhMcMass[1]->GetSumOfWeights()); 
  fhMcMass[2]->Scale(fhMcMass[0]->GetSumOfWeights()/fhMcMass[2]->GetSumOfWeights()); 

  fhMcRatio[1]->Divide(fhMcBDT[1], fhMcBDT[0], 1., 1.); 
  fhMcRatio[2]->Divide(fhMcBDT[2], fhMcBDT[0], 1., 1.); 

  c0->cd(1);
  double ymax = fhMcMass[0]->GetMaximum();
  if (fhMcMass[1]->GetMaximum() > ymax) ymax = fhMcMass[1]->GetMaximum();
  if (fhMcMass[2]->GetMaximum() > ymax) ymax = fhMcMass[2]->GetMaximum();
  fhMcMass[0]->SetMaximum(1.2*ymax); 
  fhMcMass[0]->Draw();
  fhMcMass[1]->Draw("same");
  fhMcMass[2]->Draw("same");
  tl->SetTextColor(kRed);   tl->DrawLatex(0.77, 0.6, "1");
  tl->SetTextColor(kBlack); tl->DrawLatex(0.60, 0.6, "0");
  tl->SetTextColor(kBlue);  tl->DrawLatex(0.48, 0.6, "2");

  c0->cd(2);
  fhMcBDT[0]->Draw("e");
  fhMcBDT[1]->Draw("samehist");
  fhMcBDT[2]->Draw("samehist");

  c0->cd(3);
  fhMcRatio[1]->Draw("e");
  fhMcRatio[1]->SetMaximum(2.);  fhMcRatio[1]->SetMinimum(0.);
  fhMcRatio[1]->Fit(fitstring.c_str(), "", "same"); 
  fhMcRatio[1]->GetFunction(fitstring.c_str())->SetLineColor(kRed); 
  fhMcRatio[1]->GetFunction(fitstring.c_str())->SetLineWidth(2); 
  double chi2 = fhMcRatio[1]->GetFunction(fitstring.c_str())->GetChisquare(); 
  int    ndof = fhMcRatio[1]->GetFunction(fitstring.c_str())->GetNDF();
  tl->SetTextColor(kBlack);  
  tl->DrawLatex(0.25, 0.92, Form("#chi^{2}/dof = %3.1f/%d", chi2, ndof)); 

  c0->cd(4);
  fhMcRatio[2]->Draw("e");
  fhMcRatio[2]->SetMaximum(2.);  fhMcRatio[2]->SetMinimum(0.);
  fhMcRatio[2]->Fit(fitstring.c_str(), "", "same"); 
  fhMcRatio[2]->GetFunction(fitstring.c_str())->SetLineColor(kBlue); 
  fhMcRatio[2]->GetFunction(fitstring.c_str())->SetLineWidth(2); 
  chi2 = fhMcRatio[2]->GetFunction(fitstring.c_str())->GetChisquare(); 
  ndof = fhMcRatio[2]->GetFunction(fitstring.c_str())->GetNDF();
  tl->SetTextColor(kBlack);  
  tl->DrawLatex(0.25, 0.92, Form("#chi^{2}/dof = %3.1f/%d", chi2, ndof)); 

  c0->SaveAs(Form("%s/hackedMC-bdt-for-shiftedSignalMC-%s-chan%d.pdf", fDirectory.c_str(), fitstring.c_str(), chan));


  // ----------------------------------------------------------------------
  // -- hacked MC2 analysis/display
  // ----------------------------------------------------------------------
  for (int i = 0; i < 3; ++i) {
    fhMcMass1[i]->Scale(fhMcMass0[i]->GetSumOfWeights()/fhMcMass1[i]->GetSumOfWeights());
    fhMcMass2[i]->Scale(fhMcMass0[i]->GetSumOfWeights()/fhMcMass2[i]->GetSumOfWeights());
    fhMcRatio1[i]->Divide(fhMcMass1[i], fhMcMass0[i], 1., 1.); 
    fhMcRatio2[i]->Divide(fhMcMass2[i], fhMcMass0[i], 1., 1.); 
  }

  // -- Bs
  string func("pol0"); 
  makeCanvas(1);
  
  string pdfname; 
  int color(kBlack); 
  for (int i = 0; i < 3; ++i) {
    zone(2, 1, c1);
    fhMcMass0[i]->SetMaximum(1.3*fhMcMass0[i]->GetMaximum());
    fhMcMass0[i]->SetMinimum(0.1);
    setTitles(fhMcMass0[i], "m [GeV]", "candidates", 0.05, 1.1, 1.5); 
    fhMcMass0[i]->Draw("hist");
    fhMcMass1[i]->Draw("samee");
    fhMcMass2[i]->Draw("samee0");
    
    newLegend(0.25, 0.7, 0.45, 0.85); 
    legg->AddEntry(fhMcMass0[i], Form("b>%3.1f", 0.0), "l"); 
    legg->AddEntry(fhMcMass1[i], Form("b>%3.1f", 0.1), "p"); 
    legg->AddEntry(fhMcMass2[i], Form("b>%3.1f", 0.2), "p"); 
    legg->Draw();
    
    if (1 == i ) color = kRed;
    if (2 == i ) color = kBlue;

    c1->cd(2);  
    if (0 == i) fhMcRatio1[i]->SetAxisRange(4.9, 5.6, "X"); 
    if (1 == i) fhMcRatio1[i]->SetAxisRange(5.4, 5.9, "X"); 
    if (2 == i) fhMcRatio1[i]->SetAxisRange(4.9, 5.4, "X"); 

    setTitles(fhMcRatio1[i], "m [GeV]", "ratio"); 
    fhMcRatio1[i]->SetMinimum(0.);    fhMcRatio1[i]->SetMaximum(2.);
    fhMcRatio1[i]->Fit(func.c_str(), "", "e1");
    fhMcRatio2[i]->Fit(func.c_str(), "", "same");
    fhMcRatio1[i]->GetFunction(func.c_str())->SetLineColor(color); 
    fhMcRatio2[i]->GetFunction(func.c_str())->SetLineColor(color);   fhMcRatio2[i]->GetFunction(func.c_str())->SetLineStyle(kDashed); 
    
    tl->SetTextSize(0.04);
    tl->DrawLatex(0.15, 0.92, Form("b>0.2: #chi^{2}/dof = %3.1f/%d", 
				   fhMcRatio1[i]->GetFunction(func.c_str())->GetChisquare(), 
				   fhMcRatio1[i]->GetFunction(func.c_str())->GetNDF()));
    
    tl->DrawLatex(0.55, 0.92, Form("b>0.3: #chi^{2}/dof = %3.1f/%d", 
				   fhMcRatio2[i]->GetFunction(func.c_str())->GetChisquare(), 
				   fhMcRatio2[i]->GetFunction(func.c_str())->GetNDF()));
    
    pdfname = "bs"; 
    if (1 == i) pdfname = "bx"; 
    if (2 == i) pdfname = "by"; 
    c1->SaveAs(Form("%s/hackedMC-massratio-%s-chan%d.pdf", fDirectory.c_str(), pdfname.c_str(), chan));
  }



  // ----------------------------------------------------------------------
  // -- hacked MC3 analysis/display
  // ----------------------------------------------------------------------
  gStyle->SetOptTitle(0); 
  zone(1);
  fhMcBDT5[1]->Scale(fhMcBDT5[0]->GetSumOfWeights()/fhMcBDT5[1]->GetSumOfWeights());
  fhMcBDT5[2]->Scale(fhMcBDT5[0]->GetSumOfWeights()/fhMcBDT5[2]->GetSumOfWeights());
  fhMcBDT5[3]->Scale(fhMcBDT5[0]->GetSumOfWeights()/fhMcBDT5[3]->GetSumOfWeights());
  fhMcBDT5[4]->Scale(fhMcBDT5[0]->GetSumOfWeights()/fhMcBDT5[4]->GetSumOfWeights());

  zone();
  fhMcBDT5[0]->SetMaximum(1.2*fhMcBDT5[0]->GetMaximum());
  fhMcBDT5[0]->Draw("");
  fhMcBDT5[1]->Draw("samehist");
  fhMcBDT5[2]->Draw("samehist");
  fhMcBDT5[3]->Draw("samehist");
  fhMcBDT5[4]->Draw("samehist");

  newLegend(0.25, 0.6, 0.45, 0.85); 
  legg->AddEntry(fhMcBDT5[0], "B_{s} #rightarrow #mu #mu", "p"); 
  legg->AddEntry(fhMcBDT5[1], "B_{x} #rightarrow #mu #mu", "l"); 
  legg->AddEntry(fhMcBDT5[2], "B_{y} #rightarrow #mu #mu", "l"); 
  legg->AddEntry(fhMcBDT5[3], "B_{s} #rightarrow J/#psi #phi", "f"); 
  legg->AddEntry(fhMcBDT5[4], "B^{+} #rightarrow J/#psi K", "f"); 
  legg->Draw();

  c0->SaveAs(Form("%s/hackedMC-bdt-overlays-chan%d.pdf", fDirectory.c_str(), chan));

}


// ----------------------------------------------------------------------
void plotBDT::loopFunction2(int mode) {
  if (fChan < 0) return;

  if (fSetup == "B" && fChan != 0) return;
  if (fSetup == "E" && fChan != 1) return;
  
  if (!fGoodAcceptance) return;
  if (!fGoodTracks) return;
  if (!fGoodTracksPt) return;
  if (!fGoodTracksEta) return;
  if (!fGoodMuonsPt) return;
  if (!fGoodMuonsEta) return;
  if (!fGoodJpsiCuts) return;
  if (!fGoodMuonsID) return;
  if (!fGoodHLT) return;

  if (fBDT < -1.) return; 

  if (fMode < 3) {
    fhMcMass[fMode]->Fill(fb.m); 
    if (0 == fMode) {
      if (fb.m > 5.30) fhMcBDT[fMode]->Fill(fBDT); 
    }
    
    if (1 == fMode) {
      if (fb.m > 5.63) fhMcBDT[fMode]->Fill(fBDT); 
    }
    
    if (2 == fMode) {
      if (fb.m > 5.03) fhMcBDT[fMode]->Fill(fBDT); 
    }
    
    if (fBDT>0) fhMcMass0[fMode]->Fill(fb.m);
    if (fBDT>0.1) fhMcMass1[fMode]->Fill(fb.m);
    if (fBDT>0.2) fhMcMass2[fMode]->Fill(fb.m);
  }
  
  // -- fill BDT response for hacked MC3 analysis (overlay of three signals with Jpsi X)
  if (fb.flsxy < 3) return;
  if (fb.pchi2dof < 0.1) return;
  if (fb.pt < 6.9) return;
  if (fb.maxdoca > 0.5) return;
  if (1) {
    fhMcBDT5[fMode]->Fill(fBDT);
  }
}




// ----------------------------------------------------------------------
void plotBDT::hackedMC3(int chan) {

  gStyle->SetOptTitle(0);

  string defaultCuts = "hlt&&gmuid&&m1pt>4&&m2pt>4&&flsxy>3&&pt>6.9&&alpha<0.45&&prob>0.1&&maxdoca<0.5";
  string chanCuts = (0 == chan? "(abs(m1eta)<1.4&&abs(m2eta)<1.4)" : "!(abs(m1eta)<1.4&&abs(m2eta)<1.4)");

  fF["SgBs2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hs0 = new TH1D("hs0", "B_{s} #rightarrow #mu #mu", 100, -0.6, 0.6); hs0->Sumw2(); setHist(hs0, kBlack); 
  setTitles(hs0, "b", "# candidates", 0.05, 1.1, 1.5);
  TTree *t0 = (TTree*)gDirectory->Get("events"); 
  t0->Draw("bdt>>hs0", Form("%s&&%s", defaultCuts.c_str(), chanCuts.c_str()));

  fF["SgBx2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hx0 = new TH1D("hx0", "B_{x} #rightarrow #mu #mu", 100, -0.6, 0.6); hx0->Sumw2(); setHist(hx0, kBlue); 
  t0 = (TTree*)gDirectory->Get("events"); 
  t0->Draw("bdt>>hx0", Form("%s&&%s", defaultCuts.c_str(), chanCuts.c_str()));

  fF["SgBy2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hy0 = new TH1D("hy0", "B_{y} #rightarrow #mu #mu", 100, -0.6, 0.6); hy0->Sumw2(); setHist(hy0, kRed); 
  t0 = (TTree*)gDirectory->Get("events"); 
  t0->Draw("bdt>>hy0", Form("%s&&%s", defaultCuts.c_str(), chanCuts.c_str()));

  fF["CsMc"]->cd("candAnaBs2JpsiPhi"); 
  TH1D *ha0 = new TH1D("ha0", "B_{s} #rightarrow J/#psi #phi", 100, -0.6, 0.6); ha0->Sumw2(); setFilledHist(ha0, 35, 35, 3365); 
  t0 = (TTree*)gDirectory->Get("events"); 
  t0->Draw("bdt>>ha0", Form("%s&&%s", defaultCuts.c_str(), chanCuts.c_str()));

  fF["NoMc"]->cd("candAnaBu2JpsiK"); 
  TH1D *hb0 = new TH1D("hb0", "B^{+} #rightarrow J/#psi K^{+}", 100, -0.6, 0.6); hb0->Sumw2(); setFilledHist(hb0, 45, 45, 3356); 
  t0 = (TTree*)gDirectory->Get("events"); 
  t0->Draw("bdt>>hb0", Form("%s&&%s", defaultCuts.c_str(), chanCuts.c_str()));

  hb0->Scale(hs0->GetSumOfWeights()/hb0->GetSumOfWeights());
  ha0->Scale(hs0->GetSumOfWeights()/ha0->GetSumOfWeights());
  hx0->Scale(hs0->GetSumOfWeights()/hx0->GetSumOfWeights());
  hy0->Scale(hs0->GetSumOfWeights()/hy0->GetSumOfWeights());

  zone();
  hs0->SetMaximum(1.2*hs0->GetMaximum());
  hs0->Draw("");
  hx0->Draw("samehist");
  hy0->Draw("samehist");
  ha0->Draw("samehist");
  hb0->Draw("samehist");

  newLegend(0.25, 0.6, 0.45, 0.85); 
  legg->AddEntry(hs0, "B_{s} #rightarrow #mu #mu", "p"); 
  legg->AddEntry(hx0, "B_{x} #rightarrow #mu #mu", "l"); 
  legg->AddEntry(hy0, "B_{y} #rightarrow #mu #mu", "l"); 
  legg->AddEntry(ha0, "B_{s} #rightarrow J/#psi #phi", "f"); 
  legg->AddEntry(hb0, "B^{+} #rightarrow J/#psi K", "f"); 
  legg->Draw();

  c0->SaveAs(Form("%s/hackedMC-bdt-overlays-%d.pdf", fDirectory.c_str(), chan));
  

}



// ----------------------------------------------------------------------
void plotBDT::hackedMC2(double bdtCut0, double bdtCut1, std::string func) {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  string defaultCuts = "hlt&&gmuid&&m1pt>4&&m2pt>4&&abs(m1eta)<2.4&&abs(m2eta)<2.4";

  TH1D *hrs1 = new TH1D("hrs1", Form("ratio"), 100, 4.5, 6.0); hrs1->Sumw2(); setHist(hrs1, kBlack, 20); 
  TH1D *hrs2 = new TH1D("hrs2", Form("ratio"), 100, 4.5, 6.0); hrs2->Sumw2(); setHist(hrs2, kBlack, 24);  

  TH1D *hrx1 = new TH1D("hrx1", Form("ratio"), 100, 4.5, 6.0); hrx1->Sumw2(); setHist(hrx1, kRed, 20); 
  TH1D *hrx2 = new TH1D("hrx2", Form("ratio"), 100, 4.5, 6.0); hrx2->Sumw2(); setHist(hrx2, kRed, 24);  

  TH1D *hry1 = new TH1D("hry1", Form("ratio"), 100, 4.5, 6.0); hry1->Sumw2(); setHist(hry1, kBlue, 20);  
  TH1D *hry2 = new TH1D("hry2", Form("ratio"), 100, 4.5, 6.0); hry2->Sumw2(); setHist(hry2, kBlue, 24);  

  fF["SgBs2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hs0 = new TH1D("hs0", "B_{s} #rightarrow #mu #mu", 100, 4.5, 6.0); hs0->Sumw2(); setHist(hs0, kBlack); 
  TH1D *hs1 = new TH1D("hs1", "B_{s} #rightarrow #mu #mu", 100, 4.5, 6.0); hs1->Sumw2(); setHist(hs1, kBlack); 
  TH1D *hs2 = new TH1D("hs2", "B_{s} #rightarrow #mu #mu", 100, 4.5, 6.0); hs2->Sumw2(); setHist(hs2, kBlack, 24); 
  TTree *t0 = (TTree*)gDirectory->Get("events"); 
  t0->Draw("m>>hs0", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.));
  t0->Draw("m>>hs1", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.2));
  t0->Draw("m>>hs2", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.3));
	  
  fF["SgBx2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hx0 = new TH1D("hx0", "B_{x} #rightarrow #mu #mu ", 100, 4.5, 6.0); hx0->Sumw2(); setHist(hx0, kRed); 
  TH1D *hx1 = new TH1D("hx1", "B_{x} #rightarrow #mu #mu ", 100, 4.5, 6.0); hx1->Sumw2(); setHist(hx1, kRed); 
  TH1D *hx2 = new TH1D("hx2", "B_{x} #rightarrow #mu #mu ", 100, 4.5, 6.0); hx2->Sumw2(); setHist(hx2, kRed, 24); 
  TTree *t1 = (TTree*)gDirectory->Get("events"); 
  t1->Draw("m>>hx0", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.));
  t1->Draw("m>>hx1", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.2));
  t1->Draw("m>>hx2", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.3));

  fF["SgBy2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hy0 = new TH1D("hy0", "B_{y} #rightarrow #mu #mu ", 100, 4.5, 6.0); hy0->Sumw2(); setHist(hy0, kBlue); 
  TH1D *hy1 = new TH1D("hy1", "B_{y} #rightarrow #mu #mu ", 100, 4.5, 6.0); hy1->Sumw2(); setHist(hy1, kBlue); 
  TH1D *hy2 = new TH1D("hy2", "B_{y} #rightarrow #mu #mu ", 100, 4.5, 6.0); hy2->Sumw2(); setHist(hy2, kBlue, 24); 
  TTree *t2 = (TTree*)gDirectory->Get("events"); 
  t2->Draw("m>>hy0", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.));
  t2->Draw("m>>hy1", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.2));
  t2->Draw("m>>hy2", Form("%s&&bdt>%f", defaultCuts.c_str(), 0.3));

  hs1->Scale(hs0->GetSumOfWeights()/hs1->GetSumOfWeights());
  hs2->Scale(hs0->GetSumOfWeights()/hs2->GetSumOfWeights());
  hrs1->Divide(hs1, hs0, 1., 1.);
  hrs2->Divide(hs2, hs0, 1., 1.);

  hx1->Scale(hx0->GetSumOfWeights()/hx1->GetSumOfWeights());
  hx2->Scale(hx0->GetSumOfWeights()/hx2->GetSumOfWeights());
  hrx1->Divide(hx1, hx0, 1., 1.);
  hrx2->Divide(hx2, hx0, 1., 1.);

  hy1->Scale(hy0->GetSumOfWeights()/hy1->GetSumOfWeights());
  hy2->Scale(hy0->GetSumOfWeights()/hy2->GetSumOfWeights());
  hry1->Divide(hy1, hy0, 1., 1.);
  hry2->Divide(hy2, hy0, 1., 1.);

  // -- Bs
  makeCanvas(1);
  zone(2,1, c1);
  hs0->SetMaximum(1.3*hs0->GetMaximum());
  setTitles(hs0, "m [GeV]", "candidates", 0.05, 1.1, 1.5); 
  hs0->Draw();
  hs1->Draw("samee");
  hs2->Draw("samee0");

  newLegend(0.25, 0.7, 0.45, 0.85); 
  legg->AddEntry(hs0, Form("b>%3.1f", 0.0), "l"); 
  legg->AddEntry(hs1, Form("b>%3.1f", 0.2), "p"); 
  legg->AddEntry(hs2, Form("b>%3.1f", 0.3), "p"); 
  legg->Draw();

  c1->cd(2);  
  hrs1->SetAxisRange(4.9, 5.6, "X"); 
  setTitles(hrs1, "m [GeV]", "ratio"); 
  hrs1->Fit(func.c_str(), "", "e1");
  hrs2->Fit(func.c_str(), "", "same");
  hrs1->GetFunction(func.c_str())->SetLineColor(kBlack); 
  hrs2->GetFunction(func.c_str())->SetLineColor(kBlack);   hrs2->GetFunction(func.c_str())->SetLineStyle(kDashed); 
  
  tl->SetTextSize(0.04);
  tl->DrawLatex(0.15, 0.92, Form("b>0.2: #chi^{2}/dof = %3.1f/%d", 
				hrs1->GetFunction(func.c_str())->GetChisquare(), hrs1->GetFunction(func.c_str())->GetNDF()));

  tl->DrawLatex(0.55, 0.92, Form("b>0.3: #chi^{2}/dof = %3.1f/%d", 
				hrs2->GetFunction(func.c_str())->GetChisquare(), hrs2->GetFunction(func.c_str())->GetNDF()));

  c1->SaveAs(Form("%s/hackedMC-massratio-bs.pdf", fDirectory.c_str()));

  // -- Bx
  zone(2,1, c1);
  hx0->SetMaximum(1.3*hx0->GetMaximum());
  setTitles(hx0, "m [GeV]", "candidates", 0.05, 1.1, 1.5); 
  hx0->Draw();
  hx1->Draw("samee");
  hx2->Draw("samee0");

  newLegend(0.25, 0.7, 0.45, 0.85); 
  legg->AddEntry(hx0, Form("b>%3.1f", 0.0), "l"); 
  legg->AddEntry(hx1, Form("b>%3.1f", 0.2), "p"); 
  legg->AddEntry(hx2, Form("b>%3.1f", 0.3), "p"); 
  legg->Draw();

  c1->cd(2);  
  setTitles(hrx1, "m [GeV]", "ratio"); 
  hrx1->SetAxisRange(5.4, 5.9, "X"); 
  hrx1->Fit(func.c_str(), "", "e1");
  hrx2->Fit(func.c_str(), "", "same");
  hrx1->GetFunction(func.c_str())->SetLineColor(kRed); 
  hrx2->GetFunction(func.c_str())->SetLineColor(kRed);   hrx2->GetFunction(func.c_str())->SetLineStyle(kDashed); 
  
  tl->SetTextSize(0.04);
  tl->DrawLatex(0.15, 0.92, Form("b>0.2: #chi^{2}/dof = %3.1f/%d", 
				hrx1->GetFunction(func.c_str())->GetChisquare(), hrx1->GetFunction(func.c_str())->GetNDF()));

  tl->DrawLatex(0.55, 0.92, Form("b>0.3: #chi^{2}/dof = %3.1f/%d", 
				hrx2->GetFunction(func.c_str())->GetChisquare(), hrx2->GetFunction(func.c_str())->GetNDF()));

  c1->SaveAs(Form("%s/hackedMC-massratio-bx.pdf", fDirectory.c_str()));


  // -- By
  zone(2,1, c1);
  hy0->SetMaximum(1.3*hy0->GetMaximum());
  setTitles(hy0, "m [GeV]", "candidates", 0.05, 1.1, 1.5); 
  hy0->Draw();
  hy1->Draw("samee");
  hy2->Draw("samee0");

  newLegend(0.25, 0.7, 0.45, 0.85); 
  legg->AddEntry(hy0, Form("b>%3.1f", 0.0), "l"); 
  legg->AddEntry(hy1, Form("b>%3.1f", 0.2), "p"); 
  legg->AddEntry(hy2, Form("b>%3.1f", 0.3), "p"); 
  legg->Draw();

  c1->cd(2);  
  setTitles(hry1, "m [GeV]", "ratio"); 
  hry1->SetAxisRange(4.9, 5.4, "X"); 
  hry1->Fit(func.c_str(), "", "e1");
  hry2->Fit(func.c_str(), "", "same");
  hry1->GetFunction(func.c_str())->SetLineColor(kBlue); 
  hry2->GetFunction(func.c_str())->SetLineColor(kBlue);   hry2->GetFunction(func.c_str())->SetLineStyle(kDashed); 
  
  tl->SetTextSize(0.04);
  tl->DrawLatex(0.15, 0.92, Form("b>0.2: #chi^{2}/dof = %3.1f/%d", 
				hry1->GetFunction(func.c_str())->GetChisquare(), hry1->GetFunction(func.c_str())->GetNDF()));

  tl->DrawLatex(0.55, 0.92, Form("b>0.3: #chi^{2}/dof = %3.1f/%d", 
				hry2->GetFunction(func.c_str())->GetChisquare(), hry2->GetFunction(func.c_str())->GetNDF()));

  c1->SaveAs(Form("%s/hackedMC-massratio-by.pdf", fDirectory.c_str()));
  

}

// ----------------------------------------------------------------------
void plotBDT::hackedMC1(string cuts, double xmin, double xmax, string func) {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);


  string defaultCuts = "hlt&&gmuid&&m1pt>4&&m2pt>4&&abs(m1eta)<2.4&&abs(m2eta)<2.4";

  int nbins = static_cast<int>((xmax - xmin)/0.05);

  TH1D *hsr1 = new TH1D("hsr1", "ratio 1/0", nbins, xmin, xmax); hsr1->Sumw2(); hsr1->SetMarkerColor(kRed); 
  TH1D *hsr2 = new TH1D("hsr2", "ratio 2/0", nbins, xmin, xmax); hsr2->Sumw2(); hsr2->SetMarkerColor(kBlue);
  TH1D *hsr3 = new TH1D("hsr3", "ratio 1/0", nbins, xmin, xmax); hsr3->Sumw2(); hsr3->SetMarkerColor(kRed); 
  TH1D *hsr4 = new TH1D("hsr4", "ratio 2/0", nbins, xmin, xmax); hsr4->Sumw2(); hsr4->SetMarkerColor(kBlue);

  fF["SgBs2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hs0 = new TH1D("hs0", "BDT response ", nbins, xmin, xmax); setHist(hs0, kBlack, 24); 
  TH1D *hm0 = new TH1D("hm0", "Masses", 100, 4.5, 6.0); setHist(hm0, kBlack, 24); 
  TTree *t0 = (TTree*)gDirectory->Get("events"); 
  t0->Draw("bdt>>hs0", Form("m>5.3&&%s&&%s", defaultCuts.c_str(), cuts.c_str()));
  t0->Draw("m>>hm0", Form("%s&&%s", defaultCuts.c_str(), cuts.c_str()));
	  
  fF["SgBx2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hs1 = new TH1D("hs1", "Bx -> mu mu ", nbins, xmin, xmax);  setHist(hs1, kRed, kRed); 
  TH1D *hm1 = new TH1D("hm1", "Bx -> mu mu ", 100, 4.5, 6.0); setHist(hm1, kRed, kRed); 
  TTree *t1 = (TTree*)gDirectory->Get("events"); 
  t1->Draw("bdt>>hs1", Form("m>5.63&&%s&&%s", defaultCuts.c_str(), cuts.c_str()));
  t1->Draw("m>>hm1", Form("%s&&%s", defaultCuts.c_str(), cuts.c_str()));
  
  fF["SgBy2MuMu"]->cd("candAnaMuMu"); 
  TH1D *hs2 = new TH1D("hs2", "By -> mu mu ", nbins, xmin, xmax);  setHist(hs2, kBlue, kBlue); 
  TH1D *hm2 = new TH1D("hm2", "By -> mu mu ", 100, 4.5, 6.0); setHist(hm2, kBlue, kBlue); 
  TTree *t2 = (TTree*)gDirectory->Get("events"); 
  t2->Draw("bdt>>hs2", Form("m>5.03&&%s&&%s", defaultCuts.c_str(), cuts.c_str()));
  t2->Draw("m>>hm2", Form("%s&&%s", defaultCuts.c_str(), cuts.c_str()));

  //  gStyle->SetOptStat(0);
  //  gStyle->SetOptFit(1111);
  
  zone(2,3);
  hm1->Scale(hm0->GetSumOfWeights()/hm1->GetSumOfWeights()); 
  hm2->Scale(hm0->GetSumOfWeights()/hm2->GetSumOfWeights()); 
  hm0->Draw();
  hm1->Draw("same");
  hm2->Draw("same");

  tl->SetTextColor(kRed);   tl->DrawLatex(0.77, 0.6, "1");
  tl->SetTextColor(kBlack); tl->DrawLatex(0.60, 0.6, "0");
  tl->SetTextColor(kBlue);  tl->DrawLatex(0.48, 0.6, "2");

  c0->cd(2);
  hs1->Scale(hs0->GetSumOfWeights()/hs1->GetSumOfWeights()); 
  hs2->Scale(hs0->GetSumOfWeights()/hs2->GetSumOfWeights()); 

  hs0->Draw("e");
  hs1->Draw("samehist");
  hs2->Draw("samehist");
  
  c0->cd(3);
  string lfunc = "pol0"; 
  hsr1->Divide(hs1, hs0, 1., 1.);
  hsr2->Divide(hs2, hs0, 1., 1.);
  hsr1->SetMaximum(2.); hsr1->SetMinimum(0.);
  hsr1->Draw("e");
  hsr1->Fit(lfunc.c_str()); 
  hsr1->GetFunction(lfunc.c_str())->SetLineColor(kBlack); 
  hsr1->GetFunction(lfunc.c_str())->SetLineWidth(2); 
  c0->cd(4);
  hsr2->SetMaximum(2.); hsr2->SetMinimum(0.);
  hsr2->Draw("e");
  hsr2->Fit(lfunc.c_str()); 
  hsr2->GetFunction(lfunc.c_str())->SetLineColor(kBlack); 
  hsr2->GetFunction(lfunc.c_str())->SetLineWidth(2); 

  c0->cd(5);
  lfunc = "pol1";
  hsr3->Divide(hs1, hs0, 1., 1.);
  hsr4->Divide(hs2, hs0, 1., 1.);
  hsr3->SetMaximum(2.); hsr3->SetMinimum(0.);
  hsr3->Draw("e");
  hsr3->Fit(lfunc.c_str()); 
  hsr3->GetFunction(lfunc.c_str())->SetLineColor(kBlack); 
  hsr3->GetFunction(lfunc.c_str())->SetLineWidth(2); 
  c0->cd(6);
  hsr4->SetMaximum(2.); hsr4->SetMinimum(0.);
  hsr4->Draw("e");
  hsr4->Fit(lfunc.c_str()); 
  hsr4->GetFunction(lfunc.c_str())->SetLineColor(kBlack); 
  hsr4->GetFunction(lfunc.c_str())->SetLineWidth(2); 

  c0->SaveAs(Form("%s/hackedMC-bdt-for-shiftedSignalMC-%s.pdf", fDirectory.c_str(), func.c_str()));
}




// ----------------------------------------------------------------------
void plotBDT::allCorrelationPlots(double bdtcut, std::string fname) {  
  correlationPlot(bdtcut, "pvips", 0., 5., fname);
  correlationPlot(bdtcut, "pvip", 0., 0.02, fname);

//   correlationPlot(bdtcut, "m1pt", 4., 40., fname);
//   correlationPlot(bdtcut, "m2pt", 4., 30., fname);

  correlationPlot(bdtcut, "pt", 6., 40., fname);
  correlationPlot(bdtcut, "eta", -2.4, 2.4, fname);

  correlationPlot(bdtcut, "fls3d", 0, 100, fname);
  correlationPlot(bdtcut, "alpha", 0, 0.2, fname);
  correlationPlot(bdtcut, "maxdoca", 0, 0.04, fname);
  correlationPlot(bdtcut, "chi2dof", 0, 5, fname);

  correlationPlot(bdtcut, "iso", 0.21, 1.01, fname);
  correlationPlot(bdtcut, "closetrk", 0., 21, fname);
  correlationPlot(bdtcut, "docatrk", 0., 0.03,  fname);

}


// ----------------------------------------------------------------------
void plotBDT::correlationPlot(double bdtcut, std::string var, double xmin, double xmax, std::string fname) {  

  // -- Events1
  TFile *f0 = TFile::Open(Form("weights/%s-Events0.root", fname.c_str()));

  TH1D *hs1 = new TH1D("hs1", var.c_str(), 100, xmin, xmax); 
  TH1D *hs2 = new TH1D("hs2", var.c_str(), 100, xmin, xmax); 
  TH1D *hsr0 = new TH1D("hsr0", Form("%s (BDT > %3.2f)/(all)", var.c_str(), bdtcut), 100, xmin, xmax); hsr0->Sumw2();

  TH1D *hb1 = new TH1D("hb1", var.c_str(), 100, xmin, xmax); 
  TH1D *hb2 = new TH1D("hb2", var.c_str(), 100, xmin, xmax); 
  TH1D *hbr0 = new TH1D("hbr0", Form("%s (BDT > %3.2f)/(all)", var.c_str(), bdtcut), 100, xmin, xmax); hbr0->Sumw2();

  TH2D *Hs = new TH2D("Hs", Form("SIGNAL BDT vs %s", var.c_str()), 100, xmin, xmax, 100, -1., 1.);
  TH2D *Hb = new TH2D("Hb", Form("Background BDT vs %s", var.c_str()), 100, xmin,  xmax, 100, -1., 1.);

  //  TTree *t = (TTree*)f->Get("TrainTree"); 
  TTree *t = (TTree*)f0->Get("TestTree"); 
  
  t->Draw(Form("%s>>hs1", var.c_str()), Form("BDT>%f&&classID==0", bdtcut));
  t->Draw(Form("%s>>hs2", var.c_str()), "classID==0");

  t->Draw(Form("%s>>hb1", var.c_str()), Form("BDT>%f&&classID==1", bdtcut));
  t->Draw(Form("%s>>hb2", var.c_str()), "classID==1");

  t->Draw(Form("BDT:%s>>Hs", var.c_str()), "classID==0");
  t->Draw(Form("BDT:%s>>Hb", var.c_str()), "classID==1");

  hsr0->Divide(hs1, hs2, 1., 1., "");
  hbr0->Divide(hb1, hb2, 1., 1., "");

  // -- Events1
  TFile *f1 = TFile::Open(Form("weights/%s-Events1.root", fname.c_str()));

  hs1 = new TH1D("hs1", var.c_str(), 100, xmin, xmax); 
  hs2 = new TH1D("hs2", var.c_str(), 100, xmin, xmax); 
  TH1D *hsr1 = new TH1D("hsr1", Form("%s (BDT > %3.2f)/(all)", var.c_str(), bdtcut), 100, xmin, xmax); hsr1->Sumw2();

  hb1 = new TH1D("hb1", var.c_str(), 100, xmin, xmax); 
  hb2 = new TH1D("hb2", var.c_str(), 100, xmin, xmax); 
  TH1D *hbr1 = new TH1D("hbr1", Form("%s (BDT > %3.2f)/(all)", var.c_str(), bdtcut), 100, xmin, xmax); hbr1->Sumw2();

  t = (TTree*)f1->Get("TestTree"); 
  
  t->Draw(Form("%s>>hs1", var.c_str()), Form("BDT>%f&&classID==0", bdtcut));
  t->Draw(Form("%s>>hs2", var.c_str()), "classID==0");

  t->Draw(Form("%s>>hb1", var.c_str()), Form("BDT>%f&&classID==1", bdtcut));
  t->Draw(Form("%s>>hb2", var.c_str()), "classID==1");

  hsr1->Divide(hs1, hs2, 1., 1., "");
  hbr1->Divide(hb1, hb2, 1., 1., "");

  // -- Events2
  TFile *f2 = TFile::Open(Form("weights/%s-Events2.root", fname.c_str()));

  hs1 = new TH1D("hs1", var.c_str(), 100, xmin, xmax); 
  hs2 = new TH1D("hs2", var.c_str(), 100, xmin, xmax); 
  TH1D *hsr2 = new TH1D("hsr2", Form("%s (BDT > %3.2f)/(all)", var.c_str(), bdtcut), 100, xmin, xmax); hsr2->Sumw2();

  hb1 = new TH1D("hb1", var.c_str(), 100, xmin, xmax); 
  hb2 = new TH1D("hb2", var.c_str(), 100, xmin, xmax); 
  TH1D *hbr2 = new TH1D("hbr2", Form("%s (BDT > %3.2f)/(all)", var.c_str(), bdtcut), 100, xmin, xmax); hbr2->Sumw2();

  t = (TTree*)f2->Get("TestTree"); 
  
  t->Draw(Form("%s>>hs1", var.c_str()), Form("BDT>%f&&classID==0", bdtcut));
  t->Draw(Form("%s>>hs2", var.c_str()), "classID==0");

  t->Draw(Form("%s>>hb1", var.c_str()), Form("BDT>%f&&classID==1", bdtcut));
  t->Draw(Form("%s>>hb2", var.c_str()), "classID==1");


  hsr2->Divide(hs1, hs2, 1., 1., "");
  hbr2->Divide(hb1, hb2, 1., 1., "");


  zone(2,2);
  
//   c0->cd(1);
//   hs1->Draw();
//   hs2->Scale(hs1->GetSumOfWeights()/hs2->GetSumOfWeights());
//   hs2->SetMarkerColor(kRed);  hs2->SetMarkerSize(1.);
//   hs2->Draw("esame");

//   tl->SetTextColor(kBlack);
//   tl->DrawLatex(0.5, 0.5, Form("signal && BDT>%3.1f", bdtcut)); 
//   tl->SetTextColor(kRed);
//   tl->DrawLatex(0.5, 0.45, Form("signal (same area)")); 

//   c0->cd(2);
//   hb1->Draw();
//   hb2->Scale(hb1->GetSumOfWeights()/hb2->GetSumOfWeights());
//   hb2->SetMarkerColor(kRed); hb2->SetMarkerSize(1.);
//   hb2->Draw("esame");

//   tl->SetTextColor(kBlack);
//   tl->DrawLatex(0.5, 0.5, Form("bg && BDT>%3.1f", bdtcut)); 
//   tl->SetTextColor(kRed);
//   tl->DrawLatex(0.5, 0.45, Form("bg (same area)")); 

  c0->cd(1);
  hsr0->Draw("hist");
  hsr1->SetLineColor(kRed);  hsr1->Draw("histsame");
  hsr2->SetLineColor(kBlue); hsr2->Draw("histsame");

  c0->cd(2);
  hbr0->Draw("hist");
  hbr1->SetLineColor(kRed);  hbr1->Draw("histsame");
  hbr2->SetLineColor(kBlue); hbr2->Draw("histsame");

  c0->cd(3); 
  gPad->SetLogz(0);
  Hs->DrawCopy("colz");

  c0->cd(4);
  gPad->SetLogz(0);
  Hb->Draw("colz");

  c0->SaveAs(Form("%s/%s-%s-correlations-TestTree.pdf", fDirectory.c_str(), fname.c_str(), var.c_str()));
}



// ----------------------------------------------------------------------
void plotBDT::illustrateAntiMuonSample(const char *cuts) {
  gStyle->SetOptStat(0); 
  TH1D *h1 = new TH1D("h1", "", 50, 4.9, 5.9); setFilledHist(h1, kBlue, kBlue, 3356); setTitles(h1, "m [GeV]", "candidates", 0.05, 1.1, 2.0);
  TH1D *h2 = new TH1D("h2", "", 50, 4.9, 5.9); setFilledHist(h2, kRed, kRed, 3365); 
  h1->SetMinimum(0); 

  // -- dimuons
  TTree *t = (TTree*)fF["SgData"]->Get("candAnaMuMu/events");
  t->Draw("m>>h1", Form("gmuid && %s", cuts)); 
  t->Draw("m>>h2", Form("!gmuid && %s", cuts)); 

  zone(1);
  shrinkPad(0.15, 0.22);
  h1->Draw();
  h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  h2->Draw("same"); 

  newLegend(0.55, 0.7, 0.85, 0.85, "#mu^{+} #mu{^-}"); 
  legg->AddEntry(h1, "passing #mu ID", "f"); 
  legg->AddEntry(h2, "failing #mu ID", "f"); 
  legg->Draw();

  c0->SaveAs(Form("%s/illustrateAMS-sg.pdf", fDirectory.c_str())); 

  // -- normalization
  h1->Reset();
  h2->Reset();
  t = (TTree*)fF["NoData"]->Get("candAnaBu2JpsiK/events");
  t->Draw("cm>>h1", Form("gmuid && %s", cuts)); 
  t->Draw("cm>>h2", Form("!gmuid && %s", cuts)); 

  zone(1);
  shrinkPad(0.15, 0.22);
  h1->Draw();
  h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  h2->Draw("same"); 
  
  newLegend(0.55, 0.7, 0.85, 0.85, "B^{+} #rightarrow J/#psi K^{+}"); 
  legg->AddEntry(h1, "passing #mu ID", "f"); 
  legg->AddEntry(h2, "failing #mu ID", "f"); 
  legg->Draw();

  c0->SaveAs(Form("%s/illustrateAMS-no.pdf", fDirectory.c_str())); 


  // -- control sample
  h1->Reset();
  h2->Reset();
  t = (TTree*)fF["CsData"]->Get("candAnaBs2JpsiPhi/events");
  t->Draw("cm>>h1", Form("gmuid && 0.995<mkk&&mkk<1.045 && %s", cuts)); 
  t->Draw("cm>>h2", Form("!gmuid && 0.995<mkk&&mkk<1.045 && %s", cuts)); 

  zone(1);
  shrinkPad(0.15, 0.22);
  h1->Draw();
  h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  h2->Draw("same"); 

  newLegend(0.55, 0.7, 0.85, 0.85, "B^{0}_{s} #rightarrow J/#psi #phi"); 
  legg->AddEntry(h1, "passing #mu ID", "f"); 
  legg->AddEntry(h2, "failing #mu ID", "f"); 
  legg->Draw();

  c0->SaveAs(Form("%s/illustrateAMS-cs.pdf", fDirectory.c_str())); 

}



// ----------------------------------------------------------------------
// -- stuff from tmvaglob.C
void plotBDT::GetMethodName( TString & name, TKey * mkey ) {
  if (mkey==0) return;
  name = mkey->GetName();
  name.ReplaceAll("Method_","");
}

void plotBDT::GetMethodTitle( TString & name, TKey * ikey ) {
  if (ikey==0) return;
  name = ikey->GetName();
}

void plotBDT::GetMethodTitle( TString & name, TDirectory * idir ) {
  if (idir==0) return;
  name = idir->GetName();
}


void plotBDT::SetSignalAndBackgroundStyle( TH1* sig, TH1* bkg, TH1* all )    {
      //signal
      // const Int_t FillColor__S = 38 + 150; // change of Color Scheme in ROOT-5.16.
      // convince yourself with gROOT->GetListOfColors()->Print()
  static Int_t c_SignalLine     = TColor::GetColor( "#0000ee" );
  static Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
  static Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
  static Int_t c_BackgroundFill = TColor::GetColor( "#ff0000" );

      Int_t FillColor__S = c_SignalFill;
      Int_t FillStyle__S = 3356; //1001;
      Int_t LineColor__S = c_SignalLine;
      Int_t LineWidth__S = 2;

      // background
      //Int_t icolor = UsePaperStyle ? 2 + 100 : 2;
      Int_t FillColor__B = c_BackgroundFill;
      Int_t FillStyle__B = 3365; //3554;
      Int_t LineColor__B = c_BackgroundLine;
      Int_t LineWidth__B = 2;

      Int_t FillColor__A = kBlack; 
      Int_t FillStyle__A = 3350; //3554;
      Int_t LineColor__A = kBlack;
      Int_t LineWidth__A = 2;

      if (sig != NULL) {
         sig->SetMarkerColor(LineColor__S );
         sig->SetLineColor( LineColor__S );
         sig->SetLineWidth( LineWidth__S );
         sig->SetFillStyle( FillStyle__S );
         sig->SetFillColor( FillColor__S );
      }
      
      if (bkg != NULL) {
	bkg->SetMarkerColor(LineColor__B );
	bkg->SetLineColor( LineColor__B );
	bkg->SetLineWidth( LineWidth__B );
	bkg->SetFillStyle( FillStyle__B );
	bkg->SetFillColor( FillColor__B );
      }
      
      if (all != NULL) {
	all->SetMarkerColor(LineColor__A );
	all->SetLineColor( LineColor__A );
	all->SetLineWidth( LineWidth__A );
	all->SetFillStyle( FillStyle__A );
	all->SetFillColor( FillColor__A );
      }
   }

void plotBDT::NormalizeHists( TH1* sig, TH1* bkg )   {
      if (sig->GetSumw2N() == 0) sig->Sumw2();
      if (bkg && bkg->GetSumw2N() == 0) bkg->Sumw2();
      
      if(sig->GetSumOfWeights()!=0) {
         Float_t dx = (sig->GetXaxis()->GetXmax() - sig->GetXaxis()->GetXmin())/sig->GetNbinsX();
         sig->Scale( 1.0/sig->GetSumOfWeights()/dx );
      }
      if (bkg != 0 && bkg->GetSumOfWeights()!=0) {
         Float_t dx = (bkg->GetXaxis()->GetXmax() - bkg->GetXaxis()->GetXmin())/bkg->GetNbinsX();
         bkg->Scale( 1.0/bkg->GetSumOfWeights()/dx );
      }
   }

void plotBDT::SetFrameStyle( TH1* frame, Float_t scale)   {
  frame->SetLabelOffset( 0.012, "X" );// label offset on x axis
  frame->SetLabelOffset( 0.012, "Y" );// label offset on x axis
  frame->GetXaxis()->SetTitleOffset( 1.25 );
  frame->GetYaxis()->SetTitleOffset( 1.22 );
  frame->GetXaxis()->SetTitleSize( 0.045*scale );
  frame->GetYaxis()->SetTitleSize( 0.045*scale );
  Float_t labelSize = 0.04*scale;
  frame->GetXaxis()->SetLabelSize( labelSize );
  frame->GetYaxis()->SetLabelSize( labelSize );
  
  // global style settings
  gPad->SetTicks();
  gPad->SetLeftMargin  ( 0.108*scale );
  gPad->SetRightMargin ( 0.050*scale );
  gPad->SetBottomMargin( 0.120*scale  );
}


int plotBDT::GetNumberOfTargets( TDirectory *dir )   {
      TIter next(dir->GetListOfKeys());
      TKey* key    = 0;
      Int_t noTrgts = 0;
         
      while ((key = (TKey*)next())) {
         if (key->GetCycle() != 1) continue;        
         if (TString(key->GetName()).Contains("__Regression_target")) noTrgts++;
      }
      return noTrgts;
}

int plotBDT::GetNumberOfInputVariables( TDirectory *dir )   {
  TIter next(dir->GetListOfKeys());
  TKey* key    = 0;
  Int_t noVars = 0;
  
  while ((key = (TKey*)next())) {
    if (key->GetCycle() != 1) continue;
    
    // count number of variables (signal is sufficient), exclude target(s)
    if (TString(key->GetName()).Contains("__Signal") || (TString(key->GetName()).Contains("__Regression") && !(TString(key->GetName()).Contains("__Regression_target")))) noVars++;
  }
  
  return noVars;
}



