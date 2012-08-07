#include "plotBDT.hh"

#include <TMath.h>
#include <TGaxis.h>

using namespace std; 

ClassImp(plotBDT)

// ----------------------------------------------------------------------
plotBDT::plotBDT(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 
  fNumbersFileName = fDirectory + "/anaBmm.plotBDT." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  string hfname  = fDirectory + "/anaBmm.plotBDT." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  //  if (fHistFile) fHistFile->Close();
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

  fIsMC = false;

  fDoUseBDT = true; 

  int NBINS = (fMassHi - fMassLo)/0.025;
  TH1D *h(0); 
  TProfile *p(0); 
  for (unsigned int i = 0; i < fNchan; ++i) {
    h = new TH1D(Form("hMassNoCuts%d", i), Form("hMassNoCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassNoCuts.push_back(h); 

    h = new TH1D(Form("hMass%d", i), Form("hMass%d", i), NBINS, fMassLo, fMassHi);
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

    p = new TProfile(Form("pNpvBDT%d", i), Form("pNpvBDT%d", i), 30, 0., 30.);
    fpNpvBDT.push_back(p); 
    p = new TProfile(Form("pNpvAcBDT%d", i), Form("pNpvAcBDT%d", i), 30, 0., 30.);
    fpNpvAcBDT.push_back(p); 
  }

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

  validateAllOddEven();

  if (channels & 1) {
    tmvaControlPlots();
  }

  if (channels & 2) {
    bdtScan();
    plotSSB();
  }

  if (channels & 4) {
    bdtDependencies("SgData");
    illustrateAntiMuonSample();
  }


}



// ----------------------------------------------------------------------
void plotBDT::loopFunction(int mode) {
  if (1 == mode) loopFunction1();
  if (2 == mode) loopFunction2();
}


// ----------------------------------------------------------------------
void plotBDT::loopFunction2() {

}

// ----------------------------------------------------------------------
void plotBDT::loopFunction1() {

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
  if (fb.gmuid && fb.hlt) {
    fhBDT[fChan]->Fill(fBDT);
    if (fb.m > 5.45) fhHiBDT[fChan]->Fill(fBDT);
    if (fb.m < 5.20) fhLoBDT[fChan]->Fill(fBDT);
    if (fb.m > 5.20 && fb.m < 5.45) fhInBDT[fChan]->Fill(fBDT); 

    if (fBDT > -9.) fpMassBDT[fChan]->Fill(fb.m, fBDT);
    if (fBDT > 0.)  fpMassAcBDT[fChan]->Fill(fb.m, fBDT);
    
    if (fBDT > -9.) fpNpvBDT[fChan]->Fill(fb.pvn, fBDT);
    if (fBDT > 0.)  fpNpvAcBDT[fChan]->Fill(fb.pvn, fBDT);
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
    fpNpvBDT[i]->Reset();
    fpNpvAcBDT[i]->Reset();
  }    
}


// ----------------------------------------------------------------------
void plotBDT::testLoop(string mode) { 

}

// ----------------------------------------------------------------------
void plotBDT::overlap() {
  string mode = "SgMc";
  TTree *t = getTree(mode); 
  resetHistograms();
  setupTree(t, mode);
  loopOverTree(t, mode, 2); 


}


// ----------------------------------------------------------------------
void plotBDT::bdtDependencies(string mode) {

  TH1D *h(0), *h1(0);
  TH2D *h2(0); 
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
  for (int i = 0; i < fNchan; ++i) {
    for (int j = 1; j < fpMassAcBDT[i]->GetNbinsX(); ++j) {
      if (fpMassAcBDT[i]->GetBinEntries(j) < 10) fpMassAcBDT[i]->SetBinEntries(j, 0); 
    }

    for (int j = 1; j < fpMassBDT[i]->GetNbinsX(); ++j) {
      if (fpMassBDT[i]->GetBinEntries(j) < 10) fpMassBDT[i]->SetBinEntries(j, 0); 
    }

    for (int j = 1; j < fpNpvBDT[i]->GetNbinsX(); ++j) {
      if (fpNpvBDT[i]->GetBinEntries(j) < 10) fpNpvBDT[i]->SetBinEntries(j, 0); 
    }

    for (int j = 1; j < fpNpvAcBDT[i]->GetNbinsX(); ++j) {
      if (fpNpvAcBDT[i]->GetBinEntries(j) < 10) fpNpvAcBDT[i]->SetBinEntries(j, 0); 
    }
  }
  
  for (int i = 0; i < fNchan; ++i) {
    fpMassBDT[i]->SetAxisRange(4.9, 5.9, "X");
    fpMassBDT[i]->Fit("pol0");  
    c0->SaveAs(Form("%s/dep-bdt-mass-nocuts%d.pdf", fDirectory.c_str(), i));

    fpMassAcBDT[i]->SetAxisRange(4.9, 5.9, "X");
    fpMassAcBDT[i]->Fit("pol0");  
    c0->SaveAs(Form("%s/dep-bdt-mass-aftercuts%d.pdf", fDirectory.c_str(), i));

    fpNpvBDT[i]->Fit("pol0");  
    c0->SaveAs(Form("%s/dep-bdt-npv%d.pdf", fDirectory.c_str(), i));
  }

}


// ----------------------------------------------------------------------
void plotBDT::bdtScan() { 

  TH1D *h(0), *h1(0);
  TH2D *h2(0); 
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
  double bdtcut(0.), pass(0.), all(0.), eff(0.);
  for (int i = 0; i < fNchan; ++i) {
    cout << "channel " << i << endl;
    h = new TH1D(Form("hSgEff%d", i), Form("hSgEff%d", i), 200, -1., 1.); 
    h->SetDirectory(fHistFile); 

    all  = fhMassNoCuts[i]->GetSumOfWeights();
    for (int j = 1; j < 201; ++j) {
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
  for (int i = 0; i < fNchan; ++i) {
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

  c0->SaveAs(Form("%s/bdt-lin-Bg-HiLo-Efficiency0.pdf", fDirectory.c_str())); 
  gPad->SetLogy(1); 
  H1->Draw();
  H2->Draw("same");
  c0->SaveAs(Form("%s/bdt-log-Bg-HiLo-Efficiency0.pdf", fDirectory.c_str())); 

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

  c0->SaveAs(Form("%s/bdt-lin-Bg-HiLo-Efficiency1.pdf", fDirectory.c_str())); 
  gPad->SetLogy(1);
  H1->Draw();
  H2->Draw("same");
  c0->SaveAs(Form("%s/bdt-log-Bg-HiLo-Efficiency1.pdf", fDirectory.c_str())); 

  plotEffVsBg(0);   
  plotEffVsBg(1);   

}


// ----------------------------------------------------------------------
void plotBDT::tmvaControlPlots() {

  string type[2] = {"even", "odd"};
  string XmlName;
  for (int j = 0; j < 2; ++j) {

    for (int i = 0; i < fNchan; ++i) {
      
      XmlName = "weights/" + fCuts[i]->xmlFile + Form("-%s_BDT.weights.xml", type[j].c_str()); 
      string rootfile = XmlName; 
      replaceAll(rootfile, "_BDT.weights.xml", ".root"); 
      fRootFile = TFile::Open(rootfile.c_str());
      cout << "fRootFile: " << rootfile << endl;

      fBdtString = fRootFile->GetName(); 
      fBdtString = fBdtString.substr(0, fBdtString.find(".root"));
      fBdtString = fBdtString.substr(fBdtString.find("/")+1);
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
	cout << "fBdtString: " << fBdtString << endl;

	variableRanking();
	fRootFile->Close();
      }

    }
  }

}


// ----------------------------------------------------------------------
void plotBDT::plotSSB() {
  for (int i = 0; i < fNchan; ++i) {
    string rootfile = "weights/" + fCuts[i]->xmlFile + "-combined.root";
    fRootFile = TFile::Open(rootfile.c_str());
    cout << "fRootFile: " << rootfile << endl;
    
    fBdtString = fRootFile->GetName(); 
    fBdtString = fBdtString.substr(0, fBdtString.find(".root"));
    fBdtString = fBdtString.substr(fBdtString.find("/")+1);
    cout << "fBdtString: " << fBdtString << endl;
    
    ssb();
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
  }

}


// ----------------------------------------------------------------------
void plotBDT::validateAllOddEven() {

  validateOddEven("weights/TMVA-0-odd.root", "weights/TMVA-0-even.root", "train", 0);
  validateOddEven("weights/TMVA-0-odd.root", "weights/TMVA-0-even.root", "train", 1);

  validateOddEven("weights/TMVA-1-odd.root", "weights/TMVA-1-even.root", "train", 0);
  validateOddEven("weights/TMVA-1-odd.root", "weights/TMVA-1-even.root", "train", 1);

}



// ----------------------------------------------------------------------
void plotBDT::validateOddEven(const char *fnOdd, const char *fnEven, const char *type, int classID) {

  TFile *fOdd = TFile::Open(fnOdd); 
  TFile *fEven = TFile::Open(fnEven); 

  string sOdd = fnOdd;
  string::size_type m1 = sOdd.find("weights/"); 
  sOdd = sOdd.substr(m1+8, sOdd.length()-m1+1); 
  m1 = sOdd.find(".root"); 
  sOdd = sOdd.substr(0, m1); 

  string sEven = fnEven;
  m1 = sEven.find("weights/"); 
  sEven = sEven.substr(m1+8, sEven.length()-m1+1); 
  m1 = sEven.find(".root"); 
  sEven = sEven.substr(0, m1); 

  vector<string> varNames;
  vector<double> varMin, varMax; 
  vector<int> varNbins; 
  vector<int> varLog; 
  varNames.push_back("m1pt"); varMin.push_back(0.); varMax.push_back(40.); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("m2pt"); varMin.push_back(0.); varMax.push_back(20.); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("m1eta"); varMin.push_back(-2.5); varMax.push_back(2.5); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("m2eta"); varMin.push_back(-2.5); varMax.push_back(2.5); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("pt");    varMin.push_back(0.); varMax.push_back(40.); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("eta");   varMin.push_back(-2.5); varMax.push_back(2.5); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("fls3d"); varMin.push_back(0.); varMax.push_back(120.); varNbins.push_back(120); varLog.push_back(1); 
  varNames.push_back("alpha"); varMin.push_back(0.); varMax.push_back(2.0); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("maxdoca"); varMin.push_back(0.); varMax.push_back(0.03); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("pvip"); varMin.push_back(0.); varMax.push_back(0.05); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("pvips"); varMin.push_back(0.); varMax.push_back(5); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("iso"); varMin.push_back(0.7); varMax.push_back(1.0); varNbins.push_back(300); varLog.push_back(0); 
  varNames.push_back("closetrk"); varMin.push_back(0.); varMax.push_back(21); varNbins.push_back(21); varLog.push_back(0); 
  varNames.push_back("docatrk"); varMin.push_back(0.); varMax.push_back(0.02); varNbins.push_back(200); varLog.push_back(0); 
  varNames.push_back("chi2dof"); varMin.push_back(0.); varMax.push_back(5); varNbins.push_back(200); varLog.push_back(0); 

  
  TTree *t0Odd(0), *t0Even(0); 
  if (!strcmp("test", type)) t0Odd = (TTree*)fOdd->Get("TestTree");
  if (!strcmp("train", type)) t0Odd = (TTree*)fOdd->Get("TrainTree");
  if (!strcmp("test", type)) t0Even = (TTree*)fEven->Get("TestTree");
  if (!strcmp("train", type)) t0Even = (TTree*)fEven->Get("TrainTree");
  TH1D *h1, *h2; 
  string h1Name, h2Name; 
  
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  for (int i = 0; i < varNames.size(); ++i) {
    h1Name = Form("h0_%s", varNames[i].c_str());
    h1  = new TH1D(h1Name.c_str(), varNames[i].c_str(), varNbins[i], varMin[i], varMax[i]);
    //    setFilledHist(h1, kBlue, kBlue, 3004); 
    h2Name = Form("h1_%s", varNames[i].c_str());
    h2 = new TH1D(h2Name.c_str(), h2Name.c_str(), varNbins[i], varMin[i], varMax[i]);
    //    setFilledHist(h2, kRed, kRed, 3004); 
    SetSignalAndBackgroundStyle(h1, h2);            

    t0Odd->Draw(Form("%s>>%s", varNames[i].c_str(), h1Name.c_str()), Form("classID==%d", classID)); 
    t0Even->Draw(Form("%s>>%s", varNames[i].c_str(), h2Name.c_str()), Form("classID==%d", classID)); 
    
    h1->Draw();
    h2->Draw("same");
    double ks = h1->KolmogorovTest(h2);
    tl->DrawLatex(0.2, 0.92, Form("%s: P(KS)= %4.3f", varNames[i].c_str(), ks)); 

    if (1 == varLog[i]) {
      gPad->SetLogy(1); 
    } else {
      gPad->SetLogy(0); 
    }

    c0->SaveAs(Form("%s/ks-%s_%s_%s_%s-%s.pdf", 
		    fDirectory.c_str(), (0 == classID?"sg":"bg"), type, sOdd.c_str(), sEven.c_str(), varNames[i].c_str())); 
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
  Int_t xPad = 1; // no of plots in x
  Int_t yPad = 1; // no of plots in y
  Int_t noPad = xPad * yPad ; 
  Int_t width = 600;   // size of canvas
  
  // this defines how many canvases we need
  TCanvas *c = new TCanvas( Form("canvas%d", 1), "canvas1",  200, 20, width, (Int_t)width*0.78 ); 
  
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
  Bool_t   createNewFig = kFALSE;
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
    
    // find the corredponding backgrouns histo
    TString bgname = hname;
    bgname.ReplaceAll("__Signal","__Background");
    TH1 *bgd = (TH1*)dir->Get(bgname);
    if (bgd == NULL) {
      cout << "ERROR!!! couldn't find background histo for" << hname << endl;
      exit;
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
void plotBDT::variableRanking() {
  string which = "IdTransformation";
  TH1D *h(0); 
  vector<string> splits; 
  splits.push_back("odd"); 
  splits.push_back("even"); 

  for (unsigned int i = 0; i < splits.size(); ++i) {
    which = "IdTransformation";
    h = (TH1D*)fRootFile->Get(Form("rank_%s_%s", splits[i].c_str(), which.c_str()));
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
  double bdtCut, maxSSB(-1.), maxBDT(-1.); 
  int nEvent(0); 
  for (int ibin = 0; ibin < bdtBins; ++ibin) {
    bdtCut = bdtMin + ibin*(bdtMax-bdtMin)/bdtBins;
    //    cout << " bin " << ibin << " cutting at bdt > " << bdtCut << endl;
    nEvent = t->GetEntries();
    sm->Reset();
    dm->Reset();
    for (Long64_t ievt=0; ievt<nEvent; ievt++) {
      t->GetEntry(ievt);
      if (false == hlt) continue;
      if (false == gmuid) continue;
      if (bdt < bdtCut) continue;

      if (0 == classID) {
	sm->Fill(m, w8); 
      } else {
	dm->Fill(m, w8); 
      }
      
    }
    // FIXME: wo wird hier normiert?? Durch w8, aus dem tree gelesen!
    double s = sm->Integral(sm->FindBin(5.3), sm->FindBin(5.45));
    double d = dm->Integral(0, h->GetNbinsX());
    double b = d*(5.45-5.30)/(5.9-4.9+0.25);
    if (s+b >0) {
      double ssb = s/TMath::Sqrt(s+b);
      if (ssb > maxSSB) {
	maxSSB = ssb; 
	maxBDT = bdtCut;
      }
      h->SetBinContent(ibin, ssb); 
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
  fTEX << formatTex(maxBDT, Form("%s:%s:maxSSB:bdt",  fSuffix.c_str(), fBdtString.c_str()), 2) << endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  h->Draw();
  tl->DrawLatex(0.25, 0.75, Form("S_{max} = %4.3f", maxSSB));
  tl->DrawLatex(0.25, 0.68, Form("B_{max} = %4.3f", maxBDT));
  c0->SaveAs(Form("%s/%s-ssb.pdf", fDirectory.c_str(), fBdtString.c_str())); 

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

  
  c0->SaveAs(Form("%s/bdtSgEffVsBgEvts%d.pdf", fDirectory.c_str(), offset)); 

  // -- dump table of numbers
  h0 = (TH1D*)(gDirectory->Get(Form("hBgEvts%d", offset)))->Clone("H0"); 
  h1 = (TH1D*)(gDirectory->Get(Form("hSgEff%d", offset)))->Clone("H1"); 
  int icnt(0); 
  double eff(h1->GetBinContent(1)), step(0.0025); 
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
  static Int_t c_Canvas         = kWhite; //TColor::GetColor( "#f0f0f0" );
  static Int_t c_FrameFill      = kWhite; //TColor::GetColor( "#fffffd" );
  static Int_t c_TitleBox       = kWhite; //TColor::GetColor( "#5D6B7D" );
  static Int_t c_TitleBorder     = kWhite; //TColor::GetColor( "#7D8B9D" );
  static Int_t c_TitleText      = TColor::GetColor( "#FFFFFF" );
  static Int_t c_SignalLine     = TColor::GetColor( "#0000ee" );
  static Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
  static Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
  static Int_t c_BackgroundFill = TColor::GetColor( "#ff0000" );
  static Int_t c_NovelBlue      = TColor::GetColor( "#2244a5" );

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

      if (sig != NULL) {
         sig->SetLineColor( LineColor__S );
         sig->SetLineWidth( LineWidth__S );
         sig->SetFillStyle( FillStyle__S );
         sig->SetFillColor( FillColor__S );
      }
 
      if (bkg != NULL) {
         bkg->SetLineColor( LineColor__B );
         bkg->SetLineWidth( LineWidth__B );
         bkg->SetFillStyle( FillStyle__B );
         bkg->SetFillColor( FillColor__B );
      }

      if (all != NULL) {
         all->SetLineColor( LineColor__S );
         all->SetLineWidth( LineWidth__S );
         all->SetFillStyle( FillStyle__S );
         all->SetFillColor( FillColor__S );
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


