// ----------------------------------------------------------------------
void plotAll() {

  plotCorr("pt", 15, 0., 15., "lxy/sxy ", 30, 0., 30.);
  c0.SaveAs("pt-lxy_sxy.eps");

  plotCorr("pt ", 15, 0., 15., "eta", 10, -2., 2.);
  c0.SaveAs("pt-eta.eps");

  plotCorr("pt", 15, 0., 15., "ptl0", 10, 0., 20.);
  c0.SaveAs("pt-ptl0.eps");

  plotCorr("pt", 15, 0., 15., "iso", 20, 0., 1.);
  c0.SaveAs("pt-ptl0.eps");

  plotCorr("pt", 15, 0., 15., "cosa", 20, 0.5, 1.);
  c0.SaveAs("pt-cosa.eps");
  
}


// ----------------------------------------------------------------------
void mass(const char *cut = "") {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TH1D *h1 = new TH1D("h1", "", 40, 5., 6.);
  TTree *tS = (TTree*)gFile->Get("events");
  
  tS->Draw("mass>>h1", cut);
  
  tl.SetTextSize(0.04);
  tl->DrawLatex(0.1, 0.92, cut);
  if (h1->GetEntries() > 10) {
    h1->Fit("pol1", "l", "", 5., 5.9);
    
    tl->SetNDC(kTRUE);
    tl->SetTextSize(0.06);
    tl->DrawLatex(0.2, 0.85,
		  Form("p0: %5.3f+/-%5.3f", 
		       h1->GetFunction("pol1")->GetParameter(0), 
		       h1->GetFunction("pol1")->GetParError(0)));
    tl->DrawLatex(0.2, 0.80,
		  Form("p1: %5.3f+/-%5.3f", 
		       h1->GetFunction("pol1")->GetParameter(1), 
		       h1->GetFunction("pol1")->GetParError(1)));
    tl.SetTextSize(0.04);
    tl->DrawLatex(0.1, 0.92, cut);
  }
}


void allProcesses() {
  TFile::Open("root/csg.default-001.root");
  processes("ptl1>4 &&chi2<200 && cosa>0.95 &&lxy>0.01 &&rmm<12.", "pt", "pT [GeV]", 20, 0., 60.);
  c0.SaveAs("root/sg-pt.eps");

  TFile::Open("root/cbg.default-001.root");
  processes("ptl1>4 &&chi2<200 && cosa>0.95 &&lxy>0.01 &&rmm<12.", "pt", "pT [GeV]", 20, 0., 60.);
  c0.SaveAs("root/bg-pt.eps");


  TFile::Open("root/csg.default-001.root");
  processes("pt>10 &&ptl1>4 &&chi2<200 && cosa>0.95 &&lxy>0.01 &&rmm<12.", "iso", "iso", 33, 0., 1.1);
  c0.SaveAs("root/sg-iso.eps");

  TFile::Open("root/cbg.default-001.root");
  processes("pt>10 &&ptl1>4 &&chi2<200 && cosa>0.95 &&lxy>0.01 &&rmm<12.", "iso", "iso", 33, 0., 1.1);
  c0.SaveAs("root/bg-iso.eps");

  TFile::Open("root/csg.default-001.root");
  processes("pt>10 &&ptl1>4 &&chi2<200 && cosa>0.95 &&lxy>0.01 &&rmm<12.", "isoveto", "isoveto", 10, 0., 10.);
  c0.SaveAs("root/sg-isoveto.eps");

  TFile::Open("root/cbg.default-001.root");
  processes("pt>10 &&ptl1>4 &&chi2<200 && cosa>0.95 &&lxy>0.01 &&rmm<12.", "isoveto", "isoveto", 10, 0., 10.);
  c0.SaveAs("root/bg-isoveto.eps");

}


// ----------------------------------------------------------------------
void processes(const char *cut = "", const char *var, const char *axis, int bin, double lo, double hi) {

  TTree *tS = (TTree*)gFile->Get("events");
  TH1D *hGSP = new TH1D("hGSP", "", bin, lo, hi); hGSP->Sumw2();
  TH1D *hFEX = new TH1D("hFEX", "", bin, lo, hi); hFEX->Sumw2();
  TH1D *hGGF = new TH1D("hGGF", "", bin, lo, hi); hGGF->Sumw2();
  TH1D *hSum = new TH1D("hSum", "", bin, lo, hi); hSum->Sumw2();

  tS->Draw(Form("%s >> hGSP", var), Form("%s && process == 42", cut)); 
  tS->Draw(Form("%s >> hFEX", var), Form("%s && process == 41", cut)); 
  tS->Draw(Form("%s >> hGGF", var), Form("%s && process == 40", cut)); 

  hSum->Add(hGSP);
  hSum->Add(hFEX);
  hSum->Add(hGGF);

  hGSP->Divide(hGSP, hSum, 1., 1., "b");
  hFEX->Divide(hFEX, hSum, 1., 1., "b");
  hGGF->Divide(hGGF, hSum, 1., 1., "b");
  
  hGSP->SetMinimum(0.);
  hGSP->SetMaximum(1.);

  gStyle->SetOptStat(0);
  shrinkPad(0.15, 0.15);
  setTitles(hGSP, axis, "ratio", 0.05, 1.1, 1.3);

  setHist(hGSP, kGreen, 20, 2);
  setHist(hFEX, kRed, 20, 2);
  setHist(hGGF, kBlue, 20, 2);

  hpl(hGSP, "histgreen");
  hpl(hFEX, "histsamered");
  hpl(hGGF, "histsameblue");
}


// ----------------------------------------------------------------------
void plotVar(const char *s1, const char *cut = "", int nbin, double lo, double hi) {
  
  TTree *tS = (TTree*)bb.fS[0]->Get("events");
  TTree *tB = (TTree*)bb.fM[0]->Get("events");

  TH1D *hS = new TH1D("hS", "", nbin, lo, hi);
  TH1D *hB = new TH1D("hB", "", nbin, lo, hi);
 
  tS->Draw(Form("%s>>hS", s1), cut);
  tB->Draw(Form("%s>>hB", s1), cut);

  double max = 1.1*(hS->GetMaximum() > hB->GetMaximum()? hS->GetMaximum(): hB->GetMaximum());
  hS->SetMaximum(max);

  tl->SetNDC(kTRUE);
  
  hpl(hS, "blue"); 
  hpl(hB, "redsame");
  tl->DrawLatex(0.2, 0.92, Form("BG mean: %5.3f", hB->GetMean()));

}



// ----------------------------------------------------------------------
void plotAllCuts() {
  plotCut("isoveto", 10, 0., 10., "pt>10 && abs(eta) < 2.5 && chi2<2 &&lxy>0.015", 0); // ~OK
  plotCut("cosa", 10, 0.99, 1., "pt>10 && abs(eta) < 2.5 && chi2<2 && lxy>0.015", 1);  // BAD

  plotCut("i12", 25, 0.5, 1., "pt>12 && abs(eta) < 2.5 && chi2<3", 1);    // OK: ~10%
  plotCut("lxy", 20, 0., 0.2, "pt>12 && abs(eta) < 2.5 && chi2<3", 1);    
  plotCut("lxy/sxy", 40, 0., 40., "pt>12 && abs(eta) < 2.5 && chi2<3", 1);

  plotCut("txy", 40, 0., 0.1, "pt>12 && abs(eta) < 2.5 && chi2<3", 1);  // OK for SG, bad for BG


  plotCut("i12", 40, 0.5, 1, "pt>5 && abs(eta) < 2.5 && ptl1 > 3", 1);


  plotCut("cosa", 40, 0.99, 1, "pt>5 && abs(eta) < 2.5 && ptl1 > 3 && i12>0.8", 1);
  plotCut("i12", 40, 0.5, 1, "pt>5 && abs(eta) < 2.5 && ptl1 > 3 && cosa>0.99", 1);

}

// ----------------------------------------------------------------------
void plotCut(const char *s1, int b1, double lo1, double hi1, const char *otherCuts, int lower = 1) {

  TFile *fS = new TFile("root/csg-001.default.root");
  TTree *tS = (TTree*)fS->Get("events");
  TFile *fB = new TFile("root/cbg-001.default.root");
  TTree *tB = (TTree*)fB->Get("events");

  TH1D *hbin = new TH1D("hbin", "", b1, lo1, hi1);

  TH1D *hs1 = new TH1D("hs1", "", b1, lo1, hi1);
  TH1D *hs2 = new TH1D("hs2", "", b1, lo1, hi1);

  TH1D *hb1 = new TH1D("hb1", "", b1, lo1, hi1);
  TH1D *hb2 = new TH1D("hb2", "", b1, lo1, hi1);

  double totalS = tS->GetEntries();
  double totalB = tB->GetEntries();
  
  double cnt1(0.), cnt2(0.), cnt12(0.);
  int i1;
  for (double ic1 = lo1; ic1 < hi1; ic1 += (hi1-lo1)/b1) {
    i1 = hbin->FindBin(ic1);

    if (lower) {
      cnt1  = tS->Draw(s1, Form("%s >= %f", s1, ic1));
      cnt2  = tS->Draw(s1, Form("%s", otherCuts, ic1));
      cnt12 = tS->Draw(s1, Form("%s >= %f && %s ", s1, ic1, otherCuts));
    } else {
      cnt1  = tS->Draw(s1, Form("%s <= %f", s1, ic1));
      cnt2  = tS->Draw(s1, Form("%s", otherCuts, ic1));
      cnt12 = tS->Draw(s1, Form("%s <= %f && %s ", s1, ic1, otherCuts));
    }

    hs1->SetBinContent(i1, (cnt1/totalS)*(cnt2/totalS));
    hs2->SetBinContent(i1, (cnt12/totalS));

    if (lower) {
      cnt1  = tB->Draw(s1, Form("%s >= %f", s1, ic1));
      cnt2  = tB->Draw(s1, Form("%s", otherCuts, ic1));
      cnt12 = tB->Draw(s1, Form("%s >= %f && %s ", s1, ic1, otherCuts));
    } else {
      cnt1  = tB->Draw(s1, Form("%s <= %f", s1, ic1));
      cnt2  = tB->Draw(s1, Form("%s", otherCuts, ic1));
      cnt12 = tB->Draw(s1, Form("%s <= %f && %s ", s1, ic1, otherCuts));
    }

    hb1->SetBinContent(i1, (cnt1/totalB)*(cnt2/totalB));
    hb2->SetBinContent(i1, (cnt12/totalB));

  }  

  gStyle->SetOptStat(0);

  hs1->SetMinimum(0.);
  hb1->SetMinimum(0.);
  
  c0.Clear();
  c0.Divide(1,2);

  c0.cd(1);
  hpl(hs1, "red");
  hpl(hs2, "sameblue");

  c0.cd(2);
  hpl(hb1, "red");
  hpl(hb2, "sameblue");

}



// ----------------------------------------------------------------------
void plotCorr(const char *s1, int b1, double lo1, double hi1, 
	      const char *s2, int b2, double lo2, double hi2, int signal = 0) {
  
  TFile *fS;
  if (1 == signal) {
    fS = new TFile("root/csg-001.root");
  } else {
    fS = TFile::Open("root/cbg-001.root");
  }

  TH1D *hb1 = new TH1D("hb1", "", b1, lo1, hi1);
  TH1D *hb2 = new TH1D("hb2", "", b2, lo2, hi2);

  TH2D *h1 = new TH2D("h1", "", b1, lo1, hi1, b2, lo2, hi2);
  TH2D *h2 = new TH2D("h2", "", b1, lo1, hi1, b2, lo2, hi2);
  TH2D *h3 = new TH2D("h3", "", b1, lo1, hi1, b2, lo2, hi2);


  TTree *tS = (TTree*)fS->Get("events");
  //   TTree *tB = (TTree*)fB->Get("events");

  double total = tS->GetEntries();
  double cnt1(0.), cnt2(0.), cnt12(0.);
  int i1, i2; 

  for (double ic1 = lo1; ic1 < hi1; ic1 += (hi1-lo1)/b1) {
    for (double ic2 = lo2; ic2 < hi2; ic2 += (hi2-lo2)/b2) {
      cnt1  = tS->Draw(s1, Form("%s > %f", s1, ic1));
      cnt2  = tS->Draw(s2, Form("%s > %f", s2, ic2));
      cnt12 = tS->Draw(s2, Form("%s > %f && %s > %f", s1, ic1, s2, ic2));

      i1 = hb1->FindBin(ic1);
      i2 = hb2->FindBin(ic2);

      h1->SetBinContent(i1, i2, (cnt1/total)*(cnt2/total));
      h2->SetBinContent(i1, i2, cnt12/total);

    cout << "ic1: " << ic1 << " ic2: " << ic2
	 << " i1: " << i1 << " i2: " << i2
	 << " cnt1: " << cnt1  
	 << " cnt2: " << cnt2 
	 << " cnt12: " << cnt12 
	 << endl;

    }

  }

  gStyle->SetOptStat(0);
  
  c0.Clear();
  c0.Divide(2,2);
  
  c0.cd(1);
  h1->Draw("colz");

  c0.cd(2);
  h2->Draw("colz");

  c0.cd(3);
  h3->Divide(h1, h2);
  h3->SetMinimum(0.3);
  h3->SetMaximum(1.1);
  h3->Draw("colz");
  
  c0.cd(4);
  plotVar(Form("%s:%s", s2, s1), Form("%f < %s && %s < %f && %f < %s && %s < %f", lo2, s2, s2, hi2, lo1, s1, s1, hi1), "colz");

}
