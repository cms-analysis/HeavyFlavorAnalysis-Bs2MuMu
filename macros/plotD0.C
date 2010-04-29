// -- Fit Functions
double f_Gauss(double *x, double *par) {
  // par[0] -> area
  // par[1] -> mean
  // par[2] -> sigma

  double sqrt2pi = 2.506628275;

  if (par[2] > 0.) {
    Double_t arg = (x[0] - par[1]) / par[2];
    Double_t fitval =  (par[0]/(sqrt2pi*par[2])) * TMath::Exp(-0.5*arg*arg);
    return fitval;
  }
  else {
    return -1.;
  }
}



// ----------------------------------------------------------------------
// pol1 and Gauss
double f_p1ag(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] = par 0 of pol1
  // par[4] = par 1 of pol1

  return (par[3] + par[4]*x[0] + f_Gauss(x, &par[0]));
}


// ----------------------------------------------------------------------
// pol1 
double f_p1(double *x, double *par) {
  // par[0] = par 0 of pol1
  // par[1] = par 1 of pol1

  return (par[0] + par[1]*x[0]);
}


TF1 *f0 = new TF1("f0", f_Gauss,  1.7,  2.0, 3);
TF1 *f1 = new TF1("f1", f_p1,     1.7,  2.0, 2);
TF1 *f2 = new TF1("f2", f_p1ag,   1.7,  2.0, 5);

// ----------------------------------------------------------------------
void show(const char *dfile = "/shome/ursl/root/combined.root", 
	  const char *mfile = "/shome/ursl/root/100414.v05.mc.bmmReader.default.root", 
	  const char *hname = "h14") {


  TFile *fm = TFile::Open(mfile); 
  TFile *fd = TFile::Open(dfile); 

  tdrStyle->cd();

  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 
  gStyle->SetOptTitle(0); 

  fm->cd(); 
  TH1D *hm = (TH1D*)gDirectory->Get(hname); 
  hm->SetFillStyle(1000); hm->SetFillColor(kYellow);

  fd->cd(); 
  TH1D *hd = (TH1D*)gDirectory->Get(hname); 

  // -- normalization and maximum
  double hmax = hm->GetMaximum()> hd->GetMaximum()?hm->GetMaximum():hd->GetMaximum();
  double nd = hd->Integral(hd->FindBin(1.7), hd->FindBin(2.0)); 
  double nm = hm->Integral(hm->FindBin(1.7), hm->FindBin(2.0)); 
  cout << "nd: " << nd << " nm: " << nm << endl;

  hm->Scale(nd/nm); 
  hm->SetMaximum(1.1*hmax); 

  TLatex  *tl; tl = new TLatex(); tl->SetNDC(kTRUE); tl->SetTextSize(0.04); 
  TLegend *legg; 
  TLegendEntry *legge; 
  legg = new TLegend(0.65, 0.7, 0.90, 0.85);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04);  legg->SetFillColor(0); legg->SetTextFont(42); 
  legge = legg->AddEntry(hd, "Data", "p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hm, "MC", "f"); legge->SetTextColor(kBlack);

  double xcms(0.22), ycms(0.78); 

  gPad->SetLeftMargin(0.2);
  showHist(hm, ""); 
  showHist(hd, "esame"); 
  legg->Draw();
  tl->SetTextFont(62);  tl->SetTextSize(0.04);  tl->DrawLatex(xcms, ycms, "CMS preliminary"); 
  c0.SaveAs("d0-damc-nofit.eps"); 

  fitHist(hd, "esame"); 
  legg->Draw();
  tl->SetTextFont(62);  tl->SetTextSize(0.04);  tl->DrawLatex(xcms, ycms, "CMS preliminary"); 
  c0.SaveAs("d0-damc-fit.eps"); 

  fitHist(hm, ""); 
  f0->SetParameters(f2->GetParameter(0), f2->GetParameter(1),  f2->GetParameter(2));
  double yield = f0->Integral(f2->GetParameter(1)-2.*f2->GetParameter(2), f2->GetParameter(1)+2.*f2->GetParameter(2))/hm->GetBinWidth(2); 
  double yieldE = yield*f2->GetParError(0)/f2->GetParameter(0); 
  f1->SetParameters(f2->GetParameter(3), f2->GetParameter(4));
  double bg =  f1->Integral(f2->GetParameter(1)-2.*f2->GetParameter(2), f2->GetParameter(1)+2.*f2->GetParameter(2))/hm->GetBinWidth(2); 
  cout << "MC yield: " << yield << "  " << bg << endl;

  tl->SetTextFont(62);   tl->SetTextSize(0.04); tl->DrawLatex(xcms, ycms, "MC"); 
  tl->SetTextFont(42);  tl->SetTextSize(0.04); 
  tl->DrawLatex(xcms, 0.32, Form("#mu: %5.3f #pm %5.3f GeV", f2->GetParameter(1), f2->GetParError(1))); 
  tl->DrawLatex(xcms, 0.28, Form("#sigma: %5.3f #pm %5.3f GeV", f2->GetParameter(2), f2->GetParError(2))); 
  tl->DrawLatex(xcms, 0.24, Form("N: %4.0f #pm %3.0f", yield, yieldE)); 
  tl->DrawLatex(xcms, 0.20, Form("S/ #sqrt{S+B}: %4.1f", yield/TMath::Sqrt(yield+bg))); 
  c0.SaveAs("d0-mc-fit.eps"); 

  fitHist(hd, "e"); 
  f0->SetParameters(f2->GetParameter(0), f2->GetParameter(1),  f2->GetParameter(2));
  yield = f0->Integral(f2->GetParameter(1)-2.*f2->GetParameter(2), f2->GetParameter(1)+2.*f2->GetParameter(2))/hd->GetBinWidth(2); 
  yieldE = yield*f2->GetParError(0)/f2->GetParameter(0); 
  bg =  f1->Integral(f2->GetParameter(1)-2.*f2->GetParameter(2), f2->GetParameter(1)+2.*f2->GetParameter(2))/hd->GetBinWidth(2); 
  cout << "data yield: " << yield << "  " << bg << endl;
  tl->SetTextFont(62);   tl->SetTextSize(0.04); tl->DrawLatex(xcms, ycms, "CMS preliminary"); 
  tl->SetTextFont(42);  tl->SetTextSize(0.04); 
  tl->DrawLatex(xcms, 0.32, Form("#mu: %5.3f #pm %5.3f GeV", f2->GetParameter(1), f2->GetParError(1))); 
  tl->DrawLatex(xcms, 0.28, Form("#sigma: %5.3f #pm %5.3f GeV", f2->GetParameter(2), f2->GetParError(2))); 
  tl->DrawLatex(xcms, 0.24, Form("N: %4.0f #pm %3.0f", yield, yieldE)); 
  tl->DrawLatex(xcms, 0.20, Form("S/ #sqrt{S+B}: %4.1f", yield/TMath::Sqrt(yield+bg))); 
  c0.SaveAs("d0-da-fit.eps"); 

}


// ----------------------------------------------------------------------
void showHist(TH1D *h, const char *opt = "") {
  h->SetAxisRange(1.71, 1.999); 
  h->SetMinimum(0.1); 
  h->SetLineWidth(2); 
  h->SetMarkerSize(1.3); 
  h->SetXTitle("m_{K #pi} [GeV]");
  h->SetYTitle("Candidates/0.1GeV");

  h->Draw(opt);
}


// ----------------------------------------------------------------------
void fitHist(TH1D *h, const char *opt = "") {
  h->SetAxisRange(1.71, 1.999); 
  h->SetMinimum(0.1); 
  h->SetLineWidth(2); 
  h->SetXTitle("m_{K #pi} [GeV]");
  h->SetYTitle("Candidates/0.1GeV");

  f2->SetLineWidth(2); 
  f2->SetParameters(h->GetMaximum(), 1.86, 0.015, 
		    h->GetBinContent(h->FindBin(1.7)), 
		    h->GetBinContent(h->FindBin(2.0)) - h->GetBinContent(h->FindBin(1.7))
		    ); 

  h->Fit(f2, "", opt);

}


