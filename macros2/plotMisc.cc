#include "plotMisc.hh"

#include <algorithm>

#include "TMath.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"

using namespace std; 
using std::string; 

ClassImp(plotMisc)

// ----------------------------------------------------------------------
plotMisc::plotMisc(const char *files, const char *dir, const char *cuts, int mode) : plotClass(files, dir, cuts, mode) { 
  fNumbersFileName = fDirectory + "/anaBmm.plotMisc." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  cout << "==> plotMisc files: " << files << " dir: " << dir << " cuts: " << cuts << endl;

  string hfname  = fDirectory + "/anaBmm.plotMisc." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  //  if (fHistFile) fHistFile->Close();
  fHistFile = TFile::Open(hfname.c_str(), "UPDATE");

  fIsMC = false;

  fDoUseBDT = true; 
  fSaveSmallTree = false; 

  TH1D *h = (TH1D*)fHistFile->Get("fh0M0_0"); 
  if (0 == h) {
    fNormProcessed = false; 
    for (int i = 0; i < 10; ++i) {
      h = new TH1D(Form("fh0M0_%d", i), Form("fh0M0_%d", i), 60, 4.8, 6.0); 
      fh0M0.insert(make_pair(Form("fh0M0_%d", i), h)); 
      h = new TH1D(Form("fh0M3_%d", i), Form("fh0M3_%d", i), 60, 4.8, 6.0); 
      fh0M3.insert(make_pair(Form("fh0M3_%d", i), h)); 
      h = new TH1D(Form("fh0M4_%d", i), Form("fh0M4_%d", i), 60, 4.8, 6.0); 
      fh0M4.insert(make_pair(Form("fh0M4_%d", i), h)); 
    }
  } else {
    fNormProcessed = true; 
    for (int i = 0; i < 10; ++i) {
      h = (TH1D*)fHistFile->Get(Form("fh0M0_%d", i)); 
      fh0M0.insert(make_pair(Form("fh0M0_%d", i), h)); 
      h = (TH1D*)fHistFile->Get(Form("fh0M3_%d", i)); 
      fh0M3.insert(make_pair(Form("fh0M3_%d", i), h)); 
      h = (TH1D*)fHistFile->Get(Form("fh0M4_%d", i)); 
      fh0M4.insert(make_pair(Form("fh0M4_%d", i), h)); 
    }
  }
}


// ----------------------------------------------------------------------
plotMisc::~plotMisc() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotMisc::makeAll() {
  signalMass(); 
  pidTableDisplay();

}


// ----------------------------------------------------------------------
void plotMisc::loopFunction(int function, int mode) {
  if (1 == function) loopFunction1(mode);
  if (2 == function) loopFunction2(mode);
}


// ----------------------------------------------------------------------
void plotMisc::loopFunction1(int mode) {
  int bin = etaBin(fb.eta); 
  //  cout << "eta = " << fb.eta << " -> bin: " << bin << " with masses: " << fb.m << " " << fb.m3 << " " << fb.m4 << endl;
  if (fb.gmuid && fb.hlt && fb.m1pt>4 && fb.m2pt>4) {
    fh0M0[Form("fh0M0_%d", bin)]->Fill(fb.m); 
    fh0M3[Form("fh0M3_%d", bin)]->Fill(fb.m3); 
    fh0M4[Form("fh0M4_%d", bin)]->Fill(fb.m4); 
  }
}


// ----------------------------------------------------------------------
void plotMisc::loopFunction2(int mode) {
}




// ----------------------------------------------------------------------
void plotMisc::calcTauError() {

  const int MIN(5);

  TFile *f = TFile::Open(Form("%s/anaBmm.plotBDT.%d.root", fDirectory.c_str(), fYear)); 

  TH1D *hd(0), *ha(0);
  int BDTBINS(0); 
  hd = (TH1D*)f->Get("nor_dm_0_0");
  while (hd) {
    hd = (TH1D*)f->Get(Form("nor_dm_0_%d", BDTBINS));
    ++BDTBINS;
  }
  --BDTBINS;
  hd = (TH1D*)f->Get("nor_dm_0_0");
      
  cout << "BDTBINS: " << BDTBINS << endl;

  float bdtCut(0.); 
  string title, cutVal; 

  TH1D *h(0), *hdSlope[fNchan], *haSlope[fNchan], *hDelta[fNchan];
  int NCAT(5); 
  if (2011 == fYear) NCAT = 3; 
  TH1D *hdMassCat[fNchan][NCAT]; 
  TH1D *haMassCat[fNchan][NCAT]; 
  double slope(0.), slopeE(0.); 
  int bCat(-1), histCount(0); 
  if (1) {
    c0->Clear();
    for (int ichan = 0; ichan < fNchan; ++ichan) {
      hdSlope[ichan] = new TH1D(Form("hdSlope_c%d", ichan), Form("slope vs BDT cut C%d", ichan), BDTBINS, -1., 1.); 
      hdSlope[ichan]->Sumw2();
      haSlope[ichan] = new TH1D(Form("haSlope_c%d", ichan), Form("AMS slope vs BDT cut C%d", ichan), BDTBINS, -1.005, 1.005); 
      haSlope[ichan]->Sumw2();

      for (int i = 0; i < BDTBINS; ++i) {
	//      for (int i = 0; i < BDTBINS; i = i+20) {
	hd = (TH1D*)f->Get(Form("nor_dm_%d_%d", ichan, i));
	h = (TH1D*)hd->Clone("h"); 

	title = h->GetTitle(); 
	string::size_type n2 = title.find_last_of(">"); 
	cutVal = title.substr(n2+2); 
	bdtCut = atof(cutVal.c_str()); 

	histCount = h->Integral(h->FindBin(5.45+fEpsilon), h->FindBin(5.9-fEpsilon));
	cout << hd->GetTitle() << " histCount: " << histCount << endl;
	if (histCount < MIN) continue;
	// 	for (int j = 0; j <= h->GetNbinsX(); ++j) {
	// 	  h->SetBinContent(j, h->GetBinContent(j)/histCount); 
	// 	  h->SetBinError(j, h->GetBinError(j)/histCount); 
	// 	}

	// -- hd numbers
	bgBlind(h, 2, 5.45, 5.9);
	c0->SaveAs(Form("%s/%d-amshdxcheckFit_c%d_%d.pdf", fDirectory.c_str(), fYear, ichan, i)); 
	if (h->GetFunction("f1")) {
	  slope  = h->GetFunction("f1")->GetParameter(1)/histCount; 
	  slopeE = h->GetFunction("f1")->GetParError(1)/histCount; 
	  cout << "==> HD " << fYear << " bdtbin:  "  << i << ", b > " << bdtCut << " B0: " << fBdBgExp << " Bs: " << fBsBgExp 
	       << " slope = " << slope << " +/- " << slopeE 
	       << endl;     
	  hdSlope[ichan]->SetBinContent(i, slope); 
	  hdSlope[ichan]->SetBinError(i, slopeE); 
	}
	delete h; 
      
	// -- ha numbers
	ha = (TH1D*)f->Get(Form("ams_am_%d_%d", ichan, i));
	h = (TH1D*)ha->Clone("h"); 

	histCount = h->Integral(h->FindBin(5.45+fEpsilon), h->FindBin(5.9-fEpsilon));
	cout << ha->GetTitle() << " histCount: " << histCount << endl;
	if (histCount < MIN) continue;
	// 	for (int j = 0; j <= h->GetNbinsX(); ++j) {
	// 	  h->SetBinContent(j, h->GetBinContent(j)/histCount); 
	// 	  h->SetBinError(j, h->GetBinError(j)/histCount); 
	// 	}

	bgBlind(h, 2, 5.45, 5.9);
	c0->SaveAs(Form("%s/%d-amshaxcheckFit_c%d_%d.pdf", fDirectory.c_str(), fYear, ichan, i)); 
	if (h->GetFunction("f1")) {
	  slope  = h->GetFunction("f1")->GetParameter(1)/histCount; 
	  slopeE = h->GetFunction("f1")->GetParError(1)/histCount; 
	  cout << "==> HA " << fYear << " bdtbin: "  << i << ", b > " << bdtCut << " B0: " << fBdBgExp << " Bs: " << fBsBgExp 
	       << " slope = " << slope << " +/- " << slopeE
	       << endl;     
	  haSlope[ichan]->SetBinContent(i, slope); 
	  haSlope[ichan]->SetBinError(i, slopeE); 
	}
	delete h; 


      }
    }
  
    c0->Clear(); 
    c0->Divide(1,2);
    
    setHist(hdSlope[0], kBlack); 
    setHist(haSlope[0], kBlue); 
    
    setHist(hdSlope[1], kBlack); 
    setHist(haSlope[1], kBlue); 
    
    c0->cd(1); 
    hdSlope[0]->SetAxisRange(-0.08, 0.08, "Y"); 
    hdSlope[0]->SetAxisRange(0., 0.5-fEpsilon, "X");
    hdSlope[0]->Draw();
    haSlope[0]->Draw("same");
    pl->DrawLine(0., 0., 0.5-fEpsilon, 0.); 
    
    c0->cd(2); 
    hdSlope[1]->SetAxisRange(-0.08, 0.08, "Y"); 
    hdSlope[1]->SetAxisRange(0., 0.5, "X");
    hdSlope[1]->Draw();
    haSlope[1]->Draw("same");

    pl->DrawLine(0., 0., 0.5-fEpsilon, 0.); 
    
    c0->SaveAs(Form("%s/plotMisc-%d-slope-btdCut.pdf", fDirectory.c_str(), fYear)); 
  }

  // -- mass histograms in BDT categories
  if (1) {
    TH1D *hfirst(0), *hlast(0); 
    TH1D *aSlopes[fNchan], *nSlopes[fNchan]; 
    int oldCat(-1), hmin(0), hmax(0); 
    for (int ichan = 0; ichan < fNchan; ++ichan) {
      nSlopes[ichan] = new TH1D(Form("nSlopesCat_c%d", ichan), Form("%d Channel %d", fYear, ichan), NCAT, 0., NCAT);  nSlopes[ichan]->Sumw2();
      setHist(nSlopes[ichan], kBlack); 
      aSlopes[ichan] = new TH1D(Form("aSlopesCat_c%d", ichan), Form("%d AMS Channel %d", fYear, ichan), NCAT, 0.+0.2, NCAT+0.2); aSlopes[ichan]->Sumw2();
      setHist(aSlopes[ichan], kBlue); 
      for (int j = 0; j < NCAT; ++j) {
	hdMassCat[ichan][j] = new TH1D(Form("hdMass_b%d_c%d", j, ichan), Form("%d BDT cat %d C%d", fYear, j, ichan), 
				       hd->GetNbinsX(), hd->GetBinLowEdge(1), hd->GetBinLowEdge(hd->GetNbinsX()+1)); 
	setHist(hdMassCat[ichan][j], kBlack, 24); 
	
	haMassCat[ichan][j] = new TH1D(Form("haMass_b%d_c%d", j, ichan), Form("%d AMS BDT cat %d C%d", fYear, j, ichan), 
				       hd->GetNbinsX(), hd->GetBinLowEdge(1), hd->GetBinLowEdge(hd->GetNbinsX()+1)); 
	setHist(haMassCat[ichan][j], kBlue, 24); 

	bdtCatIdx(j, ichan, hmin, hmax);
	cout << " --> chan " << ichan << " cat " << j << " hists " << hmin << " " << hmax << endl;
      
	hfirst = (TH1D*)f->Get(Form("nor_dm_%d_%d", ichan, hmin));
	hlast = (TH1D*)f->Get(Form("nor_dm_%d_%d", ichan, hmax));

	hdMassCat[ichan][j]->Add(hfirst); 
	hdMassCat[ichan][j]->Add(hlast, -1.); 
	hdMassCat[ichan][j]->Draw();
	c0->SaveAs(Form("%s/plotMisc-hdmassCat%dc%d.pdf", fDirectory.c_str(), j, ichan)); 

	hfirst = (TH1D*)f->Get(Form("ams_am_%d_%d", ichan, hmin));
	hlast = (TH1D*)f->Get(Form("ams_am_%d_%d", ichan, hmax));

	haMassCat[ichan][j]->Add(hfirst); 
	haMassCat[ichan][j]->Add(hlast, -1.); 
	haMassCat[ichan][j]->Draw();
	c0->SaveAs(Form("%s/plotMisc-hamassCat%dc%d.pdf", fDirectory.c_str(), j, ichan)); 
	
      }
    

      c0->Clear();
      double ap1, ap1E, np1, np1E; 
      for (int j = 0; j < NCAT; ++j) {
	h = hdMassCat[ichan][j]; 
	histCount = h->Integral(h->FindBin(5.45+fEpsilon), h->FindBin(5.9-fEpsilon));
	if (histCount < 4) {
	  cout << "XXX no more counts in ha, continue" << endl;
	  continue;
	}
	// 	for (int k = 0; k <= h->GetNbinsX(); ++k) {
	// 	  h->SetBinContent(k, h->GetBinContent(k)/histCount); 
	// 	  h->SetBinError(k, h->GetBinError(k)/histCount); 
	// 	}

	h->SetMinimum(0.);
	setHist(h, kBlack); 
	bgBlind(h, 2, 5.45, 5.9);
	h->GetFunction("f1")->SetLineColor(kBlack); 	h->GetFunction("f1")->SetLineWidth(2); 
	np1  = h->GetFunction("f1")->GetParameter(1)/histCount; 
	np1E = h->GetFunction("f1")->GetParError(1)/histCount; 
	nSlopes[ichan]->SetBinContent(j+1, np1); 
	nSlopes[ichan]->SetBinError(j+1, np1E); 

	h = haMassCat[ichan][j]; 
	histCount = h->Integral(h->FindBin(5.45+fEpsilon), h->FindBin(5.9-fEpsilon));
	if (histCount < 4) {
	  cout << "XXX no more counts in ha, continue" << endl;
	  continue;
	}
	// 	for (int k = 0; k <= h->GetNbinsX(); ++k) {
	// 	  h->SetBinContent(k, h->GetBinContent(k)/histCount); 
	// 	  h->SetBinError(k, h->GetBinError(k)/histCount); 
	// 	}
	h->SetMinimum(0.);
	setHist(h, kBlue); 
	bgBlind(h, 2, 5.45, 5.9);
	h->GetFunction("f1")->SetLineColor(kBlue); h->GetFunction("f1")->SetLineWidth(2); 
	ap1  = h->GetFunction("f1")->GetParameter(1)/histCount; 
	ap1E = h->GetFunction("f1")->GetParError(1)/histCount; 
	aSlopes[ichan]->SetBinContent(j+1, ap1); 
	aSlopes[ichan]->SetBinError(j+1, ap1E); 
	
	c0->Clear();
	hdMassCat[ichan][j]->Draw("");
	haMassCat[ichan][j]->Draw("same");

	tl->SetTextSize(0.03); 
	tl->SetTextColor(kBlack); 
	tl->DrawLatex(0.2, 0.92, Form("p1: %5.4f+/-%5.4f", np1, np1E)); 
	tl->SetTextColor(kBlue); 
	tl->DrawLatex(0.6, 0.92, Form("p1: %5.4f+/-%5.4f", ap1, ap1E)); 
	c0->SaveAs(Form("%s/plotMisc-%d-slope-chan%d-btdCat%d.pdf", fDirectory.c_str(), fYear, ichan, j)); 
      }

      c0->Clear(); 
      nSlopes[ichan]->SetAxisRange(-0.08, 0.08, "Y"); 
      nSlopes[ichan]->Draw();
      aSlopes[ichan]->Draw("same");
      pl->DrawLine(0., 0., NCAT, 0.); 
      c0->SaveAs(Form("%s/plotMisc-%d-slope-chan%d-overlay.pdf", fDirectory.c_str(), fYear, ichan)); 
    }
  }
}


// ----------------------------------------------------------------------
void plotMisc::pidTableDisplay() {

  PidTable *a; 
  gStyle->SetOptTitle(0);
  tl->SetTextSize(0.07); 

  // -- fake rates
  if (1) {
    gStyle->SetPaintTextFormat("5.4f");
    double xbins[] = {0., 1.4, 2.4};
    double ybins[] = {0., 4., 5., 7., 10., 15., 20., 50.}; 
    TH2D *h2 = new TH2D("h2", "", 2, xbins, 7, ybins); 
    setTitles(h2, "|#eta|", "p_{T} [GeV]"); 
    h2->SetMinimum(0.0); 
    h2->SetMaximum(0.002); 
    h2->SetMarkerSize(1.3);
    h2->SetMarkerColor(kBlack);
    
    gStyle->SetOptStat(0); 
    
    shrinkPad(0.15, 0.15, 0.25); 
    gPad->SetLogy(0); 
    a = fptFakePosKaons;  h2->Reset(); h2->SetTitle(Form("positive kaons (%d)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-fakePosKaons.pdf", fDirectory.c_str(), fYear));
    
    a = fptFakeNegKaons;  h2->Reset(); h2->SetTitle(Form("negative kaons (%d)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-fakeNegKaons.pdf", fDirectory.c_str(), fYear));
    
    a = fptFakePosPions;  h2->Reset(); h2->SetTitle(Form("positive pions (%d)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-fakePosPions.pdf", fDirectory.c_str(), fYear));
    
    a = fptFakeNegPions;  h2->Reset(); h2->SetTitle(Form("negative pions (%d)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-fakeNegPions.pdf", fDirectory.c_str(), fYear));
    
    a = fptFakePosProtons;  h2->Reset(); h2->SetTitle(Form("positive protons (%d)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-fakePosProtons.pdf", fDirectory.c_str(), fYear));
    
    a = fptFakeNegProtons;  h2->Reset(); h2->SetTitle(Form("negative protons (%d)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-fakeNegProtons.pdf", fDirectory.c_str(), fYear));
  }

  // -- muon id and trigger
  gStyle->SetPaintTextFormat("3.2f");
  if (1) {
    double xbins[] = {0.0, 0.2, 0.3, 0.6, 0.8, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};
    double ybins[] = {4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 15., 20., 50.}; 
    TH2D *h2 = new TH2D("h2", "", 11, xbins, 10, ybins); 
    setTitles(h2, "|#eta|", "p_{T} [GeV]"); 
    h2->SetMinimum(0.0); 
    h2->SetMaximum(1.0); 
    h2->SetMarkerSize(1.3);
    //    h2->SetMarkerColor(kWhite);
    
    gStyle->SetOptStat(0); 
    
    shrinkPad(0.15, 0.15, 0.25); 
    gPad->SetLogy(1); 

    // -- muon id 
    a = fptSgM;  h2->Reset(); h2->SetTitle(Form("%d seagulls muon id (data)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-sg-da-muonid.pdf", fDirectory.c_str(), fYear));

    a = fptCbM;  h2->Reset(); h2->SetTitle(Form("%d cowboys muon id (data)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-cb-da-muonid.pdf", fDirectory.c_str(), fYear));

    a = fptSgMMC; h2->Reset(); h2->SetTitle(Form("%d seagulls muon id (MC)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-sg-mc-muonid.pdf", fDirectory.c_str(), fYear));

    a = fptCbMMC;  h2->Reset(); h2->SetTitle(Form("%d cowboys muon id (MC)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-cb-mc-muonid.pdf", fDirectory.c_str(), fYear));

    // -- L1L2 
    a = fptSgT1;  h2->Reset(); h2->SetTitle(Form("%d seagulls L1L2 (data)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-sg-da-L1L2.pdf", fDirectory.c_str(), fYear));

    a = fptCbT1;  h2->Reset(); h2->SetTitle(Form("%d cowboys L1L2 (data)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-cb-da-L1L2.pdf", fDirectory.c_str(), fYear));

    a = fptSgT1MC; h2->Reset(); h2->SetTitle(Form("%d seagulls L1L2 (MC)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-sg-mc-L1L2.pdf", fDirectory.c_str(), fYear));

    a = fptCbT1MC;  h2->Reset(); h2->SetTitle(Form("%d cowboys L1L2 (MC)", fYear)); 
    a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
    c0->SaveAs(Form("%s/%d-cb-mc-L1L2.pdf", fDirectory.c_str(), fYear));


    // -- L3
    {
      gPad->SetLogy(1); 
      double xbins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
      double ybins[] = {4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 15., 20., 50.}; 
      TH2D *h2 = new TH2D("h2", "", 6, xbins, 10, ybins); 
      setTitles(h2, "|#eta|", "p_{T} [GeV]"); 
      h2->SetMinimum(0.0); 
      h2->SetMaximum(1.0); 
      h2->SetMarkerSize(1.3);
      //      h2->SetMarkerColor(kWhite);
      
      a = fptSgT2;  h2->Reset(); h2->SetTitle(Form("%d seagulls L3 (data)", fYear)); 
      a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
      c0->SaveAs(Form("%s/%d-sg-da-L3.pdf", fDirectory.c_str(), fYear));
      
      a = fptCbT2;  h2->Reset(); h2->SetTitle(Form("%d cowboys L3 (data)", fYear)); 
      a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
      c0->SaveAs(Form("%s/%d-cb-da-L3.pdf", fDirectory.c_str(), fYear));
      
      a = fptSgT2MC; h2->Reset(); h2->SetTitle(Form("%d seagulls L3 (MC)", fYear)); 
      a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
      c0->SaveAs(Form("%s/%d-sg-mc-L3.pdf", fDirectory.c_str(), fYear));
      
      a = fptCbT2MC;  h2->Reset(); h2->SetTitle(Form("%d cowboys L3 (MC)", fYear)); 
      a->eff2d(h2);  h2->Draw("colztext");  tl->DrawLatex(0.10, 0.92, h2->GetTitle()); 
      c0->SaveAs(Form("%s/%d-cb-mc-L3.pdf", fDirectory.c_str(), fYear));
    }

  }

}



// ----------------------------------------------------------------------
void plotMisc::fakeRateOverlaysDK(string mode) {

  double ymax(2.5e-3); 

  // K, data
  const int points21 = 8;
  //                                      mu12  ht12  jet12 el12 psik12 phi12 mu11  ht11
  double v_x21[points21] = {8.04, 11.5, 11.4, 10.3, 8.28, 6.56, 7.79, 11.7};  // pt
  double v_y21[points21] = {0.88e-3,  0.1e-3, 0.12e-3, 0.69e-3, 0.93e-3, 0.87e-3, 0.46e-3,  1.1e-3};  // ratio
  double e_y21[points21] = {0.26e-3,  1.6e-3, 0.16e-3, 0.28e-3, 0.23e-3, 0.20e-3, 0.22e-3,  0.3e-3};  // error
  
  // K, MC
  const int points22 = 4;
  //                     dstar  psik  phi   bskk
  double v_x22[points22] = {7.18, 8.59, 6.80, 7.17};
  double v_y22[points22] = {1.63e-3, 0.87e-3, 0.86e-3, 1.40e-3};
  double e_y22[points22] = {0.17e-3, 0.08e-3, 0.12e-3, 0.20e-3};


  // pi, data
  const int points23 = 7;
  //                      mu12  ht12  jet12 el12  mu11  ht11  mike(Ks)
  double v_x23[points23] = {8.74, 12.2, 12.1, 11.0, 8.41, 12.5, 8.1};
  double v_y23[points23] = {1.19e-3, 0.1e-3,  0.65e-3, 1.53e-3, 0.77e-3, 1.8e-3, 0.65e-3};
  double e_y23[points23] = {0.28e-3, 1.2e-3,  0.21e-3, 0.30e-3, 0.24e-3, 0.8e-3, 0.04e-3};
  
  // pi,mc
  const int points24 = 2;
  //                    dstar-12 bdpipi
  double v_x24[points24] = {7.19,  7.16};
  double v_y24[points24] = {0.87e-3,  0.86e-3};
  double e_y24[points24] = {0.13e-3,  0.02e-3};


  TH1D *hpt = new TH1D("hpt", "", 20, 0., 20.); hpt->Sumw2();
  setHist(hpt, kBlue, 20, 0.001);
  hpt->SetFillColor(kBlue); 
  //  setFilledHist(hpt);

  TGraphErrors *gD(0), *gM(0); 

  if (mode == "kaons") {
    fptFakePosKaons->projectP(hpt, -2.4, 2.4, -180., 180., 1); 
    for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
    gD = new TGraphErrors(points21, v_x21, v_y21, 0, e_y21); 
    gM = new TGraphErrors(points22, v_x22, v_y22, 0, e_y22); 
  } 

  if (mode == "pions") {
    fptFakePosPions->projectP(hpt, -2.4, 2.4, -180., 180., 1); 
    for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
    gD = new TGraphErrors(points23, v_x23, v_y23, 0, e_y23); 
    gM = new TGraphErrors(points24, v_x24, v_y24, 0, e_y24); 
  } 

  if (mode == "protons") {
    fptFakePosProtons->projectP(hpt, -2.4, 2.4, -180., 180., 1); 
    for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
  } 

  if (gD) {
    gD->SetMarkerStyle(20); 
    gD->SetMarkerColor(kBlack); 
    gD->SetMarkerSize(1.5); 
    gM->SetMarkerStyle(20); 
    gM->SetMarkerColor(kRed); 
    gM->SetMarkerSize(1.5); 
  }

  zone(); 
  shrinkPad(0.11, 0.25); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  TH1D *hFrame = new TH1D("hFrame", "", 10, 0., 20.);
  setTitles(hFrame, "p_{T} [GeV]", "fake rate", 0.05, 1.1, 2.1);
  hFrame->SetMinimum(0.); 
  hFrame->SetMaximum(ymax); 
  hFrame->Draw(); 
  hpt->Draw("e5same][");
  hpt->Draw("histsame][");
  if (gD) {
    gD->Draw("p");
    gM->Draw("p");
  }

  newLegend(0.6, 0.65, 0.80, 0.89);
  if (mode == "pions") legg->SetHeader("Pion fake rates"); 
  if (mode == "kaons") legg->SetHeader("Kaon fake rates"); 
  if (mode == "protons") legg->SetHeader("Proton fake rates"); 
  legg->AddEntry(hpt, "MC vs p_{T}", "f"); 
  if (gD) {
    legg->AddEntry(gM, "MC vs <p_{T}>", "p"); 
    legg->AddEntry(gD, "data vs <p_{T}>", "p"); 
  }
  legg->Draw(); 

  c0->SaveAs(Form("%s/%s-fakeRateOverlaysDK-%s.pdf", fDirectory.c_str(), fSuffix.c_str(), mode.c_str())); 
  
  
}


// ----------------------------------------------------------------------
void plotMisc::fakeRateOverlaysMM(string mode) {

  if (mode == "all") {
    fakeRateOverlaysMM("pions");
    fakeRateOverlaysMM("pions");

    fakeRateOverlaysMM("protons");
    fakeRateOverlaysMM("protons");
    return;
  }


  double ymax(2.5e-3); 


// TightMVA Muons:
//   Pions:
//     Pt 4-6 GeV:  0.00066 +- 0.00005
//     Pt 6-10 GeV:  0.00066 +- 0.00008
//     Pt 10-30 GeV:  0.00048 +- 0.00009
//     |Eta| 0.0-0.9:  0.00056 +- 0.00004
//     |Eta| 0.9-1.4:  0.00080 +- 0.00008
//     |Eta| 1.4-2.4:  0.00092 +- 0.00013
//     Integrated:  0.00064 +- 0.00003
//   Protons:
//     Pt 4-6 GeV:  0.00013 +- 0.00007
//     Pt 6-10 GeV:  0.00010 +- 0.00007
//     Pt 10-30 GeV:  0.00000 +- 0.00008
//     |Eta| 0.0-0.9:  0.00011 +- 0.00004
//     |Eta| 0.9-1.4:  0.00009 +- 0.00010
//     |Eta| 1.4-2.4:  0.00000 +- 0.00010
//     Integrated:  0.00009 +- 0.00004


  // pi, data
  const int points21 = 3;
  double v_x21[points21] = {5.0,     8.0,     20.}; 
  double v_y21[points21] = {0.00066, 0.00066, 0.00048};
  double e_y21[points21] = {0.00005, 0.00008, 0.00009};

  // p, data
  const int points23 = 3;
  double v_x23[points23] = {5.0, 8.0, 20.}; 
  double v_y23[points23] = {0.00013, 0.00010, 0.00000};
  double e_y23[points23] = {0.00007, 0.00007, 0.00008};
  

  TH1D *hpt = new TH1D("hpt", "", 20, 0., 20.); hpt->Sumw2();
  setHist(hpt, kBlue, 20, 0.001);
  hpt->SetFillColor(kBlue); 
  //  setFilledHist(hpt);

  TGraphErrors *gD(0); 

  if (mode == "pions") {
    fptFakePosPions->projectP(hpt, -2.4, 2.4, -180., 180., 1); 
    for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
    gD = new TGraphErrors(points21, v_x21, v_y21, 0, e_y21); 
  } 

  if (mode == "protons") {
    fptFakePosProtons->projectP(hpt, -2.4, 2.4, -180., 180., 1); 
    for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
    gD = new TGraphErrors(points23, v_x23, v_y23, 0, e_y23); 
  } 

  if (gD) {
    gD->SetMarkerStyle(20); 
    gD->SetMarkerColor(kBlack); 
    gD->SetMarkerSize(1.5); 
  }

  zone(); 
  shrinkPad(0.11, 0.25); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  TH1D *hFrame = new TH1D("hFrame", "", 10, 0., 20.);
  setTitles(hFrame, "p_{T} [GeV]", "fake rate", 0.05, 1.1, 2.1);
  hFrame->SetMinimum(0.); 
  hFrame->SetMaximum(ymax); 
  hFrame->Draw(); 
  hpt->Draw("e5same][");
  hpt->Draw("histsame][");
  if (gD) {
    gD->Draw("p");
  }

  newLegend(0.6, 0.65, 0.80, 0.89);
  if (mode == "pions") legg->SetHeader("Pion fake rates"); 
  if (mode == "kaons") legg->SetHeader("Kaon fake rates"); 
  if (mode == "protons") legg->SetHeader("Proton fake rates"); 
  legg->AddEntry(hpt, "MC", "f"); 
  if (gD) {
    legg->AddEntry(gD, "data", "p"); 
  }
  legg->Draw(); 

  c0->SaveAs(Form("%s/%s-fakeRateOverlaysMM-%s.pdf", fDirectory.c_str(), fSuffix.c_str(), mode.c_str())); 

}


// ----------------------------------------------------------------------
void plotMisc::fakeRateOverlaysMG(string mode, string charge, string chan) {

  if (mode == "all") {
    fakeRateOverlaysMG("pions", "pos", "barrel");
    fakeRateOverlaysMG("pions", "neg", "barrel");

    fakeRateOverlaysMG("kaons", "pos", "barrel");
    fakeRateOverlaysMG("kaons", "neg", "barrel");

    fakeRateOverlaysMG("pions", "pos", "endcap");
    fakeRateOverlaysMG("pions", "neg", "endcap");

    fakeRateOverlaysMG("kaons", "pos", "endcap");
    fakeRateOverlaysMG("kaons", "neg", "endcap");
    return;
  }

  double ymax(2.5e-3); 

  TH1D *hpt = new TH1D("hpt", "", 20, 0., 20.); hpt->Sumw2();
  setHist(hpt, kBlue, 20, 0.001);
  hpt->SetFillColor(kBlue); 
  //  setFilledHist(hpt);

  TH1D *marioA(0), *marioP(0), *marioE(0); 

  PidTable *fptFakePosPions = new PidTable("../macros/pidtables/130702/2012-pionPosFakeRate-mvaMuon.dat"); 
  PidTable *fptFakeNegPions = new PidTable("../macros/pidtables/130702/2012-pionNegFakeRate-mvaMuon.dat"); 

  PidTable *fptFakePosKaons = new PidTable("../macros/pidtables/130702/2012-kaonPosFakeRate-mvaMuon.dat"); 
  PidTable *fptFakeNegKaons = new PidTable("../macros/pidtables/130702/2012-kaonNegFakeRate-mvaMuon.dat"); 

  string etaString = (chan=="barrel"? "Eta_0_1.4" : "Eta_1.4_2.4"); 

  if (mode == "pions") {
    if (charge == "pos") {
      fptFakePosPions->projectP(hpt, -1.4, 1.4, -180., 180., 1); 
      for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
      TFile *fa = TFile::Open(Form("mg/hAllGenParticlePlusRecoPt_%s_pion.root", chan.c_str())); 
      marioA = (TH1D*)fa->Get(Form("hAssocTkGenParticlePlusRecoPt_%s", etaString.c_str())); 
      cout << "marioA: " << marioA << endl;
      TFile *fp = TFile::Open(Form("mg/hAssocGenParticlePlusRecoPt_%s_pion.root", chan.c_str())); 
      marioP = (TH1D*)fp->Get(Form("hAssocGenParticlePlusRecoPt_%s", etaString.c_str())); 
      cout << "marioP: " << marioP << endl;
      marioE = (TH1D*)marioA->Clone("eff"); marioE->Reset(); 
      marioE->Divide(marioP, marioA, 1., 1., "b");

    } else {
      fptFakeNegPions->projectP(hpt, -1.4, 1.4, -180., 180., 1); 
      for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
      TFile *fa = TFile::Open(Form("mg/hAllGenParticleMinusRecoPt_%s_pion.root", chan.c_str())); 
      marioA = (TH1D*)fa->Get(Form("hAssocTkGenParticleMinusRecoPt_%s", etaString.c_str()));  
      cout << "marioA: " << marioA << endl;
      TFile *fp = TFile::Open(Form("mg/hAssocGenParticleMinusRecoPt_%s_pion.root", chan.c_str())); 
      marioP = (TH1D*)fp->Get(Form("hAssocGenParticleMinusRecoPt_%s", etaString.c_str())); 
      cout << "marioP: " << marioP << endl;
      marioE = (TH1D*)marioA->Clone("eff"); marioE->Reset(); 
      marioE->Divide(marioP, marioA, 1., 1., "b");
    }
  } 

  if (mode == "kaons") {
    if (charge == "pos") {
      fptFakePosKaons->projectP(hpt, -1.4, 1.4, -180., 180., 1); 
      for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
      TFile *fa = TFile::Open(Form("mg/hAllGenParticlePlusRecoPt_%s_kaon.root", chan.c_str())); 
      marioA = (TH1D*)fa->Get(Form("hAssocTkGenParticlePlusRecoPt_%s", etaString.c_str())); 
      cout << "marioA: " << marioA << endl;
      TFile *fp = TFile::Open(Form("mg/hAssocGenParticlePlusRecoPt_%s_kaon.root", chan.c_str())); 
      marioP = (TH1D*)fp->Get(Form("hAssocGenParticlePlusRecoPt_%s", etaString.c_str())); 
      cout << "marioP: " << marioP << endl;
    } else {
      fptFakeNegKaons->projectP(hpt, -1.4, 1.4, -180., 180., 1); 
      for (int i = 1; i < hpt->GetNbinsX(); ++i) hpt->SetBinError(i, 0.5*hpt->GetBinContent(i)); 
      TFile *fa = TFile::Open(Form("mg/hAllGenParticleMinusRecoPt_%s_kaon.root", chan.c_str())); 
      marioA = (TH1D*)fa->Get(Form("hAssocTkGenParticleMinusRecoPt_%s", etaString.c_str())); 
      cout << "marioA: " << marioA << endl;
      TFile *fp = TFile::Open(Form("mg/hAssocGenParticleMinusRecoPt_%s_kaon.root", chan.c_str())); 
      marioP = (TH1D*)fp->Get(Form("hAssocGenParticleMinusRecoPt_%s", etaString.c_str())); 
      cout << "marioP: " << marioP << endl;
    }

    marioE = (TH1D*)marioA->Clone("eff"); marioE->Reset(); 
    marioE->Divide(marioP, marioA, 1., 1., "b");

  } 


  zone(); 
  shrinkPad(0.11, 0.25); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  TH1D *hFrame = new TH1D("hFrame", "", 10, 0., 20.);
  setTitles(hFrame, "p_{T} [GeV]", "fake rate", 0.05, 1.1, 2.1);
  hFrame->SetMinimum(0.); 
  hFrame->SetMaximum(ymax); 
  hFrame->Draw(); 
  hpt->Draw("e5same][");
  hpt->Draw("histsame][");
  marioE->Draw("esame");

  tl->SetTextColor(kBlack); tl->DrawLatex(0.6, 0.82, Form("%s %s", charge.c_str(), mode.c_str())); 
  tl->SetTextColor(kBlack); tl->DrawLatex(0.6, 0.75, chan.c_str()); 
  tl->SetTextColor(kBlack); tl->DrawLatex(0.25, 0.92, "Mario"); 
  tl->SetTextColor(kBlue);  tl->DrawLatex(0.6, 0.92, "Urs"); 

  c0->SaveAs(Form("%s/%s-fakeRateOverlaysMG-%s-%s-%s.pdf", fDirectory.c_str(), fSuffix.c_str(), chan.c_str(), mode.c_str(), charge.c_str())); 

}


// ----------------------------------------------------------------------
void plotMisc::massError() {

  double bdtmin(0.3); 

  TProfile *hd = new TProfile("hd", "", 100, -2.5, 2.5); 
  hd->SetMinimum(0.);
  hd->SetMaximum(0.1);

  TProfile *hm = new TProfile("hm", "", 100, -2.5, 2.5); 
  hm->SetMinimum(0.01);
  hm->SetMaximum(0.06);

  string mode("NoMc"); 
  TTree *t = getTree(mode); 
  t->Draw("me:eta>>hd", Form("bdt>%f", bdtmin), "goff");

  mode = "NoData"; 
  t = getTree(mode); 
  t->Draw("me:eta>>hm", Form("bdt>%f", bdtmin), "goff");

  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  setTitles(hd, "#eta(B)", "mass error [GeV]", 0.05, 1.1, 1.5); 
  hd->Draw();
  hm->Draw("samehist");
  
  newLegend(0.3, 0.6, 0.5, 0.7);
  legg->SetHeader("Dimuon (J/#psi) mass error"); 
  legg->AddEntry(hd, "data", "p"); 
  legg->AddEntry(hm, "MC", "f"); 
  legg->Draw(); 

  int i = 0; 
  c0->SaveAs(Form("%s/%s-me-eta-%d.pdf", fDirectory.c_str(), fSuffix.c_str(), i)); 
  

}

// ----------------------------------------------------------------------
void plotMisc::signalMass() {

  TVirtualFitter::SetMaxIterations(20000);

  fHistFile->cd(); 
  if (!fNormProcessed) {
    fF["SgMc"]->cd("candAnaMuMu"); 
    string mode("SgMc"); 
    TTree *t = getTree(mode); 
    if (0 == t) {
      cout << "plotBDT: no tree found for mode  " << mode << endl;
      return;
    } 
    int nevts2use(-1); 
    setupTree(t, mode);
    if (t->GetEntries() > 5000000) nevts2use = 5000000; 
    loopOverTree(t, fSetup, 1, nevts2use);
  }

  double eps(0.1), eps2(2*eps); 
  TH1D *hM0 = new TH1D("hM0", "vtx mass", 10, 0., 10.); hM0->Sumw2(); setHist(hM0, kBlack, 20, 2); 
  TH1D *hM3 = new TH1D("hM3", "vtx mass, corrected by mu scale", 10, 0.+eps, 10.+eps); hM3->Sumw2(); setHist(hM3, kBlue, 25, 2); 
  TH1D *hM4 = new TH1D("hM4", "Mu scale mass", 10, 0.+eps2, 10.+eps2); hM4->Sumw2(); setHist(hM4, kRed, 22, 2); 
  setTitles(hM0, "candidate #eta bin", "mass [GeV]", 0.05, 1.1, 1.5); 
  for (int i = 0; i < 10; ++i) {
    hM0->GetXaxis()->SetBinLabel(i+1, Form("%3.2f <|#eta|<%3.2f", i*0.25, (i+1)*0.25)); 
  }
  hM0->SetLabelSize(0.04, "x"); 

  TF1 *f(0); 
  double mean(0.), sigma(0.04); 
  // -- fit histograms: m0
  int ibin(0); 
  for (map<string, TH1D*>::iterator imap = fh0M0.begin(); imap != fh0M0.end(); ++imap) {  
    TH1D *h = imap->second;
    if (h->GetSumOfWeights() < 100) continue;
    sscanf(h->GetName(), "fh0M0_%d", &ibin); 
    f = fpFunc->crystalBall(h, 5.37, sigma, 1., 1.); 
    h->Fit(f); 
    c0->SaveAs(Form("%s/misc-signalMass-m0-%s.pdf", fDirectory.c_str(), imap->first.c_str())); 
    mean  = f->GetParameter(0); 
    sigma = f->GetParameter(1); 
    cout << "ibin = " << ibin << " mean: " << mean << " sigma: " << sigma << endl;
    hM0->SetBinContent(ibin+1, mean); 
    hM0->SetBinError(ibin+1, sigma); 
    c0->Modified(); 
    c0->Update();
  }

  hM0->SetMinimum(5.27); 
  hM0->SetMaximum(5.47); 
  hM0->Draw();
  c0->SaveAs(Form("%s/misc-signalMass-m0-eta.pdf", fDirectory.c_str())); 

  // -- fit histograms: m3
  mean = 5.365;
  sigma = 0.04; 
  for (map<string, TH1D*>::iterator imap = fh0M3.begin(); imap != fh0M3.end(); ++imap) {  
    TH1D *h = imap->second;
    if (h->GetSumOfWeights() < 100) continue;
    sscanf(h->GetName(), "fh0M3_%d", &ibin); 
    f = fpFunc->crystalBall(h, mean, 1.1*sigma, 1., 1.); 
    h->Fit(f); 
    c0->SaveAs(Form("%s/misc-signalMass-m3-%s.pdf", fDirectory.c_str(), imap->first.c_str())); 
    mean  = f->GetParameter(0); 
    sigma = f->GetParameter(1); 
    cout << "ibin = " << ibin << " mean: " << mean << " sigma: " << sigma << endl;
    hM3->SetBinContent(ibin+1, mean); 
    hM3->SetBinError(ibin+1, sigma); 
    c0->Modified(); 
    c0->Update();
  }

  hM3->SetMinimum(5.27); 
  hM3->SetMaximum(5.47); 
  hM3->Draw();
  c0->SaveAs(Form("%s/misc-signalMass-m3-eta.pdf", fDirectory.c_str())); 

  // -- fit histograms: m4
  mean = 5.365;
  sigma = 0.04; 
  for (map<string, TH1D*>::iterator imap = fh0M4.begin(); imap != fh0M4.end(); ++imap) {  
    TH1D *h = imap->second;
    sscanf(h->GetName(), "fh0M4_%d", &ibin); 
    if (h->GetSumOfWeights() < 100) continue;
    f = fpFunc->crystalBall(h, mean, 1.05*sigma, 1., 1.); 
    h->Fit(f); 
    c0->SaveAs(Form("%s/misc-signalMass-m4-%s.pdf", fDirectory.c_str(), imap->first.c_str())); 
    //    h->Fit("gaus");
    mean  = f->GetParameter(0); 
    sigma = f->GetParameter(1); 
    cout << "ibin = " << ibin << " mean: " << mean << " sigma: " << sigma << endl;
    hM4->SetBinContent(ibin+1, mean); 
    hM4->SetBinError(ibin+1, sigma); 
    c0->Modified(); 
    c0->Update();
    delete f; 
  }

  hM4->SetMinimum(5.27); 
  hM4->SetMaximum(5.47); 
  hM4->Draw();
  c0->SaveAs(Form("%s/misc-signalMass-m4-eta.pdf", fDirectory.c_str())); 

  gStyle->SetOptTitle(0);
  hM0->Draw();
  hM3->Draw("same");
  hM4->Draw("same");

  double mmc = 5.3696; 
  pl->DrawLine(0., mmc, 10., mmc);

  double x(0.20); 
  newLegend(x, 0.75, x+0.3, 0.88);
  legg->AddEntry(hM0, "vertex-constrained mass", "p"); 
  legg->AddEntry(hM3, "#mu-scale corrector mass", "p"); 
  legg->AddEntry(hM4, "vertex-c./#mu-scale corr. mass", "p"); 
  legg->Draw();

  c0->SaveAs(Form("%s/misc-%d-signalMass-m034-eta.pdf", fDirectory.c_str(), fYear)); 

}


// ----------------------------------------------------------------------
void plotMisc::effImpactTrkHit() {

  TFile *f = TFile::Open(Form("%s/anaBmm.plotResults.2012.root", fDirectory.c_str())); 
  TH1D *hAnPass(0), *hCoPass(0), *hAnAll(0), *hCoAll(0);
  int ybin(30); 
  for (int y = 0; y < 2; ++y) {
    ybin = 30 + y*2; 
    for (int i = 0; i < fNchan; ++i) {
      hCoPass = (TH1D*)f->Get(Form("hMassWithMuonCuts_bdt_%d_chan%d", ybin+0, i));
      hAnPass = (TH1D*)f->Get(Form("hMassWithMuonCuts_bdt_%d_chan%d", ybin+1, i));

      hCoAll  = (TH1D*)f->Get(Form("hMassWithAnaCuts_bdt_%d_chan%d", ybin+0, i));
      hAnAll  = (TH1D*)f->Get(Form("hMassWithAnaCuts_bdt_%d_chan%d", ybin+1, i));

      cout << "Year: " << (y == 0?2011:2012) << " chan " << i << " analysis: " << hAnPass->Integral()/hAnAll->Integral() << endl;
      cout << "Year: " << (y == 0?2011:2012) << " chan " << i << " correct:  " << hCoPass->Integral()/hCoAll->Integral() << endl;
      cout << "relative difference: " << (hAnPass->Integral()/hAnAll->Integral()) / (hCoPass->Integral()/hCoAll->Integral()) << endl;
    }
  }
    


}




// ----------------------------------------------------------------------
int plotMisc::etaBin(double eta) {
  const int NSTEP(10); 
  double step(2.5/NSTEP); 
  double aeta = TMath::Abs(eta); 
  for (int i = 0 ; i < NSTEP; ++i) {
    if ((i*step < aeta) && (aeta < (i+1)*step)) return i; 
  }
  
  return 9; 
}

// ----------------------------------------------------------------------
int plotMisc::bdtCat(double bdt, int chan) {

  if (2011 == fYear) {
    if (0 == chan) {
      if (0.10 < bdt && bdt <= 0.31) return 0; 
      if (0.31 < bdt && bdt <= 1.00) return 1; 
    } else if (1 == chan) {
      if (0.10 < bdt && bdt <= 0.26) return 0; 
      if (0.26 < bdt && bdt <= 1.00) return 1; 
    }
  } else if (2012 == fYear) {
    if (0 == chan) {
      if (0.10 < bdt && bdt <= 0.23) return 0; 
      if (0.23 < bdt && bdt <= 0.33) return 1; 
      if (0.33 < bdt && bdt <= 0.44) return 2; 
      if (0.44 < bdt && bdt <= 1.00) return 3; 
    } else if (1 == chan) {
      if (0.10 < bdt && bdt <= 0.22) return 0; 
      if (0.22 < bdt && bdt <= 0.29) return 1; 
      if (0.29 < bdt && bdt <= 0.45) return 2; 
      if (0.45 < bdt && bdt <= 1.00) return 3; 
    }
  }
  return -1; 

}

// ----------------------------------------------------------------------
void plotMisc::bdtCatIdx(int cat, int chan, int &hmin, int &hmax) {
  
  if (2011 == fYear) {
    if (0 == chan) {
      if (0 == cat) {hmin = 120; hmax = 199;}
      if (1 == cat) {hmin = 110; hmax = 130;}
      if (2 == cat) {hmin = 131; hmax = 199;}
    } else if (1 == chan) {
      if (0 == cat) {hmin = 120; hmax = 199;}
      if (1 == cat) {hmin = 110; hmax = 125;}
      if (2 == cat) {hmin = 126; hmax = 199;}
    }
  } else if (2012 == fYear) {
    if (0 == chan) {
      if (0 == cat) {hmin = 120; hmax = 199;}
      if (1 == cat) {hmin = 110; hmax = 122;}
      if (2 == cat) {hmin = 123; hmax = 132;}
      if (3 == cat) {hmin = 133; hmax = 143;}
      if (4 == cat) {hmin = 144; hmax = 199;}
    } else if (1 == chan) {
      if (0 == cat) {hmin = 120; hmax = 199;}
      if (1 == cat) {hmin = 110; hmax = 121;}
      if (2 == cat) {hmin = 122; hmax = 128;}
      if (3 == cat) {hmin = 129; hmax = 144;}
      if (4 == cat) {hmin = 145; hmax = 199;}
    }
  }

}
