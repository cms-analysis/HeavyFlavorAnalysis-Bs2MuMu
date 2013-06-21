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
void plotMisc::fakeRateOverlays(string mode) {

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

  c0->SaveAs(Form("%s/%s-fakeRateOverlays-%s.pdf", fDirectory.c_str(), fSuffix.c_str(), mode.c_str())); 
  
  
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
int plotMisc::etaBin(double eta) {
  const int NSTEP(10); 
  double step(2.5/NSTEP); 
  double aeta = TMath::Abs(eta); 
  for (int i = 0 ; i < NSTEP; ++i) {
    if ((i*step < aeta) && (aeta < (i+1)*step)) return i; 
  }
  
  return 9; 
}
