#include "plotMisc.hh"

#include <algorithm>

#include "TMath.h"

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
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

  fIsMC = false;

  fDoUseBDT = true; 
  fSaveSmallTree = false; 

  TH1D *h(0); 
  for (int i = 0; i < 10; ++i) {
    h = new TH1D(Form("fh0M0_%d", i), Form("fh0M0_%d", i), 60, 4.8, 6.0); 
    fh0M0.insert(make_pair(Form("fh0M0_%d", i), h)); 
    h = new TH1D(Form("fh0M3_%d", i), Form("fh0M3_%d", i), 60, 4.8, 6.0); 
    fh0M3.insert(make_pair(Form("fh0M3_%d", i), h)); 
    h = new TH1D(Form("fh0M4_%d", i), Form("fh0M4_%d", i), 60, 4.8, 6.0); 
    fh0M4.insert(make_pair(Form("fh0M4_%d", i), h)); 
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

  // -- fake rates
  double xbins[] = {0., 1.4, 2.4};
  double ybins[] = {0., 4., 5., 7., 10., 15., 20., 50.}; 
  TH2D *h2 = new TH2D("h2", "", 2, xbins, 7, ybins); 
  setTitles(h2, "|#eta|", "p_{T} [GeV]"); 
  h2->SetMinimum(0.0); 
  h2->SetMaximum(0.001); 
  h2->SetMarkerSize(1.3);
  h2->SetMarkerColor(kWhite);
  
  gStyle->SetOptStat(0); 

  shrinkPad(0.15, 0.15, 0.25); 
  a = fptFakePosKaons;  h2->Reset(); h2->SetTitle("positive kaons"); 
  a->eff2d(h2);  h2->Draw("colztext");  c0->SaveAs(Form("%s/%d-fakePosKaons.pdf", fDirectory.c_str(), fYear));

  a = fptFakeNegKaons;  h2->Reset(); h2->SetTitle("negative kaons"); 
  a->eff2d(h2);  h2->Draw("colztext");  c0->SaveAs(Form("%s/%d-fakeNegKaons.pdf", fDirectory.c_str(), fYear));
																              
  a = fptFakePosPions;  h2->Reset(); h2->SetTitle("positive pions"); 
  a->eff2d(h2);  h2->Draw("colztext");  c0->SaveAs(Form("%s/%d-fakePosPions.pdf", fDirectory.c_str(), fYear));
  
  a = fptFakeNegPions;  h2->Reset(); h2->SetTitle("negative pions"); 
  a->eff2d(h2);  h2->Draw("colztext");  c0->SaveAs(Form("%s/%d-fakeNegPions.pdf", fDirectory.c_str(), fYear));

  a = fptFakePosProtons;  h2->Reset(); h2->SetTitle("positive protons"); 
  a->eff2d(h2);  h2->Draw("colztext");  c0->SaveAs(Form("%s/%d-fakePosProtons.pdf", fDirectory.c_str(), fYear));

  a = fptFakeNegProtons;  h2->Reset(); h2->SetTitle("negative protons"); 
  a->eff2d(h2);  h2->Draw("colztext");  c0->SaveAs(Form("%s/%d-fakeNegProtons.pdf", fDirectory.c_str(), fYear));


//   PidTable *fptT1;
//   PidTable *fptT2;
//   PidTable *fptM; 

//   PidTable *fptT1MC;
//   PidTable *fptT2MC;
//   PidTable *fptMMC; 

//   // -- split into seagull and cowboys
//   PidTable *fptSgT1;
//   PidTable *fptSgT2;
//   PidTable *fptSgM; 

//   PidTable *fptSgT1MC;
//   PidTable *fptSgT2MC;
//   PidTable *fptSgMMC; 

//   PidTable *fptCbT1;
//   PidTable *fptCbT2;
//   PidTable *fptCbM; 

//   PidTable *fptCbT1MC;
//   PidTable *fptCbT2MC;
//   PidTable *fptCbMMC; 
  
//   PidTable *fptFakePosKaons, *fptFakePosPions, *fptFakePosProtons;
//   PidTable *fptFakeNegKaons, *fptFakeNegPions, *fptFakeNegProtons;

  

}



// ----------------------------------------------------------------------
void plotMisc::signalMass() {
  fF["SgMc"]->cd("candAnaMuMu"); 
  string mode("SgMc"); 
  TTree *t = getTree(mode); 
  if (0 == t) {
    cout << "plotBDT: no tree found for mode  " << mode << endl;
    return;
  } 
  setupTree(t, mode);
  loopOverTree(t, mode, 1); 

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
  // -- fit histograms: m0
  int ibin(0); 
  for (map<string, TH1D*>::iterator imap = fh0M0.begin(); imap != fh0M0.end(); ++imap) {  
    TH1D *h = imap->second;
    if (h->GetSumOfWeights() < 100) continue;
    sscanf(h->GetName(), "fh0M0_%d", &ibin); 
    f = fpFunc->pol1CrystalBall(h, 5.37, 0.04, 1., 1.); 
    h->Fit(f); 
    c0->SaveAs(Form("%s/misc-signalMass-%s.pdf", fDirectory.c_str(), imap->first.c_str())); 
    double mean  = f->GetParameter(0); 
    double sigma = f->GetParameter(1); 
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
  for (map<string, TH1D*>::iterator imap = fh0M3.begin(); imap != fh0M3.end(); ++imap) {  
    TH1D *h = imap->second;
    if (h->GetSumOfWeights() < 100) continue;
    sscanf(h->GetName(), "fh0M3_%d", &ibin); 
    f = fpFunc->pol1CrystalBall(h, 5.37, 0.04, 1., 1.); 
    h->Fit(f); 
    c0->SaveAs(Form("%s/misc-signalMass-%s.pdf", fDirectory.c_str(), imap->first.c_str())); 
    double mean  = f->GetParameter(0); 
    double sigma = f->GetParameter(1); 
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
  for (map<string, TH1D*>::iterator imap = fh0M4.begin(); imap != fh0M4.end(); ++imap) {  
    TH1D *h = imap->second;
    sscanf(h->GetName(), "fh0M4_%d", &ibin); 
    if (h->GetSumOfWeights() < 100) continue;
    f = fpFunc->pol1CrystalBall(h, 5.37, 0.04, 1., 1.); 
    h->Fit(f); 
    c0->SaveAs(Form("%s/misc-signaMass-%s.pdf", fDirectory.c_str(), imap->first.c_str())); 
    //    h->Fit("gaus");
    double mean  = f->GetParameter(0); 
    double sigma = f->GetParameter(1); 
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
