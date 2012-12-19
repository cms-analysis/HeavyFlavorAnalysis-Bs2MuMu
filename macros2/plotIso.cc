#include "plotIso.hh"

#include <algorithm>

#include "../macros/AnalysisDistribution.hh"
#include "../macros/HistCutEfficiency.hh"
#include "TMath.h"

using namespace std; 
using std::string; 

ClassImp(plotIso)

// ----------------------------------------------------------------------
plotIso::plotIso(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

  string hfname  = fDirectory + "/anaBmm.plotIso." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  //  if (fHistFile) fHistFile->Close();
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

  fNumbersFileName = fDirectory + "/anaBmm.plotIso." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
plotIso::~plotIso() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotIso::makeAll(double isoCut) {
  
  fIsoCut = isoCut; 
  fMode = 3; 
  fPreco = 5.1;
  dataVsMc("NoData", "candAnaBu2JpsiK",   "NoMc", "candAnaBu2JpsiK",   "Ao"); 

  fMode = 2; 
  dataVsMc("CsData", "candAnaBs2JpsiPhi", "CsMc", "candAnaBs2JpsiPhi", "Ao"); 

  fMode = 0; 
  dataVsMc("SgData", "candAnaMuMu", "SgMc", "candAnaMuMu", "Ao"); 

  anaTextFiles(static_cast<int>(100*isoCut)); 

}

// ----------------------------------------------------------------------
bool largerRatio(const data a, const data b) {
  return (a.ratio > b.ratio); 
}

// ----------------------------------------------------------------------
bool closerToOne(const data a, const data b) {
  return (TMath::Abs(a.ratio - 1) < TMath::Abs(b.ratio - 1)); 
}

// ----------------------------------------------------------------------
void plotIso::anaTextFiles(int cutval) {

  vector<data> mm, no, cs; 
  
  readFile(Form("plotIso-%d-SgData-SgMc.txt", cutval), mm);
  readFile(Form("plotIso-%d-NoData-NoMc.txt", cutval), no);
  readFile(Form("plotIso-%d-CsData-CsMc.txt", cutval), cs);

  sort(mm.begin(), mm.end(), largerRatio); 
  double maxDev(1.03); 

  TH2D *mr3 = new TH2D("mr3", "MM MC/data, doca = 0.3 mm", 6, 0., 6., 5, 0., 5.); relabel(mr3);
  TH2D *mr4 = new TH2D("mr4", "MM MC/data, doca = 0.4 mm", 6, 0., 6., 5, 0., 5.); relabel(mr4);
  TH2D *mr5 = new TH2D("mr5", "MM MC/data, doca = 0.5 mm", 6, 0., 6., 5, 0., 5.); relabel(mr5);

  TH2D *mm3 = new TH2D("mm3", "MM data, doca = 0.3 mm", 6, 0., 6., 5, 0., 5.); relabel(mm3);
  TH2D *mm4 = new TH2D("mm4", "MM data, doca = 0.4 mm", 6, 0., 6., 5, 0., 5.); relabel(mm4);
  TH2D *mm5 = new TH2D("mm5", "MM data, doca = 0.5 mm", 6, 0., 6., 5, 0., 5.); relabel(mm5);

  TH2D *no3 = new TH2D("no3", "NO MC/data, doca = 0.3 mm", 6, 0., 6., 5, 0., 5.); relabel(no3); no3->SetMinimum(1.0); no3->SetMaximum(1.05);
  TH2D *no4 = new TH2D("no4", "NO MC/data, doca = 0.4 mm", 6, 0., 6., 5, 0., 5.); relabel(no4); no4->SetMinimum(1.0); no4->SetMaximum(1.05);
  TH2D *no5 = new TH2D("no5", "NO MC/data, doca = 0.5 mm", 6, 0., 6., 5, 0., 5.); relabel(no5); no5->SetMinimum(1.0); no5->SetMaximum(1.05);

  TH2D *cs3 = new TH2D("cs3", "CS MC/data, doca = 0.3 mm", 6, 0., 6., 5, 0., 5.); relabel(cs3); cs3->SetMinimum(1.0); cs3->SetMaximum(1.05);
  TH2D *cs4 = new TH2D("cs4", "CS MC/data, doca = 0.4 mm", 6, 0., 6., 5, 0., 5.); relabel(cs4); cs4->SetMinimum(1.0); cs4->SetMaximum(1.05);
  TH2D *cs5 = new TH2D("cs5", "CS MC/data, doca = 0.5 mm", 6, 0., 6., 5, 0., 5.); relabel(cs5); cs5->SetMinimum(1.0); cs5->SetMaximum(1.05);
  double doca, pt, r, x; 
  int ipt, ir; 
  for (unsigned int i = 0; i < mm.size(); ++i) {
    doca = mm[i].doca; 
    pt   = mm[i].pt; 
    r    = mm[i].r;
    ipt = ptBin(pt); 
    ir  = rBin(r); 
    x    = mm[i].eff1;
    if (doca > 0.025 && doca < 0.035)  mm3->SetBinContent(ir, ipt, x); 
    if (doca > 0.035 && doca < 0.045)  mm4->SetBinContent(ir, ipt, x); 
    if (doca > 0.045 && doca < 0.055)  mm5->SetBinContent(ir, ipt, x); 
    x    = mm[i].ratio;
    if (doca > 0.025 && doca < 0.035)  mr3->SetBinContent(ir, ipt, x); 
    if (doca > 0.035 && doca < 0.045)  mr4->SetBinContent(ir, ipt, x); 
    if (doca > 0.045 && doca < 0.055)  mr5->SetBinContent(ir, ipt, x); 
  }

  for (unsigned int i = 0; i < no.size(); ++i) {
    doca = no[i].doca; 
    pt   = no[i].pt; 
    r    = no[i].r;
    x    = no[i].ratio;

    ipt = ptBin(pt); 
    ir  = rBin(r); 
    if (doca > 0.025 && doca < 0.035)  no3->SetBinContent(ir, ipt, x); 
    if (doca > 0.035 && doca < 0.045)  no4->SetBinContent(ir, ipt, x); 
    if (doca > 0.045 && doca < 0.055)  no5->SetBinContent(ir, ipt, x); 
  }

  for (unsigned int i = 0; i < cs.size(); ++i) {
    doca = cs[i].doca; 
    pt   = cs[i].pt; 
    r    = cs[i].r;
    x    = cs[i].ratio;

    ipt = ptBin(pt); 
    ir  = rBin(r); 
    if (doca > 0.025 && doca < 0.035)  cs3->SetBinContent(ir, ipt, x); 
    if (doca > 0.035 && doca < 0.045)  cs4->SetBinContent(ir, ipt, x); 
    if (doca > 0.045 && doca < 0.055)  cs5->SetBinContent(ir, ipt, x); 
  }

  gStyle->SetOptStat(0); 
  makeCanvas(1); 
  c1->cd(); 
  c1->Clear();
  c1->Divide(4,1);
  
  c1->cd(1); 
  gPad->SetRightMargin(0.15); 
  mm3->Draw("colztext"); 

  c1->cd(2); 
  gPad->SetRightMargin(0.15); 
  mr3->Draw("colztext"); 

  c1->cd(3); 
  gPad->SetRightMargin(0.15); 
  no3->Draw("colztext"); 

  c1->cd(4); 
  gPad->SetRightMargin(0.15); 
  cs3->Draw("colztext"); 
  c1->SaveAs(Form("doca3-%d.pdf", cutval)); 


  c1->cd(); 
  c1->Clear();
  c1->Divide(4,1);
  
  c1->cd(1); 
  gPad->SetRightMargin(0.15); 
  mm4->Draw("colztext"); 

  c1->cd(2); 
  gPad->SetRightMargin(0.15); 
  mr4->Draw("colztext"); 

  c1->cd(3); 
  gPad->SetRightMargin(0.15); 
  no4->Draw("colztext"); 

  c1->cd(4); 
  gPad->SetRightMargin(0.15); 
  cs4->Draw("colztext"); 
  c1->SaveAs(Form("doca4-%d.pdf", cutval)); 

  c1->cd(); 
  c1->Clear();
  c1->Divide(4,1);
  
  c1->cd(1); 
  gPad->SetRightMargin(0.15); 
  mm5->Draw("colztext"); 

  c1->cd(2); 
  gPad->SetRightMargin(0.15); 
  mr5->Draw("colztext"); 

  c1->cd(3); 
  gPad->SetRightMargin(0.15); 
  no5->Draw("colztext"); 

  c1->cd(4); 
  gPad->SetRightMargin(0.15); 
  cs5->Draw("colztext"); 
  c1->SaveAs(Form("doca5-%d.pdf", cutval)); 
  
  for (unsigned int i = 0; i < mm.size(); ++i) {
    for (unsigned int j = 0; j < no.size(); ++j) {
      for (unsigned int k = 0; k < cs.size(); ++k) {
	if (no[j].name == mm[i].name && (cs[k].name == mm[i].name)) {
	  if (no[j].ratio < maxDev &&cs[k].ratio < maxDev) 
	    cout << "MM: " << mm[i].name << ": " << mm[i].ratio << " r = " << mm[i].r << " pt = " << mm[i].pt << " doca = " << mm[i].doca 
		 << " NO: " << no[j].name << " ratio = " << no[j].ratio 
		 << " CS: " << cs[k].name << " ratio = " << cs[k].ratio 
		 << endl;
	}
      }
    }
  }

  for (unsigned int i = 0; i < mm.size(); ++i) {
    for (unsigned int j = 0; j < no.size(); ++j) {
      for (unsigned int k = 0; k < cs.size(); ++k) {
	if ((mm[i].name == "Doca05isor010pt09") && (no[j].name == mm[i].name) && (cs[k].name == mm[i].name)) {
	  cout << "*MM: " << mm[i].name << ": " << mm[i].ratio << " r = " << mm[i].r << " pt = " << mm[i].pt << " doca = " << mm[i].doca 
	       << " NO: " << no[j].name << " ratio = " << no[j].ratio 
	       << " CS: " << cs[k].name << " ratio = " << cs[k].ratio 
	       << endl;
	}
      }
    }
  }
}


// ----------------------------------------------------------------------
void plotIso::readFile(string fname, vector<data> &v) {

  int i0, i1, i2; 
  float f0, f1, f2; 
  char buffer[1000]; 
  ifstream is(fname.c_str()); 
  while (is.getline(buffer, 1000, '\n')) {
    sscanf(buffer, "Doca%disor%dpt%d: %f %f %f", &i0, &i1, &i2, &f0, &f1, &f2); 
    //    cout << "doca: " << i0*0.01 << " r = " << i1*0.1 << " pt = " << i2*0.1 << endl;
    data a; 
    a.name = Form("Doca0%disor0%dpt0%d", i0, i1, i2);
    a.doca = i0*0.01; 
    a.r    = i1*0.1; 
    a.pt   = i2*0.1;
    a.eff1 = f0; 
    a.eff2 = f1; 
    a.ratio= f2; 
    v.push_back(a); 
  }
  is.close(); 

}


// ----------------------------------------------------------------------
void plotIso::dataVsMc(string file1, string dir1, string file2, string dir2, string selection) {

  double cutval = fIsoCut; 
  ofstream OUT(Form("plotIso-%d-%s-%s.txt", static_cast<int>(cutval*100), file1.c_str(), file2.c_str()));

  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 

  char option1[100], option2[100],  loption1[100], loption2[100]; 
  int color1(1), fill1(1000), color2(1), fill2(1000), marker1(20), marker2(25); 

  if (string::npos != file1.find("No")) {
    color1 = kBlack; 
    fill1  = 3004; 
    fill1  = 3356; 
  }
  if (string::npos != file2.find("No")) {
    color2 = kBlue; 
    fill2  = 3004; 
    fill2  = 3356; 
  }
  if (string::npos != file1.find("Cs")) {
    color1 = kBlack; 
    fill1  = 3005; 
    fill1  = 3356; 
  }
  if (string::npos != file2.find("Cs")) {
    color2 = kRed; 
    fill2  = 3005; 
    fill2  = 3356; 
  }

  if (string::npos != file1.find("Sg")) {
    // -- data is first, then signal MC
    color1 = kBlack; 
    fill1  = 0; 
    sprintf(option1, "e"); 
    sprintf(loption1, "p");
    color2 = kBlue; 
    fill2  = 3005; 
    fill2  = 3356; 
    sprintf(option2, "hist"); 
    sprintf(loption2, "f");
  }
  
  // -- now fix overlays of two data files
  sprintf(option1, "e");
  sprintf(loption1, "p");
  
  sprintf(option2,  "hist");
  sprintf(loption2, "f");
  
  if (string::npos != file2.find("Data")) {
    sprintf(option1, "hist"); 
    sprintf(loption1, "f");
    sprintf(option2, "e"); 
    sprintf(loption2, "p");
    fill2   = 0; 
    color2  = kBlack; 
    marker2 = 21; 
  }
  
  if (string::npos != file1.find("Mc")) {
    sprintf(option1, "hist"); 
    sprintf(loption1, "f");
    sprintf(loption2, "p");
    fill2   = 0; 
    color2  = kBlack; 
    marker2 = 25; 
  }
 
  //  ofstream OUT("testUL.txt", ios::app);

  fF[file1]->cd(Form("%s", dir1.c_str())); 
  TH1D *h1(0), *h2(0);

  AnalysisDistribution a("A_pvz"); 

  fF[file1]->cd(Form("%s/%s", dir1.c_str(), "Isolation")); 

  TCanvas *c1;
  string cut, pdfname; 
  vector<string> doList; 
  TString hname; 

  // -- Suck in all AD's
  TObject *key;
  TIter next(gDirectory->GetListOfKeys());                           
  while ((key = (TObject*)next())) {
    hname = key->GetName();
    if (!hname.Contains("MassSi")) continue;
    hname.ReplaceAll("MassSi", ""); 
    doList.push_back(hname.Data()); 
  }



  for (unsigned int i = 0; i < doList.size(); ++i) {
    cut = Form("%s", doList[i].c_str()); 

    pdfname = Form("%s/%s_%s-%s_sbs-%d_%s_%s.pdf", fDirectory.c_str(), fSuffix.c_str(),
		   file1.c_str(), file2.c_str(), static_cast<int>(100*cutval), cut.c_str(), selection.c_str());
    
    cout << pdfname << endl;
    
    fF[file1]->cd(Form("%s/%s", dir1.c_str(), "Isolation")); 
    cout << "==> File1: "; gDirectory->pwd(); 
    cout << "==> pdf: " << pdfname << endl;
    // -- separate mm from rest because the former don't need sbs
    if (0 == fMode) {
      h1 = (TH1D*)gDirectory->Get(Form("%s%s1", cut.c_str(), selection.c_str()));
    } else {
      if (string::npos != file1.find("Mc")) {
	// -- For MC use the signal distributions directly
	h1 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
      } else {
	if (1 == fMode) h1 = a.sbsDistribution(cut.c_str(), selection.c_str());
	if (2 == fMode) h1 = a.sbsDistributionExpoGauss(cut.c_str(), selection.c_str());
	if (3 == fMode) h1 = a.sbsDistributionPol1ErrGauss(cut.c_str(), selection.c_str(), fPreco);
      }
    }

    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint){
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file1.c_str(), cut.c_str(), selection.c_str()));
    }

    fF[file2]->cd(Form("%s/%s", dir2.c_str(), "Isolation")); 
    cout << "==> File2: pwd() = "; gDirectory->pwd(); 
    if (0 == fMode) {
      h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
    } else {
      if (string::npos != file1.find("Mc")) {
	// -- For MC use the signal distributions directly
	h2 = (TH1D*)gDirectory->Get(Form("%s%s0", cut.c_str(), selection.c_str()));
      } else {
	if (1 == fMode) h2 = a.sbsDistribution(cut.c_str(), selection.c_str());
	if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection.c_str());
	if (3 == fMode) h2 = a.sbsDistributionPol1ErrGauss(cut.c_str(), selection.c_str(), fPreco);
      }
      //      if (2 == fMode) h2 = a.sbsDistributionExpoGauss(cut.c_str(), selection);
    }
    c1 = (TCanvas*)gROOT->FindObject("c1"); 
    if (fDoPrint) {
      if (c1) c1->SaveAs(Form("%s/tmp/%s_sbs_%s_%s.pdf", fDirectory.c_str(), file2.c_str(), cut.c_str(), selection.c_str()));
    }

    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    c0->cd(); 
    c0->Clear();
    double dmax = (h1->GetMaximum() > h2->GetMaximum()? 1.1*h1->GetMaximum(): 1.1*h2->GetMaximum()); 

    h1->SetMinimum(0.1);
    h1->SetMaximum(dmax); 
    gPad->SetLogy(0); 

    setHist(h1, color1, marker1, 1.5); 
    setFilledHist(h1, color1, color1, fill1); 
    h1->SetTitle("");
    h1->Draw(option1);
    
    setHist(h2, color2, marker2, 1.5); 
    setFilledHist(h2, color2, color2, fill2); 
    h2->Draw(Form("same%s", option2));

    if (string::npos != cut.find("iso")) {

      HistCutEfficiency h(h1); 
      //      h.fVerbose = 1; 
      h.eff(h1, cutval); 

      double eps1   = h.hiEff; 

      h.eff(h2, cutval); 
      double eps2   = h.hiEff; 
      
      tl->DrawLatex(0.2, 0.92, cut.c_str());
      tl->DrawLatex(0.2, 0.80, Form("#epsilon_{D}  = %4.3f", eps1));
      tl->DrawLatex(0.2, 0.72, Form("#epsilon_{M}  = %4.3f", eps2));
      tl->DrawLatex(0.2, 0.64, Form("r_{M/D} = %4.3f", eps2/eps1));
      
      OUT << Form("%s: %4.3f %4.3f %4.3f", cut.c_str(), eps1, eps2, eps2/eps1) << endl;
    }

    if (fDoPrint) c0->SaveAs(pdfname.c_str()); 
    //    break;
  }

} 

// ----------------------------------------------------------------------
int plotIso::ptBin(double pt) {

  if (pt > 0.25 && pt < 0.35) return 1; 
  if (pt > 0.45 && pt < 0.55) return 2; 
  if (pt > 0.65 && pt < 0.75) return 3; 
  if (pt > 0.85 && pt < 0.95) return 4; 
  if (pt > 1.05 && pt < 1.15) return 5; 
  cout << "don't know what to do with pt = " << pt << endl;
  return -1;
}


// ----------------------------------------------------------------------
int plotIso::rBin(double r) {
  if (r > 0.25 && r < 0.35) return 1; 
  if (r > 0.45 && r < 0.55) return 2; 
  if (r > 0.65 && r < 0.75) return 3; 
  if (r > 0.85 && r < 0.95) return 4; 
  if (r > 0.95 && r < 1.05) return 5; 
  if (r > 1.05 && r < 1.15) return 6; 
  cout << "don't know what to do with r = " << r << endl;
  return -1;
}


// ----------------------------------------------------------------------
void plotIso::relabel(TH2 *h) {
  h->GetYaxis()->SetBinLabel(1, "pT<0.3"); 
  h->GetYaxis()->SetBinLabel(2, "pT<0.5"); 
  h->GetYaxis()->SetBinLabel(3, "pT<0.7"); 
  h->GetYaxis()->SetBinLabel(4, "pT<0.9"); 
  h->GetYaxis()->SetBinLabel(5, "pT<1.1"); 

  h->GetXaxis()->SetBinLabel(1, "r<0.3"); 
  h->GetXaxis()->SetBinLabel(2, "r<0.5"); 
  h->GetXaxis()->SetBinLabel(3, "r<0.7"); 
  h->GetXaxis()->SetBinLabel(4, "r<0.9"); 
  h->GetXaxis()->SetBinLabel(5, "r<1.0"); 
  h->GetXaxis()->SetBinLabel(6, "r<1.1"); 
  
  h->SetMarkerSize(2);
}
