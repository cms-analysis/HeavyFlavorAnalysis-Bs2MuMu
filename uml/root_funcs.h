#include "TH1.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"

TObjArray Create_MassRes(TTree* tree, Double_t bdt) {
  TH2D* MassRes_hh = new TH2D("MassRes_hh", "MassRes_hh", 100, 4.9, 5.9, 40, -2.4, 2.4);
  Double_t m_t, eta_t, bdt_t;
  tree->SetBranchAddress("m",     &m_t);
  tree->SetBranchAddress("eta",   &eta_t);
  tree->SetBranchAddress("bdt",   &bdt_t);
  for (int j = 0; j < tree->GetEntries(); j++) {
    tree->GetEntry(j);
    if (bdt_t > bdt) MassRes_hh->Fill(m_t, eta_t);
  }
  TObjArray aSlices;
  MassRes_hh->FitSlicesX(0, 0, -1, 0, "Q", &aSlices);
  return aSlices;
}

TF1* Fit_MassRes(std::string file, Double_t bdt) {
  std::cout << "Mass resolution fit" << std::endl;
  TFile* MassRes_f = new TFile(file.c_str(), "UPDATE");
  TTree* MassRes_t = (TTree*)MassRes_f->Get("SgMc_bdt");
  TObjArray MassRes_toa = Create_MassRes(MassRes_t, bdt);
  TH1D* MassRes_h = (TH1D*)MassRes_toa[2];
  MassRes_h->Fit("pol6");
  TF1* MassRes_fit = MassRes_h->GetFunction("pol6");
  TCanvas c("c", "c", 600, 600);
  MassRes_h->SetYTitle("Resolution [GeV]");
  MassRes_h->SetXTitle("Candidate #eta");
  MassRes_h->SetStats(0);
  MassRes_h->Draw();
  c.Print("fig/width_fit.gif");
  c.Print("fig/width_fit.pdf");
  TF1* MassRes2_fit = (TF1*)MassRes_fit->Clone();
  return MassRes2_fit;
}

