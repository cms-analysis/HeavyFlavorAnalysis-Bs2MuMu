#include "TH1.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"

TObjArray Create_MassRes(TTree* tree) {
  TH2D* MassRes_hh = new TH2D("MassRes_hh", "MassRes_hh", 100, 4.9, 5.9, 50, -2.4, 2.4);
  Double_t m_t, eta_t;
  tree->SetBranchAddress("m",     &m_t);
  tree->SetBranchAddress("eta",   &eta_t);
  for (int j = 0; j < tree->GetEntries(); j++) {
    tree->GetEntry(j);
    MassRes_hh->Fill(m_t, eta_t);
  }
  TObjArray aSlices;
  MassRes_hh->FitSlicesX(0, 0, -1, 0, "Q", &aSlices);
  TCanvas c("c", "c", 600, 600);
  aSlices[2]->Draw();
//  c.Print("fig/slices.gif");exit(1);
  return aSlices;
}

void Fit_MassRes(std::string file, Double_t &p0, Double_t &p1, Double_t &p2) {
  std::cout << "Mass resolution fit" << std::endl;
  TFile* MassRes_f = new TFile(file.c_str(), "UPDATE");
  TTree* MassRes_t = (TTree*)MassRes_f->Get("SgMc_bdt");
  TObjArray MassRes_toa = Create_MassRes(MassRes_t);
  TH1D* MassRes_h = (TH1D*)MassRes_toa[2];
  MassRes_h->Fit("pol2");
  TF1* MassRes_fit = MassRes_h->GetFunction("pol2");
  p0 = MassRes_fit->GetParameter(0);
  p1 = MassRes_fit->GetParameter(1);
  p2 = MassRes_fit->GetParameter(2);
  return;
}

