#include "TH1.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"

TObjArray Create_MassRes(TTree* tree) {
  TH2D* MassRes_hh = new TH2D("MassRes_hh", "MassRes_hh", 100, 4.9, 5.9, 30, -2.4, 2.4);
  Double_t m_t, eta_t;
  tree->SetBranchAddress("m",     &m_t);
  tree->SetBranchAddress("eta",   &eta_t);
  for (int j = 0; j < tree->GetEntries(); j++) {
    tree->GetEntry(j);
    MassRes_hh->Fill(m_t, eta_t);
  }
  TObjArray aSlices;
  MassRes_hh->FitSlicesX(0, 0, -1, 0, "Q", &aSlices);
//  TCanvas c("c", "c", 600, 600);
//  aSlices[2]->Draw();
//  c.Print("fig/slices.gif");
  return aSlices;
}
