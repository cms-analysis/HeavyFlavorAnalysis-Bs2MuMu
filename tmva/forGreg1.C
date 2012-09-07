// ----------------------------------------------------------------------
void bla(double bdtcut = 0., double xmax, std::string var, std::string fname) {

  TFile *f = TFile::Open(Form("%s.root", fname.c_str()));

  TH1D *hs1 = new TH1D("hs1", var.c_str(), 100, 0., xmax); 
  TH1D *hs2 = new TH1D("hs2", var.c_str(), 100, 0., xmax); 
  TH1D *hsr = new TH1D("hsr", "Ratio", 100, 0., xmax); hsr->Sumw2();

  TH1D *hb1 = new TH1D("hb1", var.c_str(), 100, 0., xmax); 
  TH1D *hb2 = new TH1D("hb2", var.c_str(), 100, 0., xmax); 
  TH1D *hbr = new TH1D("hbr", "Ratio", 100, 0., xmax); hbr->Sumw2();

  TH2D *Hs = new TH2D("Hs", Form("SIGNAL BDT vs %s", var.c_str()), 100, 0., xmax, 100, -1., 1.);
  TH2D *Hb = new TH2D("Hb", Form("Background BDT vs %s", var.c_str()), 100, 0., xmax, 100, -1., 1.);

  //  TTree *t = (TTree*)f->Get("TrainTree"); 
  TTree *t = (TTree*)f->Get("TestTree"); 
  
  t->Draw(Form("%s>>hs1", var.c_str()), Form("BDT>%f&&classID==0", bdtcut));
  t->Draw(Form("%s>>hs2", var.c_str()), "classID==0");

  t->Draw(Form("%s>>hb1", var.c_str()), Form("BDT>%f&&classID==1", bdtcut));
  t->Draw(Form("%s>>hb2", var.c_str()), "classID==1");

  hsr->Divide(hs1, hs2, 1., 1., "");
  hbr->Divide(hb1, hb2, 1., 1., "");

  t->Draw(Form("BDT:%s>>Hs", var.c_str()), "classID==0");
  t->Draw(Form("BDT:%s>>Hb", var.c_str()), "classID==1");

  zone(2,3);
  
  c0.cd(1);
  hpl(hs1);
  hpl(hs2, "redesamenorm");

  tl->SetTextColor(kBlack);
  tl->DrawLatex(0.5, 0.5, Form("signal && BDT>%3.1f", bdtcut)); 
  tl->SetTextColor(kRed);
  tl->DrawLatex(0.5, 0.45, Form("signal (same area)", bdtcut)); 

  c0.cd(2);
  hpl(hb1);
  hpl(hb2, "redesamenorm");
  tl->SetTextColor(kBlack);
  tl->DrawLatex(0.5, 0.5, Form("bg && BDT>%3.1f", bdtcut)); 
  tl->SetTextColor(kRed);
  tl->DrawLatex(0.5, 0.45, Form("bg (same area)", bdtcut)); 

  c0.cd(3);
  hpl(hsr);

  c0.cd(4);
  hpl(hbr);

  c0.cd(5); 
  gPad->SetLogz(0);
  Hs->DrawCopy("colz");
//   c0.cd(6);
//   gPad->SetLogz(1);
//   Hs->DrawCopy("colz");

  c0.cd(6);
  gPad->SetLogz(0);
  Hb->Draw("colz");
//   c0.cd(4);
//   gPad->SetLogz(1);
//   Hb->Draw("colz");
 
  c0.SaveAs(Form("%s-%s-test.pdf", fname.c_str(), var.c_str()));


}
