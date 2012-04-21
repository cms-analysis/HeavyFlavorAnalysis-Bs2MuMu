#include "pdf_analysis.h"

pdf_analysis::pdf_analysis( RooWorkspace *ws, bool print, string source_name , string meth, string ch_s) {
  set_ws(ws);
  print_ = print;
  pdf_primitives_.resize(5);
  pdf_primitives_[0] = "bs";
  pdf_primitives_[1] = "bd";
  pdf_primitives_[2] = "peaking";
  pdf_primitives_[3] = "nonpeaking";
  pdf_primitives_[4] = "comb";

  source_name_ = source_name;
  pdf_name = "pdf_" + source_name_;
  rdh_name = "rdh_" + source_name_;
  meth_ = meth;
  ch_s_ = ch_s;
}

void pdf_analysis::fill_inputs(string input) {
  input_f = input;
  TFile* file = new TFile(input.c_str());
  string parse_sig;
  if (source_name_ == "bs") parse_sig = "Bs";
  if (source_name_ == "bd") parse_sig = "Bd";
  string histo_name = parse_sig + "_" + meth_ + "_chan" + ch_s_;
  cout << "input histo name: " << histo_name << endl;
  TH1D* h = (TH1D*)file->Get(histo_name.c_str());

  //// FIXMEEEEE ERRORS!!! <<<<<<<<<<<<<<<
  for (int i = 1; i < h->GetNbinsX(); i++) {
    if (h->GetBinCenter(i) < 4.90 || h->GetBinCenter(i) > 5.90) {
      h->SetBinContent(i, 0.);
      h->SetBinError(i, 0.);
    }
    if (source_name_ == "bs" || source_name_ == "bd") {
      h->SetBinError(i, h->GetBinContent(i)*0.15);
    }
    if ( source_name_ == "rare") {
      h->SetBinError(i, h->GetBinContent(i)*0.20);
    }
  }

  RooDataHist* rdh = new RooDataHist( rdh_name.c_str(), "dataset", *ws_->var("Mass"), h);
  ws_->import(*rdh);
}

void pdf_analysis::define_pdf (string pdf_name) {

  define_bs();
  define_bd();
  define_peaking();
  define_nonpeaking();
  define_comb();
}

void pdf_analysis::fit_pdf () {
  ws_->pdf( pdf_name.c_str())->fitTo(*ws_->data(rdh_name.c_str()), Extended(1), SumW2Error(1));
  if (print_) print();
}

void pdf_analysis::define_bs() {
  ws_->factory("NBs[20, 0, 200]");
  ws_->factory("BsMean[5.35, 5.32, 5.4]");
  ws_->factory("BsSigma[0.02, 0.005, 0.2]");
  ws_->factory("BsSigma2[0.04, 0.005, 0.2]");
  ws_->factory("BsAlpha[2.8, 0.1, 3.0]");
  ws_->factory("BsEnne[1., 0., 1.5]");
  ws_->factory("BscoeffGauss[0.5,0.,1.]");

  ws_->factory("Gaussian::bsGau(Mass, BsMean, BsSigma)");
  ws_->factory("CBShape::bsCB(Mass, BsMean, BsSigma2, BsAlpha, BsEnne)");
  ws_->factory("SUM::bsGauCB(BscoeffGauss*bsGau, bsCB)");
  RooExtendPdf* bsGauCBExt = new RooExtendPdf("pdf_bs", "bsGauCBExt", *ws_->pdf("bsGauCB"), *ws_->var("NBs"));
  ws_->import(*bsGauCBExt);

  return;
}

void pdf_analysis::define_bd() {

  ws_->factory("NBd[20, 0, 200]");
  ws_->factory("BdMean[5.25, 5.2, 5.29]");
  ws_->factory("BdSigma[0.1, 0.005, 0.2]");
  ws_->factory("BdSigma2[0.04, 0.005, 0.2]");
  ws_->factory("BdAlpha[0.5, 0., 20.]");
  ws_->factory("BdEnne[1., 0., 50.]");
  ws_->factory("BdcoeffGauss[0.5,0.,1.]");

  ws_->factory("Gaussian::bdGau(Mass, BdMean, BdSigma)");
  ws_->factory("CBShape::bdCB(Mass, BdMean, BdSigma2, BdAlpha, BdEnne)");
  ws_->factory("SUM::bdGauCB(BdcoeffGauss*bdGau, bdCB)");
  RooExtendPdf* bdGauCBExt = new RooExtendPdf("pdf_bd", "bdGauCBExt", *ws_->pdf("bdGauCB"), *ws_->var("NBd"));
  ws_->import(*bdGauCBExt);

  return;
}

void pdf_analysis::define_peaking() {

  ws_->factory("PeakingMean[5.25, 5.20, 5.3]");
  ws_->factory("PeakingSigma[0.050, 0.01, 0.20]");

  ws_->factory("Gaussian::PeakingGauss(Mass,PeakingMean,PeakingSigma)");

  return;
}

void pdf_analysis::define_nonpeaking() {
  ws_->factory("nonpeaking_m0[5.5]");
  ws_->factory("nonpeaking_c[-1., -10., 10]");
  ws_->factory("nonpeaking_p[0.5, -10., 10.]");

  ws_->factory("ArgusBG::NonpeakingArgus(Mass,nonpeaking_m0,nonpeaking_c,nonpeaking_p)");
  return;
}

void pdf_analysis::define_comb() {
  ws_->factory("Uniform::comb(Mass)");
  return;
}

void pdf_analysis::print() {
  RooPlot *rp = ws_->var("Mass")->frame();
  ws_->data((rdh_name).c_str())->plotOn(rp);
  ws_->pdf(pdf_name.c_str())->plotOn(rp);
  ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  string address = "fig/" + source_name_ + "_" + meth_ + "_" + ch_s_;
  c->Print( (address + ".gif").c_str());
  c->Print( (address + ".pdf").c_str());
  delete rp;
  delete c;
  return;
}
