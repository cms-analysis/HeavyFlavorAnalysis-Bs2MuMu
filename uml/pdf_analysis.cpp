#include "pdf_analysis.h"

pdf_analysis::pdf_analysis(bool print, string meth, string ch_s, string range) {
  ws_ = new RooWorkspace("ws", "ws");
  print_ = print;
  meth_ = meth;
  ch_s_ = ch_s;
  range_ = range;
  RooRealVar* Mass = new RooRealVar("Mass", "Candidate invariant mass", 4.90, 5.90, "GeV/c^{2}");
  ws_->import(*Mass);
  ws_->var("Mass")->setRange("sb_lo", 4.90, 5.20);  
  ws_->var("Mass")->setRange("blind", 5.20, 5.45);
  ws_->var("Mass")->setRange("sb_hi", 5.45, 5.90);
  ws_->var("Mass")->setRange("overall", 4.90, 5.90);
}

void pdf_analysis::define_pdfs () {
  
  define_bs();
  define_bd();
  define_peaking();
  define_nonpeaking();
  define_comb();
  
  
  define_signals();
  define_rare();
  define_bkg();
  define_signalsrare();
  
  define_all();
  define_total();
}

void pdf_analysis::fit_pdf (string pdf, RooAbsData* data, bool extended) {
  rds_ = data;
  pdf_name = "pdf_" + pdf;
  if (extended) pdf_name = "pdf_ext_" + pdf;
  rdh_name = rds_->GetName();
  cout << "fitting " << rdh_name << " in range " << range_ << " with " << pdf_name << ":" << endl;
  ws_->pdf( pdf_name.c_str())->Print();
  
  ws_->pdf( pdf_name.c_str())->fitTo(*rds_, Extended(extended), SumW2Error(1), Range(range_.c_str()));
  if (print_) print();
  set_pdf_constant(pdf_name);
}

void pdf_analysis::set_pdf_constant(string name) {
  RooArgSet * set = ws_->pdf(name.c_str())->getVariables();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (!(name == "Mass")) {
      size_t found;
      found = name.find("N");
      if (found == string::npos) ws_->var(var_Obj->GetName())->setConstant(1);
    }
  }
}

void pdf_analysis::define_bs() {
  ws_->factory("N_bs[10, 0, 100]");
  ws_->factory("Mean_bs[5.35, 5.32, 5.4]");
  ws_->factory("Sigma_bs[0.02, 0.005, 0.2]");
  ws_->factory("Sigma2_bs[0.04, 0.005, 0.2]");
  ws_->factory("Alpha_bs[2.8, 0.1, 3.0]");
  ws_->factory("Enne_bs[1., 0., 1.5]");
  ws_->factory("CoeffGauss_bs[0.5,0.,1.]");

  ws_->factory("Gaussian::Gau_bs(Mass, Mean_bs, Sigma_bs)");
  ws_->factory("CBShape::CB_bs(Mass, Mean_bs, Sigma2_bs, Alpha_bs, Enne_bs)");
  ws_->factory("SUM::pdf_bs(CoeffGauss_bs*Gau_bs, CB_bs)");
  
  RooExtendPdf* bsGauCBExt = new RooExtendPdf("pdf_ext_bs", "GauCBExt_bs", *ws_->pdf("pdf_bs"), *ws_->var("N_bs"), range_.c_str());
  ws_->import(*bsGauCBExt);

  return;
}

void pdf_analysis::define_bd() {

  ws_->factory("N_bd[10, 0, 100]");
  ws_->factory("Mean_bd[5.25, 5.20, 5.29]");
  ws_->factory("Sigma_bd[0.02, 0.005, 0.2]");
  ws_->factory("Sigma2_bd[0.04, 0.005, 0.2]");
  ws_->factory("Alpha_bd[2.8, 0.1, 3.0]");
  ws_->factory("Enne_bd[1., 0., 1.5]");
  ws_->factory("CoeffGauss_bd[0.5,0.,1.]");

  ws_->factory("Gaussian::Gau_bd(Mass, Mean_bd, Sigma_bd)");
  ws_->factory("CBShape::CB_bd(Mass, Mean_bd, Sigma2_bd, Alpha_bd, Enne_bd)");
  ws_->factory("SUM::pdf_bd(CoeffGauss_bd*Gau_bd, CB_bd)");
  
  RooExtendPdf* bdGauCBExt = new RooExtendPdf("pdf_ext_bd", "GauCBExt_bd", *ws_->pdf("pdf_bd"), *ws_->var("N_bd"), range_.c_str());
  ws_->import(*bdGauCBExt);

  return;
}

void pdf_analysis::define_peaking() {
  
  ws_->factory("N_peaking[10,0,100]");
  ws_->factory("Mean_peaking[5.25, 5.20, 5.3]");
  ws_->factory("Sigma_peaking[0.050, 0.01, 0.20]");
  
  ws_->factory("Gaussian::pdf_peaking(Mass,Mean_peaking,Sigma_peaking)");
  
  RooExtendPdf* peakingGaussExt = new RooExtendPdf("pdf_ext_peaking", "peakingGaussExt", *ws_->pdf("pdf_peaking"), *ws_->var("N_peaking"), range_.c_str());
  ws_->import(*peakingGaussExt);
  
  return;
}

void pdf_analysis::define_nonpeaking() {
  
  ws_->factory("N_nonpeaking[10, 0, 100]"); 
  ws_->factory("m0_nonpeaking[5.5]");
  ws_->factory("c_nonpeaking[-1., -10., 20]");
  ws_->factory("p_nonpeaking[0.5, -10., 10.]");

  ws_->factory("ArgusBG::pdf_nonpeaking(Mass,m0_nonpeaking,c_nonpeaking,p_nonpeaking)");
  
  RooExtendPdf* nonpeakingGaussExt = new RooExtendPdf("pdf_ext_nonpeaking", "nonpeakingGaussExt", *ws_->pdf("pdf_nonpeaking"), *ws_->var("N_nonpeaking"), range_.c_str());
  ws_->import(*nonpeakingGaussExt);
  return;
}

void pdf_analysis::define_comb() {
  
  ws_->factory("N_comb[10, 0, 100]");
  
  ws_->factory("Uniform::pdf_comb(Mass)");
  
  RooExtendPdf* combExt = new RooExtendPdf("pdf_ext_comb", "combExt", *ws_->pdf("pdf_comb"), *ws_->var("N_comb"), range_.c_str());
  ws_->import(*combExt);
  return;
}

void pdf_analysis::define_signals() {
  
  ws_->factory("N_signals[10, 0, 100]");
  ws_->factory("bsfraction_signals[0.5, 0.0, 1.0]"); 
  
  ws_->factory("SUM::pdf_signals(bsfraction_signals*pdf_bs, pdf_bd)");
  
  RooExtendPdf* signalsExt = new RooExtendPdf("pdf_ext_signals", "signalsExt", *ws_->pdf("pdf_signals"), *ws_->var("N_signals"), range_.c_str());
  ws_->import(*signalsExt);
  return; 
}

void pdf_analysis::define_rare() {
  
  ws_->factory("N_rare[10, 0, 100]");
  ws_->factory("peakingfraction_rare[0.5, 0.0, 1.0]"); 
  
  ws_->factory("SUM::pdf_rare(peakingfraction_rare*pdf_peaking, pdf_nonpeaking)");
  
  RooExtendPdf* rareExt = new RooExtendPdf("pdf_ext_rare", "rareExt", *ws_->pdf("pdf_rare"), *ws_->var("N_rare"), range_.c_str());
  ws_->import(*rareExt);
  return; 
}

void pdf_analysis::define_signalsrare() {
  
  ws_->factory("signalfraction_signalsrare[0.5, 0.0, 1.0]"); 
  ws_->factory("SUM::pdf_signalsrare(signalfraction_signalsrare*pdf_signals, pdf_rare)");
  return;
}

void pdf_analysis::define_bkg() {
  
  ws_->factory("N_bkg[10, 0, 100]");
  ws_->factory("rarefraction_bkg[0.5, 0.0, 1.0]"); 
  
  ws_->factory("SUM::pdf_bkg(rarefraction_bkg*pdf_rare, pdf_comb)");
  
  RooExtendPdf* bkgExt = new RooExtendPdf("pdf_ext_bkg", "bkgExt", *ws_->pdf("pdf_bkg"), *ws_->var("N_bkg"), range_.c_str());
  ws_->import(*bkgExt);
  return; 
}

void pdf_analysis::define_all() {
  
  ws_->factory("N_all[10, 0, 100]");
  ws_->factory("signalsfraction_all[0.5, 0.0, 1.0]"); 
  
  ws_->factory("SUM::pdf_all(signalsfraction_all*pdf_signals, pdf_bkg)");
  
  RooExtendPdf* allExt = new RooExtendPdf("pdf_ext_all", "allExt", *ws_->pdf("pdf_all"), *ws_->var("N_all"), range_.c_str());
  ws_->import(*allExt);
  return; 
}

void pdf_analysis::define_total() {
  
  
  ws_->factory("SUM::pdf_ext_total(N_bs*pdf_bs, N_bd*pdf_bd, N_rare*pdf_rare, N_comb*pdf_comb)");
  
  return; 
}

void pdf_analysis::print() {
  RooPlot *rp = ws_->var("Mass")->frame();
  rds_->plotOn(rp, Binning(20));
  ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kRed), Range(range_.c_str()));
  ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
  
  RooArgSet * set = ws_->pdf(pdf_name.c_str())->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  int i = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (!(name == pdf_name)) {
      i++;  
      ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kRed + 50*i),  LineStyle(1), Range(range_.c_str()));
    }
  }
  
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  string address = "fig/" + pdf_name + "_" + meth_ + "_" + ch_s_;
  c->Print( (address + ".gif").c_str());
  //c->Print( (address + ".pdf").c_str());
  delete rp;
  delete c;
  return;
}
