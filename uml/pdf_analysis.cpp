#include "pdf_analysis.h"

pdf_analysis::pdf_analysis(bool print, string meth, string ch_s, string range, bool SM, bool bd_constr) {
  cout << "analysis constructor" << endl;
  ws_ = new RooWorkspace("ws", "ws");
  print_ = print;
  meth_ = meth;
  ch_s_ = ch_s;
  range_ = range;
  SM_ = SM;
  bd_constr_ = bd_constr;
  channels = 0;  
}

void pdf_analysis::define_pdfs () {
  
  Mass = new RooRealVar("Mass", "Candidate invariant mass", 4.90, 5.90, "GeV/c^{2}");
  ws_->import(*Mass);
  ws_->var("Mass")->setRange("sb_lo", 4.90, 5.20);
  ws_->var("Mass")->setRange("blind", 5.20, 5.45);
  ws_->var("Mass")->setRange("sb_hi", 5.45, 5.90);
  ws_->var("Mass")->setRange("overall", 4.90, 5.90);

  define_bs();
  define_bd();
  define_peaking();
  define_nonpeaking();
  define_comb();
  
  define_signals();
  define_rare();
  define_bkg_fractional();
  define_bkg_extended();
  define_signalsrare();
  //define_bscomb();

  define_total_fractional();
  define_total_extended();

}

void pdf_analysis::fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error, bool hesse) {
  rds_ = data;
  pdf_name = "pdf_" + pdf;
  if (extended) pdf_name = "pdf_ext_" + pdf;
  string rdh_name = rds_->GetName();
  cout << "fitting " << rdh_name << " in range " << range_ << " with " << pdf_name << ":" << endl;
  ws_->pdf( pdf_name.c_str())->Print();
  
  RFR = ws_->pdf( pdf_name.c_str())->fitTo(*rds_, Extended(extended), SumW2Error(sumw2error), Range(range_.c_str())/*, SumCoefRange(range_.c_str())*/, NumCPU(2), Hesse(hesse));
  if (print_) print();
  set_pdf_constant(pdf_name);
}

void pdf_analysis::set_pdf_constant(string name) {
  RooArgSet * set = ws_->pdf(name.c_str())->getVariables();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (!(name == "Mass") && !(name == "Bd_over_Bs")) {
      size_t found;
      found = name.find("N");
      if (found == string::npos) ws_->var(var_Obj->GetName())->setConstant(1);
    }
  }
}

void pdf_analysis::define_bs() {
  ws_->factory("N_bs[0, 100]");
  ws_->factory("Mean_bs[5.35, 5.32, 5.4]");
  ws_->factory("Sigma_bs[0.02, 0.005, 0.2]");
  ws_->factory("Sigma2_bs[0.04, 0.005, 0.2]");
  ws_->factory("Alpha_bs[2.8, 0.1, 3.0]");
  ws_->factory("Enne_bs[1., 0., 1.5]");
  ws_->factory("CoeffGauss_bs[0.5,0.,1.]");
  ws_->factory("Gaussian::Gau_bs(Mass, Mean_bs, Sigma_bs)");
  ws_->factory("CBShape::CB_bs(Mass, Mean_bs, Sigma2_bs, Alpha_bs, Enne_bs)");

  ws_->factory("SUM::pdf_bs(CoeffGauss_bs*Gau_bs, CB_bs)");
  
}

void pdf_analysis::define_bd() {

  ws_->factory("N_bd[0, 100]");
  ws_->factory("Mean_bd[5.25, 5.20, 5.29]");
  ws_->factory("Sigma_bd[0.02, 0.005, 0.2]");
  ws_->factory("Sigma2_bd[0.04, 0.005, 0.2]");
  ws_->factory("Alpha_bd[2.8, 0.1, 3.0]");
  ws_->factory("Enne_bd[1., 0., 1.5]");
  ws_->factory("CoeffGauss_bd[0.5,0.,1.]");
  ws_->factory("Gaussian::Gau_bd(Mass, Mean_bd, Sigma_bd)");
  ws_->factory("CBShape::CB_bd(Mass, Mean_bd, Sigma2_bd, Alpha_bd, Enne_bd)");

  ws_->factory("SUM::pdf_bd(CoeffGauss_bd*Gau_bd, CB_bd)");
  
  if (SM_) {
    RooConstVar *SM_Bd_over_Bs = new RooConstVar("SM_Bd_over_Bs", "SM_Bd_over_Bs", ratio_);
    RooFormulaVar *N_bd_constr = new RooFormulaVar("N_bd_constr", "@0*@1", RooArgList(*ws_->var("N_bs"), *SM_Bd_over_Bs));
    ws_->import(*N_bd_constr);
  }
  if (bd_constr_) {
    RooRealVar * Bd_over_Bs = new RooRealVar("Bd_over_Bs", "Bd_over_Bs", 0., 2., "");
    RooFormulaVar *N_bd_constr = new RooFormulaVar("N_bd_constr", "@0*@1", RooArgList(*ws_->var("N_bs"), *Bd_over_Bs));
    ws_->import(*N_bd_constr);
  }
}

void pdf_analysis::define_peaking() {
  
  ws_->factory("N_peaking[0, 100]");
  ws_->factory("Mean_peaking[5.25, 5.20, 5.3]");
  ws_->factory("Sigma_peaking[0.050, 0.01, 0.20]");
  
  ws_->factory("Gaussian::pdf_peaking(Mass,Mean_peaking,Sigma_peaking)");
  
}

void pdf_analysis::define_nonpeaking() {
  
  ws_->factory("N_nonpeaking[0, 100]");
  ws_->factory("m0_nonpeaking[5.5]");
  ws_->factory("c_nonpeaking[1., 0.1, 20]");
  ws_->factory("p_nonpeaking[0.5, 0.1, 5.]");

  ws_->factory("ArgusBG::pdf_nonpeaking(Mass,m0_nonpeaking,c_nonpeaking,p_nonpeaking)");
  
}

void pdf_analysis::define_comb() {
  
  ws_->factory("N_comb[0, 100]");
  
  ws_->factory("Uniform::pdf_comb(Mass)");
  
}

void pdf_analysis::define_signals() {
  
  ws_->factory("N_signals[0, 100]");
  ws_->factory("bsfraction_signals[0.5, 0.0, 1.0]"); 
  
  ws_->factory("SUM::pdf_signals(bsfraction_signals*pdf_bs, pdf_bd)");
  
  ws_->factory("SUM::pdf_ext_signals(N_bs*pdf_bs, N_bd*pdf_bd)");

}

void pdf_analysis::define_rare() {
  
  ws_->factory("N_rare[0, 100]");
  ws_->factory("peakingfraction_rare[0.5, 0.0, 1.0]"); 
  ws_->factory("SUM::pdf_rare(peakingfraction_rare*pdf_peaking, pdf_nonpeaking)");
}

void pdf_analysis::define_rare2(RooDataHist* data) {
  ws_->factory("N_hist[0, 100]");
  RooHistPdf* pdf_rare = new RooHistPdf("pdf_hist", "pdf_hist", *ws_->var("Mass"), *data, 4);
  ws_->import(*pdf_rare);
}

void pdf_analysis::define_rare3() {
  ws_->factory("N_expo[0, 100]");
  ws_->factory("Exponential::pdf_expo(Mass,Alpha1[-10.,10.])");
  //ws_->factory("SUM::pdf_expo(expo1fraction_rare[0.,1.]*Exponential::expo1(Mass,Alpha1[-10.,10.]),Exponential::expo2(Mass,Alpha2[-10.,10.]))");
}

void pdf_analysis::define_signalsrare() {
  
  ws_->factory("signalfraction_signalsrare[0.5, 0.0, 1.0]"); 
  ws_->factory("SUM::pdf_signalsrare(signalfraction_signalsrare*pdf_signals, pdf_rare)");
  ws_->factory("SUM::pdf_ext_signalsrare(N_bs*pdf_bs, N_bd*pdf_bd, N_rare*pdf_rare)");
}

void pdf_analysis::define_bkg_fractional() {
  
  ws_->factory("N_bkg[0, 100]");
  ws_->factory("rarefraction_bkg[0.5, 0.0, 1.0]"); 
  
  ws_->factory("SUM::pdf_bkg(rarefraction_bkg*pdf_rare, pdf_comb)");
  
  RooExtendPdf* bkgExt = new RooExtendPdf("pdf_fract_ext_bkg", "bkgExt", *ws_->pdf("pdf_bkg"), *ws_->var("N_bkg"), /* range_.c_str()*/ "sb_lo,sb_hi"); /// WARNING
  ws_->import(*bkgExt);
}

void pdf_analysis::define_bkg_extended() {

  RooExtendPdf rare_ext("pdf_ext_rare", "rare_ext", *ws_->pdf("pdf_rare"), *ws_->var("N_rare"), "sb_lo,sb_hi");
  RooExtendPdf comb_ext("pdf_ext_comb", "comb_ext", *ws_->pdf("pdf_comb"), *ws_->var("N_comb"), "sb_lo,sb_hi");
  RooAddPdf bkg_ext("pdf_ext_bkg", "pdf_ext_bkg", RooArgList(rare_ext, comb_ext));
  ws_->import(bkg_ext);
//  ws_->factory("SUM::pdf_ext_bkg(N_rare*pdf_rare, N_comb*pdf_comb)");

}


void pdf_analysis::define_total_fractional() {
  
  ws_->factory("N_all[0, 100]");
  ws_->factory("signalsfraction_total[0.5, 0.0, 1.0]");
  
  ws_->factory("SUM::pdf_frac_total(signalsfraction_total*pdf_signals, pdf_bkg)");
  
  RooExtendPdf* allExt = new RooExtendPdf("pdf_frac_ext_total", "allExt", *ws_->pdf("pdf_frac_total"), *ws_->var("N_all"), range_.c_str());
  ws_->import(*allExt);
}

void pdf_analysis::define_total_extended() {
  
  if (SM_ || bd_constr_) ws_->factory("SUM::pdf_ext_total(N_bs*pdf_bs, N_bd_constr*pdf_bd, N_rare*pdf_rare, N_comb*pdf_comb)");
  else ws_->factory("SUM::pdf_ext_total(N_bs*pdf_bs, N_bd*pdf_bd, N_rare*pdf_rare, N_comb*pdf_comb)");
  
  return; 
}

string pdf_analysis::define_pdf_sum(string name) {

  vector <string> pdfs;
  size_t found;
  found = name.find("bs");
  if (found != string::npos) pdfs.push_back("bs");
  found = name.find("bd");
  if (found != string::npos) pdfs.push_back("bd");
  found = name.find("rare");
  if (found != string::npos) pdfs.push_back("rare");
  found = name.find("hist");
  if (found != string::npos) pdfs.push_back("hist");
  found = name.find("expo");
  if (found != string::npos) pdfs.push_back("expo");
  found = name.find("comb");
  if (found != string::npos) pdfs.push_back("comb");

  string pdf_sum = "SUM::pdf_ext_";
  string title;
  for (unsigned int i = 0; i < pdfs.size(); i++) {
    pdf_sum += pdfs[i];
    title += pdfs[i];
  }
  pdf_sum += "(";
  for (unsigned int i = 0; i < pdfs.size(); i++) {
    pdf_sum += "N_";
    pdf_sum += pdfs[i];
    pdf_sum += "*pdf_";
    pdf_sum += pdfs[i];
    if (i != pdfs.size() -1) pdf_sum += ",";
  }
  pdf_sum += ")";
  cout << "formed pdf: " << pdf_sum << endl;
  ws_->factory(pdf_sum.c_str());
  return (title);
}

void pdf_analysis::print(string output, RooWorkspace* ws) {
  if (ws == 0)  ws = ws_;
  int colors[11] = {632, 400, 616, 432, 800, 416, 820, 840, 860, 880, 900};
  
  RooPlot *rp = ws->var("Mass")->frame();
  rds_->plotOn(rp, Binning(20), RefreshNorm());
  ws->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue), Range(range_.c_str()));
  ws->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
  
  RooArgSet * set = ws->pdf(pdf_name.c_str())->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  int i = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (i > 11) i = 0;
      size_t found1 = pdf_name.find("total");
      if (found1 == string::npos) ws->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
      else {
        if (name=="pdf_bs" || name=="pdf_bd" || name=="pdf_rare" || name=="pdf_comb") ws->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
      }
      i++;
    }
  }
  
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  string address = "fig/" + pdf_name + "_" + meth_ + "_" + ch_s_ + output;
  if (SM_) address += "_SM";
  if (bd_constr_) address += "_BdConst";
  c->Print( (address + ".gif").c_str());
  c->Print( (address + ".pdf").c_str());
  delete rp;
  delete c;
  return;
}


