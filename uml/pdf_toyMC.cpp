/*
 * File:   pdf_toyMC.cpp
 * Author: lucamartini
 * 
 * Created on 26 aprile 2012, 11.03
 */

#include "pdf_toyMC.h"

pdf_toyMC::pdf_toyMC(string input_estimates, bool print, string meth, string ch_s, string range): pdf_analysis(print, meth, ch_s, range) {
  
  input_estimates_ = input_estimates;
  pull_h_bs = new TH1D("pull_h_bs", "pull_h_bs", 100, -5., 5.);
  pull_h_bd = new TH1D("pull_h_bd", "pull_h_bd", 100, -5., 5.);
  
}

void pdf_toyMC::parse_estimate(){
  char buffer[1024];
  char cutName[128];
  float cut;
  int ok;
  FILE *estimate_file = fopen(input_estimates_.c_str(), "r");
  cout << "event estimates:" << endl;
  while (fgets(buffer, sizeof(buffer), estimate_file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    // skip comment line
    if (buffer[0] == '#') continue;
    ok = sscanf(buffer, "%s %f", cutName, &cut);
    if (!parse(cutName, cut)) {
      cout << "==> lmtreeReader: Error parsing variable " << cutName << ". abort..." << endl;
      exit(1);
    }
  }
  if (estimate_file) fclose(estimate_file);
  return;
}

bool pdf_toyMC::parse(char *cutName, float cut) {
  // parse the options...
  if (!strcmp(cutName,"bs")) {
    estimate_bs = (int)cut;
    cout << "bs: " << estimate_bs << endl;
    return true;
  }
  if (!strcmp(cutName,"bd")) {
    estimate_bd = (int)cut;
    cout << "bd: " << estimate_bd << endl;
    return true;
  }
  if (!strcmp(cutName,"rare")) {
    estimate_rare = (int)cut;
    cout << "rare: " << estimate_rare << endl;
    return true;
  }
  if (!strcmp(cutName,"comb")) {
    estimate_comb = (int)cut;
    cout << "comb: " << estimate_comb << endl;
    return true;
  }
  return false;
}

void pdf_toyMC::generate(int NExp, string pdf_toy) {
  
  pull_bs = new RooRealVar("pull_bs", "pull_bs", -5., 5.);
  pull_bd = new RooRealVar("pull_bd", "pull_bd", -5., 5.);
  pull_rds_bs = new RooDataSet("pull_rds_bs", "pull_rds_bs", *pull_bs);
  pull_rds_bd = new RooDataSet("pull_rds_bd", "pull_rds_bd", *pull_bd);
  
  parse_estimate();
  
  RooRandom::randomGenerator()->SetSeed(0);
  double p_bs = 0., p_bd = 0.;
  for (int i = 1; i <= NExp; i++) {
    pdf_toy_ = pdf_toy;
    RooDataSet* data_bs   = ws_->pdf("pdf_bs")->generate(*ws_->var("Mass"), NumEvents((int)estimate_bs), Extended());
    RooDataSet* data_bd   = ws_->pdf("pdf_bd")->generate(*ws_->var("Mass"), NumEvents((int)estimate_bd), Extended());
    RooDataSet* data_rare = ws_->pdf("pdf_rare")->generate(*ws_->var("Mass"), NumEvents((int)estimate_rare), Extended());
    RooDataSet* data_comb = ws_->pdf("pdf_comb")->generate(*ws_->var("Mass"), NumEvents((int)estimate_comb), Extended());
    
    RooDataSet* data = new RooDataSet("data", "data", *ws_->var("Mass"));
    data->append(*data_bs);
    if (pdf_toy_ != "bs") {
      data->append(*data_bd);
      if (pdf_toy_ != "signals") {
        data->append(*data_rare);
        if (pdf_toy_ != "signalsrare") {
          data->append(*data_comb);
        }
      }
    }
    
    double printlevel = -1;
    if (i == 1) printlevel = 1;
    //////
    fit_pdf(pdf_toy_, data, printlevel);
    //////
    if (i == 1) {
      rds_ = data;
      pdf_name = pdf_toy_;
      print("_toyMC");
    }
      
    /// test statistics
    // PULL
    double bs_pull = (ws_->var("N_bs")->getVal() - data_bs->numEntries()) / ws_->var("N_bs")->getError();
    double bd_pull = (ws_->var("N_bd")->getVal() - data_bd->numEntries()) / ws_->var("N_bd")->getError();
    pull_bs->setVal(bs_pull);
    pull_bd->setVal(bd_pull);
    pull_rds_bs->add(*pull_bs);
    pull_rds_bd->add(*pull_bd);
    
    // H_0
    if (pdf_toy_ == "pdf_ext_total") {
      RooDataSet* bkg = new RooDataSet("bkg", "bkg", *ws_->var("Mass"));
      bkg->append(*data_rare);
      bkg->append(*data_comb);
      fit_pdf("total", bkg, printlevel);
      if (ws_->var("N_bs")->getVal() >= estimate_bs) p_bs++;
      if (ws_->var("N_bd")->getVal() >= estimate_bd) p_bd++;
    }
  }

  pdf_toy_ = pdf_toy;
  fit_pulls();
  
  double p_value_bs = p_bs / NExp;
  double p_value_bd = p_bd / NExp;
  cout << "p-value bs = " << p_value_bs << endl;
  cout << "p-value bd = " << p_value_bd << endl;

  return;
}

void pdf_toyMC::fit_pdf (string pdf, RooAbsData* data, int printlevel) {
  pdf_toy_ = "pdf_ext_" + pdf;
  ws_->pdf( pdf_toy_.c_str())->fitTo(*data, Extended(true), SumW2Error(0), Range(range_.c_str()), PrintLevel(printlevel));
}

void pdf_toyMC::fit_pulls() {
  RooRealVar* mean_bs = new RooRealVar("mean_bs", "mean_bs", -5., 5.);
  RooRealVar* sigma_bs = new RooRealVar("sigma_bs", "sigma_bs", 0.001, 5.);
  RooGaussian* gauss_bs = new RooGaussian("gauss_bs", "gauss_bs", *pull_bs, *mean_bs, *sigma_bs);
  gauss_bs->fitTo(*pull_rds_bs);
  
  RooPlot *rp_bs = pull_bs->frame();
  pull_rds_bs->plotOn(rp_bs, Binning(40));
  gauss_bs->plotOn(rp_bs, LineColor(kBlue));
  gauss_bs->paramOn(rp_bs, Layout(0.66, 0.9, 0.9));
  TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 600, 600); 
  rp_bs->Draw();
  string address = "fig/bs_pull_" + meth_ + "_" + ch_s_ + "_" + pdf_toy_;
  canvas_bs->Print( (address + ".gif").c_str());
  delete rp_bs;
  delete canvas_bs;
  
  RooRealVar* mean_bd = new RooRealVar("mean_bd", "mean_bd", -5., 5.);
  RooRealVar* sigma_bd = new RooRealVar("sigma_bd", "sigma_bd", 0.001, 5.);
  RooGaussian* gauss_bd = new RooGaussian("gauss_bd", "gauss_bd", *pull_bd, *mean_bd, *sigma_bd);
  gauss_bd->fitTo(*pull_rds_bd);
  
  RooPlot *rp_bd = pull_bd->frame();
  pull_rds_bd->plotOn(rp_bd, Binning(40));
  gauss_bd->plotOn(rp_bd, LineColor(kBlue));
  gauss_bd->paramOn(rp_bd, Layout(0.66, 0.9, 0.9));
  TCanvas* canvas_bd = new TCanvas("canvas_bd", "canvas_bd", 600, 600); 
  rp_bd->Draw();
  address = "fig/bd_pull_" + meth_ + "_" + ch_s_ + "_" + pdf_toy_;
  canvas_bd->Print( (address + ".gif").c_str());
  delete rp_bd;
  delete canvas_bd;
  
}

void pdf_toyMC::unset_constant() {
  cout << "unset constants" << endl;
  RooArgSet set = ws_->allVars();
  TIterator* it = set.createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    ws_->var(var_Obj->GetName())->setConstant(0);
  }
}
