/*
 * File:   pdf_toyMC.cpp
 * Author: lucamartini
 * 
 * Created on 26 aprile 2012, 11.03
 */

#include "pdf_toyMC.h"

pdf_toyMC::pdf_toyMC(bool print, int inputs, string input_estimates, string input_cuts, string meth, string range, bool SM, bool bd_constr, TTree *input_tree, string bias, bool simul, bool pee_, bool bdt_fit, string ch_s, int sig): pdf_fitData( print,  inputs,  input_estimates, input_cuts,  meth,  range,  SM,  bd_constr,  input_tree, simul, pee_, bdt_fit, ch_s, sig) {
  cout << "pdf_toyMC constructor" << endl;
  bias_ = bias;
}

pdf_toyMC::~pdf_toyMC() {
  cout << "pdf_toyMC destructor" << endl;
}

void pdf_toyMC::generate(int NExp, string pdf_toy, string test_pdf) {

  if (!simul_) channels = 1;
  pull_bs.resize(channels);
  pull_bd.resize(channels);
  pull_rare.resize(channels);
  pull_comb.resize(channels);
  pull_rds_bs.resize(channels);
  pull_rds_bd.resize(channels);
  pull_rds_rare.resize(channels);
  pull_rds_comb.resize(channels);
//  vector <double> p_bs(channels);
//  vector <double> p_bd(channels);
  vector <TH1D*> correlation_h(channels);
  vector <TH1D*> bs_mean_h(channels);
  vector <TH1D*> bd_mean_h(channels);
  vector <TH1D*> rare_mean_h(channels);
  vector <TH1D*> comb_mean_h(channels);
  vector <TH1D*> bs_sigma_h(channels);
  vector <TH1D*> bd_sigma_h(channels);
  vector <TH1D*> rare_sigma_h(channels);
  vector <TH1D*> comb_sigma_h(channels);
  vector <TH2D*> corr_Nbs_Nbd_vs_N_bs_h(channels);

  TH1D * sign_h = new TH1D("sign_h", "sign_h", 100, 0, 10);

  for (int i = 0; i < channels; i++) {
    pull_bs[i] = new RooRealVar("pull_bs", "pull_bs", -8., 8.);
    pull_bd[i] = new RooRealVar("pull_bd", "pull_bd", -8., 8.);
    pull_rare[i] = new RooRealVar("pull_rare", "pull_rare", -8., 8.);
    pull_comb[i] = new RooRealVar("pull_comb", "pull_comb", -8., 8.);
    pull_rds_bs[i] = new RooDataSet("pull_rds_bs", "pull_rds_bs", *pull_bs[i]);
    pull_rds_bd[i] = new RooDataSet("pull_rds_bd", "pull_rds_bd", *pull_bd[i]);
    pull_rds_rare[i] = new RooDataSet("pull_rds_rare", "pull_rds_rare", *pull_rare[i]);
    pull_rds_comb[i] = new RooDataSet("pull_rds_comb", "pull_rds_comb", *pull_comb[i]);
//    p_bs[i] = 0;
//    p_bd[i] = 0;
    correlation_h[i] = new TH1D(Form("correlation_%d_h", i), "correlation_h", 100, -1., 1.);
    bs_mean_h[i] = new TH1D(name("bs_mean_h", i), name("bs_mean_h", i), 100, 0., max(100., estimate_bs[i]*2));
    bd_mean_h[i] = new TH1D(name("bd_mean_h", i), name("bd_mean_h", i), 100, 0., max(100., estimate_bd[i]*2));
    rare_mean_h[i] = new TH1D(name("rare_mean_h", i), name("rare_mean_h", i), 100, 0., max(100., estimate_rare[i]*2));
    comb_mean_h[i] = new TH1D(name("comb_mean_h", i), name("comb_mean_h", i), 100, 0., max(100., estimate_comb[i]*2));
    bs_sigma_h[i] = new TH1D(name("bs_sigma_h", i), name("bs_sigma_h", i), 100, 0., max(100., estimate_bs[i]*2));
    bd_sigma_h[i] = new TH1D(name("bd_sigma_h", i), name("bd_sigma_h", i), 100, 0., max(100., estimate_bd[i]*2));
    rare_sigma_h[i] = new TH1D(name("rare_sigma_h", i), name("rare_sigma_h", i), 100, 0., max(100., estimate_rare[i]*2));
    comb_sigma_h[i] = new TH1D(name("comb_sigma_h", i), name("comb_sigma_h", i), 100, 0., max(100., estimate_comb[i]*2));
    corr_Nbs_Nbd_vs_N_bs_h[i] = new TH2D(Form("corr_Nbs_Nbd_vs_N_bs_%d_h", i), "corr_Nbs_Nbd_vs_N_bs_h", 30, 0., 30., 100, -1., 1.);
  }
  if (!simul_) {
    pdf_name = "total";
    if (pdf_toy != "total") {
      pdf_name = define_pdf_sum(pdf_toy);
    }
    pdf_test_= "total";
    if (test_pdf != "total") {
      pdf_test_ = define_pdf_sum(test_pdf);
    }
  }
  else {
    pdf_name = "simul";
    pdf_test_= "simul";
  }
  
  bool bd_b = false, corr0 = false, corrneg = false;

  for (int i = 1; i <= NExp; i++) {
    if (i%100 == 0) cout << "Exp # " << i << " / " << NExp << endl;
    pdf_toy_ = pdf_toy;
    RooWorkspace* ws_temp = (RooWorkspace*)ws_->Clone("ws_temp");
    RooDataSet* data = new RooDataSet("data", "data", RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes")));
    double printlevel = -1;
    if (i == 1) printlevel = 1;

    for (int j = 0; j < channels; j++) {
      ws_temp->var(name("N_bs", j))->setVal((int)estimate_bs[j]);
      if (!SM_ && !bd_constr_) ws_temp->var(name("N_bd", j))->setVal((int)estimate_bd[j]);
      ws_temp->var(name("N_rare", j))->setVal((int)estimate_rare[j]);
      ws_temp->var(name("N_comb", j))->setVal((int)estimate_comb[j]);

      ws_temp->var(name("N_bs", j))->setConstant(kFALSE);
      if (!SM_ && !bd_constr_) ws_temp->var(name("N_bd", j))->setConstant(kFALSE);
      ws_temp->var(name("N_rare", j))->setConstant(kFALSE);
      ws_temp->var(name("N_comb", j))->setConstant(kFALSE);
    }
    if (channels == 1) {
      if (pdf_toy_ == "total") {
        bd_b = true;
        data = ws_temp->pdf("pdf_ext_total")->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes")), Extended(1));
      }
      else {
        size_t found;
        found = pdf_toy_.find("bd");
        if (found != string::npos) bd_b = true;
        data = ws_temp->pdf( ("pdf_ext_" + pdf_toy_).c_str())->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes")), Extended(1));
      }
    }
    else {
      bd_b = true;
      RooCategory* cat =  (RooCategory*)ws_temp->obj("channels");
      data = ws_temp->pdf("pdf_ext_simul")->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes"), *cat), Extended(1));
    }
    do_bias(ws_temp);
    //////
    RFR = pdf_toyMC::fit_pdf(pdf_test_, data, printlevel, ws_temp);
    //////
    pdf_name = pdf_toy_;
    rds_ = data;

    if (!simul_) {
      if (i == 2) {
        ws_ = ws_temp;
        print("_first", ws_temp);
      }
      if (i == NExp) {
        ws_ = ws_temp;
        print("_last", ws_temp);
      }
    }
    /// pull
    for (int j = 0; j < channels; j++) {
      bs_mean_h[j]->Fill(ws_temp->var(name("N_bs", j))->getVal());
      if (!(SM_ || bd_constr_)) bd_mean_h[j]->Fill(ws_temp->var(name("N_bd", j))->getVal());
      rare_mean_h[j]->Fill(ws_temp->var(name("N_rare", j))->getVal());
      comb_mean_h[j]->Fill(ws_temp->var(name("N_comb", j))->getVal());
      bs_sigma_h[j]->Fill(ws_temp->var(name("N_bs", j))->getError());
      if (!(SM_ || bd_constr_)) bd_sigma_h[j]->Fill(ws_temp->var(name("N_bd", j))->getError());
      rare_sigma_h[j]->Fill(ws_temp->var(name("N_rare", j))->getError());
      comb_sigma_h[j]->Fill(ws_temp->var(name("N_comb", j))->getError());
      // PULL
      double bs_pull = (ws_temp->var(name("N_bs", j))->getVal() - estimate_bs[j]) / ws_temp->var(name("N_bs", j))->getError();
      double bd_pull;
      if (!(SM_ || bd_constr_)) bd_pull = (ws_temp->var(name("N_bd", j))->getVal() - estimate_bd[j]) / ws_temp->var(name("N_bd", j))->getError();
      double rare_pull = (ws_temp->var(name("N_rare", j))->getVal() - estimate_rare[j]) / ws_temp->var(name("N_rare", j))->getError();
      double comb_pull = (ws_temp->var(name("N_comb", j))->getVal() - estimate_comb[j]) / ws_temp->var(name("N_comb", j))->getError();
      pull_bs[j]->setVal(bs_pull);
      pull_rds_bs[j]->add(*pull_bs[j]);
      if (!(SM_ || bd_constr_)) {
        pull_bd[j]->setVal(bd_pull);
        pull_rds_bd[j]->add(*pull_bd[j]);
      }
      pull_rare[j]->setVal(rare_pull);
      pull_rds_rare[j]->add(*pull_rare[j]);
      pull_comb[j]->setVal(comb_pull);
      pull_rds_comb[j]->add(*pull_comb[j]);
    }
    if (i == NExp) ws_temp->pdf(pdf_toy_.c_str())->Print();

    if (sign == 0) {
      sign_h->Fill(sig_hand(data, printlevel, ws_temp));
    }

    delete data;
    delete RFR;
  }

  pdf_toy_ = pdf_toy;
  fit_pulls(pull_bs, pull_rds_bs);
  if (!(SM_ || bd_constr_)) fit_pulls(pull_bd, pull_rds_bd);
  fit_pulls(pull_rare, pull_rds_rare);
  fit_pulls(pull_comb, pull_rds_comb);

  print_histos(bs_mean_h);
  if (!(SM_ || bd_constr_)) print_histos(bd_mean_h);
  print_histos(rare_mean_h);
  print_histos(comb_mean_h);

  print_histos(bs_sigma_h);
  if (!(SM_ || bd_constr_)) print_histos(bd_sigma_h);
  print_histos(rare_sigma_h);
  print_histos(comb_sigma_h);

  if (!SM_ && bd_b && !simul_ && false) {
    TCanvas* corr_c = new TCanvas("corr_c", "corr_c", 1200, 600);
    corr_c->Divide(2);
    corr_c->cd(1);
    correlation_h[0]->Draw();
    corr_c->cd(2);
    corr_Nbs_Nbd_vs_N_bs_h[0]->Draw();
    string address = "fig/corr_" + meth_ + "_" + ch_s_ + "_" + pdf_toy;
    corr_c->Print((address + ".gif").c_str());
    corr_c->Print((address + ".pdf").c_str());
    delete corr_c;
  }

  if (sign == 0) {
    TCanvas* sign_c = new TCanvas("sign_c", "sign_c", 600, 600);
    sign_h->Draw("e");
    ostringstream address_oss;
    if (!simul_) address_oss << "fig/sign_" << meth_ << "_" << ch_s_ << "_" << pdf_toy_;
    else address_oss << "fig/sign" << meth_ << "_simul_" << pdf_toy_;
    string address = address_oss.str();
    if (SM_) address += "_SM";
    if (bd_constr_) address += "_bd_const";
    if (pee) address += "_pee";
    if (bdt_fit_) address += "_2D";
    sign_c->Print( (address + ".gif").c_str());
    sign_c->Print( (address + ".pdf").c_str());
    delete sign_c;
  }
  return;
}

void pdf_toyMC::do_bias(RooWorkspace* ws) {

  if (bias_.compare("no")) {
    for (int k = 0; k < channels; k++) {
      size_t found;
      found = bias_.find("tau");
      if (found != string::npos) {
        double value =  ws->var(name("tau", ch_i_))->getVal();
        found = bias_.find("+");
        if (found != string::npos) {
          double error = ws->var(name("tau", ch_i_))->getErrorHi();
          ws->var(name("C0", ch_i_))->setVal(value + error);
        }
        else {
          double error = ws->var(name("tau", ch_i_))->getErrorLo();
          ws->var(name("C0", ch_i_))->setVal(value - error);
        }
      }
      found = bias_.find("c");
      if (found != string::npos) {
        double value =  ws->var(name("C0", ch_i_))->getVal();
        found = bias_.find("+");
        if (found != string::npos) {
          double error = ws->var(name("C0", ch_i_))->getErrorHi();
          ws->var(name("C0", ch_i_))->setVal(value + error);
        }
        else {
          double error = ws->var(name("C0", ch_i_))->getErrorLo();
          ws->var(name("C0", ch_i_))->setVal(value - error);
        }
      }
      found = bias_.find("sig");
      if (found != string::npos) {
        found = bias_.find("+");
        if (found != string::npos) {
          ws->var(name("Mean_bs", ch_i_))->setVal(ws->var(name("Mean_bs", ch_i_))->getVal() + ws->var(name("Mean_bs", ch_i_))->getErrorHi());
          ws->var(name("Sigma_bs", ch_i_))->setVal(ws->var(name("Sigma_bs", ch_i_))->getVal() + ws->var(name("Sigma_bs", ch_i_))->getErrorHi());
          ws->var(name("Alpha_bs", ch_i_))->setVal(ws->var(name("Alpha_bs", ch_i_))->getVal() + ws->var(name("Alpha_bs", ch_i_))->getErrorHi());
          ws->var(name("Enne_bs", ch_i_))->setVal(ws->var(name("Enne_bs", ch_i_))->getVal() + ws->var(name("Enne_bs", ch_i_))->getErrorHi());
          ws->var(name("Sigma2_bs", ch_i_))->setVal(ws->var(name("Sigma2_bs", ch_i_))->getVal() + ws->var(name("Sigma2_bs", ch_i_))->getErrorHi());
          ws->var(name("CoeffGauss_bs", ch_i_))->setVal(ws->var(name("CoeffGauss_bs", ch_i_))->getVal() + ws->var(name("CoeffGauss_bs", ch_i_))->getErrorHi());
        }
        else {
          ws->var(name("Mean_bs", ch_i_))->setVal(ws->var(name("Mean_bs", ch_i_))->getVal() - ws->var(name("Mean_bs", ch_i_))->getErrorLo());
          ws->var(name("Sigma_bs", ch_i_))->setVal(ws->var(name("Sigma_bs", ch_i_))->getVal() - ws->var(name("Sigma_bs", ch_i_))->getErrorLo());
          ws->var(name("Alpha_bs", ch_i_))->setVal(ws->var(name("Alpha_bs", ch_i_))->getVal() - ws->var(name("Alpha_bs", ch_i_))->getErrorLo());
          ws->var(name("Enne_bs", ch_i_))->setVal(ws->var(name("Enne_bs", ch_i_))->getVal() - ws->var(name("Enne_bs", ch_i_))->getErrorLo());
          ws->var(name("Sigma2_bs", ch_i_))->setVal(ws->var(name("Sigma2_bs", ch_i_))->getVal() - ws->var(name("Sigma2_bs", ch_i_))->getErrorLo());
          ws->var(name("CoeffGauss_bs", ch_i_))->setVal(ws->var(name("CoeffGauss_bs", ch_i_))->getVal() - ws->var(name("CoeffGauss_bs", ch_i_))->getErrorLo());
        }
      }
    }
  }

}

RooFitResult* pdf_toyMC::fit_pdf(string pdf, RooAbsData* data, int printlevel, RooWorkspace* ws) {

  pdf_toy_ = "pdf_ext_" + pdf;
  RooFitResult* result;
  if (!pee) result = ws->pdf(pdf_toy_.c_str())->fitTo(*data, Extended(true), SumW2Error(0), PrintLevel(printlevel), PrintEvalErrors(-1), Save(kTRUE), NumCPU(2));
  else result = ws->pdf(pdf_toy_.c_str())->fitTo(*data, ConditionalObservables(*ws->var("MassRes")), Extended(true), SumW2Error(0), PrintLevel(printlevel), PrintEvalErrors(-1), Save(kTRUE), NumCPU(2));
  return result;
}

void pdf_toyMC::fit_pulls(vector <RooRealVar*> pull,  vector <RooDataSet*> rds) {
  for (int i = 0; i < channels; i++) {
    RooRealVar* mean_bs = new RooRealVar("mean_bs", "mean_bs", -5., 5.);
    RooRealVar* sigma_bs = new RooRealVar("sigma_bs", "sigma_bs", 0.001, 5.);
    RooGaussian* gauss_bs = new RooGaussian("gauss_bs", "gauss_bs", *pull[i], *mean_bs, *sigma_bs);
    gauss_bs->fitTo(*rds[i]);
  
    RooPlot *rp_bs = pull[i]->frame();
    rds[i]->plotOn(rp_bs, Binning(40));
    gauss_bs->plotOn(rp_bs, LineColor(kBlue));
    gauss_bs->paramOn(rp_bs, Layout(0.66, 0.9, 0.9));
    TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 600, 600);
    rp_bs->Draw();
    ostringstream address_oss;
    if (!simul_) address_oss << "fig/" << pull[0]->GetName() << "_" << meth_ << "_" << ch_s_ << "_" << pdf_toy_;
    else address_oss << "fig/" << pull[0]->GetName() << "_" << meth_ << "_" << i << "_" << pdf_toy_;
    string address = address_oss.str();
    if (SM_) address += "_SM";
    if (bd_constr_) address += "_bd_const";
    if (pee) address += "_pee";
    if (bdt_fit_) address += "_2D";
    canvas_bs->Print( (address + ".gif").c_str());
    canvas_bs->Print( (address + ".pdf").c_str());
    delete rp_bs;
    delete canvas_bs;
  }
}

void pdf_toyMC::print_histos(vector <TH1D*> histos) {
  for (int j = 0; j < channels; j++) {
    cout << "channel " << j << " " << histos[j]->GetName() << " mean = " << histos[j]->GetMean() << endl;
    TCanvas* N_mean_c = new TCanvas("N_mean_c", "N_mean_c", 600, 600);
    histos[j]->Draw("e");
    ostringstream address;
    if (!simul_) address << "fig/" << histos[j]->GetName() << "_" << meth_ << "_" << ch_s_ << "_" + pdf_toy_;
    else address << "fig/"<< histos[j]->GetName() << "_" << meth_ << "_simul_" + pdf_toy_;
    if (SM_) address << "_SM";
    if (bd_constr_) address << "_bd_const";
    if (pee) address << "_pee";
    if (bdt_fit_) address << "_2D";
    N_mean_c->Print((address.str() + ".gif").c_str());
    N_mean_c->Print((address.str() + ".pdf").c_str());
    delete N_mean_c;
  }
}

Double_t pdf_toyMC::sig_hand(RooAbsData* data, int printlevel, RooWorkspace* ws_temp) {
  Double_t minNLL = RFR->minNll();
  for (int j = 0; j < channels; j++) {
    ws_temp->var(name("N_bs", j))->setVal(0);
    ws_temp->var(name("N_bs", j))->setConstant(1);
    if ( !(SM_ || bd_constr_)) {
      ws_temp->var(name("N_bd", j))->setVal(0);
      ws_temp->var(name("N_bd", j))->setConstant(1);
    }
    else if (bd_constr_) {
      ws_temp->var("Bd_over_Bs")->setVal(0);
      ws_temp->var("Bd_over_Bs")->setConstant(1);
    }
  }
  RooFitResult * rfr_H0 = pdf_toyMC::fit_pdf(pdf_test_, data, printlevel, ws_temp);
  for (int j = 0; j < channels; j++) {
    ws_temp->var(name("N_bs", j))->setVal(estimate_bs[j]);
    ws_temp->var(name("N_bs", j))->setConstant(0);
    if ( !(SM_ || bd_constr_)) {
      ws_temp->var(name("N_bd", j))->setVal(estimate_bd[j]);
      ws_temp->var(name("N_bd", j))->setConstant(0);
    }
    else if (bd_constr_) {
      ws_temp->var("Bd_over_Bs")->setVal(estimate_bd[j]/estimate_bs[j]);
      ws_temp->var("Bd_over_Bs")->setConstant(0);
    }
  }
  Double_t newNLL = rfr_H0->minNll();
  Double_t deltaLL = newNLL - minNLL;
  Double_t signif = deltaLL>0 ? sqrt(2*deltaLL) : -sqrt(-2*deltaLL) ;
  return signif;
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

void pdf_toyMC::mcstudy(int NExp, string pdf_toy) {

  for (int i = 0; i < channels; i++) {
    ws_->var(name("N_bs", i))->setVal(estimate_bs[i]);
    ws_->var(name("N_rare", i))->setVal(estimate_rare[i]);
    ws_->var(name("N_comb", i))->setVal(estimate_comb[i]);
    if (!SM_ && !bd_constr_) ws_->var(name("N_bd", i))->setVal(estimate_bd[i]);
  }
  if (bd_constr_) {
    double ratio = (double) estimate_bd[0] / estimate_bs[0]; // it's the same in every channel
    ws_->var("Bd_over_Bs")->setVal(ratio);
  }

  pdf_name = "pdf_ext_total";
  if (pdf_toy != "total") {
    pdf_name = "pdf_ext_" + define_pdf_sum(pdf_toy);
  }
  if (simul_) pdf_name = "pdf_ext_simul";

  RooArgSet obsv(*ws_->var("Mass"), *ws_->var("bdt"), "obsv");
  if (simul_) obsv.add(*ws_->cat("channels"));

  RooMCStudy * mcstudy;
  if (bias_.compare("no")) {
    RooWorkspace* fitws = (RooWorkspace*)ws_->Clone("fitws");
    do_bias(fitws);
    if (!pee) mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv,  Binned(kFALSE), Extended(kTRUE), FitOptions(Save(kTRUE)), Silence(), FitModel(*fitws->pdf(pdf_name.c_str())));
    else mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv, ConditionalObservables(*ws_->var("MassRes")), Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(Save(kTRUE)), FitModel(*fitws->pdf(pdf_name.c_str())));
  }
  else {
    if (!pee) mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv,  Binned(kFALSE), Extended(kTRUE), FitOptions(Save(kTRUE)), Silence());
    else mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv, ConditionalObservables(*ws_->var("MassRes")), Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(PrintLevel(-1), Save(kTRUE), PrintEvalErrors(-1)));
  }
  vector <RooDLLSignificanceMCSModule*> sigModule;
  sigModule.resize(channels);
  for (int i = 0; i < channels; i++) {
    sigModule[i] = new RooDLLSignificanceMCSModule(*ws_->var(name("N_bs", i)), 0);
    mcstudy->addModule(*sigModule[i]) ;
  }

  mcstudy->generateAndFit(NExp, 0, kTRUE);

  for (int j = 0; j < 4; j++) {
    if ((SM_ || bd_constr_) && j == 1) continue;
    for (int i = 0; i < channels; i++) {
      RooPlot* frame1_bs = mcstudy->plotParam(*ws_->var(name("N_" + source[j], i)), Bins(20)) ;
      RooPlot* frame2_bs = mcstudy->plotError(*ws_->var(name("N_" + source[j], i)), Bins(20)) ;
      RooPlot* frame3_bs = mcstudy->plotPull(*ws_->var(name("N_" + source[j], i)), Bins(20), FitGauss(kTRUE)) ;
      TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 600, 600);
      frame1_bs->Draw();
      TPaveText* stats1 = new TPaveText(0.57, 0.66, 0.9, 0.9, "NDCR");
      ostringstream entries;
      entries << "entries = " << mcstudy->fitParDataSet().sumEntries();
      ostringstream mean;
      mean << "mean = " << mcstudy->fitParDataSet().mean(*ws_->var(name("N_" + source[j], i)));
      ostringstream sigma;
      sigma << "sigma = " << mcstudy->fitParDataSet().sigma(*ws_->var(name("N_" + source[j], i)));
      stats1->AddText(entries.str().c_str());
      stats1->AddText(mean.str().c_str());
      stats1->AddText(sigma.str().c_str());
      stats1->Draw();
      ostringstream index;
      index << i;
      string address;
      if (simul_) address = "fig/RooMCStudy_mean_simul_" + source[j] + "_" + meth_;
      else address = "fig/RooMCStudy_mean_" + source[j] + "_" + index.str() + "_" + meth_;
      if (SM_) address += "_SM";
      if (bd_constr_) address += "_bdConstr";
      canvas1->Print((address + ".gif").c_str());
      canvas1->Print((address + ".pdf").c_str());
      delete frame1_bs;
      delete canvas1;
      TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 600, 600);
      frame2_bs->Draw();
      TPaveText* stats2 = new TPaveText(0.57, 0.66, 0.9, 0.9, "NDCR");
      entries.str("");
      entries << "entries = " << mcstudy->fitParDataSet().sumEntries();
      mean.str("");
      mean << "mean = " << mcstudy->fitParDataSet().mean(*ws_->var(name("N_" + source[j], i)));
      sigma.str("");
      sigma << "sigma = " << mcstudy->fitParDataSet().sigma(*ws_->var(name("N_" + source[j], i)));
      stats2->AddText(entries.str().c_str());
      stats2->AddText(mean.str().c_str());
      stats2->AddText(sigma.str().c_str());
      stats2->Draw();
      if (simul_) address = "fig/RooMCStudy_sigma_simul_" + source[j] + "_" + meth_;
      else address = "fig/RooMCStudy_sigma_" + source[j] + "_" + index.str() + "_" + meth_;
      if (SM_) address += "_SM";
      if (bd_constr_) address += "_bdConstr";
      canvas2->Print((address + ".gif").c_str());
      canvas2->Print((address + ".pdf").c_str());
      delete frame2_bs;
      delete canvas2;
      TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 600, 600);
      frame3_bs->Draw();
      if (simul_) address = "fig/RooMCStudy_pull_simul_" + source[j] + "_" + meth_;
      else address = "fig/RooMCStudy_pull_" + source[j] + "_" + index.str() + "_" + meth_;
      if (SM_) address += "_SM";
      if (bd_constr_) address += "_bdConstr";
      canvas3->Print((address + ".gif").c_str());
      canvas3->Print((address + ".pdf").c_str());
      delete frame3_bs;
      delete canvas3;
    }
  }

  TCanvas* sig_c = new TCanvas("sig_c", "sig_c", 600*channels, 600);
  sig_c->Divide(channels);
  for (int i = 0; i < channels; i++) {
    sig_c->cd(i+1);
    TH1* sig_h = mcstudy->fitParDataSet().createHistogram(name("significance_nullhypo_N_bs", i));
    sig_h->Draw();
  }
  string address = "fig/RooMCStudy_simul_sig_bs_" + meth_;
  if (SM_) address += "_SM";
  sig_c->Print((address + ".gif").c_str());
  sig_c->Print((address + ".pdf").c_str());
  delete sig_c;
}

void pdf_toyMC::print(string output, RooWorkspace* ws) {
  int colors[] = {632, 400, 616, 432, 800, 416, 820, 840, 860, 880, 900};
  RooPlot *rp = ws_->var("Mass")->frame();
  RooPlot *rp_bdt = ws_->var("bdt")->frame();
  rds_->plotOn(rp, Binning(20));
  if (bdt_fit_) rds_->plotOn(rp_bdt, Binning(100));
  if (!pee) {
    ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue));
    if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, LineColor(kBlue));
  }
  else {
    ws_->pdf(pdf_name.c_str())->Print();
    TH1* mass_eta_h;
    mass_eta_h = ws_->pdf(pdf_name.c_str())->createHistogram("fit", *ws_->var("Mass"), Binning(50), YVar(*ws_->var("MassRes"), Binning(50))) ;
    mass_eta_h->SetLineColor(kBlue) ;
    mass_eta_h->GetXaxis()->SetTitleOffset(2.) ;
    mass_eta_h->GetYaxis()->SetTitleOffset(2.) ;
    mass_eta_h->GetZaxis()->SetTitleOffset(2.5) ;
    RooPlot* frame = Mass->frame() ;
    for (Int_t ibin=0 ; ibin<100; ibin+=20) {
      ws_->var("eta")->setBin(ibin) ;
      ws_->pdf(pdf_name.c_str())->plotOn(frame, Normalization(5)) ;
    }
    TCanvas* cetad = new TCanvas("cetad", "cetad", 1200, 600);
    cetad->Divide(3);
    cetad->cd(1);
    mass_eta_h->Draw("surf") ;
    cetad->cd(2);
    ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue));
    if(bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, LineColor(kBlue));
    if(!no_legend) {
      ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
      if(bdt_fit_) ws_->pdf(pdf_name.c_str())->paramOn(rp_bdt, Layout(0.50, 0.9, 0.9));
    }
    rp->Draw();
    if(bdt_fit_) {
      cetad->cd(3);
      rp_bdt->Draw();
    }
    string address;
    if (simul_) address = "fig/" + pdf_name + "_MassEta_" + meth_ + "_simul" + output;
    else address = "fig/" + pdf_name + "_MassEta_" + meth_ + "_" + ch_s_ + output;
    if (SM_) address += "_SM";
    if (bd_constr_) address += "_BdConst";
    if (bdt_fit_) address += "_2D";
    cetad->Print( (address + ".gif").c_str());
    cetad->Print( (address + ".pdf").c_str());
    delete cetad;
    delete frame;
    delete mass_eta_h;
    cout << pdf_name << endl;
  }
  if(!no_legend) {
    ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
    if(bdt_fit_) ws_->pdf(pdf_name.c_str())->paramOn(rp_bdt, Layout(0.50, 0.9, 0.9));
  }

  RooArgSet * set = ws_->pdf(pdf_name.c_str())->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  int i = 0;
  while((var_Obj = it->Next()) && !pee){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (i > 11) i = 0;
      size_t found1 = pdf_name.find("total");
      if (found1 == string::npos) {
        ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2));
        if(bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2));
      }
      else {
        if (name=="pdf_bs" || name=="pdf_bd" || name=="pdf_rare" || name=="pdf_comb") {
          ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2));
          if(bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2));
        }
      }
      i++;
    }
  }

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  string address;
  if (simul_) address = "fig/" + pdf_name + "_" + meth_ + "_simul" + output;
  else address = "fig/" + pdf_name + "_" + meth_ + "_" + ch_s_ + output;
  if (SM_) address += "_SM";
  if (bd_constr_) address += "_BdConst";
  if (bdt_fit_) address += "_2D";
  if (pee) address += "_PEE";
  c->Print( (address + ".gif").c_str());
  c->Print( (address + ".pdf").c_str());
  delete rp;
  delete c;

  if (bdt_fit_) {
    TCanvas* d = new TCanvas("d", "d", 600, 600);
    rp_bdt->Draw();
    string address_b;
    if (simul_) address_b = "fig/BDT_" + pdf_name + "_" + meth_ + "_simul" + output;
    else address_b = "fig/BDT_" + pdf_name + "_" + meth_ + "_" + ch_s_ + output;
    if (SM_) address_b += "_SM";
    if (bd_constr_) address_b += "_BdConst";
    if (bdt_fit_) address_b += "_2D";
    if (pee) address_b += "_PEE";
    d->Print( (address_b + ".gif").c_str());
    d->Print( (address_b + ".pdf").c_str());
    delete d;
  }
  delete rp_bdt;
  return;
}
