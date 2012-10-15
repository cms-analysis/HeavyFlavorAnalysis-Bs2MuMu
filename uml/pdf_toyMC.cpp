/*
 * File:   pdf_toyMC.cpp
 * Author: lucamartini
 * 
 * Created on 26 aprile 2012, 11.03
 */

#include "pdf_toyMC.h"

pdf_toyMC::pdf_toyMC(bool print, int inputs, int inputs_bdt, string input_estimates, string meth, string range, int BF, bool SM, bool bd_constr, TTree* input_tree, bool simul, bool simulbdt, bool pee_, bool bdt_fit, string ch_s, int sig, string bias): pdf_fitData( print,  inputs,  inputs_bdt,  input_estimates,  meth,  range, BF, SM,  bd_constr,  input_tree,  simul,  simulbdt,  pee_,  bdt_fit,  ch_s,  sig) {
  cout << "pdf_toyMC constructor" << endl;
  bias_ = bias;
}

pdf_toyMC::~pdf_toyMC() {
  cout << "pdf_toyMC destructor" << endl;
}

void pdf_toyMC::generate(int NExp, string pdf_toy, string test_pdf) {

  if (!simul_) channels = 1;
  if (!simul_bdt_) channels_bdt = 1;
  residual_bs.resize(channels);
  residual_bd.resize(channels);
  pull_bs.resize(channels);
  pull_bd.resize(channels);
  pull_rare.resize(channels);
  pull_comb.resize(channels);
  residual_rds_bs.resize(channels);
  residual_rds_bd.resize(channels);
  pull_rds_bs.resize(channels);
  pull_rds_bd.resize(channels);
  pull_rds_rare.resize(channels);
  pull_rds_comb.resize(channels);

  for (int k = 0; k < channels; k++) {
    residual_bs[k].resize(channels);
    residual_bd[k].resize(channels);
    pull_bs[k].resize(channels_bdt);
    pull_bd[k].resize(channels_bdt);
    pull_rare[k].resize(channels_bdt);
    pull_comb[k].resize(channels_bdt);
    residual_rds_bs[k].resize(channels);
    residual_rds_bd[k].resize(channels);
    pull_rds_bs[k].resize(channels_bdt);
    pull_rds_bd[k].resize(channels_bdt);
    pull_rds_rare[k].resize(channels_bdt);
    pull_rds_comb[k].resize(channels_bdt);
  }

  vector < vector  <TH1D*> > correlation_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > bs_mean_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > bd_mean_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > rare_mean_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > comb_mean_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > bs_sigma_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > bd_sigma_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > rare_sigma_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > comb_sigma_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH2D*> > corr_Nbs_Nbd_vs_N_bs_h(channels, vector <TH2D*> (channels_bdt));

  TH1D * sign_h = new TH1D("sign_h", "sign_h", 100, 0, 10);
  TH1D * BF_bs_mean_h = new TH1D("BF_bs_mean_h", "BF_bs_mean_h", 100, 1.e-10, 1.e-8);
  TH1D * BF_bs_sigma_h = new TH1D("BF_bs_sigma_h", "BF_bs_sigma_h", 100, 1.e-10, 1.e-8);
  pull_BF_bs = new RooRealVar("pull_BF_bs", "pull_BF_bs", -8., 8.);
  pull_rds_BF_bs = new RooDataSet("pull_rds_BF_bs", "pull_rds_BF_bs", *pull_BF_bs);
  TH1D * BF_bd_mean_h = new TH1D("BF_bd_mean_h", "BF_bd_mean_h", 100, 1.e-10, 1.e-8);
  TH1D * BF_bd_sigma_h = new TH1D("BF_bd_sigma_h", "BF_bd_sigma_h", 100, 1.e-10, 1.e-8);
  pull_BF_bd = new RooRealVar("pull_BF_bd", "pull_BF_bd", -8., 8.);
  pull_rds_BF_bd = new RooDataSet("pull_rds_BF_bd", "pull_rds_BF_bd", *pull_BF_bd);
  for (int i = 0; i < channels; i++) { cout << i << endl;
    for (int j = 0; j < channels_bdt; j++) { cout << j << endl;
      residual_bs[i][j] = new RooRealVar(name("residual_bs", i, j), "residual_bs", -20., 20.);
      residual_bd[i][j] = new RooRealVar(name("residual_bd", i, j), "residual_bd", -20., 20.);
      pull_bs[i][j] = new RooRealVar(name("pull_bs", i, j), "pull_bs", -8., 8.);
      pull_bd[i][j] = new RooRealVar(name("pull_bd", i, j), "pull_bd", -8., 8.);
      pull_rare[i][j] = new RooRealVar(name("pull_rare", i, j), "pull_rare", -8., 8.);
      pull_comb[i][j] = new RooRealVar(name("pull_comb", i, j), "pull_comb", -8., 8.);
      residual_rds_bs[i][j] = new RooDataSet(name("residual_rds_bs", i, j), "residual_rds_bs", *residual_bs[i][j]);
      residual_rds_bd[i][j] = new RooDataSet(name("residual_rds_bd", i, j), "residual_rds_bd", *residual_bd[i][j]);
      pull_rds_bs[i][j] = new RooDataSet(name("pull_rds_bs", i, j), "pull_rds_bs", *pull_bs[i][j]);
      pull_rds_bd[i][j] = new RooDataSet(name("pull_rds_bd", i, j), "pull_rds_bd", *pull_bd[i][j]);
      pull_rds_rare[i][j] = new RooDataSet(name("pull_rds_rare", i, j), "pull_rds_rare", *pull_rare[i][j]);
      pull_rds_comb[i][j] = new RooDataSet(name("pull_rds_comb", i, j), "pull_rds_comb", *pull_comb[i][j]);
      correlation_h[i][j] = new TH1D(name("correlation_h", i, j), name("correlation_h", i, j), 100, -1., 1.);
      bs_mean_h[i][j] = new TH1D(name("bs_mean_h", i, j), name("bs_mean_h", i, j), 100, 0., 100);
      bs_sigma_h[i][j] = new TH1D(name("bs_sigma_h", i, j), name("bs_sigma_h", i, j), 100, 0., 100);
      bd_mean_h[i][j] = new TH1D(name("bd_mean_h", i, j), name("bd_mean_h", i, j), 100, 0., 100);
      bd_sigma_h[i][j] = new TH1D(name("bd_sigma_h", i, j), name("bd_sigma_h", i, j), 100, 0., 100);
      rare_mean_h[i][j] = new TH1D(name("rare_mean_h", i, j), name("rare_mean_h", i, j), 100, 0., 100);
      comb_mean_h[i][j] = new TH1D(name("comb_mean_h", i, j), name("comb_mean_h", i, j), 100, 0., 100);
      rare_sigma_h[i][j] = new TH1D(name("rare_sigma_h", i, j), name("rare_sigma_h", i, j), 100, 0., 100);
      comb_sigma_h[i][j] = new TH1D(name("comb_sigma_h", i, j), name("comb_sigma_h", i, j), 100, 0., 100);
      corr_Nbs_Nbd_vs_N_bs_h[i][j] = new TH2D(Form("corr_Nbs_Nbd_vs_N_bs_%d_%d_h", i, j), "corr_Nbs_Nbd_vs_N_bs_h", 30, 0., 30., 100, -1., 1.);
    }
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
  vector < vector < Double_t > > estimate_bs_formula(channels, vector < Double_t > (channels_bdt, 0.));
  vector < vector < Double_t > > estimate_bd_formula(channels, vector < Double_t > (channels_bdt, 0.));
  if (BF_ > 0) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < channels_bdt; j++) {
        estimate_bs_formula[i][j] = ws_->function(name("N_bs_formula", i, j))->getVal();
        if (BF_ > 1) estimate_bd_formula[i][j] = ws_->function(name("N_bd_formula", i, j))->getVal();
      }
    }
  }

  for (int k = 1; k <= NExp; k++) {
    if (k%100 == 0) cout << "Exp # " << k << " / " << NExp << endl;
    pdf_toy_ = pdf_toy;
    RooWorkspace* ws_temp = (RooWorkspace*)ws_->Clone("ws_temp");
    RooDataSet* data = new RooDataSet("data", "data", RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes")));
    double printlevel = -1;
    if (k == 1) printlevel = 1;

    if (!simul_bdt_) { /// 2D
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
    }
    else { /// 1D
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < channels_bdt; j++) {
          ws_temp->var(name("N_bs", i, j))->setVal((int)estimate2D_bs[i][j]);
          if (!SM_ && !bd_constr_) ws_temp->var(name("N_bd", i, j))->setVal((int)estimate2D_bd[i][j]);
          ws_temp->var(name("N_rare", i, j))->setVal((int)estimate2D_rare[i][j]);
          ws_temp->var(name("N_comb", i, j))->setVal((int)estimate2D_comb[i][j]);
          ws_temp->var(name("N_bs", i, j))->setConstant(kFALSE);
          if (!SM_ && !bd_constr_) ws_temp->var(name("N_bd", i, j))->setConstant(kFALSE);
          ws_temp->var(name("N_rare", i, j))->setConstant(kFALSE);
          ws_temp->var(name("N_comb", i, j))->setConstant(kFALSE);
        }
      }
    }
    if (channels == 1) { /// 1D non simul 1 channel
      if (pdf_toy_ == "total") { /// 1D non simul total channel
        bd_b = true;
        /*if (BF_==0) */data = ws_temp->pdf("pdf_ext_total")->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes")), Extended(1));
//        data = ws_temp->pdf("pdf_ext_total_test")->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes")), Extended(1));
      }
      else { /// /// 1D non simul non-total channel
        size_t found;
        found = pdf_toy_.find("bd");
        if (found != string::npos) bd_b = true;
        data = ws_temp->pdf( ("pdf_ext_" + pdf_toy_).c_str())->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes")), Extended(1));
      }
    }
    else if (!simul_bdt_) { /// 1D simul
      bd_b = true;
      RooCategory* cat = (RooCategory*)ws_temp->obj("etacat");
      data->addColumn(*cat);
      RooArgSet vars(*ws_temp->var("Mass"), *ws_temp->var("bdt"), *ws_temp->var("MassRes"), *cat);
      /*if (BF_==0) */data = ws_temp->pdf("pdf_ext_simul")->generate(vars, Extended(1));
//      else data = ws_temp->pdf("pdf_ext_simul_test")->generate(vars, Extended(1));
    }
    else { /// 2D
      data->addColumn(*channels_cat);
      data->addColumn(*bdt_cat);
      vector < vector <RooDataSet*> > data_i(channels, vector <RooDataSet* > (channels_bdt));
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < channels_bdt; j++) {
          ws_->var(name("N_bs", i, j))->setVal(estimate2D_bs[i][j]);
          if (!SM_ && !bd_constr_) ws_->var(name("N_bd", i, j))->setVal(estimate2D_bd[i][j]);
          else if (bd_constr_) ws_->var("bd_over_bs")->setVal(estimate2D_bd[i][j]/estimate2D_bd[i][j]);
          ws_->var(name("N_rare", i, j))->setVal(estimate2D_rare[i][j]);
          ws_->var(name("N_comb", i, j))->setVal(estimate2D_comb[i][j]);
          data_i[i][j] = ws_->pdf(name("pdf_ext_total", i, j))->generate(RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->var("bdt")));
          channels_cat->setIndex(i);
          bdt_cat->setIndex(j);
          data_i[i][j]->addColumn(*channels_cat);
          data_i[i][j]->addColumn(*bdt_cat);
          data->append(*data_i[i][j]);
        }
      }
      data->SetName("global_data");
    }
    do_bias(ws_temp);
    ///////////
    //////
    RFR = pdf_toyMC::fit_pdf(pdf_test_, data, printlevel, ws_temp);
    //////
    //////////
    pdf_name = pdf_toy_;
    rds_ = data;

    if (!simul_) {
      if (k == 1) {
        ws_ = ws_temp;
        print("_first", ws_temp);
      }
      if (k == NExp) {
        ws_ = ws_temp;
        print("_last", ws_temp);
      }
    }
    /// pull
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < channels_bdt; j++) {
        if (BF_==0) bs_mean_h[i][j]->Fill(ws_temp->var(name("N_bs", i, j))->getVal());
        else bs_mean_h[i][j]->Fill(ws_temp->function(name("N_bs_formula", i, j))->getVal());
        if (!(SM_ || bd_constr_ || BF_ == 2)) bd_mean_h[i][j]->Fill(ws_temp->var(name("N_bd", i, j))->getVal());
        else if (BF_ == 2) bd_mean_h[i][j]->Fill(ws_temp->function(name("N_bd_formula", i, j))->getVal());
        rare_mean_h[i][j]->Fill(ws_temp->var(name("N_rare", i, j))->getVal());
        comb_mean_h[i][j]->Fill(ws_temp->var(name("N_comb", i, j))->getVal());
        if (BF_==0) bs_sigma_h[i][j]->Fill(ws_temp->var(name("N_bs", i, j))->getError());
        if (!(SM_ || bd_constr_ || BF_ == 2)) bd_sigma_h[i][j]->Fill(ws_temp->var(name("N_bd", i, j))->getError());
        rare_sigma_h[i][j]->Fill(ws_temp->var(name("N_rare", i, j))->getError());
        comb_sigma_h[i][j]->Fill(ws_temp->var(name("N_comb", i, j))->getError());
        // PULL
        double bs_pull, bd_pull, rare_pull, comb_pull, bs_residual, bd_residual;
        if (!simul_bdt_) {
          if (BF_==0) bs_pull = (ws_temp->var(name("N_bs", i, j))->getVal() - estimate_bs[i]) / ws_temp->var(name("N_bs", i, j))->getError();
          else {
            bs_residual = (ws_temp->function(name("N_bs_formula", i, j))->getVal() - estimate_bs_formula[i][0]);
//            cout << ws_temp->function(name("N_bs_formula", i, j))->getVal()  << "  " << estimate_bs[i] << endl;
          }
          if (!(SM_ || bd_constr_ || BF_ == 2)) bd_pull = (ws_temp->var(name("N_bd", i, j))->getVal() - estimate_bd[i]) / ws_temp->var(name("N_bd", i, j))->getError();
          else if (BF_ == 2) bd_residual = (ws_temp->function(name("N_bd_formula", i, j))->getVal() - estimate_bd_formula[i][0]);
          rare_pull = (ws_temp->var(name("N_rare", i, j))->getVal() - estimate_rare[i]) / ws_temp->var(name("N_rare", i, j))->getError();
          comb_pull = (ws_temp->var(name("N_comb", i, j))->getVal() - estimate_comb[i]) / ws_temp->var(name("N_comb", i, j))->getError();
        }
        else {
          if (BF_==0) bs_pull = (ws_temp->var(name("N_bs", i, j))->getVal() - estimate2D_bs[i][j]) / ws_temp->var(name("N_bs", i, j))->getError();
          else bs_residual = (ws_temp->function(name("N_bs_formula", i, j))->getVal() - estimate_bs_formula[i][j]);
          if (!(SM_ || bd_constr_|| BF_ == 2)) bd_pull = (ws_temp->var(name("N_bd", i, j))->getVal() - estimate2D_bd[i][j]) / ws_temp->var(name("N_bd", i, j))->getError();
          else if (BF_ == 2) bd_residual = (ws_temp->function(name("N_bd_formula", i, j))->getVal() - estimate_bd_formula[i][j]);
          rare_pull = (ws_temp->var(name("N_rare", i, j))->getVal() - estimate2D_rare[i][j]) / ws_temp->var(name("N_rare", i, j))->getError();
          comb_pull = (ws_temp->var(name("N_comb", i, j))->getVal() - estimate2D_comb[i][j]) / ws_temp->var(name("N_comb", i, j))->getError();
        }

        if (BF_==0) {
          pull_bs[i][j]->setVal(bs_pull);
          pull_rds_bs[i][j]->add(*pull_bs[i][j]);
        }
        else {
          residual_bs[i][j]->setVal(bs_residual);
          residual_rds_bs[i][j]->add(*residual_bs[i][j]);
        }
        if (!(SM_ || bd_constr_ || BF_==2)) {
          pull_bd[i][j]->setVal(bd_pull);
          pull_rds_bd[i][j]->add(*pull_bd[i][j]);
        }
        else if (BF_==2) {
          residual_bd[i][j]->setVal(bd_residual);
          residual_rds_bd[i][j]->add(*residual_bd[i][j]);
        }
        pull_rare[i][j]->setVal(rare_pull);
        pull_rds_rare[i][j]->add(*pull_rare[i][j]);
        pull_comb[i][j]->setVal(comb_pull);
        pull_rds_comb[i][j]->add(*pull_comb[i][j]);
        if (BF_==0) correlation_h[i][j]->Fill(RFR->correlation(*ws_temp->var(name("N_bs", i, j)), *ws_temp->var(name("N_bd", i, j))));
        else if (BF_==1) correlation_h[i][j]->Fill(RFR->correlation(*ws_temp->var("BF_bs"), *ws_temp->var(name("N_bd", i, j))));
        else correlation_h[i][j]->Fill(RFR->correlation(*ws_temp->var("BF_bs"), *ws_temp->var("BF_bd")));
      }
    }
    if (BF_>0) {
      BF_bs_mean_h->Fill(ws_temp->var("BF_bs")->getVal());
      BF_bs_sigma_h->Fill(ws_temp->var("BF_bs")->getError());
      double BF_pull = (ws_temp->var("BF_bs")->getVal() - Bs2MuMu_SM_BF_val) / ws_temp->var("BF_bs")->getError();
      pull_BF_bs->setVal(BF_pull);
      pull_rds_BF_bs->add(*pull_BF_bs);
    }
    if (BF_>1) {
      BF_bd_mean_h->Fill(ws_temp->var("BF_bd")->getVal());
      BF_bd_sigma_h->Fill(ws_temp->var("BF_bd")->getError());
      double BF_pull = (ws_temp->var("BF_bd")->getVal() - Bd2MuMu_SM_BF_val) / ws_temp->var("BF_bd")->getError();
      pull_BF_bd->setVal(BF_pull);
      pull_rds_BF_bd->add(*pull_BF_bd);
    }
    if (k == NExp) ws_temp->pdf(pdf_toy_.c_str())->Print();

    if (sign == 0) {
      sign_h->Fill(pdf_toyMC::sig_hand(data, printlevel, ws_temp));
    }

    delete data;
    delete RFR;
  }

  pdf_toy_ = pdf_toy;

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < channels_bdt; j++) {
      if (BF_==0) fit_pulls(pull_bs[i][j], pull_rds_bs[i][j], i, j);
      if (BF_>0) fit_pulls(residual_bs[i][j], residual_rds_bs[i][j], i, j);
      if (!(SM_ || bd_constr_ || BF_==2)) fit_pulls(pull_bd[i][j], pull_rds_bd[i][j], i, j);
      else if (BF_==2) fit_pulls(residual_bd[i][j], residual_rds_bd[i][j], i, j);
      fit_pulls(pull_rare[i][j], pull_rds_rare[i][j], i, j);
      fit_pulls(pull_comb[i][j], pull_rds_comb[i][j], i, j);

      print_histos(bs_mean_h[i][j], i, j);
      if (!(SM_ || bd_constr_)) print_histos(bd_mean_h[i][j], i, j);
      print_histos(rare_mean_h[i][j], i, j);
      print_histos(comb_mean_h[i][j], i, j);

      if (BF_==0) print_histos(bs_sigma_h[i][j], i, j);
      if (!(SM_ || bd_constr_)) print_histos(bd_sigma_h[i][j], i, j);
      print_histos(rare_sigma_h[i][j], i, j);
      print_histos(comb_sigma_h[i][j], i, j);
    }
  }
  if (BF_ > 0) {
    fit_pulls(pull_BF_bs, pull_rds_BF_bs, 0, 0);
    print_histos(BF_bs_mean_h, 0, 0);
    print_histos(BF_bs_sigma_h, 0, 0);
  }
  if (BF_ > 1) {
    fit_pulls(pull_BF_bd, pull_rds_BF_bd, 0, 0);
    print_histos(BF_bd_mean_h, 0, 0);
    print_histos(BF_bd_sigma_h, 0, 0);
  }
  if (!SM_ && bd_b && !simul_) {
    TCanvas* corr_c = new TCanvas("corr_c", "corr_c", 600*channels, 600*channels_bdt);
    corr_c->Divide(channels, channels_bdt);
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < channels_bdt; j++) {
        corr_c->cd(j + 1);
        correlation_h[i][j]->Draw();
      }
    }
    corr_c->Print((get_address("corr") + ".gif").c_str());
    corr_c->Print((get_address("corr") + ".pdf").c_str());
    delete corr_c;
  }

  if (sign == 0) {
    TCanvas* sign_c = new TCanvas("sign_c", "sign_c", 600, 600);
    sign_h->Draw("e");
    sign_c->Print((get_address("sign", "", false) + ".gif").c_str());
    sign_c->Print((get_address("sign", "", false) + ".pdf").c_str());
    delete sign_c;
  }
  return;
}

void pdf_toyMC::do_bias(RooWorkspace* ws) {

  if (bias_.compare("no")) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < channels_bdt; j++) {
        if (!simul_ && !simul_bdt_) i = ch_i_;
        size_t found;
        found = bias_.find("tau");
        if (found != string::npos) {
          double value =  ws->var(name("tau", i, j))->getVal();
          found = bias_.find("+");
          if (found != string::npos) {
            double error = ws->var(name("tau", i, j))->getErrorHi();
            ws->var(name("tau", i, j))->setVal(value + error);
          }
          else {
            double error = ws->var(name("tau", i, j))->getErrorLo();
            ws->var(name("tau", i, j))->setVal(value - error);
          }
        }
        found = bias_.find("c");
        if (found != string::npos) {
          double value =  ws->var(name("C0", i, j))->getVal();
          found = bias_.find("+");
          if (found != string::npos) {
            double error = ws->var(name("C0", i, j))->getErrorHi();
            ws->var(name("C0", i, j))->setVal(value + error);
          }
          else {
            double error = ws->var(name("C0", i, j))->getErrorLo();
            ws->var(name("C0", i, j))->setVal(value - error);
          }
        }
        found = bias_.find("sig");
        if (found != string::npos) {
          found = bias_.find("+");
          if (found != string::npos) {
            ws->var(name("Mean_bs", i, j))->setVal(ws->var(name("Mean_bs", i, j))->getVal() + ws->var(name("Mean_bs", i, j))->getErrorHi());
            ws->var(name("Sigma_bs", i, j))->setVal(ws->var(name("Sigma_bs", i, j))->getVal() + ws->var(name("Sigma_bs", i, j))->getErrorHi());
            ws->var(name("Alpha_bs", i, j))->setVal(ws->var(name("Alpha_bs", i, j))->getVal() + ws->var(name("Alpha_bs", i, j))->getErrorHi());
            ws->var(name("Enne_bs", i, j))->setVal(ws->var(name("Enne_bs", i, j))->getVal() + ws->var(name("Enne_bs", i, j))->getErrorHi());
            if (!pee) {
              ws->var(name("Sigma2_bs", i, j))->setVal(ws->var(name("Sigma2_bs", i, j))->getVal() + ws->var(name("Sigma2_bs", i, j))->getErrorHi());
              ws->var(name("CoeffGauss_bs", i, j))->setVal(ws->var(name("CoeffGauss_bs", i, j))->getVal() + ws->var(name("CoeffGauss_bs", i, j))->getErrorHi());
            }
          }
          else {
            ws->var(name("Mean_bs", i, j))->setVal(ws->var(name("Mean_bs", i, j))->getVal() - ws->var(name("Mean_bs", i, j))->getErrorLo());
            ws->var(name("Sigma_bs", i, j))->setVal(ws->var(name("Sigma_bs", i, j))->getVal() - ws->var(name("Sigma_bs", i, j))->getErrorLo());
            ws->var(name("Alpha_bs", i, j))->setVal(ws->var(name("Alpha_bs", i, j))->getVal() - ws->var(name("Alpha_bs", i, j))->getErrorLo());
            ws->var(name("Enne_bs", i, j))->setVal(ws->var(name("Enne_bs", i, j))->getVal() - ws->var(name("Enne_bs", i, j))->getErrorLo());
            if (!pee) {
              ws->var(name("Sigma2_bs", i, j))->setVal(ws->var(name("Sigma2_bs", i, j))->getVal() - ws->var(name("Sigma2_bs", i, j))->getErrorLo());
              ws->var(name("CoeffGauss_bs", i, j))->setVal(ws->var(name("CoeffGauss_bs", i, j))->getVal() - ws->var(name("CoeffGauss_bs", i, j))->getErrorLo());
            }
          }
        }
      }
      if (!simul_ && !simul_bdt_) break;
    }
  }
}

RooFitResult* pdf_toyMC::fit_pdf(string pdf, RooAbsData* data, int printlevel, RooWorkspace* ws) {

  pdf_toy_ = "pdf_ext_" + pdf;
  RooFitResult* result;
  if (printlevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  if (!pee) result = ws->pdf(pdf_toy_.c_str())->fitTo(*data, Extended(true), SumW2Error(0), PrintLevel(printlevel), PrintEvalErrors(-1), Save(kTRUE), NumCPU(2));
  else result = ws->pdf(pdf_toy_.c_str())->fitTo(*data, ConditionalObservables(*ws->var("MassRes")), Extended(true), SumW2Error(0), PrintLevel(printlevel), PrintEvalErrors(-1), Save(kTRUE), NumCPU(2));
  if (printlevel < 0) RooMsgService::instance().cleanup();
  return result;
}

void pdf_toyMC::fit_pulls(RooRealVar* pull, RooDataSet* rds, int i, int j) {

  RooRealVar* mean_bs = new RooRealVar("mean_bs", "mean_bs", -5., 5.);
  RooRealVar* sigma_bs = new RooRealVar("sigma_bs", "sigma_bs", 0.001, 5.);
  RooGaussian* gauss_bs = new RooGaussian("gauss_bs", "gauss_bs", *pull, *mean_bs, *sigma_bs);
  gauss_bs->fitTo(*rds);
  
  RooPlot *rp_bs = pull->frame();
  rds->plotOn(rp_bs, Binning(40));
  gauss_bs->plotOn(rp_bs, LineColor(kBlue));
  gauss_bs->paramOn(rp_bs, Layout(0.66, 0.9, 0.9));
  TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 600, 600);
  rp_bs->Draw();
  channel = simul_ ? i : ch_i_;
  channel_bdt = simul_bdt_ ? j : 0;
  canvas_bs->Print((get_address(pull->GetName(), pdf_toy_) + ".gif").c_str());
  canvas_bs->Print((get_address(pull->GetName(), pdf_toy_) + ".pdf").c_str());
  delete rp_bs;
  delete canvas_bs;
}

void pdf_toyMC::print_histos(TH1D* histos, int i, int j) {

  cout << "etacat " << i << "  bdtcat " << j << " " << histos->GetName() << " mean = " << histos->GetMean() << endl;
  TCanvas* N_mean_c = new TCanvas("N_mean_c", "N_mean_c", 600, 600);
  histos->Draw("e");
  channel = simul_ ? i : ch_i_;
  channel_bdt = simul_bdt_ ? j : 0;
  N_mean_c->Print((get_address(histos->GetName(), pdf_toy_) + ".gif").c_str());
  N_mean_c->Print((get_address(histos->GetName(), pdf_toy_) + ".pdf").c_str());
  delete N_mean_c;
}

Double_t pdf_toyMC::sig_hand(RooAbsData* data, int printlevel, RooWorkspace* ws_temp) {
  Double_t minNLL = RFR->minNll();
  if (BF_==0) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < channels_bdt; j++) {
        ws_temp->var(name("N_bs", i, j))->setVal(0);
        ws_temp->var(name("N_bs", i, j))->setConstant(1);
        /*      if ( !(SM_ || bd_constr_)) {
        ws_temp->var(name("N_bd", i, j))->setVal(0);
        ws_temp->var(name("N_bd", i, j))->setConstant(1);
      }
      else */if (bd_constr_) {
          ws_temp->var("Bd_over_Bs")->setVal(0);
          ws_temp->var("Bd_over_Bs")->setConstant(1);
        }
      }
    }
  }
  else {
    ws_temp->var("BF_bs")->setVal(0);
    ws_temp->var("BF_bs")->setConstant(1);
  }

  RooFitResult * rfr_H0 = pdf_toyMC::fit_pdf(pdf_test_, data, printlevel, ws_temp);

  if (BF_==0) {
    if (!simul_bdt_) {
      for (int j = 0; j < channels; j++) {
        ws_temp->var(name("N_bs", j))->setVal(estimate_bs[j]);
        ws_temp->var(name("N_bs", j))->setConstant(0);
        /*      if ( !(SM_ || bd_constr_)) {
        ws_temp->var(name("N_bd", j))->setVal(estimate_bd[j]);
        ws_temp->var(name("N_bd", j))->setConstant(0);
      }
      else */if (bd_constr_) {
          ws_temp->var("Bd_over_Bs")->setVal(estimate_bd[j]/estimate_bs[j]);
          ws_temp->var("Bd_over_Bs")->setConstant(0);
        }
      }
    }
    else {
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < channels_bdt; j++) {
          ws_temp->var(name("N_bs", i, j))->setVal(estimate2D_bs[i][j]);
          ws_temp->var(name("N_bs", i, j))->setConstant(0);
          /*        if ( !(SM_ || bd_constr_)) {
          ws_temp->var(name("N_bd", i, j))->setVal(estimate2D_bd[i][j]);
          ws_temp->var(name("N_bd", i, j))->setConstant(0);
        }
        else */if (bd_constr_) {
            ws_temp->var("Bd_over_Bs")->setVal(estimate2D_bd[i][j]/estimate2D_bs[i][j]);
            ws_temp->var("Bd_over_Bs")->setConstant(0);
          }
        }
      }
    }
  }
  else {
    ws_temp->var("BF_bs")->setVal(Bs2MuMu_SM_BF_val);
    ws_temp->var("BF_bs")->setConstant(0);
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

void pdf_toyMC::mcstudy(int NExp, string pdf_toy, string test_pdf) {

  if (!simul_bdt_) {
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
  }
  else {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < channels_bdt; j++) {
        ws_->var(name("N_bs", i, j))->setVal(estimate2D_bs[i][j]);
        ws_->var(name("N_rare", i, j))->setVal(estimate2D_rare[i][j]);
        ws_->var(name("N_comb", i, j))->setVal(estimate2D_comb[i][j]);
        if (!SM_ && !bd_constr_) ws_->var(name("N_bd", i, j))->setVal(estimate2D_bd[i][j]);
      }
    }
    if (bd_constr_) {
      double ratio = (double) estimate2D_bd[0][0] / estimate2D_bs[0][0]; // it's the same in every channel
      ws_->var("Bd_over_Bs")->setVal(ratio);
    }
  }

  pdf_name = "pdf_ext_total";
  if (pdf_toy != "total") {
    pdf_name = "pdf_ext_" + define_pdf_sum(pdf_toy);
  }
  if (simul_) pdf_name = "pdf_ext_simul";

  pdf_test_= "pdf_ext_total";
  if (test_pdf != "total") {
    pdf_test_ = "pdf_ext_" + define_pdf_sum(test_pdf);
  }
  if (simul_) pdf_test_ = "pdf_ext_simul";

  RooArgSet obsv(*ws_->var("Mass"), *ws_->var("bdt"), "obsv");
  if (simul_ && !simul_bdt_) obsv.add(*ws_->cat("etacat"));
  if (simul_bdt_) {
    cout << "RooMCStudy seems not to work with RooSuperCategory" << endl;
    obsv.add(*ws_->cat("bdtcat"));
    exit(0);
  }
  if (pee) obsv.add(*ws_->var("MassRes"));

  RooMCStudy * mcstudy;
  if (bias_.compare("no")) { /// bias
    RooWorkspace* fitws = (RooWorkspace*)ws_->Clone("fitws");
    do_bias(fitws);
    if (!pee) mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv, FitModel(*fitws->pdf(pdf_name.c_str())), Binned(kFALSE), Extended(kTRUE), FitOptions(Save(kTRUE)), Silence());
    else mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv, FitModel(*fitws->pdf(pdf_name.c_str())), Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(Save(kTRUE), ConditionalObservables(*ws_->var("MassRes"))));
  }
  else if (pdf_test_ == pdf_name) { /// same pdf
    if (!pee) mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv, Binned(kFALSE), Extended(kTRUE), FitOptions(Save(kTRUE)), Silence());
    else mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv, Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(Save(kTRUE), ConditionalObservables(*ws_->var("MassRes"))));
  }
  else { /// different pdf
    if (!pee) mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv, Binned(kFALSE), FitModel(*ws_->pdf(pdf_test_.c_str())), Extended(kTRUE), FitOptions(Save(kTRUE)), Silence());
    else mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), obsv, FitModel(*ws_->pdf(pdf_test_.c_str())), Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(Save(kTRUE), ConditionalObservables(*ws_->var("MassRes"))));
  }

  vector <vector <RooDLLSignificanceMCSModule*> > sigModule(channels, vector <RooDLLSignificanceMCSModule*> (channels_bdt));
  if (BF_==0) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < channels_bdt; j++) {
        sigModule[i][j] = new RooDLLSignificanceMCSModule(*ws_->var(name("N_bs", i, j)), 0);
        mcstudy->addModule(*sigModule[i][j]);
      }
    }
  }
  else {
    sigModule[0][0] = new RooDLLSignificanceMCSModule(*ws_->var("BF_bs"), 0);
    mcstudy->addModule(*sigModule[0][0]);
  }

  /////////
  ////
  mcstudy->generateAndFit(NExp, 0, kTRUE);
  ////
  ////////

  for (int k = 0; k < 4; k++) {
    if ((SM_ || bd_constr_) && k == 1) continue;
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < channels_bdt; j++) {
        if (k == 0 && BF_ > 0 && (i > 0 || j > 0)) continue;
        if (k == 1 && BF_ > 1 && (i > 0 || j > 0)) continue;
        string name_;
        if (k == 0 && BF_ == 1) {
          name_ = "BF_bs";
          source[0] = "BF_bs";
        }
        else if (k == 1 && BF_ == 2) {
          name_ = "BF_bd";
          source[0] = "BF_bd";
        }
        else name_ = name("N_" + source[k], i, j);
        RooPlot* frame1_bs = mcstudy->plotParam(*ws_->var(name_.c_str()), Bins(20)) ;
        RooPlot* frame2_bs = mcstudy->plotError(*ws_->var(name_.c_str()), Bins(20)) ;
        RooPlot* frame3_bs = mcstudy->plotPull(*ws_->var(name_.c_str()), Bins(20), FitGauss(kTRUE)) ;
        TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 600, 600);
        frame1_bs->Draw();
        TPaveText* stats1 = new TPaveText(0.57, 0.66, 0.9, 0.9, "NDCR");
        ostringstream entries;
        entries << "entries = " << mcstudy->fitParDataSet().sumEntries();
        ostringstream mean;
        mean << "mean = " << mcstudy->fitParDataSet().mean(*ws_->var(name_.c_str()));
        ostringstream sigma;
        sigma << "sigma = " << mcstudy->fitParDataSet().sigma(*ws_->var(name_.c_str()));
        stats1->AddText(entries.str().c_str());
        stats1->AddText(mean.str().c_str());
        stats1->AddText(sigma.str().c_str());
        stats1->Draw();
        channel = simul_ ? i : ch_i_;
        channel_bdt = simul_bdt_ ? j : 0;
        canvas1->Print((get_address("RooMCStudy_mean", source[k], true) + ".gif").c_str());
        canvas1->Print((get_address("RooMCStudy_mean", source[k], true) + ".pdf").c_str());
        delete frame1_bs;
        delete canvas1;
        TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 600, 600);
        frame2_bs->Draw();
        TPaveText* stats2 = new TPaveText(0.57, 0.66, 0.9, 0.9, "NDCR");
        entries.str("");
        entries << "entries = " << mcstudy->fitParDataSet().sumEntries();
        mean.str("");
        mean << "mean = " << mcstudy->fitParDataSet().mean(*ws_->var(name_.c_str()));
        sigma.str("");
        sigma << "sigma = " << mcstudy->fitParDataSet().sigma(*ws_->var(name_.c_str()));
        stats2->AddText(entries.str().c_str());
        stats2->AddText(mean.str().c_str());
        stats2->AddText(sigma.str().c_str());
        stats2->Draw();
        canvas2->Print((get_address("RooMCStudy_sigma", source[k], true) + ".gif").c_str());
        canvas2->Print((get_address("RooMCStudy_sigma", source[k], true) + ".pdf").c_str());
        delete frame2_bs;
        delete canvas2;
        TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 600, 600);
        frame3_bs->Draw();
        canvas3->Print((get_address("RooMCStudy_pull", source[k], true) + ".gif").c_str());
        canvas3->Print((get_address("RooMCStudy_pull", source[k], true) + ".pdf").c_str());
        delete frame3_bs;
        delete canvas3;
      }
    }
  }
  TCanvas* sig_c = new TCanvas("sig_c", "sig_c", 600*channels, 600*channels_bdt);
  if (BF_) sig_c->SetCanvasSize(600, 600);
  else sig_c->Divide(channels, channels_bdt);
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < channels_bdt; j++) {
      if (BF_ > 0 && (i > 0 || j > 0)) continue;
      string name_(name("significance_nullhypo_N_bs", i, j));
      if (BF_) name_ = "significance_nullhypo_BF_bs";
      sig_c->cd(i+1);
      TH1* sig_h = mcstudy->fitParDataSet().createHistogram(name_.c_str());
      sig_h->Rebin();
      sig_h->GetYaxis()->SetTitleOffset(1.3);
      sig_h->GetXaxis()->SetTitle("Gaussian significance of #Delta(-ln L) w.r.t. bkg only");
      sig_h->SetTitle(0);
      cout << ">3 sigma fraction is " << sig_h->Integral(sig_h->FindBin(3.), sig_h->FindBin(sig_h->GetNbinsX())) / sig_h->Integral() << endl;
      sig_h->Draw();
      sig_c->Update();
      sig_h->SetStats(0);
      TLatex* t = new TLatex();
      t->SetNDC();
      t->SetTextFont(sig_h->GetXaxis()->GetTitleFont());
      t->SetTextSize(sig_h->GetXaxis()->GetTitleSize());
      char leg[256];
      sprintf(leg, "CMS simulation #sqrt{s} = 7 TeV, L = 30 fb^{-1}, assumed BF = 3.5 #times 10^{-9}");
      //t->DrawLatex(0.05, 0.93, leg);
      sprintf(leg, "mean = %.2f", sig_h->GetMean());
      t->DrawLatex(0.7, 0.85, leg);
      sprintf(leg, "RMS = %.2f", sig_h->GetRMS());
      t->DrawLatex(0.7, 0.80, leg);
      sprintf(leg, "|#eta_{#mu}|<1.4");
      t->DrawLatex(0.7, 0.75, leg);
    }
  }
  sig_c->Print((get_address("RooMCStudy_sig", "bs") + ".gif").c_str());
  sig_c->Print((get_address("RooMCStudy_sig", "bs") + ".pdf").c_str());
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
    TCanvas* cetad = new TCanvas("cetad", "cetad", 1200, 600);
    cetad->Divide(bdt_fit_ ? 3 : 2);
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
    cetad->Print((get_address(pdf_name, (string)("MassEta"+output)) + ".gif").c_str());
    cetad->Print((get_address(pdf_name, (string)("MassEta"+output)) + ".pdf").c_str());
    delete cetad;
    delete mass_eta_h;
  }
  if(!no_legend) {
    ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
    if(bdt_fit_) ws_->pdf(pdf_name.c_str())->paramOn(rp_bdt, Layout(0.50, 0.9, 0.9));
  }

  /// components
  RooArgSet * set = ws_->pdf(pdf_name.c_str())->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  int i = 0;
  while((var_Obj = it->Next())/* && !pee*/){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (i > 11) i = 0;
      size_t found1 = pdf_name.find("total");
      if (found1 == string::npos) {
        ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2));
        if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2));
      }
      else {
        if (name=="pdf_bs" || name=="pdf_bd" || name=="pdf_rare" || name=="pdf_comb") {
          ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2));
          if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2));
        }
      }
      i++;
    }
  }

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  c->Print((get_address(pdf_name, output) + ".gif").c_str());
  c->Print((get_address(pdf_name, output) + ".pdf").c_str());
  delete rp;
  delete c;

  if (bdt_fit_) {
    TCanvas* d = new TCanvas("d", "d", 600, 600);
    rp_bdt->Draw();
    d->Print((get_address("BDT", (string)(pdf_name+output)) + ".gif").c_str());
    d->Print((get_address("BDT", (string)(pdf_name+output)) + ".pdf").c_str());
    delete d;
  }
  delete rp_bdt;
  return;
}
