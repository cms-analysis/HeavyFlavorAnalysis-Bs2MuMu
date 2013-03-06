/*
 * File:   pdf_toyMC.cpp
 * Author: lucamartini
 * 
 * Created on 26 aprile 2012, 11.03
 */

#include "pdf_toyMC.h"

pdf_toyMC::pdf_toyMC(bool print, string input_estimates, string range, int BF, bool SM, bool bd_constr, int simul, int simulbdt, int simulall, bool pee_, bool bdt_fit, string ch_s, int sig, bool asimov, bool syste, bool randomsyste, bool rare_constr_, int nexp, bool bd, string years, string bias): pdf_fitData( print, input_estimates, range, BF, SM, bd_constr, simul, simulbdt, simulall, pee_, bdt_fit, ch_s, sig, asimov, syste, randomsyste, rare_constr_, nexp, bd, years) {
  cout << "pdf_toyMC constructor" << endl;
  bias_ = bias;

  residual_bs.resize(channels);
  residual_bd.resize(channels);
  pull_bs.resize(channels);
  pull_bd.resize(channels);
  pull_semi.resize(channels);
  pull_comb.resize(channels);
  residual_rds_bs.resize(channels);
  residual_rds_bd.resize(channels);
  pull_rds_bs.resize(channels);
  pull_rds_bd.resize(channels);
  pull_rds_semi.resize(channels);
  pull_rds_comb.resize(channels);

  for (int k = 0; k < channels; k++) {
    residual_bs[k].resize(channels_bdt);
    residual_bd[k].resize(channels_bdt);
    pull_bs[k].resize(channels_bdt);
    pull_bd[k].resize(channels_bdt);
    pull_semi[k].resize(channels_bdt);
    pull_comb[k].resize(channels_bdt);
    residual_rds_bs[k].resize(channels_bdt);
    residual_rds_bd[k].resize(channels_bdt);
    pull_rds_bs[k].resize(channels_bdt);
    pull_rds_bd[k].resize(channels_bdt);
    pull_rds_semi[k].resize(channels_bdt);
    pull_rds_comb[k].resize(channels_bdt);
  }

  pull_BF_bs = new RooRealVar("pull_BF_bs", "pull_BF_bs", -8., 8.);
  pull_rds_BF_bs = new RooDataSet("pull_rds_BF_bs", "pull_rds_BF_bs", *pull_BF_bs);
  pull_BF_bd = new RooRealVar("pull_BF_bd", "pull_BF_bd", -8., 8.);
  pull_rds_BF_bd = new RooDataSet("pull_rds_BF_bd", "pull_rds_BF_bd", *pull_BF_bd);

}

pdf_toyMC::~pdf_toyMC() {
  cout << "pdf_toyMC destructor" << endl;
}

void pdf_toyMC::generate(string pdf_toy, string pdf_test) {

  vector < vector  <TH1D*> > correlation_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > bs_mean_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > bd_mean_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > semi_mean_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > comb_mean_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > bs_sigma_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > bd_sigma_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > semi_sigma_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH1D*> > comb_sigma_h(channels, vector <TH1D*> (channels_bdt));
  vector < vector  <TH2D*> > corr_Nbs_Nbd_vs_N_bs_h(channels, vector <TH2D*> (channels_bdt));

  TH1D * sign_h = new TH1D("sign_h", "sign_h", 100, 0, 10);
  TH1D * BF_bs_mean_h = new TH1D("BF_bs_mean_h", "BF_bs_mean_h", 100, 1.e-10, 1.e-8);
  TH1D * BF_bs_sigma_h = new TH1D("BF_bs_sigma_h", "BF_bs_sigma_h", 100, 1.e-10, 1.e-8);
  TH1D * BF_bd_mean_h = new TH1D("BF_bd_mean_h", "BF_bd_mean_h", 100, 1.e-10, 1.e-8);
  TH1D * BF_bd_sigma_h = new TH1D("BF_bd_sigma_h", "BF_bd_sigma_h", 100, 1.e-10, 1.e-8);
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      residual_bs[i][j] = new RooRealVar(name("residual_bs", i, j), "residual_bs", -20., 20.);
      residual_bd[i][j] = new RooRealVar(name("residual_bd", i, j), "residual_bd", -20., 20.);
      pull_bs[i][j] = new RooRealVar(name("pull_bs", i, j), "pull_bs", -8., 8.);
      pull_bd[i][j] = new RooRealVar(name("pull_bd", i, j), "pull_bd", -8., 8.);
      pull_semi[i][j] = new RooRealVar(name("pull_semi", i, j), "pull_semi", -8., 8.);
      pull_comb[i][j] = new RooRealVar(name("pull_comb", i, j), "pull_comb", -8., 8.);
      residual_rds_bs[i][j] = new RooDataSet(name("residual_rds_bs", i, j), "residual_rds_bs", *residual_bs[i][j]);
      residual_rds_bd[i][j] = new RooDataSet(name("residual_rds_bd", i, j), "residual_rds_bd", *residual_bd[i][j]);
      pull_rds_bs[i][j] = new RooDataSet(name("pull_rds_bs", i, j), "pull_rds_bs", *pull_bs[i][j]);
      pull_rds_bd[i][j] = new RooDataSet(name("pull_rds_bd", i, j), "pull_rds_bd", *pull_bd[i][j]);
      pull_rds_semi[i][j] = new RooDataSet(name("pull_rds_semi", i, j), "pull_rds_semi", *pull_semi[i][j]);
      pull_rds_comb[i][j] = new RooDataSet(name("pull_rds_comb", i, j), "pull_rds_comb", *pull_comb[i][j]);
      correlation_h[i][j] = new TH1D(name("correlation_h", i, j), name("correlation_h", i, j), 100, -1., 1.);
      bs_mean_h[i][j] = new TH1D(name("bs_mean_h", i, j), name("bs_mean_h", i, j), 100, 0., 100);
      bs_sigma_h[i][j] = new TH1D(name("bs_sigma_h", i, j), name("bs_sigma_h", i, j), 100, 0., 100);
      bd_mean_h[i][j] = new TH1D(name("bd_mean_h", i, j), name("bd_mean_h", i, j), 100, 0., 100);
      bd_sigma_h[i][j] = new TH1D(name("bd_sigma_h", i, j), name("bd_sigma_h", i, j), 100, 0., 100);
      semi_mean_h[i][j] = new TH1D(name("semi_mean_h", i, j), name("semi_mean_h", i, j), 100, 0., 100);
      comb_mean_h[i][j] = new TH1D(name("comb_mean_h", i, j), name("comb_mean_h", i, j), 100, 0., 100);
      semi_sigma_h[i][j] = new TH1D(name("semi_sigma_h", i, j), name("semi_sigma_h", i, j), 100, 0., 100);
      comb_sigma_h[i][j] = new TH1D(name("comb_sigma_h", i, j), name("comb_sigma_h", i, j), 100, 0., 100);
      corr_Nbs_Nbd_vs_N_bs_h[i][j] = new TH2D(Form("corr_Nbs_Nbd_vs_N_bs_%d_%d_h", i, j), "corr_Nbs_Nbd_vs_N_bs_h", 30, 0., 30., 100, -1., 1.);
    }
  }
  bool bd_b = false, corr0 = false, corrneg = false;
  vector < vector < Double_t > > estimate_bs_formula(channels, vector < Double_t > (channels_bdt, 0.));
  vector < vector < Double_t > > estimate_bd_formula(channels, vector < Double_t > (channels_bdt, 0.));
  if (BF_ > 0) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < bdt_index_max(i); j++) {
        estimate_bs_formula[i][j] = ws_->function(name("N_bs_formula", i, j))->getVal();
        if (BF_ > 1) estimate_bd_formula[i][j] = ws_->function(name("N_bd_formula", i, j))->getVal();
      }
    }
  }
  cout << red_color_bold << "START OF EXPERIMENTS!" << default_console_color << endl;
  for (int k = 1; k <= NExp; k++) {
    if (k % 100 == 0) cout << "Exp # " << k << " / " << NExp << endl;
    double printlevel = -1;
    if (k == 1) printlevel = 1;
    RooWorkspace* ws_temp = (RooWorkspace*)ws_->Clone("ws_temp");
/// vars
    RooArgSet vars(*ws_temp->var("Mass"), *ws_temp->var("MassRes"), "vars");
    if (simul_ && !simul_bdt_ && !simul_all_) vars.add(*ws_temp->cat("etacat"));
    if (simul_ && simul_bdt_ && !simul_all_) {
    	vars.add(*ws_temp->cat("bdtcat"));
    	vars.add(*ws_temp->cat("etacat"));
    }
    if (simul_ && !simul_bdt_ && simul_all_) vars.add(*ws_temp->cat("allcat"));
    if (bdt_fit_) vars.add(*ws_temp->var("bdt"));
    RooDataSet* data = new RooDataSet("data", "data", vars);
/// setup
    if (!simul_bdt_ && !simul_all_) { /// simple 1D or 2D
      for (int j = 0; j < channels; j++) {
        ws_temp->var(name("N_bs", j))->setVal((int)estimate_bs[j]);
        if (!SM_ && !bd_constr_) ws_temp->var(name("N_bd", j))->setVal((int)estimate_bd[j]);
        ws_temp->var(name("N_comb", j))->setVal((int)estimate_comb[j]);
        ws_temp->var(name("N_bs", j))->setConstant(kFALSE);
        if (!SM_ && !bd_constr_) ws_temp->var(name("N_bd", j))->setConstant(kFALSE);
        ws_temp->var(name("N_comb", j))->setConstant(kFALSE);
        if (!rare_constr_) {
        	ws_temp->var(name("N_semi", j))->setVal((int)estimate_semi[j]);
        	ws_temp->var(name("N_semi", j))->setConstant(kFALSE);
        }
      }
    }
    else { /// 1D with 2 categories
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < bdt_index_max(i); j++) {
          ws_temp->var(name("N_bs", i, j))->setVal((int)estimate2D_bs[i][j]);
          if (!SM_ && !bd_constr_) ws_temp->var(name("N_bd", i, j))->setVal((int)estimate2D_bd[i][j]);
          ws_temp->var(name("N_comb", i, j))->setVal((int)estimate2D_comb[i][j]);
          ws_temp->var(name("N_bs", i, j))->setConstant(kFALSE);
          if (!SM_ && !bd_constr_) ws_temp->var(name("N_bd", i, j))->setConstant(kFALSE);
          ws_temp->var(name("N_comb", i, j))->setConstant(kFALSE);
          if (!rare_constr_) {
          	ws_temp->var(name("N_semi", i, j))->setConstant(kFALSE);
          	ws_temp->var(name("N_semi", i, j))->setVal((int)estimate2D_semi[i][j]);
          }
        }
      }
    }
    if (BF_ > 0) {
      ws_temp->var("BF_bs")->setVal(Bs2MuMu_SM_BF_val);
      if (BF_ > 1) {
        ws_temp->var("BF_bd")->setVal(Bd2MuMu_SM_BF_val);
      }
    }

/// generation
//    if (syst) randomize_constraints(ws_temp);
    if (!simul_bdt_ && !simul_all_) { /// simple 1D or 2D fit
      bd_b = true;
      data = ws_temp->pdf(pdfname.c_str())->generate(vars, Extended());
    }
    else {
      RooArgSet set(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->cat("etacat"), *ws_->cat("bdtcat"));
      if (simul_all_) set.add(*ws_->cat("allcat"));
      for (int i = 0; i < channels; i++) {
      	if (years_=="0" && i > 1) continue;
      	if (years_=="1" && i < 2) continue;
      	for (int j = 0; j < bdt_index_max(i); j++) {
          RooDataSet* data_i = ws_->pdf(name("pdf_ext_total", i, j))->generate(RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes")), Extended());
          channels_cat->setIndex(i);
          bdt_cat->setIndex(j);
          data_i->addColumn(*channels_cat);
          data_i->addColumn(*bdt_cat);
          if (simul_all_) {
          	all_cat->setIndex(super_index(i, j));
          	data_i->addColumn(*all_cat);
          }
          data->append(*data_i);
      	}
      }
    }
    data->SetName("global_data");

    do_bias(ws_temp);
    ///////////
    //////
    RFR = pdf_toyMC::fit_pdf(data, printlevel, ws_temp);
    //////
    //////////

    if (simul_) {
      if (k == 1) {
        print_each_channel("Mass", "_first", ws_temp, data);
      }
      if (k == NExp) {
        print_each_channel("Mass", "_last", ws_temp, data);
      }
    }
    /// pull
    for (int i = 0; i < channels; i++) {
    	if (years_=="0" && i > 1) continue;
    	if (years_=="1" && i < 2) continue;
      for (int j = 0; j < bdt_index_max(i); j++) {
        if (BF_ == 0) bs_mean_h[i][j]->Fill(ws_temp->var(name("N_bs", i, j))->getVal());
        else bs_mean_h[i][j]->Fill(ws_temp->function(name("N_bs_formula", i, j))->getVal());
        if (!(SM_ || bd_constr_ || BF_ == 2)) bd_mean_h[i][j]->Fill(ws_temp->var(name("N_bd", i, j))->getVal());
        else if (BF_ == 2) bd_mean_h[i][j]->Fill(ws_temp->function(name("N_bd_formula", i, j))->getVal());
        if (!rare_constr_) semi_mean_h[i][j]->Fill(ws_temp->var(name("N_semi", i, j))->getVal());
        comb_mean_h[i][j]->Fill(ws_temp->var(name("N_comb", i, j))->getVal());

        if (BF_ == 0) bs_sigma_h[i][j]->Fill(ws_temp->var(name("N_bs", i, j))->getError());
        if (!(SM_ || bd_constr_ || BF_ == 2)) bd_sigma_h[i][j]->Fill(ws_temp->var(name("N_bd", i, j))->getError());
        if (!rare_constr_) semi_sigma_h[i][j]->Fill(ws_temp->var(name("N_semi", i, j))->getError());
        comb_sigma_h[i][j]->Fill(ws_temp->var(name("N_comb", i, j))->getError());
        // PULL
        double bs_pull, bd_pull, semi_pull, comb_pull, bs_residual, bd_residual;
        if (!simul_bdt_ && !simul_all_) {
          if (BF_ == 0) bs_pull = (ws_temp->var(name("N_bs", i, j))->getVal() - estimate_bs[i]) / ws_temp->var(name("N_bs", i, j))->getError();
          bs_residual = (ws_temp->function(name("N_bs_formula", i, j))->getVal() - estimate_bs_formula[i][0]);
          if (!(SM_ || bd_constr_ || BF_ == 2)) bd_pull = (ws_temp->var(name("N_bd", i, j))->getVal() - estimate_bd[i]) / ws_temp->var(name("N_bd", i, j))->getError();
          else if (BF_ == 2) bd_residual = (ws_temp->function(name("N_bd_formula", i, j))->getVal() - estimate_bd_formula[i][0]);
          if (!rare_constr_) semi_pull = (ws_temp->var(name("N_semi", i, j))->getVal() - estimate_semi[i]) / ws_temp->var(name("N_semi", i, j))->getError();
          comb_pull = (ws_temp->var(name("N_comb", i, j))->getVal() - estimate_comb[i]) / ws_temp->var(name("N_comb", i, j))->getError();
        }
        else {
          if (BF_ == 0) bs_pull = (ws_temp->var(name("N_bs", i, j))->getVal() - estimate2D_bs[i][j]) / ws_temp->var(name("N_bs", i, j))->getError();
          else bs_residual = (ws_temp->function(name("N_bs_formula", i, j))->getVal() - estimate_bs_formula[i][j]);
          if (!(SM_ || bd_constr_|| BF_ == 2)) bd_pull = (ws_temp->var(name("N_bd", i, j))->getVal() - estimate2D_bd[i][j]) / ws_temp->var(name("N_bd", i, j))->getError();
          else if (BF_ == 2) bd_residual = (ws_temp->function(name("N_bd_formula", i, j))->getVal() - estimate_bd_formula[i][j]);
          if (!rare_constr_) semi_pull = (ws_temp->var(name("N_semi", i, j))->getVal() - estimate2D_semi[i][j]) / ws_temp->var(name("N_semi", i, j))->getError();
          comb_pull = (ws_temp->var(name("N_comb", i, j))->getVal() - estimate2D_comb[i][j]) / ws_temp->var(name("N_comb", i, j))->getError();
        }

        if (BF_ == 0) {
          pull_bs[i][j]->setVal(bs_pull);
          pull_rds_bs[i][j]->add(*pull_bs[i][j]);
        }
        else {
          residual_bs[i][j]->setVal(bs_residual);
          residual_rds_bs[i][j]->add(*residual_bs[i][j]);
        }
        if (!(SM_ || bd_constr_ || BF_ == 2)) {
          pull_bd[i][j]->setVal(bd_pull);
          pull_rds_bd[i][j]->add(*pull_bd[i][j]);
        }
        else if (BF_ == 2) {
          residual_bd[i][j]->setVal(bd_residual);
          residual_rds_bd[i][j]->add(*residual_bd[i][j]);
        }
        if (!rare_constr_) {
        	pull_semi[i][j]->setVal(semi_pull);
        	pull_rds_semi[i][j]->add(*pull_semi[i][j]);
        }
        pull_comb[i][j]->setVal(comb_pull);
        pull_rds_comb[i][j]->add(*pull_comb[i][j]);
        if (BF_ == 0) correlation_h[i][j]->Fill(RFR->correlation(*ws_temp->var(name("N_bs", i, j)), *ws_temp->var(name("N_bd", i, j))));
        else if (BF_ == 1) correlation_h[i][j]->Fill(RFR->correlation(*ws_temp->var("BF_bs"), *ws_temp->var(name("N_bd", i, j))));
        else correlation_h[i][j]->Fill(RFR->correlation(*ws_temp->var("BF_bs"), *ws_temp->var("BF_bd")));
      }
    }
    if (BF_ > 0) {
      BF_bs_mean_h->Fill(ws_temp->var("BF_bs")->getVal());
      BF_bs_sigma_h->Fill(ws_temp->var("BF_bs")->getError());
      double BF_pull = (ws_temp->var("BF_bs")->getVal() - Bs2MuMu_SM_BF_val) / ws_temp->var("BF_bs")->getError();
      pull_BF_bs->setVal(BF_pull);
      pull_rds_BF_bs->add(*pull_BF_bs);
    }
    if (BF_ > 1) {
      BF_bd_mean_h->Fill(ws_temp->var("BF_bd")->getVal());
      BF_bd_sigma_h->Fill(ws_temp->var("BF_bd")->getError());
      double BF_pull = (ws_temp->var("BF_bd")->getVal() - Bd2MuMu_SM_BF_val) / ws_temp->var("BF_bd")->getError();
      pull_BF_bd->setVal(BF_pull);
      pull_rds_BF_bd->add(*pull_BF_bd);
    }
    if (k == NExp) ws_temp->pdf(pdfname.c_str())->Print();

    if (sign == 0) {
      sign_h->Fill(pdf_toyMC::sig_hand(data, printlevel, ws_temp));
    }

    delete data;
//    delete RFR;
    delete ws_temp;
  }
  cout << red_color_bold << "END OF EXPERIMENTS!" << default_console_color << endl;
  for (int i = 0; i < channels; i++) {
  	if (years_=="0" && i > 1) continue;
  	if (years_=="1" && i < 2) continue;
    for (int j = 0; j < bdt_index_max(i); j++) {
      if (BF_ == 0) fit_pulls(pull_bs[i][j], pull_rds_bs[i][j], i, j);
      if (BF_ > 0) fit_pulls(residual_bs[i][j], residual_rds_bs[i][j], i, j);
      if (!(SM_ || bd_constr_ || BF_ == 2)) fit_pulls(pull_bd[i][j], pull_rds_bd[i][j], i, j);
      else if (BF_ == 2) fit_pulls(residual_bd[i][j], residual_rds_bd[i][j], i, j);
      if (!rare_constr_) fit_pulls(pull_semi[i][j], pull_rds_semi[i][j], i, j);
      fit_pulls(pull_comb[i][j], pull_rds_comb[i][j], i, j);

      print_histos(bs_mean_h[i][j], i, j);
      if (!(SM_ || bd_constr_)) print_histos(bd_mean_h[i][j], i, j);
      if (!rare_constr_) print_histos(semi_mean_h[i][j], i, j);
      print_histos(comb_mean_h[i][j], i, j);

      if (BF_ == 0) print_histos(bs_sigma_h[i][j], i, j);
      if (!(SM_ || bd_constr_)) print_histos(bd_sigma_h[i][j], i, j);
      if (!rare_constr_) print_histos(semi_sigma_h[i][j], i, j);
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
    	if (years_=="0" && i > 1) continue;
    	if (years_=="1" && i < 2) continue;
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
    sign_c->Print((get_address("sign", Bd ? "Bd" : "Bs", false) + ".gif").c_str());
    sign_c->Print((get_address("sign", Bd ? "Bd" : "Bs", false) + ".pdf").c_str());
    sign_c->Print((get_address("sign", Bd ? "Bd" : "Bs", false) + ".root").c_str());
    sign_c->Print((get_address("sign", Bd ? "Bd" : "Bs", false) + ".C").c_str());
    delete sign_c;
  }
  return;
}

void pdf_toyMC::do_bias(RooWorkspace* ws) {

  if (bias_.compare("no")) {
    for (int i = 0; i < channels; i++) {
    	if (years_=="0" && i > 1) continue;
    	if (years_=="1" && i < 2) continue;
      for (int j = 0; j < bdt_index_max(i); j++) {
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

RooFitResult* pdf_toyMC::fit_pdf(RooAbsData* data, int printlevel, RooWorkspace* ws) {
  RooFitResult* result;
  if (printlevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  result = ws->pdf(pdfname.c_str())->fitTo(*data, pee ? ConditionalObservables(*ws->var("MassRes")) : RooCmdArg::none(), Extended(true), SumW2Error(0), PrintLevel(printlevel), PrintEvalErrors(1)/*, Verbose(10)*/, Save(kTRUE), NumCPU( simul_all_? 1 : 2));
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
  rds->statOn(rp_bs, Layout(0.55, 0.9, 0.9));
  gauss_bs->plotOn(rp_bs, LineColor(kBlue));
  gauss_bs->paramOn(rp_bs, Layout(0.55, 0.9, 0.7));
  TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 600, 600);
  rp_bs->Draw();
  channel = simul_ ? i : ch_i_;
  channel_bdt = simul_bdt_ ? j : 0;
  canvas_bs->Print((get_address(pull->GetName(), pdfname) + ".gif").c_str());
  canvas_bs->Print((get_address(pull->GetName(), pdfname) + ".pdf").c_str());
  delete rp_bs;
  delete canvas_bs;
}

void pdf_toyMC::print_histos(TH1D* histos, int i, int j) {

  cout << "etacat " << i << "  bdtcat " << j << " " << histos->GetName() << " mean = " << histos->GetMean() << endl;
  TCanvas* N_mean_c = new TCanvas("N_mean_c", "N_mean_c", 600, 600);
  histos->Draw("e");
  channel = simul_ ? i : ch_i_;
  channel_bdt = simul_bdt_ ? j : 0;
  N_mean_c->Print((get_address(histos->GetName(), pdfname) + ".gif").c_str());
  N_mean_c->Print((get_address(histos->GetName(), pdfname) + ".pdf").c_str());
  delete N_mean_c;
}

Double_t pdf_toyMC::sig_hand(RooAbsData* data, int printlevel, RooWorkspace* ws_temp) {
  Double_t minNLL = RFR->minNll();
  string alt_name("N_bs");
  if (BF_ > 0) alt_name = "BF_bs";
  if (Bd) alt_name = "N_bd";
  if (Bd && BF_ > 0) alt_name = "BF_bd";
  if (BF_ == 0) {
    for (int i = 0; i < channels; i++) {
    	if (years_=="0" && i > 1) continue;
    	if (years_=="1" && i < 2) continue;
      for (int j = 0; j < bdt_index_max(i); j++) {
        ws_temp->var(name(alt_name.c_str(), i, j))->setVal(0);
        ws_temp->var(name(alt_name.c_str(), i, j))->setConstant(1);
        if (bd_constr_) {
          ws_temp->var("Bd_over_Bs")->setVal(0);
          ws_temp->var("Bd_over_Bs")->setConstant(1);
        }
      }
    }
  }
  else {
    ws_temp->var(alt_name.c_str())->setVal(0);
    ws_temp->var(alt_name.c_str())->setConstant(1);
  }

  RooFitResult * rfr_H0 = pdf_toyMC::fit_pdf(data, printlevel, ws_temp);

  if (BF_ == 0) {
    if (!simul_bdt_) {
      for (int i = 0; i < channels; i++) {
      	if (years_=="0" && i > 1) continue;
      	if (years_=="1" && i < 2) continue;
        ws_temp->var(name(alt_name.c_str(), i))->setVal(estimate_bs[i]);
        ws_temp->var(name(alt_name.c_str(), i))->setConstant(0);
        if (bd_constr_) {
          ws_temp->var("Bd_over_Bs")->setVal(estimate_bd[i]/estimate_bs[i]);
          ws_temp->var("Bd_over_Bs")->setConstant(0);
        }
      }
    }
    else {
      for (int i = 0; i < channels; i++) {
      	if (years_=="0" && i > 1) continue;
      	if (years_=="1" && i < 2) continue;
        for (int j = 0; j < bdt_index_max(i); j++) {
          ws_temp->var(name(alt_name.c_str(), i, j))->setVal(estimate2D_bs[i][j]);
          ws_temp->var(name(alt_name.c_str(), i, j))->setConstant(0);
          if (bd_constr_) {
            ws_temp->var("Bd_over_Bs")->setVal(estimate2D_bd[i][j]/estimate2D_bs[i][j]);
            ws_temp->var("Bd_over_Bs")->setConstant(0);
          }
        }
      }
    }
  }
  else {
    if (!Bd) ws_temp->var("BF_bs")->setVal(Bs2MuMu_SM_BF_val);
    else ws_temp->var("BF_bd")->setVal(Bd2MuMu_SM_BF_val);
    ws_temp->var(alt_name.c_str())->setConstant(0);
  }
  Double_t newNLL = rfr_H0->minNll();
  Double_t deltaLL = newNLL - minNLL;
//  Double_t signif = deltaLL>0 ? sqrt(2*deltaLL) : -sqrt(-2*deltaLL) ;
  Double_t signif = deltaLL>0 ? sqrt(2*deltaLL) : 0;
  return signif;
}

void pdf_toyMC::mcstudy(string pdf_toy, string pdf_test) {
  if (!simul_bdt_ && !simul_all_) {
    for (int i = 0; i < channels; i++) {
      ws_->var(name("N_bs", i))->setVal(estimate_bs[i]);
      if (!rare_constr_) ws_->var(name("N_semi", i))->setVal(estimate_semi[i]);
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
      for (int j = 0; j < bdt_index_max(i); j++) {
        ws_->var(name("N_bs", i, j))->setVal(estimate2D_bs[i][j]);
        if (!rare_constr_) ws_->var(name("N_semi", i, j))->setVal(estimate2D_semi[i][j]);
        ws_->var(name("N_comb", i, j))->setVal(estimate2D_comb[i][j]);
        if (!SM_ && !bd_constr_) ws_->var(name("N_bd", i, j))->setVal(estimate2D_bd[i][j]);
      }
    }
    if (bd_constr_) {
      double ratio = (double) estimate2D_bd[0][0] / estimate2D_bs[0][0]; // it's the same in every channel
      ws_->var("Bd_over_Bs")->setVal(ratio);
    }
  }
  RooArgSet obsv(*ws_->var("Mass"), "obsv");
  if (bdt_fit_) obsv.add(*ws_->var("bdt"));
  if (simul_ && !simul_bdt_ && !simul_all_) obsv.add(*ws_->cat("etacat"));
  if (simul_bdt_) {
    cout << "RooMCStudy seems not to work with RooSuperCategory" << endl;
    obsv.add(*ws_->cat("bdtcat"));
    exit(0);
  }
  if (simul_all_) obsv.add(*ws_->cat("allcat"));
  if (pee) obsv.add(*ws_->var("MassRes"));

  RooMCStudy * mcstudy;
  if (bias_.compare("no")) { /// bias
    RooWorkspace* fitws = (RooWorkspace*)ws_->Clone("fitws");
    do_bias(fitws);
    mcstudy = new RooMCStudy( *ws_->pdf(pdfname.c_str()), obsv, FitModel(*fitws->pdf(pdfname.c_str())), Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(pee ? ConditionalObservables(*ws_->var("MassRes")) : RooCmdArg::none(), NumCPU(2)));
  }
  else if (pdf_test == pdfname) { /// same pdf
    mcstudy = new RooMCStudy( *ws_->pdf(pdfname.c_str()), obsv, Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(pee ? ConditionalObservables(*ws_->var("MassRes")) : RooCmdArg::none(), NumCPU(2)));
  }
  else { /// different pdf
    mcstudy = new RooMCStudy( *ws_->pdf(pdfname.c_str()), obsv, FitModel(*ws_->pdf(pdf_test.c_str())), Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(pee ? ConditionalObservables(*ws_->var("MassRes")) : RooCmdArg::none(), NumCPU(2)));
  }

  vector <vector <RooDLLSignificanceMCSModule*> > sigModule(channels, vector <RooDLLSignificanceMCSModule*> (channels_bdt));
  if (BF_ == 0) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < bdt_index_max(i); j++) {
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
  mcstudy->generateAndFit(NExp, 0, 0);
  ////
  ////////

  for (int k = 0; k < 4; k++) {
    if ((SM_ || bd_constr_) && k == 1) continue;
    if (rare_constr_ && k == 3) continue;
    for (int i = 0; i < channels; i++) {
    	if (years_=="0" && i > 1) continue;
    	if (years_=="1" && i < 2) continue;
      for (int j = 0; j < bdt_index_max(i); j++) {
        if (k == 0 && BF_ > 0 && (i > 0 || j > 0)) continue;
        if (k == 1 && BF_ > 1 && (i > 0 || j > 0)) continue;
        string name_;
        if (k == 0 && BF_ > 0) {
          name_ = "BF_bs";
          source[0] = "BF_bs";
        }
        else if (k == 1 && BF_ > 1) {
          name_ = "BF_bd";
          source[1] = "BF_bd";
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
  	if (years_=="0" && i > 1) continue;
  	if (years_=="1" && i < 2) continue;
    for (int j = 0; j < bdt_index_max(i); j++) {
      if (BF_ > 0 && (i > 0 || j > 0)) continue;
      string name_(name("significance_nullhypo_N_bs", i, j));
      if (BF_ > 0) name_ = "significance_nullhypo_BF_bs";
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
