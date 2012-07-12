/*
 * File:   pdf_toyMC.cpp
 * Author: lucamartini
 * 
 * Created on 26 aprile 2012, 11.03
 */

#include "pdf_toyMC.h"

pdf_toyMC::pdf_toyMC(bool print, int inputs, string input_estimates, string meth, string range, bool SM, bool bd_constr, TTree *input_tree, string bias, bool simul, string ch_s): pdf_fitData( print,  inputs,  input_estimates,  meth,  range,  SM,  bd_constr,  input_tree,   simul,  ch_s) {
  cout << "pdf_toyMC constructor" << endl;
  bias_ = bias;

  pull_h_bs = new TH1D("pull_h_bs", "pull_h_bs", 20, -5., 5.);
  pull_h_bd = new TH1D("pull_h_bd", "pull_h_bd", 20, -5., 5.);
  RooRandom::randomGenerator()->SetSeed(0);
}

void pdf_toyMC::generate(int NExp, string pdf_toy, string test_pdf) {
  
  pull_bs = new RooRealVar("pull_bs", "pull_bs", -5., 5.);
  pull_bd = new RooRealVar("pull_bd", "pull_bd", -5., 5.);
  pull_rds_bs = new RooDataSet("pull_rds_bs", "pull_rds_bs", *pull_bs);
  pull_rds_bd = new RooDataSet("pull_rds_bd", "pull_rds_bd", *pull_bd);
  
  pdf_name = "total";
  if (pdf_toy != "total") {
    pdf_name = define_pdf_sum(pdf_toy);
  }
  pdf_test_= "total";
  if (test_pdf != "total") {
    pdf_test_ = define_pdf_sum(test_pdf);
  }
  
  double p_bs = 0., p_bd = 0.;
  bool bd_b = false, corr0 = false, corrneg = false;
  TH1D* correlation_h = new TH1D("correlation_h", "correlation_h", 100, -1., 1.);
  TH1D* bs_mean_h = new TH1D("bs_mean_h", "bs_mean_h", 100, 0., 20.);
  TH2D* corr_Nbs_Nbd_vs_N_bs_h = new TH2D("corr_Nbs_Nbd_vs_N_bs_h", "corr_Nbs_Nbd_vs_N_bs_h", 30, 0., 30., 100, -1., 1.);

  for (int i = 1; i <= NExp; i++) {
    if (i%100 == 0) cout << "Exp # " << i << " / " << NExp << endl;
    pdf_toy_ = pdf_toy;

    RooWorkspace* ws_temp = (RooWorkspace*)ws_->Clone("ws_temp");
    RooDataSet* data = new RooDataSet("data", "data", *ws_temp->var("Mass"));
    RooDataSet* data_bs;
    RooDataSet* data_bd;
    RooDataSet* data_rare;
    RooDataSet* data_hist;
    RooDataSet* data_expo;
    RooDataSet* data_comb;
    if (pdf_toy_ == "total") {
      data_bs   = ws_temp->pdf("pdf_bs")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_bs[0]), Extended(1));
      data_bd   = ws_temp->pdf("pdf_bd")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_bd[0]), Extended(1));
      data_rare = ws_temp->pdf("pdf_rare")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_rare[0]), Extended(1));
      data_comb = ws_temp->pdf("pdf_comb")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_comb[0]), Extended(1));
      data->append(*data_bs);
      data->append(*data_bd);
      data->append(*data_rare);
      data->append(*data_comb);
      bd_b = true;
    }
    else {
      size_t found;
      found = pdf_toy_.find("bs");
      if (found != string::npos) {
        data_bs   = ws_temp->pdf("pdf_bs")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_bs[0]), Extended(1));
        data->append(*data_bs);
      }
      found = pdf_toy_.find("bd");
      if (found != string::npos) {
        data_bd   = ws_temp->pdf("pdf_bd")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_bd[0]), Extended(1));
        data->append(*data_bd);
        bd_b = true;
      }
      found = pdf_toy_.find("rare");
      if (found != string::npos) {
        data_rare = ws_temp->pdf("pdf_rare")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_rare[0]), Extended(1));
        data->append(*data_rare);
      }
      found = pdf_toy_.find("comb");
      if (found != string::npos) {
        data_comb = ws_temp->pdf("pdf_comb")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_comb[0]), Extended(1));
        data->append(*data_comb);
      }
      found = pdf_toy_.find("hist");
      if (found != string::npos) {
        data_hist = ws_temp->pdf("pdf_hist")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_rare[0]), Extended(1));
        data->append(*data_hist);
      }
      found = pdf_toy_.find("expo");
      if (found != string::npos) {
        data_expo = ws_temp->pdf("pdf_expo")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_rare[0]), Extended(1));
        data->append(*data_expo);
      }
    }

    ws_temp->var("N_bs")->setVal(estimate_bs[0]);
    if (!SM_ && !bd_constr_) ws_temp->var("N_bd")->setVal(estimate_bd[0]);
    RooRealVar* temp_var = ws_temp->var("N_rare");
    if (temp_var) ws_temp->var("N_rare")->setVal(estimate_rare[0]);
    temp_var = ws_temp->var("N_hist");
    if (temp_var) ws_temp->var("N_hist")->setVal(estimate_rare[0]);
    temp_var = ws_temp->var("N_expo");
    if (temp_var) ws_temp->var("N_expo")->setVal(estimate_rare[0]);
    ws_temp->var("N_comb")->setVal(estimate_comb[0]);

    double printlevel = -1;
    if (i == 1) printlevel = 1;
    //////
    RFR = pdf_toyMC::fit_pdf(pdf_test_, data, printlevel, ws_temp);
    //////
    pdf_name = pdf_toy_;
    rds_ = data;
    if (i == 1) {
      pdf_analysis::print("_first", ws_temp);
    }
    if (i == NExp) {
      pdf_analysis::print("_last", ws_temp);
    }
    /// test statistics
    double N_bs = ws_temp->var("N_bs")->getVal();
    bs_mean_h->Fill(N_bs);
    // PULL
    double bs_pull = (ws_temp->var("N_bs")->getVal() - estimate_bs[0]) / ws_temp->var("N_bs")->getError();
    pull_bs->setVal(bs_pull);
    pull_rds_bs->add(*pull_bs);
    if (!SM_ && bd_b) {
      double bd_pull = (ws_temp->var("N_bd")->getVal() - estimate_bd[0]) / ws_temp->var("N_bd")->getError();
      pull_bd->setVal(bd_pull);
      pull_rds_bd->add(*pull_bd);

      double correlation = RFR->correlation(*ws_temp->var("N_bs"), *ws_temp->var("N_bd"));
      correlation_h->Fill(correlation);
      corr_Nbs_Nbd_vs_N_bs_h->Fill(data_bs->numEntries() , correlation);
      if (correlation > 0. && corr0 == false) {
        pdf_analysis::print("_corrzero", ws_temp);
        corr0 = true;
      }
      if (correlation < -0.3 && corrneg == false) {
        pdf_analysis::print("_corrneg", ws_temp);
        corrneg = true;
      }
    }
    
    // H_0
    if (pdf_toy_ == "pdf_ext_total") {
      RooDataSet* bkg = new RooDataSet("bkg", "bkg", *ws_temp->var("Mass"));
      bkg->append(*data_rare);
      bkg->append(*data_comb);
      fit_pdf("total", bkg, printlevel, ws_temp);
      if (ws_temp->var("N_bs")->getVal() >= estimate_bs[0]) p_bs++;
      if (!SM_) if (ws_temp->var("N_bd")->getVal() >= estimate_bd[0]) p_bd++;
      delete bkg;
    }
    if (i == NExp) ws_temp->pdf(pdf_toy_.c_str())->Print();
    delete data;
    delete data_bs;
    delete data_bd;
    delete data_comb;
    delete data_rare;
    delete ws_temp;
    delete RFR;
  }

  pdf_toy_ = pdf_toy;
  fit_pulls();
  
  cout << "N_bs mean = " << bs_mean_h->GetMean() << endl;
  double p_value_bs = p_bs / NExp;
  cout << "p-value bs = " << p_value_bs << endl;
  if (!SM_) {
    double p_value_bd = p_bd / NExp;
    cout << "p-value bd = " << p_value_bd << endl;
  }

  if (!SM_ && bd_b) {
    TCanvas* corr_c = new TCanvas("corr_c", "corr_c", 1200, 600);
    corr_c->Divide(3);
    corr_c->cd(1);
    correlation_h->Draw();
    corr_c->cd(2);
    corr_Nbs_Nbd_vs_N_bs_h->Draw();
    string address = "fig/corr_" + meth_ + "_" + ch_s_ + "_" + pdf_toy;
    corr_c->Print((address + ".gif").c_str());
    corr_c->Print((address + ".pdf").c_str());
    delete corr_c;
  }

  return;
}

RooFitResult* pdf_toyMC::fit_pdf (string pdf, RooAbsData* data, int printlevel, RooWorkspace* ws) {
  pdf_toy_ = "pdf_ext_" + pdf;
  if (bias_.compare("no")) {
    size_t found;
    found = bias_.find("c");
    if (found != string::npos) {
      double value =  ws->var("c_nonpeaking")->getVal();
      found = bias_.find("+");
      if (found != string::npos) {
        double error = ws->var("c_nonpeaking")->getErrorHi();
        ws->var("c_nonpeaking")->setVal(value + error);
      }
      else {
        double error = ws->var("c_nonpeaking")->getErrorLo();
        ws->var("c_nonpeaking")->setVal(value - error);
      }
    }
    found = bias_.find("p");
    if (found != string::npos) {
      double value =  ws->var("p_nonpeaking")->getVal();
      found = bias_.find("+");
      if (found != string::npos) {
        double error = ws->var("p_nonpeaking")->getErrorHi();
        ws->var("p_nonpeaking")->setVal(value + error);
      }
      else {
        double error = ws->var("p_nonpeaking")->getErrorLo();
        ws->var("p_nonpeaking")->setVal(value - error);
      }
    }
  }
  RooFitResult* result = ws->pdf(pdf_toy_.c_str())->fitTo(*data, Extended(true), SumW2Error(0), Range(range_.c_str()), PrintLevel(printlevel), PrintEvalErrors(-1), Save(kTRUE));
  return result;
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
  if (SM_) address += "_SM";
  canvas_bs->Print( (address + ".gif").c_str());
  canvas_bs->Print( (address + ".pdf").c_str());
  delete rp_bs;
  delete canvas_bs;
  
  if (!SM_) {
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
    canvas_bd->Print( (address + ".pdf").c_str());
    delete rp_bd;
    delete canvas_bd;
  }
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

  if (!simul_) {
    ws_->var("N_bs")->setVal(estimate_bs[0]);
    ws_->var("N_bd")->setVal(estimate_bd[0]);
    ws_->var("N_rare")->setVal(estimate_rare[0]);
    ws_->var("N_comb")->setVal(estimate_comb[0]);

    pdf_name = "pdf_ext_total";
    if (pdf_toy != "total") {
      pdf_name = "pdf_ext_" + define_pdf_sum(pdf_toy);
    }

    RooMCStudy* mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), *ws_->var("Mass"), Binned(kFALSE), Silence(), Extended(kTRUE), FitOptions(PrintLevel(-1), Save(kTRUE), PrintEvalErrors(-1)));

    cout << pdf_toy << endl;
    RooDLLSignificanceMCSModule sigModule(*ws_->var("N_bs"), 0);
    mcstudy->addModule(sigModule) ;
    cout << "toy MC on pdf: " << endl;
    ws_->pdf(pdf_name.c_str())->Print();
    cout << endl;
    mcstudy->generateAndFit(NExp, 0, kTRUE);
    cout << "mcstudy:" << endl;
    mcstudy->Print();
    cout << "fitParDataSet:"<< endl;
    mcstudy->fitParDataSet().Print();
    cout << "fitResult(0):"<< endl;
    mcstudy->fitResult(0)->Print();
    cout << "genData(0):"<< endl;
    mcstudy->genData(0)->Print();
    cout << "fitParams(0):"<< endl;
    mcstudy->fitParams(0)->Print();

    // Make plots of the distributions of mean, the error on mean and the pull of mean
    RooPlot* frame1_bs = mcstudy->plotParam(*ws_->var("N_bs"), Bins(20)) ;
    RooPlot* frame2_bs = mcstudy->plotError(*ws_->var("N_bs"), Bins(20), Range(0., 8.)) ;
    RooPlot* frame3_bs = mcstudy->plotPull(*ws_->var("N_bs"), Bins(20), FitGauss(kTRUE)) ;
    TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 1200, 600) ;
    canvas_bs->Divide(3,1) ;
    canvas_bs->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1_bs->GetYaxis()->SetTitleOffset(1.4) ; frame1_bs->Draw() ;
    canvas_bs->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_bs->GetYaxis()->SetTitleOffset(1.4) ; frame2_bs->Draw() ;
    canvas_bs->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3_bs->GetYaxis()->SetTitleOffset(1.4) ; frame3_bs->Draw() ;
    string address = "fig/RooMCStudy_bs_" + meth_ + "_" + ch_s_ + "_" + pdf_toy;
    if (SM_) address += "_SM";
    canvas_bs->Print((address + ".gif").c_str());
    canvas_bs->Print((address + ".pdf").c_str());
    delete frame1_bs;
    delete frame2_bs;
    delete frame3_bs;
    delete canvas_bs;
    TCanvas* sig_c = new TCanvas("sig_c", "sig_c", 600, 600);
    TH1* sig_h = mcstudy->fitParDataSet().createHistogram("significance_nullhypo_N_bs");
    sig_h->Draw();
    address = "fig/RooMCStudy_sig_bs_" + meth_ + "_" + ch_s_ + "_" + pdf_toy;
    if (SM_) address += "_SM";
    sig_c->Print((address + ".gif").c_str());
    sig_c->Print((address + ".pdf").c_str());
    delete sig_c;

    if (!SM_) {
      size_t found;
      found = pdf_name.find("bd");
      if (found != string::npos || pdf_toy == "total") {
        TH1* corr_Nbs_Nbd  = mcstudy->fitParDataSet().createHistogram("corr_Nbs_Nbd", *ws_->var("N_bs"), YVar(*ws_->var("N_bd")));
        corr_Nbs_Nbd->GetXaxis()->SetRangeUser(0., 30.);
        corr_Nbs_Nbd->GetYaxis()->SetRangeUser(0., 10.);
        TCanvas* corr_c = new TCanvas("corr_c", "corr_c", 1200, 600);
        corr_c->Divide(3);
        corr_c->cd(1);
        corr_Nbs_Nbd->Draw("box");
        TH1D* corr_Nbs_Nbd_h = new TH1D("corr_Nbs_Nbd_h", "corr_Nbs_Nbd_h", 100, -1., 1.);
        TH2D* corr_Nbs_Nbd_vs_N_bs_h = new TH2D("corr_Nbs_Nbd_vs_N_bs_h", "corr_Nbs_Nbd_vs_N_bs_h", 30, 0., 30., 100, -1., 1.);
        for (int i = 0; i < NExp; i++) {
          double correlation = mcstudy->fitResult(i)->correlation(*ws_->var("N_bs"), *ws_->var("N_bd"));
          corr_Nbs_Nbd_h->Fill(correlation);
          double bs_pull = mcstudy->fitParams(i)->getRealValue("N_bspull");
          double bs_fit = mcstudy->fitParams(i)->getRealValue("N_bs");
          double bs_err = mcstudy->fitParams(i)->getRealValue("N_bserr");
          double bs_gen = bs_fit - bs_err*bs_pull;
          corr_Nbs_Nbd_vs_N_bs_h->Fill(bs_gen, correlation);
        }
        corr_c->cd(2);
        corr_Nbs_Nbd_h->Draw();
        corr_c->cd(3);
        corr_Nbs_Nbd_vs_N_bs_h->Draw();
        address = "fig/RooMCStudy_corr_" + meth_ + "_" + ch_s_ + "_" + pdf_toy;
        corr_c->Print((address + ".gif").c_str());
        corr_c->Print((address + ".pdf").c_str());
        delete corr_c;

        RooPlot* frame1_bd = mcstudy->plotParam(*ws_->var("N_bd"), Bins(20)) ;
        RooPlot* frame2_bd = mcstudy->plotError(*ws_->var("N_bd"), Bins(20)) ;
        RooPlot* frame3_bd = mcstudy->plotPull(*ws_->var("N_bd"), Bins(20), FitGauss(kTRUE)) ;
        TCanvas* canvas_bd = new TCanvas("canvas_bd", "canvas_bd", 900, 900) ;
        canvas_bd->Divide(3,1) ;
        canvas_bd->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1_bd->GetYaxis()->SetTitleOffset(1.4) ; frame1_bd->Draw() ;
        canvas_bd->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_bd->GetYaxis()->SetTitleOffset(1.4) ; frame2_bd->Draw() ;
        canvas_bd->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3_bd->GetYaxis()->SetTitleOffset(1.4) ; frame3_bd->Draw() ;
        address = "fig/RooMCStudy_bd_" + meth_ + "_" + ch_s_ + "_" + pdf_toy;
        canvas_bd->Print((address + ".gif").c_str());
        canvas_bd->Print((address + ".pdf").c_str());
        delete frame1_bd;
        delete frame2_bd;
        delete frame3_bd;
        delete canvas_bd;
      }
    }
  }
  else {

    RooArgSet* all_vars = simul_pdf->getVariables();
    TIterator* vars_it = all_vars->createIterator();
    RooRealVar *arg_var = 0;
    while ( (arg_var = (RooRealVar*)vars_it->Next())) {
      cout << arg_var->GetName() << " : " << arg_var->getVal() << endl;
    }

    RooCategory* channel__ = (RooCategory*)ws_->obj("channel");
    RooRealVar* mass__ = ws_->var("Mass");
    RooMCStudy * mcstudy = new RooMCStudy( *simul_pdf, RooArgSet(*mass__, *channel__),  Binned(kFALSE), Extended(kTRUE), FitOptions(Save(kTRUE)), Silence());
    mcstudy->generateAndFit(NExp, 0, kTRUE);
    RooPlot* frame0_bs = mcstudy->plotPull(*ws_->var("N_bs_0"), Bins(20), FitGauss(kTRUE), Range(-5., 5.)) ;
    RooPlot* frame1_bs = mcstudy->plotPull(*ws_->var("N_bs_1"), Bins(20), FitGauss(kTRUE), Range(-5., 5.)) ;
    TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 800, 600) ;
    canvas_bs->Divide(2,1) ;
    canvas_bs->cd(1) ; gPad->SetLeftMargin(0.15) ; frame0_bs->GetYaxis()->SetTitleOffset(1.4) ; frame0_bs->Draw() ;
    canvas_bs->cd(2) ; gPad->SetLeftMargin(0.15) ; frame1_bs->GetYaxis()->SetTitleOffset(1.4) ; frame1_bs->Draw() ;
    string address = "fig/RooMCStudy_simul_bs_" + meth_ + "_" + pdf_toy;
    if (SM_) address += "_SM";
    if (bd_constr_) address += "_bdConstr";
    canvas_bs->Print((address + ".gif").c_str());
    canvas_bs->Print((address + ".pdf").c_str());
    delete frame0_bs;
    delete frame1_bs;
    delete canvas_bs;
  }
  return;
}

void pdf_toyMC::set_ws(RooWorkspace *ws) {
  ws_ = ws;
  RooConstVar * SM_constr = (RooConstVar *) ws->obj("SM_Bd_over_Bs");
  if (SM_constr) {
    SM_ = true;
  }
  else SM_ = false;
}
