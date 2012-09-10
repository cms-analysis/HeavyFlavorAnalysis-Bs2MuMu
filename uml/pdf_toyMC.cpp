/*
 * File:   pdf_toyMC.cpp
 * Author: lucamartini
 * 
 * Created on 26 aprile 2012, 11.03
 */

#include "pdf_toyMC.h"

pdf_toyMC::pdf_toyMC(bool print, int inputs, string input_estimates, string input_cuts, string meth, string range, bool SM, bool bd_constr, TTree *input_tree, string bias, bool simul, string ch_s): pdf_fitData( print,  inputs,  input_estimates, input_cuts,  meth,  range,  SM,  bd_constr,  input_tree,   simul,  ch_s) {
  cout << "pdf_toyMC constructor" << endl;
  bias_ = bias;

  RooRandom::randomGenerator()->SetSeed(0);
}

void pdf_toyMC::generate(int NExp, string pdf_toy, string test_pdf) {

  if (!simul_) channels = 1;
  pull_bs.resize(channels);
  pull_bd.resize(channels);
  pull_rds_bs.resize(channels);
  pull_rds_bd.resize(channels);
  vector <double> p_bs(channels);
  vector <double> p_bd(channels);
  vector <TH1D*> correlation_h(channels);
  vector <TH1D*> bs_mean_h(channels);
  vector <TH2D*> corr_Nbs_Nbd_vs_N_bs_h(channels);
  for (int i = 0; i < channels; i++) {
    pull_bs[i] = new RooRealVar("pull_bs", "pull_bs", -5., 5.);
    pull_bd[i] = new RooRealVar("pull_bd", "pull_bd", -5., 5.);
    pull_rds_bs[i] = new RooDataSet("pull_rds_bs", "pull_rds_bs", *pull_bs[i]);
    pull_rds_bd[i] = new RooDataSet("pull_rds_bd", "pull_rds_bd", *pull_bd[i]);
    p_bs[i] = 0;
    p_bd[i] = 0;
    correlation_h[i] = new TH1D(Form("correlation_%d_h", i), "correlation_h", 100, -1., 1.);
    bs_mean_h[i] = new TH1D(Form("bs_mean_%d_h", i), "bs_mean_h", 100, 0., 40.);
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
    cout << ws_->GetName() << endl;
    RooWorkspace* ws_temp = (RooWorkspace*)ws_->Clone("ws_temp");
    RooDataSet* data = new RooDataSet("data", "data", *ws_temp->var("Mass"));
    double printlevel = -1;
    if (i == 1) printlevel = 1;

    RooDataSet* data_bs;
    RooDataSet* data_bd;
    RooDataSet* data_rare;
    RooDataSet* data_hist;
    RooDataSet* data_expo;
    RooDataSet* data_comb;
    if (channels == 1) {
      if (pdf_toy_ == "total") {
        if (!pee) {
          data_bs   = ws_temp->pdf("pdf_bs")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_bs[0]), Extended(1));
          data_bd   = ws_temp->pdf("pdf_bd")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_bd[0]), Extended(1));
          data_rare = ws_temp->pdf("pdf_rare")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_rare[0]), Extended(1));
          data_comb = ws_temp->pdf("pdf_comb")->generate(*ws_temp->var("Mass"), NumEvents((int)estimate_comb[0]), Extended(1));
        }
        else {
          data_bs   = ws_temp->pdf("pdf_bs")->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("MassRes")), NumEvents((int)estimate_bs[0]), Extended(1)); cout << "81" << endl;
          data_bd   = ws_temp->pdf("pdf_bd")->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("MassRes")), NumEvents((int)estimate_bd[0]), Extended(1));
          data_rare = ws_temp->pdf("pdf_rare")->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("MassRes")), NumEvents((int)estimate_rare[0]), Extended(1));
          data_comb = ws_temp->pdf("pdf_comb")->generate(RooArgSet(*ws_temp->var("Mass"), *ws_temp->var("MassRes")), NumEvents((int)estimate_comb[0]), Extended(1));
        }
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
    }
    else {
      bd_b = true;
      for (int i = 0; i < channels; i++) {
        ws_temp->var(Form("N_bs_%d",i))->setVal(estimate_bs[i]);
        if (!SM_ && !bd_constr_) ws_temp->var(Form("N_bd_%d",i))->setVal(estimate_bd[i]);
        RooRealVar* temp_var = ws_temp->var(Form("N_rare_%d",i));
        if (temp_var) ws_temp->var(Form("N_rare_%d",i))->setVal(estimate_rare[i]);
        temp_var = ws_temp->var(Form("N_hist_%d",i));
        if (temp_var) ws_temp->var(Form("N_hist_%d",i))->setVal(estimate_rare[i]);
        temp_var = ws_temp->var(Form("N_expo_%d",i));
        if (temp_var) ws_temp->var(Form("N_expo_%d",i))->setVal(estimate_rare[i]);
        ws_temp->var(Form("N_comb_%d",i))->setVal(estimate_comb[i]);
      }
      RooCategory* cat =  (RooCategory*)ws_temp->obj("channels");
      data = ws_temp->pdf("pdf_ext_simul")->generate(RooArgSet(*ws_temp->var("Mass"),*cat), Extended(1));
    }
    //////
    RFR = pdf_toyMC::fit_pdf(pdf_test_, data, printlevel, ws_temp);
    //////
    pdf_name = pdf_toy_;
    rds_ = data;
    if (!simul_) {
      if (i == 2) {
        ws_=ws_temp;
        print("_first", ws_temp);
      }
      if (i == NExp) {
        ws_=ws_temp;
        print("_last", ws_temp);
      }
      /// test statistics
      double N_bs = ws_temp->var("N_bs")->getVal();
      bs_mean_h[0]->Fill(N_bs);
      // PULL
      double bs_pull = (ws_temp->var("N_bs")->getVal() - estimate_bs[0]) / ws_temp->var("N_bs")->getError();
      pull_bs[0]->setVal(bs_pull);
      pull_rds_bs[0]->add(*pull_bs[0]);
      if (!SM_ && bd_b) {
        double bd_pull = (ws_temp->var("N_bd")->getVal() - estimate_bd[0]) / ws_temp->var("N_bd")->getError();
        pull_bd[0]->setVal(bd_pull);
        pull_rds_bd[0]->add(*pull_bd[0]);

        double correlation = RFR->correlation(*ws_temp->var("N_bs"), *ws_temp->var("N_bd"));
        correlation_h[0]->Fill(correlation);
        corr_Nbs_Nbd_vs_N_bs_h[0]->Fill(data_bs->numEntries() , correlation);
        if (correlation > 0. && corr0 == false) {
          ws_=ws_temp;
          print("_corrzero", ws_temp);
          corr0 = true;
        }
        if (correlation < -0.3 && corrneg == false) {
          ws_=ws_temp;
          print("_corrneg", ws_temp);
          corrneg = true;
        }
      }
    }
    else {
      for (int j = 0; j < channels; j++) {
      /// test statistics
        bs_mean_h[j]->Fill(ws_temp->var(Form("N_bs_%d",j))->getVal());
      // PULL
        double bs_pull = (ws_temp->var(Form("N_bs_%d",j))->getVal() - estimate_bs[j]) / ws_temp->var(Form("N_bs_%d",j))->getError();
        pull_bs[j]->setVal(bs_pull);
        pull_rds_bs[j]->add(*pull_bs[j]);
      }
    }
    // H_0
    if (pdf_toy_ == "pdf_ext_total") {
      RooDataSet* bkg = new RooDataSet("bkg", "bkg", *ws_temp->var("Mass"));
      bkg->append(*data_rare);
      bkg->append(*data_comb);
      fit_pdf("total", bkg, printlevel, ws_temp);
      if (ws_temp->var("N_bs")->getVal() >= estimate_bs[0]) p_bs[0]++;
      if (!SM_) if (ws_temp->var("N_bd")->getVal() >= estimate_bd[0]) p_bd[0]++;
      delete bkg;
    }
    if (i == NExp) ws_temp->pdf(pdf_toy_.c_str())->Print();
    delete data;
    if (!simul_) {
      delete data_bs;
      delete data_bd;
      delete data_comb;
      delete data_rare;
//      delete ws_temp;
    }
    delete RFR;
  }

  pdf_toy_ = pdf_toy;
  fit_pulls();

  for (int j = 0; j < channels; j++) {
    cout << "channel " << j << "  N_bs mean = " << bs_mean_h[j]->GetMean() << endl;
    TCanvas* N_mean_c = new TCanvas("N_mean_c", "N_mean_c", 600, 600);
    bs_mean_h[j]->Draw();
    ostringstream address;
    if (!simul_) address << "fig/bs_mean_" << meth_ << "_" << ch_s_ << "_" + pdf_toy;
    else address << "fig/bs_mean_" << meth_ << "_" << j << "_" + pdf_toy;
    if (SM_) address << "_SM";
    if (bd_constr_) address << "_bd_const";
    if (pee) address << "_pee";
    N_mean_c->Print((address.str() + ".gif").c_str());
    N_mean_c->Print((address.str() + ".pdf").c_str());
    delete N_mean_c;
  }
//  double p_value_bs = p_bs / NExp;
//  cout << "p-value bs = " << p_value_bs << endl;
//  if (!SM_) {
//    double p_value_bd = p_bd / NExp;
//    cout << "p-value bd = " << p_value_bd << endl;
//  }

  if (!SM_ && bd_b && !simul_) {
    TCanvas* corr_c = new TCanvas("corr_c", "corr_c", 1200, 600);
    corr_c->Divide(3);
    corr_c->cd(1);
    correlation_h[0]->Draw();
    corr_c->cd(2);
    corr_Nbs_Nbd_vs_N_bs_h[0]->Draw();
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
  RooFitResult* result;
  if (!pee) result = ws->pdf(pdf_toy_.c_str())->fitTo(*data, Extended(true), SumW2Error(0), Range(range_.c_str()), PrintLevel(printlevel), PrintEvalErrors(-1), Save(kTRUE));
  else result = ws->pdf(pdf_toy_.c_str())->fitTo(*data, ConditionalObservables(*ws->var("MassRes")), Extended(true), SumW2Error(0), Range(range_.c_str()), PrintLevel(printlevel), PrintEvalErrors(-1), Save(kTRUE));
  return result;
}

void pdf_toyMC::fit_pulls() {
  for (int i = 0; i < channels; i++) {
    RooRealVar* mean_bs = new RooRealVar("mean_bs", "mean_bs", -5., 5.);
    RooRealVar* sigma_bs = new RooRealVar("sigma_bs", "sigma_bs", 0.001, 5.);
    RooGaussian* gauss_bs = new RooGaussian("gauss_bs", "gauss_bs", *pull_bs[i], *mean_bs, *sigma_bs);
    gauss_bs->fitTo(*pull_rds_bs[i]);
  
    RooPlot *rp_bs = pull_bs[i]->frame();
    pull_rds_bs[i]->plotOn(rp_bs, Binning(40));
    gauss_bs->plotOn(rp_bs, LineColor(kBlue));
    gauss_bs->paramOn(rp_bs, Layout(0.66, 0.9, 0.9));
    TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 600, 600);
    rp_bs->Draw();
    ostringstream address_oss;
    if (!simul_) address_oss << "fig/bs_pull_" << meth_ << "_" << ch_s_ << "_" << pdf_toy_;
    else address_oss << "fig/bs_pull_" << meth_ << "_" << i << "_" << pdf_toy_;
    string address = address_oss.str();
    if (SM_) address += "_SM";
    canvas_bs->Print( (address + ".gif").c_str());
    canvas_bs->Print( (address + ".pdf").c_str());
    delete rp_bs;
    delete canvas_bs;

  }
//  if (!SM_ && !simul_) {
//    RooRealVar* mean_bd = new RooRealVar("mean_bd", "mean_bd", -5., 5.);
//    RooRealVar* sigma_bd = new RooRealVar("sigma_bd", "sigma_bd", 0.001, 5.);
//    RooGaussian* gauss_bd = new RooGaussian("gauss_bd", "gauss_bd", *pull_bd, *mean_bd, *sigma_bd);
//    gauss_bd->fitTo(*pull_rds_bd);
  
//    RooPlot *rp_bd = pull_bd->frame();
//    pull_rds_bd->plotOn(rp_bd, Binning(40));
//    gauss_bd->plotOn(rp_bd, LineColor(kBlue));
//    gauss_bd->paramOn(rp_bd, Layout(0.66, 0.9, 0.9));
//    TCanvas* canvas_bd = new TCanvas("canvas_bd", "canvas_bd", 600, 600);
//    rp_bd->Draw();
//    address = "fig/bd_pull_" + meth_ + "_" + ch_s_ + "_" + pdf_toy_;
//    canvas_bd->Print( (address + ".gif").c_str());
//    canvas_bd->Print( (address + ".pdf").c_str());
//    delete rp_bd;
//    delete canvas_bd;
//  }
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
    if (!SM_ && !bd_constr_) ws_->var("N_bd")->setVal(estimate_bd[0]);
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
    mcstudy->generateAndFit(NExp, 0);

    // Make plots of the distributions of mean, the error on mean and the pull of mean
    RooPlot* frame1_bs = mcstudy->plotParam(*ws_->var("N_bs"), Bins(20)) ;
    RooPlot* frame2_bs = mcstudy->plotError(*ws_->var("N_bs"), Bins(20), Range(0., 15.)) ;
    RooPlot* frame3_bs = mcstudy->plotPull(*ws_->var("N_bs"), Bins(20), FitGauss(kTRUE), Range(-8., 8.)) ;
    TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 600, 600) ;
//    canvas_bs->Divide(3,1) ;
//    canvas_bs->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1_bs->GetYaxis()->SetTitleOffset(1.4) ; frame1_bs->Draw() ;
//    canvas_bs->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_bs->GetYaxis()->SetTitleOffset(1.4) ; frame2_bs->Draw() ;
//    canvas_bs->cd(3) ;
    gPad->SetLeftMargin(0.15) ; frame3_bs->GetYaxis()->SetTitleOffset(1.4) ; frame3_bs->Draw() ;
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
    sprintf(leg, "|#eta_{#mu}|<1.6");
    t->DrawLatex(0.7, 0.75, leg);
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
        TCanvas* corr_c = new TCanvas("corr_c", "corr_c", 600, 600);
        //corr_c->Divide(3);
//        corr_c->cd(1);
//        corr_Nbs_Nbd->Draw("box");
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
//        corr_c->cd(2);
        corr_Nbs_Nbd_h->Draw();
//        corr_c->cd(3);
//        corr_Nbs_Nbd_vs_N_bs_h->Draw();
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
  else { // simultaneous
    for (int i = 0; i < channels; i++) {
      ws_->var(Form("N_bs_%d", i))->setVal(estimate_bs[i]);
      ws_->var(Form("N_rare_%d", i))->setVal(estimate_rare[i]);
      ws_->var(Form("N_comb_%d", i))->setVal(estimate_comb[i]);
      if (!SM_ && !bd_constr_) ws_->var(Form("N_bd_%d", i))->setVal(estimate_bd[i]);
    }
    if (bd_constr_) {
       double ratio = (double) estimate_bd[0] / estimate_bs[0]; // it's the same
       ws_->var("Bd_over_Bs")->setVal(ratio);
     }

    simul_pdf = (RooSimultaneous*)ws_->obj("pdf_ext_simul");
    RooCategory* channel__ = (RooCategory*)ws_->obj("channels");
    RooRealVar* mass__ = ws_->var("Mass");
    RooMCStudy * mcstudy = new RooMCStudy( *simul_pdf, RooArgSet(*mass__, *channel__),  Binned(kFALSE), Extended(kTRUE), FitOptions(Save(kTRUE)), Silence());
    vector <RooDLLSignificanceMCSModule*> sigModule;
    sigModule.resize(channels);
    for (int i = 0; i < channels; i++) {
      ostringstream varname;
      varname << "N_bs_" << i;
      sigModule[i] = new  RooDLLSignificanceMCSModule(*ws_->var(varname.str().c_str()), 0);
      mcstudy->addModule(*sigModule[i]) ;
    }

    mcstudy->generateAndFit(NExp, 0, kTRUE);

  //  RooPlot* frame0_bs = mcstudy->plotPull(*ws_->var("N_bs_0"), Bins(20), FitGauss(kTRUE)/*, Range(-5., 5.)*/) ;
  //  RooPlot* frame1_bs = mcstudy->plotPull(*ws_->var("N_bs_1"), Bins(20), FitGauss(kTRUE)/*, Range(-5., 5.)*/) ;

    for (int i = 0; i < channels; i++) {
      RooPlot* frame1_bs = mcstudy->plotParam(*ws_->var(Form("N_bs_%d", i)), Bins(20)) ;
      RooPlot* frame2_bs = mcstudy->plotError(*ws_->var(Form("N_bs_%d", i)), Bins(20), Range(0., 15.)) ;
      RooPlot* frame3_bs = mcstudy->plotPull(*ws_->var(Form("N_bs_%d", i)), Bins(20), FitGauss(kTRUE), Range(-5., 5.)) ;
      TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 1200, 600) ;
      canvas_bs->Divide(3,1) ;
      canvas_bs->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1_bs->GetYaxis()->SetTitleOffset(1.4) ; frame1_bs->Draw() ;
      canvas_bs->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_bs->GetYaxis()->SetTitleOffset(1.4) ; frame2_bs->Draw() ;
      canvas_bs->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3_bs->GetYaxis()->SetTitleOffset(1.4) ; frame3_bs->Draw() ;
      ostringstream index;
      index << i;
      string address = "fig/RooMCStudy_simul_bs_" + index.str() + "_" + meth_;
      if (SM_) address += "_SM";
      if (bd_constr_) address += "_bdConstr";
      canvas_bs->Print((address + ".gif").c_str());
      canvas_bs->Print((address + ".pdf").c_str());
      delete frame1_bs;
      delete frame2_bs;
      delete frame3_bs;
      delete canvas_bs;
    }
    TCanvas* sig_c = new TCanvas("sig_c", "sig_c", 400*channels, 600);
    sig_c->Divide(channels);
    for (int i = 0; i < channels; i++) {
      sig_c->cd(i+1);
      ostringstream varname;
      varname << "significance_nullhypo_N_bs_" << i;
      TH1* sig_h = mcstudy->fitParDataSet().createHistogram(varname.str().c_str());
      sig_h->Draw();
    }
    string address = "fig/RooMCStudy_simul_sig_bs_" + meth_;
    if (SM_) address += "_SM";
    sig_c->Print((address + ".gif").c_str());
    sig_c->Print((address + ".pdf").c_str());
    delete sig_c;
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

void pdf_toyMC::print(string output, RooWorkspace* ws) {
  int colors[11] = {632, 400, 616, 432, 800, 416, 820, 840, 860, 880, 900};
  RooPlot *rp = ws_->var("Mass")->frame();
  rds_->plotOn(rp, Binning(20));
  if (!pee) ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue), Range(range_.c_str())/*, ProjectionRange("eta_all"), Normalization((rds_->sumEntries(), RooAbsReal::NumEvent))*/);
  else {
    ws_->pdf(pdf_name.c_str())->Print();
    TH1* mass_eta_h;
    /*if (simul_) mass_eta_h = ws_->pdf(pdf_name.c_str())->createHistogram("fit", *ws_->var("Mass"), Binning(50), YVar(*ws_->var(Form("eta_%d", channel)), Binning(50))) ;
    else*/ mass_eta_h = ws_->pdf(pdf_name.c_str())->createHistogram("fit", *ws_->var("Mass"), Binning(50), YVar(*ws_->var("eta"), Binning(50))) ;
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
    cetad->Divide(2);
    cetad->cd(1);
    mass_eta_h->Draw("surf") ;
    cetad->cd(2);
    ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue), Range(range_.c_str()), Normalization(rds_->sumEntries(), RooAbsReal::NumEvent)/*, ProjWData(*ws_->data(Form("etardh_%s", pdf_name.c_str())))*/, ProjectionRange("eta_all"));
    if(!no_legend) ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
    rp->Draw();
    string address;
    if (simul_) address = "fig/" + pdf_name + "_MassEta_" + meth_ + "_simul" + output;
    else address = "fig/" + pdf_name + "_MassEta_" + meth_ + "_" + ch_s_ + output;
    if (SM_) address += "_SM";
    if (bd_constr_) address += "_BdConst";
    cetad->Print( (address + ".gif").c_str());
    cetad->Print( (address + ".pdf").c_str());
    delete cetad;
    delete frame;
    delete mass_eta_h;
    cout << pdf_name << endl;
  }
  if(!no_legend) ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));

  RooArgSet * set = ws_->pdf(pdf_name.c_str())->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  int i = 0;
  while((var_Obj = it->Next()) && !pee){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (i > 11) i = 0;
      size_t found1 = pdf_name.find("total");
      if (found1 == string::npos) ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
      else {
        if (name=="pdf_bs" || name=="pdf_bd" || name=="pdf_rare" || name=="pdf_comb") ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
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
  if (pee) address += "_PEE";
  c->Print( (address + ".gif").c_str());
  c->Print( (address + ".pdf").c_str());
  delete rp;
  delete c;

  return;
}
