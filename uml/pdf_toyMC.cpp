/*
 * File:   pdf_toyMC.cpp
 * Author: lucamartini
 * 
 * Created on 26 aprile 2012, 11.03
 */

#include "pdf_toyMC.h"

pdf_toyMC::pdf_toyMC(string input_estimates, bool print, string meth, string ch_s, string range): pdf_analysis(print, meth, ch_s, range) {
  
  input_estimates_ = input_estimates;
  pull_h_bs = new TH1D("pull_h_bs", "pull_h_bs", 20, -5., 5.);
  pull_h_bd = new TH1D("pull_h_bd", "pull_h_bd", 20, -5., 5.);
  RooRandom::randomGenerator()->SetSeed(0);
  parse_estimate();
  
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
      cout << "==> Error parsing variable " << cutName << ". abort..." << endl;
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
  
  pdf_name = "total";
  if (pdf_toy != "total") {
    pdf_name = define_pdf_sum(pdf_toy);
  }
  
  double p_bs = 0., p_bd = 0.;
  for (int i = 1; i <= NExp; i++) {
    pdf_toy_ = pdf_toy;
    RooDataSet* data_bs   = ws_->pdf("pdf_bs")->generate(*ws_->var("Mass"), NumEvents((int)estimate_bs), Extended(1));
    RooDataSet* data_bd   = ws_->pdf("pdf_bd")->generate(*ws_->var("Mass"), NumEvents((int)estimate_bd), Extended(1));
    RooDataSet* data_rare = ws_->pdf("pdf_rare")->generate(*ws_->var("Mass"), NumEvents((int)estimate_rare), Extended(1));
    RooDataSet* data_comb = ws_->pdf("pdf_comb")->generate(*ws_->var("Mass"), NumEvents((int)estimate_comb), Extended(1));
    
    RooDataSet* data = new RooDataSet("data", "data", *ws_->var("Mass"));
    if (pdf_toy_ == "total") {
      data->append(*data_bs);
      data->append(*data_bd);
      data->append(*data_rare);
      data->append(*data_comb);
    }
    else {
      size_t found;
      found = pdf_toy_.find("bs");
      if (found != string::npos) data->append(*data_bs);
      found = pdf_toy_.find("bd");
      if (found != string::npos) data->append(*data_bd);
      found = pdf_toy_.find("rare");
      if (found != string::npos) data->append(*data_rare);
      found = pdf_toy_.find("comb");
      if (found != string::npos) data->append(*data_comb);
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
    double bs_pull = (ws_->var("N_bs")->getVal() - estimate_bs) / ws_->var("N_bs")->getError();
    pull_bs->setVal(bs_pull);
    pull_rds_bs->add(*pull_bs);
    if (!SM_) {
      double bd_pull = (ws_->var("N_bd")->getVal() - estimate_bd) / ws_->var("N_bd")->getError();
      pull_bd->setVal(bd_pull);
      pull_rds_bd->add(*pull_bd);
    }
    
    // H_0
    if (pdf_toy_ == "pdf_ext_total") {
      RooDataSet* bkg = new RooDataSet("bkg", "bkg", *ws_->var("Mass"));
      bkg->append(*data_rare);
      bkg->append(*data_comb);
      fit_pdf("total", bkg, printlevel);
      if (ws_->var("N_bs")->getVal() >= estimate_bs) p_bs++;
      if (!SM_) if (ws_->var("N_bs")->getVal() >= estimate_bd) p_bd++;
    }
  }

  pdf_toy_ = pdf_toy;
  fit_pulls();
  
  double p_value_bs = p_bs / NExp;
  cout << "p-value bs = " << p_value_bs << endl;
  if (!SM_) {
    double p_value_bd = p_bd / NExp;
    cout << "p-value bd = " << p_value_bd << endl;
  }
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

  parse_estimate();
  ws_->var("N_bs")->setVal(estimate_bs);
  ws_->var("N_bd")->setVal(estimate_bd);
  ws_->var("N_rare")->setVal(estimate_rare);
  ws_->var("N_comb")->setVal(estimate_comb);

  pdf_name = "pdf_ext_total";
  if (pdf_toy != "total") {
    pdf_name = "pdf_ext_" + define_pdf_sum(pdf_toy);
  }

  RooMCStudy* mcstudy = new RooMCStudy( *ws_->pdf(pdf_name.c_str()), *ws_->var("Mass"), Binned(kFALSE), Silence(), Extended(1), FitOptions(Save(kTRUE), PrintEvalErrors(0))) ;
  cout << pdf_toy << endl;
//  if(!pdf_toy.compare("bs")) {
//    RooDLLSignificanceMCSModule sigModule(*ws_->var("N_bs"), 0);
//    mcstudy->addModule(sigModule) ;
//  }
  cout << "toy MC on pdf: " << endl;
  ws_->pdf(pdf_name.c_str())->Print();
  mcstudy->generateAndFit(NExp) ;
  mcstudy->Print();
  mcstudy->fitParDataSet().Print();

  // Make plots of the distributions of mean, the error on mean and the pull of mean
  RooPlot* frame1_bs = mcstudy->plotParam(*ws_->var("N_bs"), Bins(20)) ;
  RooPlot* frame2_bs = mcstudy->plotError(*ws_->var("N_bs"), Bins(20), Range(0., 8.)) ;
  RooPlot* frame3_bs = mcstudy->plotPull(*ws_->var("N_bs"), Bins(20), FitGauss(kTRUE)) ;
  TCanvas* canvas_bs = new TCanvas("canvas_bs", "canvas_bs", 900, 500) ;
  canvas_bs->Divide(3,1) ;
  canvas_bs->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1_bs->GetYaxis()->SetTitleOffset(1.4) ; frame1_bs->Draw() ;
  canvas_bs->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_bs->GetYaxis()->SetTitleOffset(1.4) ; frame2_bs->Draw() ;
  canvas_bs->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3_bs->GetYaxis()->SetTitleOffset(1.4) ; frame3_bs->Draw() ;
  string address = "fig/RooMCStudy_bs_" + pdf_toy + ".gif";
  canvas_bs->Print(address.c_str());
  address = "fig/RooMCStudy_bs_" + pdf_toy + ".pdf";
  canvas_bs->Print(address.c_str());
  delete frame1_bs;
  delete frame2_bs;
  delete frame3_bs;
  delete canvas_bs;
//  if(!pdf_toy.compare("bs")) {
//    TCanvas* sig_c = new TCanvas("sig_c", "sig_c", 600, 600);
//    TH1* sig_h = mcstudy->fitParDataSet().createHistogram("significance_nullhypo_N_bs");
//    cout << ">>>>>>>" << sig_h->GetXaxis()->GetTitle() << endl;
//    sig_h->Draw();
//    sig_c->Print("./fig/RooMCStudy_sig.gif");
//    sig_c->Print("./fig/RooMCStudy_sig.pdf");
//    delete sig_c;
//  }

  if (!SM_) {
    size_t found;
    found = pdf_name.find("bd");
    if (found != string::npos) {
      RooPlot* frame1_bd = mcstudy->plotParam(*ws_->var("N_bd"), Bins(20)) ;
      RooPlot* frame2_bd = mcstudy->plotError(*ws_->var("N_bd"), Bins(20)) ;
      RooPlot* frame3_bd = mcstudy->plotPull(*ws_->var("N_bd"), Bins(20), FitGauss(kTRUE)) ;
      TCanvas* canvas_bd = new TCanvas("canvas_bd", "canvas_bd", 900, 900) ;
      canvas_bd->Divide(3,1) ;
      canvas_bd->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1_bd->GetYaxis()->SetTitleOffset(1.4) ; frame1_bd->Draw() ;
      canvas_bd->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_bd->GetYaxis()->SetTitleOffset(1.4) ; frame2_bd->Draw() ;
      canvas_bd->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3_bd->GetYaxis()->SetTitleOffset(1.4) ; frame3_bd->Draw() ;
      address = "fig/RooMCStudy_bd_" + pdf_toy + ".gif";
      canvas_bd->Print(address.c_str());
      address = "fig/RooMCStudy_bd_" + pdf_toy + ".pdf";
      canvas_bd->Print(address.c_str());
      delete frame1_bd;
      delete frame2_bd;
      delete frame3_bd;
      delete canvas_bd;
    }
  }
  return;
}

void pdf_toyMC::pvalue(int NExp) {

  ws_->var("N_bs")->setVal((double)estimate_bs);
  ws_->var("N_bd")->setVal((double)estimate_bd);
  ws_->var("N_rare")->setVal((double)estimate_rare);
  ws_->var("N_comb")->setVal((double)estimate_comb);

  ws_->defineSet("obs", "Mass");
  ws_->defineSet("poi", "N_bs");

  /// RooRandom::randomGenerator()->SetSeed(1231);

  RooDataSet *dataset = new RooDataSet("dataset", "dataset", *ws_->set("obs"));
  RooDataSet* data_bs   = ws_->pdf("pdf_bs")->generate(*ws_->var("Mass"), NumEvents((int)estimate_bs), Extended(1));
  RooDataSet* data_bd   = ws_->pdf("pdf_bd")->generate(*ws_->var("Mass"), NumEvents((int)estimate_bd), Extended(1));
  RooDataSet* data_rare = ws_->pdf("pdf_rare")->generate(*ws_->var("Mass"), NumEvents((int)estimate_rare), Extended(1));
  RooDataSet* data_comb = ws_->pdf("pdf_comb")->generate(*ws_->var("Mass"), NumEvents((int)estimate_comb), Extended(1));

  dataset->append(*data_bs);
  dataset->append(*data_bd);
  dataset->append(*data_rare);
  dataset->append(*data_comb);
  dataset->Print();

  RooPlot* plot = ws_->var("Mass")->frame();
  dataset->plotOn(plot, Binning(20));
  TCanvas* plot_c = new TCanvas("plot_c", "plot_c", 600, 600);
  plot->Draw();
  plot_c->Print("./fig/MCdataset.gif");
  plot_c->Print("./fig/MCdataset.pdf");
  delete plot;
  delete plot_c;

  using namespace RooStats;

  ModelConfig model;
  model.SetWorkspace(*ws_);
  model.SetPdf(*ws_->pdf("pdf_ext_total"));
  //RooDataSet* data = pdf_ext_total->generate(*Mass, Extended(1));
  ProfileLikelihoodCalculator plc;
  plc.SetData( *dataset);
  plc.SetModel(model);
  RooArgSet poi(*ws_->var("N_bs"));
  poi.setRealValue("N_bs",0);
  plc.SetNullParameters( poi);
  HypoTestResult* htr = plc.GetHypoTest();
  cout << "The p-value for the null is " << htr->NullPValue() << endl;
  cout << "The significance for the null is " << htr->Significance() << endl;
  //return;
  //dataset = ws_->pdf("pdf_ext_total")->generate(*ws_->var("Mass"), Extended(1));
  //dataset->add(*ws_->set("obs"));

//  ProofConfig* pc = NULL;
//  pc = new ProofConfig(*ws_, 2, "workers=2", kTRUE); // machine with 2 cores

  ModelConfig* H0 = new ModelConfig("H0", "background only hypothesis", ws_);
  H0->SetPdf(*ws_->pdf("pdf_ext_total"));
  H0->SetParametersOfInterest(*ws_->set("poi"));
  H0->SetObservables(*ws_->set("obs"));
  ws_->var("N_bs")->setVal(0.0);
  H0->SetSnapshot(*ws_->set("poi"));

  ModelConfig* H1 = new ModelConfig("H1", "background + signal hypothesis", ws_);
  H1->SetPdf(*ws_->pdf("pdf_ext_total"));
  H1->SetParametersOfInterest(*ws_->set("poi"));
  H1->SetObservables(*ws_->set("obs"));
  ws_->var("N_bs")->setVal((double)estimate_bs);
  H1->SetSnapshot(*ws_->set("poi"));

  ProfileLikelihoodTestStat pl_ts(*ws_->pdf("pdf_ext_total"));
  pl_ts.SetOneSidedDiscovery(true);
  ToyMCSampler *mcSampler_pl = new ToyMCSampler(pl_ts, NExp);
//  if(pc) mcSampler_pl->SetProofConfig(pc);

  FrequentistCalculator frequCalc(*dataset, *H1,*H0, mcSampler_pl); // null = bModel interpreted as signal, alt = s+b interpreted as bkg
  //frequCalc.SetToys(NExp, 1);

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  HypoTestResult *htr_pl = frequCalc.GetHypoTest();
  RooMsgService::instance().cleanup();

  cout << endl << "ProfileLikelihoodTestStat + frequentist" << endl;
  htr_pl->Print();
  cout << "The p-value for the alternative is " << htr_pl->AlternatePValue() << endl;

  TCanvas c("canvas_roostats", "canvas_roostats", 600, 600);
  HypoTestPlot *htp = new HypoTestPlot(*htr_pl, 30); // 30 bins
  htp->SetLogYaxis(true);
  htp->Draw();
  string address = "fig/RooStats.gif";
  c.Print(address.c_str());
  address = "fig/RooStats.pdf";
  c.Print(address.c_str());

//  double data_significance = htr_pl->Significance();

////  HybridPlot* myHybridPlot = htr_pl->GetPlot("myHybridPlot","Plot of results with HybridCalculatorOriginal",100);

////  /// compute the mean expected significance from toys
////  double mean_sb_toys_test_stat = htp->GetSBmean();
////  htr_pl->SetDataTestStatistics(mean_sb_toys_test_stat);
////  double toys_significance = htr_pl->Significance();

////  std::cout << " - significance of data  = " << data_significance << std::endl;
////  std::cout << " - mean significance of toys  = " << toys_significance << std::endl;

  return;
}

void pdf_toyMC::set_ws(RooWorkspace *ws) {
  ws_ = ws;
  RooConstVar * SM_constr = (RooConstVar *) ws->obj("SM_Bs_over_Bd");
  if (SM_constr) {
    SM_ = true;
  }
  else SM_ = false;
}
