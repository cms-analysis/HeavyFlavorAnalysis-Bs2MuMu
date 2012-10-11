#include "pdf_analysis.h"

pdf_analysis::pdf_analysis(bool print, string meth, string ch_s, string range, bool BF, bool SM, bool bd_constr, bool simul, bool simulbdt, bool pee_, bool bdt_fit) {
  cout << "analysis constructor" << endl;
  print_ = print;
  meth_ = meth;
  ch_s_ = ch_s;
  ch_i_ = atoi(ch_s_.c_str());
  range_ = range;
  SM_ = SM;
  bd_constr_ = bd_constr;
  channels = 1;
  channels_bdt = 1;
  verbosity = 1;
  old_tree = false;

  no_legend = false;
  channel = atoi(ch_s.c_str());

  pee = pee_;
  simul_ = simul;
  simul_bdt_ = simulbdt;
  bdt_fit_ = bdt_fit;
  BF_ = BF;

}

void pdf_analysis::initialize () {
  cout << "inizialization" << endl;

  source.resize(4);
  source[0] = "bs";
  source[1] = "bd";
  source[2] = "rare";
  source[3] = "comb";

  ws_ = new RooWorkspace("ws", "ws");
  Mass = new RooRealVar("Mass", "Candidate invariant mass", 4.90, 5.90, "GeV/c^{2}");
  ws_->import(*Mass);
  ws_->var("Mass")->setRange("sb_lo", 4.90, 5.20);
  ws_->var("Mass")->setRange("blind", 5.20, 5.45);
  ws_->var("Mass")->setRange("sb_hi", 5.45, 5.90);
  ws_->var("Mass")->setRange("overall", 4.90, 5.90);

  eta = new RooRealVar("eta", "Candidate pseudorapidity", -2.4, 2.4, "");
  ws_->import(*eta);
  ws_->var("eta")->setRange("eta_all", -2.4, 2.4);
  ws_->var("eta")->setRange("eta_barrel", -1.4, 1.4);
  ws_->var("eta")->setRange("eta_neg_endcap", -2.4, -1.4);
  ws_->var("eta")->setRange("eta_pos_endcap", 1.4, 2.4);

  m1eta = new RooRealVar("m1eta", "first muon pseudorapidity", -2.4, 2.4, "");
  ws_->import(*m1eta);
  m2eta = new RooRealVar("m2eta", "second muon pseudorapidity", -2.4, 2.4, "");
  ws_->import(*m2eta);
  weight = new RooRealVar("weight", "event weight", 0., 1000000.);
  ws_->import(*weight);
  bdt  = new RooRealVar("bdt", "bdt cut", -1., 1.);
  ws_->import(*bdt);

  channels_cat = new RooCategory("etacat", "eta channels");
  for (int i = 0; i < 4; i++) {
    channels_cat->defineType(Form("etacat_%d", i), i);
  }
  ws_->import(*channels_cat);
  bdt_cat = new RooCategory("bdtcat", "bdt channels");
  for (int i = 0; i < 3; i++) {
    bdt_cat->defineType(Form("bdtcat_%d", i), i);
  }
  ws_->import(*bdt_cat);
  RooArgSet cat_set(*ws_->cat("etacat"), *ws_->cat("bdtcat"));
  super_cat = new RooSuperCategory("super_cat", "super_cat", cat_set);
  ws_->import(*super_cat);

  MassRes = new RooRealVar("MassRes", "mass resolution", 0.02, 0.15, "GeV/c^{2}");
  ws_->import(*MassRes);

  obs = new RooArgSet(*ws_->var("Mass"), *ws_->var("bdt"), "obs");
  //ws_->import(*obs, RecycleConflictNodes());

  getBFnumbers("anaBmm.plotResults.default-11.tex", "input/external_numbers.txt");

}

void pdf_analysis::define_pdfs () {

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < channels_bdt; j++) {
      define_bs(i, j);
      define_bd(i, j);
      define_peaking(i, j);
      define_nonpeaking(i, j);
      define_rare(i, j);
      define_comb(i, j);

      define_bf(i, j);

      define_total_extended(i, j);
    }
  }
}

void pdf_analysis::fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error, bool hesse) {

  pdf_name = "pdf_" + pdf;
  if (extended) pdf_name = "pdf_ext_" + pdf;
  
  RooAbsData* subdata;
  if (!simul_bdt_) subdata = data->reduce(Form("etacat==etacat::etacat_%d", channel));
  else subdata = data->reduce(Form("etacat==etacat::etacat_%d&&bdtcat==bdtcat::bdtcat_%d", channel, channel_bdt));
  rds_ = subdata;
  string rdh_name = subdata->GetName();
  cout << "fitting " << rdh_name << endl;
  subdata->Print();
  cout << " in range " << range_ << " with " << pdf_name << ":" << endl;
  ws_->pdf( pdf_name.c_str())->Print();
  if (!pee) {
    RFR = ws_->pdf( pdf_name.c_str())->fitTo(*subdata, Extended(extended), SumW2Error(sumw2error), NumCPU(2), Hesse(hesse), Save());
    if (print_) print(subdata);
  }
  else {
    cout << "WARNING: range option does not work with pee" << endl;
    RFR = ws_->pdf( pdf_name.c_str())->fitTo(*subdata, ConditionalObservables(*ws_->var("MassRes")), Extended(extended), SumW2Error(sumw2error), NumCPU(2), Hesse(hesse), Save());
    if (print_) print(subdata, pdf);
  }
  set_pdf_constant(pdf_name);
}

void pdf_analysis::set_pdf_constant(string name) {
  RooArgSet * set = ws_->pdf(name.c_str())->getVariables();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (!(name == "Mass") && !(name == "Bd_over_Bs") && !(name == "MassRes") && !(name == "etacat")) {
      size_t found;
      found = name.find("N");
      if (found == string::npos) ws_->var(var_Obj->GetName())->setConstant(1);
    }
  }
}

void pdf_analysis::define_bs(int i, int j) {

  RooRealVar N_bs(name("N_bs", i, j), "N_bs", 0, 10000);
  ws_->import(N_bs);

  RooRealVar Mean_bs(name("Mean_bs", i, j), "Mean_bs", 5.35, 5.32, 5.4);
  RooRealVar Sigma_bs(name("Sigma_bs", i, j), "Sigma_bs", 0.02, 0.005, 0.2);
  RooRealVar Alpha_bs(name("Alpha_bs", i, j), "Alpha_bs", 2.8, 0.1, 3.0);
  RooRealVar Enne_bs(name("Enne_bs", i, j), "Enne_bs", 1., 0., 10.);

  if (!pee) {
    RooRealVar Sigma2_bs(name("Sigma2_bs", i, j), "Sigma2_bs", 0.04, 0.005, 0.2);
    RooRealVar CoeffGauss_bs(name("CoeffGauss_bs", i, j), "CoeffGauss_bs", 0.5, 0., 1.);

    RooGaussian Gau_bs(name("Gau_bs", i, j), "Gau_bs", *ws_->var("Mass"), Mean_bs, Sigma_bs);
    RooCBShape CB_bs(name("CB_bs", i, j), "CB_bs", *ws_->var("Mass"), Mean_bs, Sigma2_bs, Alpha_bs, Enne_bs);
    if (!bdt_fit_) {
      RooAddPdf pdf_bs(name("pdf_bs", i, j), "pdf_bs", RooArgList(Gau_bs, CB_bs),  CoeffGauss_bs);
      ws_->import(pdf_bs);
    }
    else {
      RooAddPdf pdf_bs_mass(name("pdf_bs_mass", i, j), "pdf_bs_mass", RooArgList(Gau_bs, CB_bs),  CoeffGauss_bs);
      RooProdPdf pdf_bs(name("pdf_bs", i, j), "pdf_bs", pdf_bs_mass, *ws_->pdf(name("bdt_pdf_bs", i, j)));
      ws_->import(pdf_bs);
    }
  }
  else {
    RooRealVar PeeK_bs(name("PeeK_bs", i, j), "PeeK_bs", 1., 0.1, 10.);
    RooFormulaVar SigmaRes_bs(name("SigmaRes_bs", i, j), "@0*@1", RooArgList(*ws_->var("MassRes"), PeeK_bs));
    ws_->import(SigmaRes_bs);
    RooCBShape CB_bs(name("CB_bs", i, j), "CB_bs", *ws_->var("Mass"), Mean_bs, *ws_->function(name("SigmaRes_bs", i, j)), Alpha_bs, Enne_bs);
    if (!bdt_fit_) {
      RooProdPdf pdf_bs (name("pdf_bs", i, j), "pdf_bs", *ws_->pdf(name("MassRes_pdf_bs", i, j)), Conditional(CB_bs, *ws_->var("Mass")));
      //    RooProdPdf pdf_bs (name("pdf_bs", i, j), "pdf_bs", *ws_->var("MassRes"), Conditional(CB_bs, *ws_->var("Mass")));
      ws_->import(pdf_bs);
    }
    else {
      RooProdPdf pdf_bs_mass(name("pdf_bs_mass", i, j), "pdf_bs_mass", *ws_->pdf(name("MassRes_pdf_bs", i, j)), Conditional(CB_bs, *ws_->var("Mass")));
      RooProdPdf pdf_bs(name("pdf_bs", i, j), "pdf_bs", pdf_bs_mass, *ws_->pdf(name("bdt_pdf_bs", i, j)));
      ws_->import(pdf_bs);
    }
  }
}

void pdf_analysis::define_bd(int i, int j) {

  RooRealVar Mean_bd(name("Mean_bd", i, j), "Mean_bd", 5.25, 5.20, 5.29);
  RooRealVar Sigma_bd(name("Sigma_bd", i, j), "Sigma_bd", 0.02, 0.005, 0.2);
  RooRealVar Alpha_bd(name("Alpha_bd", i, j), "Alpha_bd", 2.8, 0.1, 3.0);
  RooRealVar Enne_bd(name("Enne_bd", i, j), "Enne_bd", 1., 0., 10.);

  if (!pee) {
    RooRealVar Sigma2_bd(name("Sigma2_bd", i, j), "Sigma2_bd", 0.04, 0.005, 0.2);
    RooRealVar CoeffGauss_bd(name("CoeffGauss_bd", i, j), "CoeffGauss_bd", 0.5, 0., 1.);

    RooGaussian Gau_bd(name("Gau_bd", i, j), "Gau_bd", *ws_->var("Mass"), Mean_bd, Sigma_bd);
    RooCBShape CB_bd(name("CB_bd", i, j), "CB_bd", *ws_->var("Mass"), Mean_bd, Sigma2_bd, Alpha_bd, Enne_bd);
    if (!bdt_fit_) {
      RooAddPdf pdf_bd(name("pdf_bd", i, j), "pdf_bd", RooArgList(Gau_bd, CB_bd),  CoeffGauss_bd);
      ws_->import(pdf_bd);
    }
    else {
      RooAddPdf pdf_bd_mass(name("pdf_bd_mass", i, j), "pdf_bd_mass", RooArgList(Gau_bd, CB_bd),  CoeffGauss_bd);
      RooProdPdf pdf_bd(name("pdf_bd", i, j), "pdf_bd", pdf_bd_mass, *ws_->pdf(name("bdt_pdf_bd", i, j)));
      ws_->import(pdf_bd);
    }
  }
  else {
    RooRealVar PeeK_bd(name("PeeK_bd", i, j), "PeeK_bd", 1., 0.1, 10.);
    RooFormulaVar SigmaRes_bd(name("SigmaRes_bd", i, j), "@0*@1", RooArgList(*ws_->var("MassRes"), PeeK_bd));
    ws_->import(SigmaRes_bd);
    RooCBShape CB_bd(name("CB_bd", i, j), "CB_bd", *ws_->var("Mass"), Mean_bd, *ws_->function(name("SigmaRes_bd", i, j)), Alpha_bd, Enne_bd);
    if (!bdt_fit_) {
      RooProdPdf pdf_bd (name("pdf_bd", i, j), "pdf_bd", *ws_->pdf(name("MassRes_pdf_bd", i, j)), Conditional(CB_bd, *ws_->var("Mass")));
      //    RooProdPdf pdf_bd (name("pdf_bd", i, j), "pdf_bd", *ws_->var("MassRes"), Conditional(CB_bd, *ws_->var("Mass")));
      ws_->import(pdf_bd);
    }
    else {
      RooProdPdf pdf_bd_mass(name("pdf_bd_mass", i, j), "pdf_bd_mass", *ws_->pdf(name("MassRes_pdf_bd", i, j)), Conditional(CB_bd, *ws_->var("Mass")));
      RooProdPdf pdf_bd(name("pdf_bd", i, j), "pdf_bd", pdf_bd_mass,*ws_->pdf(name("bdt_pdf_bd", i, j)));
      ws_->import(pdf_bd);
    }
  }

  if (SM_) {
    RooConstVar SM_Bd_over_Bs("SM_Bd_over_Bs", "SM_Bd_over_Bs", ratio_);
    ws_->import(SM_Bd_over_Bs);
    RooFormulaVar N_bd_constr(name("N_bd_constr", i, j), "N_bd_constr", "@0*@1", RooArgList(*ws_->var(name("N_bs", i, j)), *(RooConstVar*)ws_->obj("SM_Bd_over_Bs")));
    ws_->import(N_bd_constr);
  }
  else if (bd_constr_) {
    RooRealVar Bd_over_Bs("Bd_over_Bs", "Bd_over_Bs", 0., 10., "");
    ws_->import(Bd_over_Bs);
    RooFormulaVar N_bd_constr(name("N_bd_constr", i, j), "N_bd_constr", "@0*@1", RooArgList(*ws_->var(name("N_bs", i, j)), *ws_->var("Bd_over_Bs")));
    ws_->import(N_bd_constr);
  }
  else {
    RooRealVar N_bd(name("N_bd", i, j), "N_bd", 0, 10000);
    ws_->import(N_bd);
  }
}

void pdf_analysis::define_peaking(int i, int j) {
  RooRealVar N_peak(name("N_peak", i, j), "N_peak", 0, 10000);
  ws_->import(N_peak);

  if (old_tree) {
    RooRealVar Mean_peak(name("Mean_peak", i, j), "Mean_peak", 5.25, 5.20, 5.3);
    RooRealVar Sigma_peak(name("Sigma_peak", i, j), "Sigma_peak", 0.050, 0.01, 0.20);
    RooGaussian pdf_peak(name("pdf_peak", i, j), "pdf_peak", *ws_->var("Mass"), Mean_peak, Sigma_peak);
    ws_->import(pdf_peak);
  }
  else {
    RooRealVar Mean_peak(name("Mean_peak", i, j), "Mean_peak", 5.1, 4.9, 5.4);
    RooRealVar Sigma_peak(name("Sigma_peak", i, j), "Sigma_peak", 0.02, 0.005, 0.2);
    RooRealVar Sigma2_peak(name("Sigma2_peak", i, j), "Sigma2_peak", 0.04, 0.005, 0.2);
    RooRealVar Alpha_peak(name("Alpha_peak", i, j), "Alpha_peak", 2.8, 0., 100.0);
    RooRealVar Enne_peak(name("Enne_peak", i, j), "Enne_peak", 1., 0., 10.);
    RooRealVar CoeffGauss_peak(name("CoeffGauss_peak", i, j), "CoeffGauss_peak", 0.5, 0., 1.);
    RooGaussian Gau_peak(name("Gau_peak", i, j), "Gau_peak", *ws_->var("Mass"), Mean_peak, Sigma_peak);
    RooCBShape CB_peak(name("CB_peak", i, j), "CB_peak", *ws_->var("Mass"), Mean_peak, Sigma2_peak, Alpha_peak, Enne_peak);

    if (!pee) {
      if (!bdt_fit_) {
        RooAddPdf pdf_peak(name("pdf_peak", i, j), "pdf_peak", RooArgList(Gau_peak, CB_peak),  CoeffGauss_peak);
        ws_->import(pdf_peak);
      }
      else {
        RooAddPdf pdf_peak_mass(name("pdf_peak_mass", i, j), "pdf_peak_mass", RooArgList(Gau_peak, CB_peak),  CoeffGauss_peak);
        RooProdPdf pdf_peak(name("pdf_peak", i, j),"pdf_peak",pdf_peak_mass,*ws_->pdf(name("bdt_pdf_peak", i, j)));
        ws_->import(pdf_peak);
      }
    }
    else {
      RooAddPdf mass_peak(name("mass_peak", i, j), "mass_peak", RooArgList(Gau_peak, CB_peak),  CoeffGauss_peak);
      //RooProdPdf pdf_peak (name("pdf_peak", i, j), "pdf_peak", *ws_->var("MassRes"), Conditional(mass_peak, *ws_->var("Mass")));
      if (!bdt_fit_) {
        RooProdPdf pdf_peak (name("pdf_peak", i, j), "pdf_peak", *ws_->pdf(name("MassRes_pdf_peak", i, j)), Conditional(mass_peak, *ws_->var("Mass")));
        ws_->import(pdf_peak);
      }
      else {
        RooProdPdf pdf_peak_mass (name("pdf_peak_mass", i, j), "pdf_peak_mass", *ws_->pdf(name("MassRes_pdf_peak", i, j)), Conditional(mass_peak, *ws_->var("Mass")));
        RooProdPdf pdf_peak(name("pdf_peak", i, j),"pdf_peak",pdf_peak_mass,*ws_->pdf(name("bdt_pdf_peak", i, j)));
        ws_->import(pdf_peak);
      }
    }
  }
}

void pdf_analysis::define_nonpeaking(int i, int j) {

  RooRealVar N_semi(name("N_semi", i, j), "N_semi", 0, 10000);
  ws_->import(N_semi);

  if (old_tree) {
    RooRealVar m0_semi(name("m0_semi", i, j), "m0_semi", 5., 6.);
    RooRealVar c_semi(name("c_semi", i, j), "c_semi", 1., 0.1, 20);
    RooRealVar p_semi(name("p_semi", i, j), "p_semi", 0.5, 0.1, 5.);
    RooArgusBG ArgusBG(name("pdf_semi", i, j), "pdf_semi", *ws_->var("Mass"), m0_semi, c_semi, p_semi);
    ws_->import(ArgusBG);
  }
  else {
    RooRealVar C0(name("C0", i, j), "C0", 0.1, 0.001, 10., "");
    RooRealVar C1(name("C1", i, j), "C1", 0.1, 0.001, 10., "");
    RooRealVar C2(name("C2", i, j), "C2", 0., -2., 2., "");
    RooRealVar C3(name("C3", i, j), "C3", 0., -2., 2., "");
    RooRealVar C4(name("C4", i, j), "C4", 0., -2., 2., "");
    RooRealVar tau(name("tau", i, j), "tau", -5,-20.,-0.01, "");
    RooArgList poly_coeffs(C0, C1, C2, C3);
    RooChebychev poly(name("poly", i, j), "poly", *ws_->var("Mass"), poly_coeffs);
    RooExponential expo(name("expo", i, j), "expo", *ws_->var("Mass"), tau);
    if (!pee) {
      if (!bdt_fit_) {
        RooProdPdf pdf_semi(name("pdf_semi", i, j), "pdf_semi", expo, poly);
        ws_->import(pdf_semi);
      }
      else {
        RooProdPdf pdf_semi_mass(name("pdf_semi_mass", i, j), "pdf_semi_mass", expo, poly);
        RooProdPdf pdf_semi(name("pdf_semi", i, j),"pdf_semi",pdf_semi_mass,*ws_->pdf(name("bdt_pdf_semi", i, j)));
        ws_->import(pdf_semi);
      }
    }
    else {
      RooProdPdf mass_semi(name("mass_semi", i, j), "pdf_semi", expo, poly);
      //RooProdPdf pdf_semi (name("pdf_semi", i, j), "pdf_semi", *ws_->var("MassRes"), Conditional(mass_semi, *ws_->var("Mass")));
      if (!bdt_fit_) {
        RooProdPdf pdf_semi (name("pdf_semi", i, j), "pdf_semi", *ws_->pdf(name("MassRes_pdf_semi", i, j)), Conditional(mass_semi, *ws_->var("Mass")));
        ws_->import(pdf_semi);
      }
      else {
        RooProdPdf pdf_semi_mass (name("pdf_semi_mass", i, j), "pdf_semi_mass", *ws_->pdf(name("MassRes_pdf_semi", i, j)), Conditional(mass_semi, *ws_->var("Mass")));
        RooProdPdf pdf_semi(name("pdf_semi", i, j),"pdf_semi",pdf_semi_mass,*ws_->pdf(name("bdt_pdf_semi", i, j)));
        ws_->import(pdf_semi);
      }
    }
  }
}

void pdf_analysis::define_comb(int i, int j) {

  RooRealVar N_comb(name("N_comb", i, j), "N_comb", 0, 10000);
  ws_->import(N_comb);

  if (!pee) {
    if (!bdt_fit_) {
      RooUniform pdf_comb(name("pdf_comb", i, j), "N_comb", *ws_->var("Mass"));
      ws_->import(pdf_comb);
    }
    else {
      RooUniform pdf_comb_mass(name("pdf_comb_mass", i, j), "N_comb", *ws_->var("Mass"));
      RooProdPdf pdf_comb(name("pdf_comb", i, j),"pdf_comb",pdf_comb_mass,*ws_->pdf(name("bdt_pdf_comb", i, j)));
      ws_->import(pdf_comb);
    }
  }
  else  {
    RooUniform mass_comb(name("mass_comb", i, j), "N_comb", *ws_->var("Mass"));
    //RooProdPdf pdf_comb (name("pdf_comb", i, j), "pdf_comb", *ws_->var("MassRes"), Conditional(mass_comb, *ws_->var("Mass")));
    if (!bdt_fit_) {
      RooProdPdf pdf_comb (name("pdf_comb", i, j), "pdf_comb", *ws_->pdf(name("MassRes_pdf_comb", i, j)), Conditional(mass_comb, *ws_->var("Mass")));
      ws_->import(pdf_comb);
    }
    else {
      RooProdPdf pdf_comb_mass (name("pdf_comb_mass", i, j), "pdf_comb_mass", *ws_->pdf(name("MassRes_pdf_comb", i, j)), Conditional(mass_comb, *ws_->var("Mass")));
      RooProdPdf pdf_comb(name("pdf_comb", i, j),"pdf_comb", pdf_comb_mass, *ws_->pdf(name("bdt_pdf_comb", i, j)));
      ws_->import(pdf_comb);
    }
  }
}

void pdf_analysis::define_signals(int i, int j) {

  RooRealVar N_signals(name("N_signals", i, j), "N_signals", 0., 1000);
  ws_->import(N_signals);

  RooRealVar bsfrac_signals(name("bsfrac_signals", i, j), "bsfrac_signals", 0.5, 0.0, 1.0);
  RooAddPdf pdf_signals(name("pdf_signals", i, j), "pdf_signals", RooArgSet(*ws_->pdf(name("pdf_bs", i, j)), *ws_->pdf(name("pdf_bd", i, j)) ), bsfrac_signals);
  ws_->import(pdf_signals);
}

void pdf_analysis::define_rare(int i, int j) {
  
  RooRealVar N_rare(name("N_rare", i, j), "N_rare", 0, 10000);
  ws_->import(N_rare);

  RooRealVar peakfrac_rare(name("peakfrac_rare", i, j), "peakfrac_rare", 0.5, 0.0, 1.0);
  RooAddPdf pdf_rare(name("pdf_rare", i, j), "pdf_rare", RooArgSet(*ws_->pdf(name("pdf_peak", i, j)), *ws_->pdf(name("pdf_semi", i, j)) ), peakfrac_rare);
  ws_->import(pdf_rare);
}

void pdf_analysis::define_rare2(RooDataHist* data, int i, int j) {
  ws_->factory("N_hist[0, 10000]");
  RooHistPdf* pdf_rare = new RooHistPdf("pdf_hist", "pdf_hist", *ws_->var("Mass"), *data, 4);
  ws_->import(*pdf_rare);
}

void pdf_analysis::define_rare3(int i, int j) {
  RooRealVar N_expo(name("N_expo3", i, j), "N_expo3", 0, 10000);
  ws_->import(N_expo);
  RooRealVar tau3(name("tau3", i, j), "tau3", -5, -20., -0.01);
  if (!pee) {
    RooExponential pdf_expo(name("pdf_expo3", i, j), "pdf_expo3", *ws_->var("Mass"), tau3);
    ws_->import(pdf_expo);
  }
  else {
    RooExponential expo3(name("expo3", i, j), "expo3", *ws_->var("Mass"), tau3);
    RooProdPdf pdf_expo(name("pdf_expo3", i, j), "pdf_expo3", *ws_->pdf(name("MassRes_pdf_rare", i, j)), Conditional(expo3, *ws_->var("Mass")));
    ws_->import(pdf_expo);
  }
  //ws_->factory("SUM::pdf_expo(expo1frac_rare[0.,1.]*Exponential::expo1(Mass,Alpha1[-10.,10.]),Exponential::expo2(Mass,Alpha2[-10.,10.]))");
}

void pdf_analysis::define_signalsrare(int i, int j) {
  
  ws_->factory("signalfrac_signalsrare[0.5, 0.0, 1.0]");
  ws_->factory("SUM::pdf_signalsrare(signalfrac_signalsrare*pdf_signals, pdf_rare)");

  if (!SM_ && !bd_constr_) ws_->factory("SUM::pdf_ext_signalsrare(N_bs*pdf_bs, N_bd*pdf_bd, N_rare*pdf_rare)");
}

void pdf_analysis::define_bkg_fractional(int i, int j) {
  
  ws_->factory("N_bkg[0, 10000]");
  ws_->factory("rarefrac_bkg[0.5, 0.0, 1.0]");
  
  ws_->factory("SUM::pdf_bkg(rarefrac_bkg*pdf_rare, pdf_comb)");
  
  RooExtendPdf* bkgExt = new RooExtendPdf("pdf_fract_ext_bkg", "bkgExt", *ws_->pdf("pdf_bkg"), *ws_->var("N_bkg"), /* range_.c_str()*/ "sb_lo,sb_hi"); /// WARNING
  ws_->import(*bkgExt);
}

void pdf_analysis::define_bkg_extended(int i, int j) {

//  RooExtendPdf rare_ext("pdf_ext_rare", "rare_ext", *ws_->pdf("pdf_rare"), *ws_->var("N_rare"), "sb_lo,sb_hi");
//  RooExtendPdf comb_ext("pdf_ext_comb", "comb_ext", *ws_->pdf("pdf_comb"), *ws_->var("N_comb"), "sb_lo,sb_hi");
//  RooAddPdf bkg_ext("pdf_ext_bkg", "pdf_ext_bkg", RooArgList(rare_ext, comb_ext));
//  ws_->import(bkg_ext);
  ws_->factory("SUM::pdf_ext_bkg(N_rare*pdf_rare, N_comb*pdf_comb)");

}


void pdf_analysis::define_total_fractional(int i, int j) {
  
  ws_->factory("N_all[0, 10000]");
  ws_->factory("signalsfrac_total[0.5, 0.0, 1.0]");
  
  ws_->factory("SUM::pdf_frac_total(signalsfrac_total*pdf_signals, pdf_bkg)");
  
  RooExtendPdf* allExt = new RooExtendPdf("pdf_frac_ext_total", "allExt", *ws_->pdf("pdf_frac_total"), *ws_->var("N_all"), range_.c_str());
  ws_->import(*allExt);
}

void pdf_analysis::define_total_extended(int i, int j) {
  RooArgList pdf_list(*ws_->pdf(name("pdf_bs", i, j)), *ws_->pdf(name("pdf_bd", i, j)), *ws_->pdf(name("pdf_rare", i, j)), *ws_->pdf(name("pdf_comb", i, j)));
  if (!BF_) {
    if (SM_ || bd_constr_) {
      RooArgList N_list(*ws_->var(name("N_bs", i, j)), *ws_->function(name("N_bd_constr", i, j)), *ws_->var(name("N_rare", i, j)), *ws_->var(name("N_comb", i, j)));
      RooAddPdf pdf_ext_total(name("pdf_ext_total", i, j), "pdf_ext_total", pdf_list, N_list);
      ws_->import(pdf_ext_total);
    }
    else {
      RooArgList N_list(*ws_->var(name("N_bs", i, j)), *ws_->var(name("N_bd", i, j)), *ws_->var(name("N_rare", i, j)), *ws_->var(name("N_comb", i, j)));
      RooAddPdf pdf_ext_total(name("pdf_ext_total", i, j), "pdf_ext_total", pdf_list, N_list);
      ws_->import(pdf_ext_total);
    }
  }
  else {
    if (SM_ || bd_constr_) {
      RooArgList N_list(*ws_->function(name("N_bs_formula", i, j)), *ws_->function(name("N_bd_constr", i, j)), *ws_->var(name("N_rare", i, j)), *ws_->var(name("N_comb", i, j)));
      RooAddPdf pdf_ext_total(name("pdf_ext_total", i, j), "pdf_ext_total", pdf_list, N_list);
      ws_->import(pdf_ext_total);
      RooArgList N_list_test(*ws_->function(name("N_bs", i, j)), *ws_->var(name("N_bd_constr", i, j)), *ws_->var(name("N_rare", i, j)), *ws_->var(name("N_comb", i, j)));
      RooAddPdf pdf_ext_total_test(name("pdf_ext_total_test", i, j), "pdf_ext_total", pdf_list, N_list_test);
      ws_->import(pdf_ext_total_test);
    }
    else {
      RooArgList N_list(*ws_->function(name("N_bs_formula", i, j)), *ws_->var(name("N_bd", i, j)), *ws_->var(name("N_rare", i, j)), *ws_->var(name("N_comb", i, j)));
      RooAddPdf pdf_ext_total(name("pdf_ext_total", i, j), "pdf_ext_total", pdf_list, N_list);
      ws_->import(pdf_ext_total);
      RooArgList N_list_test(*ws_->function(name("N_bs", i, j)), *ws_->var(name("N_bd", i, j)), *ws_->var(name("N_rare", i, j)), *ws_->var(name("N_comb", i, j)));
      RooAddPdf pdf_ext_total_test(name("pdf_ext_total_test", i, j), "pdf_ext_total", pdf_list, N_list_test);
      ws_->import(pdf_ext_total_test);
    }
  }
  return; 
}

void pdf_analysis::define_bf(int i, int j) {

  RooRealVar N_bu_var(name("N_bu", i, j), "N_bu", N_bu_val[i], 1, 10000000, "");

  RooRealVar eff_bu_var(name("eff_bu", i, j), "eff_bu", eff_bu_val[i], eff_bu_val[i]-eff_bu_err[i], eff_bu_val[i]+eff_bu_err[i], "");
  RooRealVar eff_bs_var(name("eff_bs", i, j), "eff_bs", eff_bs_val[i], eff_bs_val[i]-eff_bs_err[i], eff_bs_val[i]+eff_bs_err[i], "");

  RooRealVar fs_over_fu_var("fs_over_fu", "fs_over_fu", fs_over_fu_val, fs_over_fu_val-fs_over_fu_err, fs_over_fu_val+fs_over_fu_err, "");
  RooRealVar Bu2JpsiK_BF_var("Bu2JpsiK_BF", "Bu2JpsiK_BF", Bu2JpsiK_BF_val, Bu2JpsiK_BF_val-Bu2JpsiK_BF_err, Bu2JpsiK_BF_val+Bu2JpsiK_BF_err, "");
  RooRealVar Jpsi2MuMu_BF_var("Jpsi2MuMu_BF", "Jpsi2MuMu_BF", Jpsi2MuMu_BF_val, Jpsi2MuMu_BF_val-Jpsi2MuMu_BF_err, Jpsi2MuMu_BF_val+Jpsi2MuMu_BF_err, "");

  N_bu_var.setConstant(1);
  eff_bu_var.setConstant(1);
  eff_bs_var.setConstant(1);
  fs_over_fu_var.setConstant(1);
  Bu2JpsiK_BF_var.setConstant(1);
  Jpsi2MuMu_BF_var.setConstant(1);

  RooRealVar BF("BF", "Bs2MuMu branching fraction", Bs2MuMu_SM_BF_val, 1.0e-11, 1e-8);
  if (i == 0 && j == 0) ws_->import(BF);

  RooFormulaVar K(name("K", i, j), "K(i) * N_bs(i) / N_bu(i) = BF",  "@0*@1*@2/@3/@4", RooArgList(Bu2JpsiK_BF_var, Jpsi2MuMu_BF_var, eff_bu_var, eff_bs_var, fs_over_fu_var));
  RooFormulaVar N_bs_constr(name("N_bs_formula", i, j), "N_bs(i) = BF / K(i) * N_Bu = ",  "@0/@1*@2", RooArgList( *ws_->var("BF"), K, N_bu_var));
  cout << "channel: " << i << ";" << j << " expected Bs = " << N_bs_constr.getVal() << endl;
  ws_->import(N_bs_constr);
}

void pdf_analysis::define_simul(bool simulbdt) {
  if (!simulbdt) {
    RooSimultaneous pdf_sim("pdf_ext_simul", "simultaneous pdf", *ws_->cat("etacat"));
    for (int i = 0; i < channels; i++) {
      pdf_sim.addPdf(*ws_->pdf(name("pdf_ext_total", i)), name("etacat", i));
    }
    ws_->import(pdf_sim);
    if (BF_) {
      RooSimultaneous pdf_sim_test("pdf_ext_simul_test", "simultaneous pdf", *ws_->cat("etacat"));
      for (int i = 0; i < channels; i++) {
        pdf_sim_test.addPdf(*ws_->pdf(name("pdf_ext_total_test", i)), name("etacat", i));
      }
      ws_->import(pdf_sim_test);
    }
  }
  else {
    RooSuperCategory* rsc = dynamic_cast<RooSuperCategory*> (ws_->obj("super_cat"));
    RooSimultaneous pdf_sim("pdf_ext_simul", "simultaneous pdf", *rsc);
    for (int g = 0; g< channels*channels_bdt ;g++) {
      rsc->setIndex(g);
      RooArgSet icl = rsc->inputCatList();
      RooCategory* eta_c = (RooCategory*)icl.find("etacat");
      RooCategory* bdt_c = (RooCategory*)icl.find("bdtcat");
      Int_t eta_i = eta_c->getIndex();
      Int_t bdt_i = bdt_c->getIndex();
      cout << rsc->getLabel() << " " <<  eta_i << " " <<  bdt_i<< endl;
      pdf_sim.addPdf(*ws_->pdf(Form("pdf_ext_total_%d_%d", eta_i, bdt_i)), rsc->getLabel());
    }
//    for (int i = 0; i < channels; i++) {
//      for (int j = 0; j < channels_bdt; j++) {
//        pdf_sim.addPdf(*ws_->pdf(name("pdf_ext_total", i, j)), Form("{etacat_%d;bdtcat_%d}", i, j));
//      }
//    }
    ws_->import(pdf_sim/*, RecycleConflictNodes()*/);
  }
}

string pdf_analysis::define_pdf_sum(string name, int i) {

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
  found = name.find("expo3");
  if (found != string::npos) pdfs.push_back("expo3");
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
    if (pdfs[i]=="hist" || pdfs[i]=="expo3") pdf_sum += "rare";
    else pdf_sum += pdfs[i];
    pdf_sum += "*pdf_";
    pdf_sum += pdfs[i];
    if (i != pdfs.size() -1) pdf_sum += ",";
  }
  pdf_sum += ")";
  cout << "formed pdf: " << pdf_sum << endl;
  ws_->factory(pdf_sum.c_str());
  return (title);
}

void pdf_analysis::print(RooAbsData* data, string output) {
  int colors[11] = {632, 400, 616, 432, 800, 416, 820, 840, 860, 880, 900};
  RooAbsData* subdata_res;
//  if (pee) subdata_res = ws_->data(Form("MassRes_rdh_%s", output.c_str()))->reduce(Form("etacat==etacat::etacat_%d", channel));
  if (pee) subdata_res = (RooAbsData*)ws_->data(Form("MassRes_rdh_%s", output.c_str()))->Clone();
  RooPlot *rp = ws_->var("Mass")->frame();
  RooPlot *rp_bdt = ws_->var("bdt")->frame();
  data->plotOn(rp, Binning(20));
  if (bdt_fit_) data->plotOn(rp_bdt, Binning(20));
  if (!pee) {
    ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue));
    if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, LineColor(kBlue));
  }
  else {
    if (!simul_) {
      ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue), ProjWData(RooArgSet(*ws_->var("MassRes")), *subdata_res, kFALSE));
      if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, LineColor(kBlue), ProjWData(RooArgSet(*ws_->var("MassRes")), *subdata_res, kFALSE));
    }
    else {
      ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue), ProjWData(RooArgSet(*ws_->var("MassRes")), *subdata_res, kFALSE));
      if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, LineColor(kBlue), ProjWData(RooArgSet(*ws_->var("MassRes")), *subdata_res, kFALSE));
    }
    TH1* mass_eta_h;
    mass_eta_h = ws_->pdf(pdf_name.c_str())->createHistogram("fit", *ws_->var("Mass"), Binning(50), YVar(*ws_->var("MassRes"), Binning(30))) ;
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
    RooPlot *rp_res = ws_->var("MassRes")->frame();
    *ws_->data(Form("MassRes_rdh_%s", output.c_str()))->plotOn(rp_res);
    rp_res->Draw();
    cetad->Print( (get_address("MassEta", pdf_name) + ".gif").c_str());
    cetad->Print( (get_address("MassEta", pdf_name) + ".pdf").c_str());

    delete cetad;
    delete frame;
    delete mass_eta_h;
    delete rp_res;
  }
  if(!no_legend) {
    ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));
    if (bdt_fit_) ws_->pdf(pdf_name.c_str())->paramOn(rp_bdt, Layout(0.50, 0.9, 0.9));
  }
  
  //components
  RooArgSet * set = ws_->pdf(pdf_name.c_str())->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  int i = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (i > 11) i = 0;
      size_t found1 = pdf_name.find("total");
      //cout << name << endl;
      if (found1 == string::npos) {
        if (!pee) {
          ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
        }
        else {
          size_t found2 = pdf_name.find("SigmaRes");
          if (found2 == string::npos) {
            ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *subdata_res, kFALSE), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
          }
        }
      }
      else {
        if (!pee) {
          if (name=="pdf_bs" || name=="pdf_bd" || name=="pdf_rare" || name=="pdf_comb") {
            ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
            if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, Components(*ws_->pdf(var_Obj->GetName())), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
          }
        }
        else {
          size_t found2 = pdf_name.find("SigmaRes");
          if (found2 == string::npos) {
            if (name=="pdf_bs" || name=="pdf_bd" || name=="pdf_rare" || name=="pdf_comb") {
              ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *subdata_res, kFALSE), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
              if (bdt_fit_) ws_->pdf(pdf_name.c_str())->plotOn(rp_bdt, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *subdata_res, kFALSE), LineColor(colors[i]),  LineStyle(1), LineWidth(2), Range(range_.c_str()));
            }
          }
        }
      }
      i++;
    }
  }
  
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  c->Print((get_address(pdf_name, "") + ".pdf").c_str());
  c->Print((get_address(pdf_name, "") + ".gif").c_str());
  delete rp;
  delete c;

  if (bdt_fit_) {
    TCanvas* c_bdt = new TCanvas("c_bdt", "c_bdt", 600, 600);
    rp_bdt->Draw();
    c->Print((get_address("BDT", pdf_name) + ".pdf").c_str());
    c->Print((get_address("BDT", pdf_name) + ".gif").c_str());
    delete c_bdt;
  }
  delete rp_bdt;

  return;
}

void pdf_analysis::simsplit() {
  cout << "simsplitting" << endl;
  simul_ = true;
  ostringstream splitter;
  ostringstream splitted;
  splitter << "SIMCLONE::pdf_ext_simul(pdf_ext_total, $SplitParam({";
  RooArgSet* allvars = ws_->pdf("pdf_ext_total")->getVariables();
  TIterator* it = allvars->createIterator();
  RooRealVar* var_Obj = 0;
  int cont = 0;
  while ( (var_Obj = (RooRealVar*)it->Next()) ) {
    string name = var_Obj->GetName();
    if ( !(name == "Mass") && !(name == "Bd_over_Bs") && !(name == "SM_Bd_over_Bs") && !(name == "bdt") && !(name == "MassRes")) {
      if (cont != 0) {
        splitter << ", ";
        splitted << ", ";
      }
      splitter << name;
      splitted << name;
      cont++;
    }
  }
  /// splitter << ",etapdf_bs"; /// warning
  splitter << "}, {etacat[";
  for (int i = 0; i < 2; i++) {
    if (i != 0) splitter << ",";
    splitter << "etacat_" << i ;
  }
  splitter << "],bdtcat[";
  for (int i = 0; i < 3; i++) {
    if (i != 0) splitter << ",";
    splitter << "bdtcat_" << i ;
  }
  splitter << "]}))";
  cout << splitter.str() << endl;
  ws_->factory(splitter.str().c_str());
  //RooSimWSTool sct(*ws_);

  ws_->Print();
  cout << splitted.str() << endl;
  //RooSimultaneous* model_sim2 = sct.build("pdf_ext_simul", "pdf_ext_total", SplitParam(splitted.str().c_str(), "etacat,bdtcat"));
  //model_sim2->Print();


}

double pdf_analysis::getErrorHigh(RooRealVar *var) {
  double value = var->getVal();
  if (value == var->getMax()) return 0;
  double error = var->getErrorHi();
  if (error == 0) return var->getError();
  else return error;
}

double pdf_analysis::getErrorLow(RooRealVar *var) {
  double value = var->getVal();
  if (value == var->getMin() || (strcmp(var->GetName(),"BF") && (value <= 0.0001 && value >= -0.0001))) return 0;
  double error = var->getErrorLo();
  if (error == 0) {
    if (value - var->getError() < 0.0) return (-1 * value);
    else return (-1*var->getError());
  }
  else return error;
}

RooHistPdf* pdf_analysis::define_MassRes_pdf(RooDataSet *rds, string name) {
  RooArgList varlist(*MassRes, *weight);
  RooDataSet* subdata_res = new RooDataSet("subdata_res", "subdata_res", varlist, "weight");
  const RooArgSet* aRow;
  TH1D histo(rds->GetTitle(), rds->GetTitle(), 40, -2.4, 2.4);
  for (Int_t j = 0; j < rds->numEntries(); j++) {
    aRow = rds->get(j);
    RooRealVar* massres = (RooRealVar*)aRow->find("MassRes");
    double Weight = rds->weight();
    RooArgSet varlist_tmp_res(*massres);
    if (aRow->getCatIndex("etacat") == channel) {
      if (!simul_bdt_ || aRow->getCatIndex("bdtcat") == channel_bdt) {
        subdata_res->add(varlist_tmp_res, Weight);
        //      histo.Fill(massres->getVal(), Weight);
      }
    }
  }
  cout << "resolution entries = " <<  subdata_res->sumEntries() << endl;
  ostringstream name_rdh;
  name_rdh << "MassRes_rdh_" << name;
  if (simul_) name_rdh << "_" << channel;
  if (simul_bdt_) name_rdh << "_" << channel_bdt;
  RooDataHist *MassRes_rdh = subdata_res->binnedClone(name_rdh.str().c_str());
  //RooDataHist *MassRes_rdh = new RooDataHist(name_rdh.str().c_str(), name_rdh.str().c_str(), *ws_->var("MassRes"), &histo);
  ostringstream name_pdf;
  name_pdf << "MassRes_pdf_" << name;
  if (simul_) name_pdf << "_" << channel;
  if (simul_bdt_) name_pdf << "_" << channel_bdt;
  RooHistPdf * MassRes_rhpdf = new RooHistPdf(name_pdf.str().c_str(), name_pdf.str().c_str(), RooArgList(*ws_->var("MassRes")), *MassRes_rdh);
  ws_->import(*MassRes_rhpdf);

  return MassRes_rhpdf;
}

RooHistPdf* pdf_analysis::define_bdt_pdf(RooDataSet *rds, string name) {
  RooArgList varlist(*bdt, *weight);
  RooDataSet* subdata_bdt = new RooDataSet("subdata_bdt", "subdata_bdt", varlist, "weight");
  const RooArgSet* aRow;
  TH1D histo(rds->GetTitle(), rds->GetTitle(), 20, -1., 1.);
  for (Int_t j = 0; j < rds->numEntries(); j++) {
    aRow = rds->get(j);
    RooRealVar* BDT = (RooRealVar*)aRow->find("bdt");
    double Weight = rds->weight();
    RooArgSet varlist_tmp_bdt(*BDT);
    if (aRow->getCatIndex("etacat") == channel) {
      if (!simul_bdt_ || aRow->getCatIndex("bdtcat") == channel_bdt) {
        //      subdata_bdt->add(varlist_tmp_bdt, Weight);
        histo.Fill(BDT->getVal(), Weight);
      }
    }
  }
//  cout << "bdt entries = " <<  subdata_bdt->sumEntries() << endl;
  cout << "bdt entries = " <<  histo.GetEntries() << endl;
  ostringstream name_rdh;
  name_rdh << "bdt_rdh_" << name;
  if (simul_) name_rdh << "_" << channel;
  if (simul_bdt_) name_rdh << "_" << channel_bdt;
  //RooDataHist *bdt_rdh = subdata_bdt->binnedClone(name_rdh.str().c_str());
  RooDataHist *bdt_rdh = new RooDataHist(name_rdh.str().c_str(), name_rdh.str().c_str(), *ws_->var("bdt"), &histo);
  ostringstream name_pdf;
  name_pdf << "bdt_pdf_" << name;
  if (simul_) name_pdf << "_" << channel;
  if (simul_bdt_) name_pdf << "_" << channel_bdt;
  RooHistPdf * bdt_rhpdf = new RooHistPdf(name_pdf.str().c_str(), name_pdf.str().c_str(), RooArgList(*ws_->var("bdt")), *bdt_rdh);
  ws_->import(*bdt_rhpdf);

  return bdt_rhpdf;
}

RooDataHist pdf_analysis::getRandom_rdh() {
  TH1D* temp_h = new TH1D("temp_h", "temp_h", 100, -2.4, 2.4);
  temp_h->FillRandom("gaus", 1000);
  RooDataHist temp_rdh("temp_rdh", "temp_rdh", *ws_->var("eta"), temp_h);
  return temp_rdh;
}

const char* pdf_analysis::name(string name, int i, int j) {
  if (!simul_) return name.c_str();
  else {
    if (!simul_bdt_) return Form("%s_%d", name.c_str(), i);
    else return Form("%s_%d_%d", name.c_str(), i, j);
  }
}

void pdf_analysis::set_rare_normalization(string input, bool extended) {
  FILE *estimate_file = fopen(input.c_str(), "r");
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < channels_bdt; j++) {
      char buffer[1024];
      char cutName[128];
      float cut;
      while (fgets(buffer, sizeof(buffer), estimate_file)) {
        if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
        if (buffer[0] == '#') continue;
        sscanf(buffer, "%s %f", cutName, &cut);
        if ( (simul_ && !strcmp(cutName, Form("peakfrac_rare_%d", i))) || (!simul_ && !strcmp(cutName, Form("peakfrac_rare_%d", channel))) ) {
          ws_->var(name("peakfrac_rare", i, j))->setConstant(kFALSE);
          ws_->var(name("peakfrac_rare", i, j))->setVal(cut);
          ws_->var(name("peakfrac_rare", i, j))->setConstant(kTRUE);
          cout << name("peakfrac_rare", i, j) << " set val to " << cut << endl;
        }
        if (extended) {
          if ( (simul_ && !strcmp(cutName, Form("N_rare_%d", i))) || (!simul_ && !strcmp(cutName, Form("N_rare_%d", channel))) ) {
            ws_->var(name("N_rare", i, j))->setConstant(kFALSE);
            ws_->var(name("N_rare", i, j))->setVal(cut);
            ws_->var(name("N_rare", i, j))->setConstant(kTRUE);
            cout << name("N_rare", i, j) << " set val to " << cut << endl;
          }
        }
      }
      rewind(estimate_file);
    }
  }
  if (estimate_file) fclose(estimate_file);
}

void pdf_analysis::gen_and_fit(string pdfname) {
  pdf_name = pdfname;
  RooAbsPdf* pdf = ws_->pdf(pdf_name.c_str());
  pdf->Print();

  RooDataSet* protodata = (RooDataSet*)ws_->data(name("MassRes_rdh_bs", channel));
  protodata->Print();
  RooPlot* rp0 = ws_->var("MassRes")->frame();
  protodata->plotOn(rp0);
  TCanvas* canvas0 = new TCanvas("canvas0", "canvas0", 600, 600);
  rp0->Draw();
  canvas0->Print((get_address("sample0") + ".gif").c_str());
  canvas0->Print((get_address("sample0") + ".pdf").c_str());

  delete rp0;
  delete canvas0;

  RooDataSet* data = pdf->generate( RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->var("bdt")), 50);
  data->Print();
  RooPlot* rp_31 = ws_->var("Mass")->frame();
  RooPlot* rp_32 = ws_->var("MassRes")->frame();
  RooPlot* rp_33 = ws_->var("bdt")->frame();
  TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 1800, 600);
  canvas3->Divide(3,1);
  data->plotOn(rp_31);
  data->plotOn(rp_32);
  data->plotOn(rp_33);
  canvas3->cd(1);
  rp_31->Draw();
  canvas3->cd(2);
  rp_32->Draw();
  canvas3->cd(3);
  rp_33->Draw();
  canvas3->Print((get_address("dataa") + ".gif").c_str());
  canvas3->Print((get_address("dataa") + ".pdf").c_str());
  delete canvas3;
  delete rp_31;
  delete rp_32;
  delete rp_33;

  if (pee) pdf->fitTo(*data, ConditionalObservables(*ws_->var("MassRes")));
  else pdf->fitTo(*data);

  RooPlot* rp = ws_->var("Mass")->frame();
  data->plotOn(rp);
  pdf->paramOn(rp);
  if (pee) {
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("MassRes")), *data, kFALSE), LineColor(kOrange));
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("MassRes")), *data, kFALSE), Components("pdf_bs"), LineColor(kRed));
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("MassRes")), *data, kFALSE), Components("pdf_bd"), LineColor(kBlue));
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("MassRes")), *data, kFALSE), Components("pdf_rare"), LineColor(kGreen));
    pdf->plotOn(rp, ProjWData(RooArgSet(*ws_->var("MassRes")), *data, kFALSE), Components("pdf_comb"), LineColor(kCyan));
  }
  else {
    pdf->plotOn(rp, LineColor(kOrange));
    pdf->plotOn(rp, Components("pdf_bs"), LineColor(kRed));
    pdf->plotOn(rp, Components("pdf_bd"), LineColor(kBlue));
    pdf->plotOn(rp, Components("pdf_rare"), LineColor(kGreen));
    pdf->plotOn(rp, Components("pdf_comb"), LineColor(kCyan));
  }
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  rp->Draw();
  canvas->Print((get_address("sample") + ".gif").c_str());
  canvas->Print((get_address("sample") + ".pdf").c_str());
  delete rp;
  delete canvas;

  if (bdt_fit_) {
    RooPlot* rp_bdt = ws_->var("bdt")->frame();
    data->plotOn(rp_bdt);
    pdf->plotOn(rp_bdt, LineColor(kOrange));
    pdf->plotOn(rp_bdt, Components("pdf_bs"), LineColor(kRed));
    pdf->plotOn(rp_bdt, Components("pdf_bd"), LineColor(kBlue));
    pdf->plotOn(rp_bdt, Components("pdf_rare"), LineColor(kGreen));
    pdf->plotOn(rp_bdt, Components("pdf_comb"), LineColor(kCyan));
    pdf->paramOn(rp_bdt);
    TCanvas* canvas_bdt = new TCanvas("canvas_bdt", "canvas_bdt", 600, 600);
    rp_bdt->Draw();
    canvas_bdt->Print((get_address("BDT_sample") + ".gif").c_str());
    canvas_bdt->Print((get_address("BDT_sample") + ".pdf").c_str());
    delete rp_bdt;
    delete canvas_bdt;
  }

  if (pee) {
    TH1* mass_eta_h = pdf->createHistogram("fit", *ws_->var("Mass"), Binning(50), YVar(*ws_->var("MassRes"), Binning(30))) ;
    TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 600, 600);
    mass_eta_h->Draw("surf");
    canvas1->Print((get_address("sample1_signals") + ".gif").c_str());
    canvas1->Print((get_address("sample1_signals") + ".pdf").c_str());
    delete canvas1;
  }
}

string pdf_analysis::get_address(string name, string pdf, bool channeling) {
  ostringstream address;
  address <<  "fig/" << name;
  if (pdf == "") address << "_" << meth_;
  else address << "_" << pdf << "_" << meth_;
  if (simul_) {
    address  << "_simul_";
    if (channeling) address << channel;
  }
  if (!simul_) address  << "_" << ch_s_;
  if (simul_ && simul_bdt_) {
    address << "_simulBdt_";
    if (channeling) address << channels_bdt;
  }
  if (BF_) address << "_BF";
  if (SM_) address << "_SM";
  if (bd_constr_) address << "_BdConst";
  if (pee) address << "_PEE";
  if (bdt_fit_) address << "_2D";
  return address.str();
}

void pdf_analysis::getBFnumbers(string eff_filename, string numbers_filename) {
  eff_bd_val.resize(channels);
  eff_bs_val.resize(channels);
  eff_bu_val.resize(channels);
  N_bu_val.resize(channels);
  eff_bd_err.resize(channels);
  eff_bs_err.resize(channels);
  eff_bu_err.resize(channels);
  N_bu_err.resize(channels);
  BF_bs_val.resize(channels);
  BF_bd_val.resize(channels);
  BF_bs_err.resize(channels);
  BF_bd_err.resize(channels);

  parse_external_numbers(numbers_filename);
  if (simul_bdt_) { cout << "I don't have efficiencies for the bdt channels..." << endl; return;}
  parse_efficiency_numbers(eff_filename);
  if (channels == 4) parse_efficiency_numbers(eff_filename, 2);
}

void pdf_analysis::parse_efficiency_numbers(string filename, int offset) {

  string file_address = "./input/";
  int channels_;
  if (offset == 0) {
    file_address += "2011/" + filename;
    channels_ = channels;
  }
  else if (offset == 2) {
    file_address += "2012/" + filename;
    channels_ = channels / 2;
  }
  else {
    cout << "not ready " << endl;
    exit(1);
  }
  cout << "parsing " << file_address << endl;
  FILE *file = fopen(file_address.c_str(), "r");
  if (!file) {cout << "file " << file_address << " does not exist"; exit(1);}

  char buffer[2048];
  char left[1024];
  double number;

  vector < pair<string, string> > end_bd(channels_);
  vector < pair<string, string> > end_bs(channels_);
  vector < pair<string, string> > end_bu(channels_);
  vector < pair<string, string> > end_Nbu(channels_);
  int ii = -1;
  for (int i = 0; i < channels_; i++) {
    if (simul_) ii = i;
    else ii = ch_i_;
    end_bd[i] = make_pair(Form("N-EFF-TOT-BDMM%d:val", ii), Form("N-EFF-TOT-BDMM%d:tot", ii));
    end_bs[i] = make_pair(Form("N-EFF-TOT-BSMM%d:val", ii), Form("N-EFF-TOT-BSMM%d:tot", ii));
    end_bu[i] = make_pair(Form("N-EFF-TOT-BPLUS%d:val", ii), Form("N-EFF-TOT-BPLUS%d:tot", ii));
    end_Nbu[i] = make_pair(Form("N-OBS-BPLUS%d:val", ii), Form("N-OBS-BPLUS%d:tot", ii));
  }

  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '\045') continue;
    if (buffer[0] == '\040') continue;
    sscanf(buffer, "%s   {\\ensuremath{{%lf } } }", left, &number);
    string left_s(left);
    for (int i = 0; i < channels_; i++) {
      size_t found;
      found = left_s.find(end_bd[i].first);
      if (found != string::npos) eff_bd_val[i+offset] = number;
      found = left_s.find(end_bd[i].second);
      if (found != string::npos)eff_bd_err[i+offset] = number;
      found = left_s.find(end_bs[i].first);
      if (found != string::npos) eff_bs_val[i+offset] = number;
      found = left_s.find(end_bs[i].second);
      if (found != string::npos) eff_bs_err[i+offset] = number;
      found = left_s.find(end_bu[i].first);
      if (found != string::npos) eff_bu_val[i+offset] = number;
      found = left_s.find(end_bu[i].second);
      if (found != string::npos) eff_bu_err[i+offset] = number;
      found = left_s.find(end_Nbu[i].first);
      if (found != string::npos) N_bu_val[i+offset] = number;
      found = left_s.find(end_Nbu[i].second);
      if (found != string::npos) N_bu_err[i+offset] = number;
    }
  }
  for (int i = 0; i < channels_ && i < 2; i++) {
    if (simul_) ii = i;
    else ii = ch_i_;
    cout << "etacat " << ii + offset << endl;
    cout << "bd eff = " << eff_bd_val[i+offset] << " \\pm " << eff_bd_err[i+offset] << endl;
    cout << "bs eff = " << eff_bs_val[i+offset] << " \\pm " << eff_bs_err[i+offset] << endl;
    cout << "bu eff = " << eff_bu_val[i+offset] << " \\pm " << eff_bu_err[i+offset] << endl;
    cout << "N bu   = " << N_bu_val[i+offset] << " \\pm " << N_bu_err[i+offset] << endl;
    cout << endl;
  }
}

void pdf_analysis::parse_external_numbers(string filename) {
  cout << "parsing " << filename << endl;
  FILE *file = fopen(filename.c_str(), "r");
  if (!file) {cout << "file " << filename << " does not exist"; exit(1);}
  char buffer[1024];
  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '\045') continue;
    sscanf(buffer, "fsfu\t%lf\t%lf", &fs_over_fu_val, &fs_over_fu_err);
    sscanf(buffer, "BuToKpsiK\t%lf\t%lf", &Bu2JpsiK_BF_val, &Bu2JpsiK_BF_err);
    sscanf(buffer, "JpsiToMuMu\t%lf\t%lf", &Jpsi2MuMu_BF_val, &Jpsi2MuMu_BF_err);
    sscanf(buffer, "Bs2MuMu_BF\t%lf\t%lf", &Bs2MuMu_SM_BF_val, &Bs2MuMu_SM_BF_err);
  }
  cout << "Jpsi2MuMu_SM_BF " <<  Bs2MuMu_SM_BF_val << " \\pm " << Bs2MuMu_SM_BF_err << endl;
  cout << "fs/fu        " <<  fs_over_fu_val << " \\pm " << fs_over_fu_err << endl;
  cout << "Bu2JpsiK_BF  " <<  Bu2JpsiK_BF_val << " \\pm " << Bu2JpsiK_BF_err << endl;
  cout << "Jpsi2MuMu_BF " <<  Jpsi2MuMu_BF_val << " \\pm " << Jpsi2MuMu_BF_err << endl;
}
