#include "pdf_analysis.h"

pdf_analysis::pdf_analysis(bool print, string meth, string ch_s, string range, bool SM, bool bd_constr) {
  cout << "analysis constructor" << endl;
  print_ = print;
  meth_ = meth;
  ch_s_ = ch_s;
  range_ = range;
  SM_ = SM;
  bd_constr_ = bd_constr;
  channels = 1;
  verbosity = 1;
  old_tree = false;
  pee = false;
  simul_ = false;
  no_legend = false;
}

void pdf_analysis::initialize () {
  cout << "inizialization" << endl;
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
}

void pdf_analysis::define_pdfs () {

  for (int i = 0; i < channels; i++) {
    define_bs(i);
    define_bd(i);
    define_peaking(i);
    define_nonpeaking(i);
    define_comb(i);

    define_signals(i);
    define_rare(i);
    define_bkg_fractional(i);
    define_bkg_extended(i);
    define_signalsrare(i);
    //define_bscomb();

    define_total_fractional(i);
    define_total_extended(i);
  }
}

void pdf_analysis::fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error, bool hesse) {
  rds_ = data;
  pdf_name = "pdf_" + pdf;
  if (extended) pdf_name = "pdf_ext_" + pdf;
  string rdh_name = rds_->GetName();
  cout << "fitting " << rdh_name << " in range " << range_ << " with " << pdf_name << ":" << endl;
  ws_->pdf( pdf_name.c_str())->Print();
  
  if (!pee) RFR = ws_->pdf( pdf_name.c_str())->fitTo(*rds_, Extended(extended), SumW2Error(sumw2error), Range(range_.c_str())/*, SumCoefRange(range_.c_str())*/, NumCPU(2), Hesse(hesse), Save());
  else RFR = ws_->pdf( pdf_name.c_str())->fitTo(*rds_, Extended(extended), SumW2Error(sumw2error), Range(range_.c_str())/*, SumCoefRange(range_.c_str())*/, NumCPU(2), Hesse(hesse), Save(), ConditionalObservables(*ws_->var("eta")));

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

void pdf_analysis::define_bs(int i = 0) {
  RooRealVar N_bs("N_bs", "N_bs", 0, 100);
  ws_->import(N_bs);


  ws_->factory("Mean_bs[5.35, 5.32, 5.4]");
  ws_->factory("Sigma_bs[0.02, 0.005, 0.2]");
  ws_->factory("Alpha_bs[2.8, 0.1, 3.0]");
  ws_->factory("Enne_bs[1., 0., 10.]");
  if (!pee) {
    ws_->factory("Sigma2_bs[0.04, 0.005, 0.2]");
    ws_->factory("CoeffGauss_bs[0.5,0.,1.]");
    ws_->factory("Gaussian::Gau_bs(Mass, Mean_bs, Sigma_bs)");
    ws_->factory("CBShape::CB_bs(Mass, Mean_bs, Sigma2_bs, Alpha_bs, Enne_bs)");

    ws_->factory("SUM::pdf_bs(CoeffGauss_bs*Gau_bs, CB_bs)");
  }
  else {
    RooFormulaVar *SigmaRes = new RooFormulaVar("SigmaRes", "0.0078*@0*@0+0.035", RooArgList(*ws_->var("eta")));
    ws_->import(*SigmaRes);
    ws_->factory("CBShape::CB_bs(Mass, Mean_bs, SigmaRes, Alpha_bs, Enne_bs)");

    if (!simul_) {
      RooProdPdf* pdf_bs = new RooProdPdf("pdf_bs", "pdf_bs", *ws_->pdf("etapdf_bs"), Conditional(*ws_->pdf("CB_bs"), *ws_->var("Mass")));
      ws_->import(*pdf_bs);
    }
    else {
      RooDataHist temp_rdh(getRandom_rdh());
      RooHistPdf* etapdf_bs = new RooHistPdf("etapdf_bs", "etapdf_bs", *ws_->var("eta"), temp_rdh);
      RooProdPdf* pdf_bs = new RooProdPdf("pdf_bs", "pdf_bs", *etapdf_bs, Conditional(*ws_->pdf("CB_bs"), *ws_->var("Mass")));
      ws_->import(*pdf_bs);
    }
  }
}

void pdf_analysis::define_bd(int i = 0) {

  ws_->factory("Mean_bd[5.25, 5.20, 5.29]");
  ws_->factory("Sigma_bd[0.02, 0.005, 0.2]");
  ws_->factory("Alpha_bd[2.8, 0.1, 3.0]");
  ws_->factory("Enne_bd[1., 0., 1.5]");
  if (!pee) {
    ws_->factory("Sigma2_bd[0.04, 0.005, 0.2]");

    ws_->factory("CoeffGauss_bd[0.5,0.,1.]");
    ws_->factory("Gaussian::Gau_bd(Mass, Mean_bd, Sigma_bd)");
    ws_->factory("CBShape::CB_bd(Mass, Mean_bd, Sigma2_bd, Alpha_bd, Enne_bd)");

    ws_->factory("SUM::pdf_bd(CoeffGauss_bd*Gau_bd, CB_bd)");
  }
  else {
    ws_->factory("CBShape::CB_bd(Mass, Mean_bd, SigmaRes, Alpha_bd, Enne_bd)");
    if (!simul_) {
      RooProdPdf* pdf_bd = new RooProdPdf("pdf_bd", "pdf_bd", *ws_->pdf("etapdf_bd"), Conditional(*ws_->pdf("CB_bd"), *ws_->var("Mass")));
      ws_->import(*pdf_bd);
    }
    else {
      RooDataHist temp_rdh(getRandom_rdh());
      RooHistPdf* etapdf_bd = new RooHistPdf("etapdf_bd", "etapdf_bd", *ws_->var("eta"), temp_rdh);
      RooProdPdf* pdf_bd = new RooProdPdf("pdf_bd", "pdf_bd", *etapdf_bd, Conditional(*ws_->pdf("CB_bd"), *ws_->var("Mass")));
      ws_->import(*pdf_bd);
    }
  }
  if (SM_) {
    RooConstVar *SM_Bd_over_Bs = new RooConstVar("SM_Bd_over_Bs", "SM_Bd_over_Bs", ratio_);
    RooFormulaVar *N_bd_constr = new RooFormulaVar("N_bd_constr", "@0*@1", RooArgList(*ws_->var("N_bs"), *SM_Bd_over_Bs));
    ws_->import(*N_bd_constr);
  }
  else if (bd_constr_) {
    RooRealVar * Bd_over_Bs = new RooRealVar("Bd_over_Bs", "Bd_over_Bs", 0., 10., "");
    RooFormulaVar *N_bd_constr = new RooFormulaVar("N_bd_constr", "@0*@1", RooArgList(*ws_->var("N_bs"), *Bd_over_Bs));
    ws_->import(*N_bd_constr);
  }
  else ws_->factory("N_bd[0, 100]");
}

void pdf_analysis::define_peaking(int i = 0) {
  
  ws_->factory("N_peaking[0, 100]");

  if (old_tree) {
    ws_->factory("Mean_peaking[5.25, 5.20, 5.3]");
    ws_->factory("Sigma_peaking[0.050, 0.01, 0.20]");
  
    ws_->factory("Gaussian::pdf_peaking(Mass,Mean_peaking,Sigma_peaking)");
  }
  else {
    ws_->factory("Mean_peaking[5.1, 4.9, 5.4]");
    ws_->factory("Sigma_peaking[0.02, 0.005, 0.2]");
    ws_->factory("Sigma2_peaking[0.04, 0.005, 0.2]");
    ws_->factory("Alpha_peaking[2.8, 0., 100.0]");
    ws_->factory("Enne_peaking[1., 0., 100.0]");
    ws_->factory("CoeffGauss_peaking[0.5,0.,1.]");
    ws_->factory("Gaussian::Gau_peaking(Mass, Mean_peaking, Sigma_peaking)");
    ws_->factory("CBShape::CB_peaking(Mass, Mean_peaking, Sigma2_peaking, Alpha_peaking, Enne_peaking)");

    if (!pee) ws_->factory("SUM::pdf_peaking(CoeffGauss_peaking*Gau_peaking, CB_peaking)");
    else {
      ws_->factory("SUM::mass_peaking(CoeffGauss_peaking*Gau_peaking, CB_peaking)");
      if (!simul_) {
        RooProdPdf* pdf_peaking = new RooProdPdf("pdf_peaking", "pdf_peaking", *ws_->pdf("etapdf_peaking"), Conditional(*ws_->pdf("mass_peaking"), *ws_->var("Mass")));
        ws_->import(*pdf_peaking);
      }
      else {
        RooDataHist temp_rdh(getRandom_rdh());
        RooHistPdf* etapdf_peaking = new RooHistPdf("etapdf_peaking", "etapdf_peaking", *ws_->var("eta"), temp_rdh);
        RooProdPdf* pdf_peaking = new RooProdPdf("pdf_peaking", "pdf_peaking", *etapdf_peaking, Conditional(*ws_->pdf("mass_peaking"), *ws_->var("Mass")));
        ws_->import(*pdf_peaking);
      }
    }
  }
}

void pdf_analysis::define_nonpeaking(int i = 0) {

  if (old_tree) {
    ws_->factory("N_nonpeaking[0, 100]");
    ws_->factory("m0_nonpeaking[5., 6.]");
    ws_->factory("c_nonpeaking[1., 0.1, 20]");
    ws_->factory("p_nonpeaking[0.5, 0.1, 5.]");

    ws_->factory("ArgusBG::pdf_nonpeaking(Mass,m0_nonpeaking,c_nonpeaking,p_nonpeaking)");
  }
  else {
    RooRealVar * C0 = new RooRealVar("C0", "C0", 0.1, 0.001, 10., "");
    RooRealVar * C1 = new RooRealVar("C1", "C1", 0.1, 0.001, 10., "");
    RooRealVar * C2 = new RooRealVar("C2", "C2", 0., -2., 2., "");
    RooRealVar * C3 = new RooRealVar("C3", "C3", 0., -2., 2., "");
    RooArgList lista(*C0, *C1, *C2, *C3);
    //RooPolynomial *poly = new RooPolynomial("poly","poly", *ws_->var("Mass"), lista);
    RooChebychev *poly = new RooChebychev("poly", "poly", *ws_->var("Mass"), lista);
    ws_->import(*poly);
    if (!pee) {
      ws_->factory("PROD::pdf_nonpeaking(Exponential::semiexpo(Mass,Alpha_semi[-5,-20.,-0.01]),poly)");
      //ws_->factory("Exponential::pdf_nonpeaking(Mass,Alpha1[-1,-20.,-0.01])");
      //ws_->factory("SUM::pdf_nonpeaking(expo1frac_rare[0.,1.]*Exponential::expo1(Mass,Alpha1[-10.,10.]),Exponential::expo2(Mass,Alpha2[-10.,10.]))");
    }
    else {
      ws_->factory("PROD::mass_nonpeaking(Exponential::semiexpo(Mass,Alpha_semi[-5,-20.,-0.01]),poly)");
      if (!simul_) {
        RooProdPdf* pdf_nonpeaking = new RooProdPdf("pdf_nonpeaking", "pdf_nonpeaking", *ws_->pdf("etapdf_nonpeaking"), Conditional(*ws_->pdf("mass_nonpeaking"), *ws_->var("Mass")));
        ws_->import(*pdf_nonpeaking);
      }
      else {
        RooDataHist temp_rdh(getRandom_rdh());
        RooHistPdf* etapdf_nonpeaking = new RooHistPdf("etapdf_nonpeaking", "etapdf_nonpeaking", *ws_->var("eta"), temp_rdh);
        RooProdPdf* pdf_nonpeaking = new RooProdPdf("pdf_nonpeaking", "pdf_nonpeaking", *etapdf_nonpeaking, Conditional(*ws_->pdf("mass_nonpeaking"), *ws_->var("Mass")));
        ws_->import(*pdf_nonpeaking);
      }
    }
  }
}

void pdf_analysis::define_comb(int i = 0) {
  
  ws_->factory("N_comb[0, 100]");
  if (!pee) ws_->factory("Uniform::pdf_comb(Mass)");
  else {
    ws_->factory("Uniform::mass_comb(Mass)");
    if (!simul_) {
      RooProdPdf* pdf_comb = new RooProdPdf("pdf_comb", "pdf_comb", *ws_->pdf("etapdf_comb"), Conditional(*ws_->pdf("mass_comb"), *ws_->var("Mass")));
      ws_->import(*pdf_comb);
    }
    else {
      RooDataHist temp_rdh(getRandom_rdh());
      RooHistPdf* etapdf_comb = new RooHistPdf("etapdf_comb", "etapdf_comb", *ws_->var("eta"), temp_rdh);
      RooProdPdf* pdf_comb = new RooProdPdf("pdf_comb", "pdf_comb", *etapdf_comb, Conditional(*ws_->pdf("mass_comb"), *ws_->var("Mass")));
      ws_->import(*pdf_comb);
    }
  }
}

void pdf_analysis::define_signals(int i = 0) {
  
  ws_->factory("N_signals[0, 100]");
  ws_->factory("bsfrac_signals[0.5, 0.0, 1.0]");

  ws_->factory("SUM::pdf_signals(bsfrac_signals*pdf_bs, pdf_bd)");
  
  if (!SM_ && !bd_constr_) ws_->factory("SUM::pdf_ext_signals(N_bs*pdf_bs, N_bd*pdf_bd)");
}

void pdf_analysis::define_rare(int i = 0) {
  
  ws_->factory("N_rare[0, 100]");
  ws_->factory("peakingfrac_rare[0.5, 0.0, 1.0]");
  ws_->factory("SUM::pdf_rare(peakingfrac_rare*pdf_peaking, pdf_nonpeaking)");
}

void pdf_analysis::define_rare2(RooDataHist* data, int i = 0) {
  ws_->factory("N_hist[0, 100]");
  RooHistPdf* pdf_rare = new RooHistPdf("pdf_hist", "pdf_hist", *ws_->var("Mass"), *data, 4);
  ws_->import(*pdf_rare);
}

void pdf_analysis::define_rare3(int i = 0) {
  ws_->factory("N_expo[0, 100]");
  ws_->factory("Exponential::pdf_expo(Mass,Alpha1[-10.,10.])");
  //ws_->factory("SUM::pdf_expo(expo1frac_rare[0.,1.]*Exponential::expo1(Mass,Alpha1[-10.,10.]),Exponential::expo2(Mass,Alpha2[-10.,10.]))");
}

void pdf_analysis::define_signalsrare(int i = 0) {
  
  ws_->factory("signalfrac_signalsrare[0.5, 0.0, 1.0]");
  ws_->factory("SUM::pdf_signalsrare(signalfrac_signalsrare*pdf_signals, pdf_rare)");

  if (!SM_ && !bd_constr_) ws_->factory("SUM::pdf_ext_signalsrare(N_bs*pdf_bs, N_bd*pdf_bd, N_rare*pdf_rare)");
}

void pdf_analysis::define_bkg_fractional(int i = 0) {
  
  ws_->factory("N_bkg[0, 100]");
  ws_->factory("rarefrac_bkg[0.5, 0.0, 1.0]");
  
  ws_->factory("SUM::pdf_bkg(rarefrac_bkg*pdf_rare, pdf_comb)");
  
  RooExtendPdf* bkgExt = new RooExtendPdf("pdf_fract_ext_bkg", "bkgExt", *ws_->pdf("pdf_bkg"), *ws_->var("N_bkg"), /* range_.c_str()*/ "sb_lo,sb_hi"); /// WARNING
  ws_->import(*bkgExt);
}

void pdf_analysis::define_bkg_extended(int i = 0) {

//  RooExtendPdf rare_ext("pdf_ext_rare", "rare_ext", *ws_->pdf("pdf_rare"), *ws_->var("N_rare"), "sb_lo,sb_hi");
//  RooExtendPdf comb_ext("pdf_ext_comb", "comb_ext", *ws_->pdf("pdf_comb"), *ws_->var("N_comb"), "sb_lo,sb_hi");
//  RooAddPdf bkg_ext("pdf_ext_bkg", "pdf_ext_bkg", RooArgList(rare_ext, comb_ext));
//  ws_->import(bkg_ext);
  ws_->factory("SUM::pdf_ext_bkg(N_rare*pdf_rare, N_comb*pdf_comb)");

}


void pdf_analysis::define_total_fractional(int i = 0) {
  
  ws_->factory("N_all[0, 100]");
  ws_->factory("signalsfrac_total[0.5, 0.0, 1.0]");
  
  ws_->factory("SUM::pdf_frac_total(signalsfrac_total*pdf_signals, pdf_bkg)");
  
  RooExtendPdf* allExt = new RooExtendPdf("pdf_frac_ext_total", "allExt", *ws_->pdf("pdf_frac_total"), *ws_->var("N_all"), range_.c_str());
  ws_->import(*allExt);
}

void pdf_analysis::define_total_extended(int i = 0) {
  
  if (SM_ || bd_constr_) ws_->factory("SUM::pdf_ext_total(N_bs*pdf_bs, N_bd_constr*pdf_bd, N_rare*pdf_rare, N_comb*pdf_comb)");
  else ws_->factory("SUM::pdf_ext_total(N_bs*pdf_bs, N_bd*pdf_bd, N_rare*pdf_rare, N_comb*pdf_comb)");
  
  return; 
}

string pdf_analysis::define_pdf_sum(string name, int i = 0) {

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
  int colors[11] = {632, 400, 616, 432, 800, 416, 820, 840, 860, 880, 900};
  RooPlot *rp = ws_->var("Mass")->frame();
  rds_->plotOn(rp, Binning(20));
  if (!pee) ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue), Range(range_.c_str())/*, ProjectionRange("eta_all"), Normalization((rds_->sumEntries(), RooAbsReal::NumEvent))*/);
  else {
    ws_->pdf(pdf_name.c_str())->Print();
    TH1* mass_eta_h;
    /*if (simul_) mass_eta_h = ws_->pdf(pdf_name.c_str())->createHistogram("fit", *ws_->var("Mass"), Binning(50), YVar(*ws_->var(Form("eta_channel_%d", channel)), Binning(50))) ;
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

void pdf_analysis::simsplit() {
  cout << "simsplitting" << endl;
  simul_ = true;
  ostringstream splitter;
  splitter << "SIMCLONE::pdf_ext_simul(pdf_ext_total, $SplitParam({";
  RooArgSet* allvars = ws_->pdf("pdf_ext_total")->getVariables();
  TIterator* it = allvars->createIterator();
  RooRealVar* var_Obj = 0;
  int cont = 0;
  while ( (var_Obj = (RooRealVar*)it->Next()) ) {
    string name = var_Obj->GetName();
    if ( !(name == "Mass") /*&& !(name == "eta")*/ && !(name == "Bd_over_Bs") && !(name == "SM_Bd_over_Bs")) {
      if (cont != 0) splitter << ", ";
      splitter << name;
      cont++;
    }
  }
  /// splitter << ",etapdf_bs"; /// warning
  splitter << "}, channel[";
  for (int i = 0; i < channels; i++) {
    if (i != 0) splitter << ",";
    splitter << "channel_" << i ;
  }
  splitter << "]))";
  cout << splitter.str() << endl;
  ws_->factory(splitter.str().c_str());

  ws_->Print();
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
  if (value == var->getMin() || (value <= 0.0001 && value >= -0.0001)) return 0;
  double error = var->getErrorLo();
  if (error == 0) {
    if (value - var->getError() < 0.0) return (-1 * value);
    else return (-1*var->getError());
  }
  else return error;
}

RooHistPdf* pdf_analysis::define_etapdf(RooDataSet *rds, string name) {
  RooRealVar* weight = new RooRealVar("weight", "weight", 0., 100000., "");
  RooDataSet* eta_rds = new RooDataSet("eta_rds", "eta_rds", RooArgSet(*ws_->var("eta"), *weight), Import(*rds), WeightVar("weight"));
  string name_ = "etardh_" + name;
  RooDataHist *eta_rdh = eta_rds->binnedClone(name_.c_str());
  name_ = "etapdf_" + name;
  RooHistPdf * eta_rhpdf = new RooHistPdf(name_.c_str(), name_.c_str(), *ws_->var("eta"), *eta_rdh);
  ws_->import(*eta_rhpdf);
  return eta_rhpdf;
}

RooDataHist pdf_analysis::getRandom_rdh() {
  TH1D* temp_h = new TH1D("temp_h", "temp_h", 100, -2.4, 2.4);
  temp_h->FillRandom("gaus", 1000);
  RooDataHist temp_rdh("temp_rdh", "temp_rdh", *ws_->var("eta"), temp_h);
  return temp_rdh;
}
