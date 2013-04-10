#include "pdf_fitData.h"

pdf_fitData::pdf_fitData(bool print, string input_estimates, string range, int BF, bool SM, bool bd_constr, int simul, int simulbdt, int simulall, bool pee_, bool bdt_fit, string ch_s, int sig, bool asimov, bool syste, bool randomsyste, bool rare_constr, int nexp, bool bd, string years): pdf_analysis(print, ch_s, range, BF, SM, bd_constr, simul, simulbdt, simulall, pee_, bdt_fit) {
  cout << "fitData constructor" << endl;
  input_estimates_ = input_estimates;
  estimate_bs.resize(channels);
  estimate_bd.resize(channels);
  estimate_semi.resize(channels);
  estimate_comb.resize(channels);

  systematics_bs.resize(channels);
  systematics_bd.resize(channels);
  systematics_semi.resize(channels);
  systematics_comb.resize(channels);

  estimate2D_channel.resize(channels, vector<double> (channels_bdt));
  estimate2D_bs.resize(channels, vector<double> (channels_bdt));
  estimate2D_bd.resize(channels, vector<double> (channels_bdt));
  estimate2D_semi.resize(channels, vector<double> (channels_bdt));
  estimate2D_comb.resize(channels, vector<double> (channels_bdt));
  systematics2D_channel.resize(channels, vector<double> (channels_bdt));
  systematics2D_bs.resize(channels, vector<double> (channels_bdt));
  systematics2D_bd.resize(channels, vector<double> (channels_bdt));
  systematics2D_semi.resize(channels, vector<double> (channels_bdt));
  systematics2D_comb.resize(channels, vector<double> (channels_bdt));

  lumi = -1;
  parse_estimate();
  random = false;

  RooRandom::randomGenerator()->SetSeed(0);

  sign = sig;

  syst = syste;
  randomsyst = randomsyste;
  rare_constr_ = rare_constr;

  /// the pdf
  pdfname = "pdf_ext_total";
  if (simul_ && !syst) pdfname = "pdf_ext_simul_noconstr";
  if (simul_ && syst) pdfname = "pdf_ext_simul";
  cout << red_color_bold << ">>>>>>>>>>>>>>> the name of the fitting pdf is " << pdfname << " <<<<<<<<<<<<<<<<<<<" << default_console_color << endl;

  NExp = nexp;
  Bd = bd;
  SMIsNull = false;
  years_ = years;

  asimov_ = asimov;
}

pdf_fitData::~pdf_fitData() {
  cout << "pdf_fitData destructor" << endl;
}

void pdf_fitData::parse_estimate() {
  char buffer[1024];
  char cutName[128];
  double cut;
  FILE *estimate_file = fopen(input_estimates_.c_str(), "r");
  cout << "event estimates in " << input_estimates_ << " :" << endl;
  while (fgets(buffer, sizeof(buffer), estimate_file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '#') continue;
    sscanf(buffer, "%s %lf", cutName, &cut);
    if (!parse(cutName, cut)) {
      cout << "==> Error parsing variable " << cutName << endl;
      exit(EXIT_FAILURE);
    }
  }
  if (estimate_file) fclose(estimate_file);
  return;
}

bool pdf_fitData::parse(char *cutName, float cut) {
  if (lumi == -1) {
    char test_cut[128];
    sprintf(test_cut, "lumi");
    if (!strcmp(cutName, test_cut)) {
      lumi = (double)cut;
      cout << "lumi is: " << lumi << " times the true luminosity" << endl;
      return true;
    }
  }
  if (simul_ && !simul_bdt_ && !simul_all_) {
    for (int i = 0; i < channels; i++) {
      char test_cut[128];
      sprintf(test_cut, "bs_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_bs[i] = (int)cut;
        cout << "bs[" << i <<"]: " << estimate_bs[i] << endl;
        return true;
      }
      sprintf(test_cut, "bd_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_bd[i] = (int)cut;
        cout << "bd[" << i <<"]: " << estimate_bd[i] << endl;
        return true;
      }
      sprintf(test_cut, "semi_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_semi[i] = (int)cut;
        cout << "semi[" << i <<"]: " << estimate_semi[i] << endl;
        return true;
      }
      sprintf(test_cut, "comb_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_comb[i] = (int)cut;
        cout << "comb[" << i <<"]: " << estimate_comb[i] << endl;
        return true;
      }
    }
  }
  else if (!simul_) {
    int i = atoi(ch_s_.c_str());
    char test_cut[128];
    sprintf(test_cut, "bs_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_bs[0] = (int)cut;
      cout << "bs[" << i <<"]: " << estimate_bs[0] << endl;
      return true;
    }
    sprintf(test_cut, "bd_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_bd[0] = (int)cut;
      cout << "bd[" << i <<"]: " << estimate_bd[0] << endl;
      return true;
    }
    sprintf(test_cut, "semi_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_semi[0] = (int)cut;
      cout << "semi[" << i <<"]: " << estimate_semi[0] << endl;
      return true;
    }
    sprintf(test_cut, "comb_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_comb[0] = (int)cut;
      cout << "comb[" << i <<"]: " << estimate_comb[0] << endl;
      return true;
    }
    return true;
  }
  else {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < bdt_index_max(i); j++) {
        char test_cut[128];
        sprintf(test_cut, "bs_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_bs[i][j] = (int)cut;
          cout << "bs[" << i <<"][" << j << "]: " << estimate2D_bs[i][j] << endl;
          return true;
        }
        sprintf(test_cut, "bd_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_bd[i][j] = (int)cut;
          cout << "bd[" << i <<"][" << j << "]: " << estimate2D_bd[i][j] << endl;
          return true;
        }
        sprintf(test_cut, "semi_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_semi[i][j] = (int)cut;
          cout << "semi[" << i <<"][" << j << "]: " << estimate2D_semi[i][j] << endl;
          return true;
        }
        sprintf(test_cut, "comb_%d_%d", i, j);
        if (!strcmp(cutName, test_cut)) {
          estimate2D_comb[i][j] = (int)cut;
          cout << "comb[" << i <<"][" << j << "]: " << estimate2D_comb[i][j] << endl;
          return true;
        }
      }
    }
  }
  return true;
}

RooFitResult* pdf_fitData::fit_pdf(bool do_not_import, string pdf_name) {
	if (pdf_name == "") pdf_name = pdfname;
  if (!simul_) {
  	pdf_name = "pdf_ext_total";
    RooAbsData* subdata = global_data->reduce(Form("etacat==etacat::etacat_%d", channel));
    global_data = (RooDataSet*)subdata;
  }
  cout << red_color_bold << ">>>>>>>>>>>>>>>>> fitting " << global_data->GetName() << " in range " << range_ << " with " << pdf_name << default_console_color << endl;
  global_data->Print();
  ws_->pdf(pdf_name.c_str())->Print();
  RFR = ws_->pdf(pdf_name.c_str())->fitTo(*global_data, Extended(), Save(1), Minos(asimov_ ? false : true), pee ? ConditionalObservables(*ws_->var("MassRes")) : RooCmdArg::none()/*, syst ? Constrain(*ws_->set("constr")) : RooCmdArg::none()*/);
  if (!do_not_import) ws_->import(*global_data);
  if (verbosity > 0) RFR->Print();

//  RooArgSet projw_set(*ws_->cat("etacat"), *ws_->var("MassRes"));
//  RooPlot* final_p = ws_->var("Mass")->frame(Bins(25));
//  global_data->plotOn(final_p, Cut("etacat==etacat::etacat_2"));
//  ws_->pdf(pdf_name.c_str())->plotOn(final_p, Slice(*ws_->cat("etacat"), "etacat_2"), ProjWData(projw_set, *global_data, false));
//  TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",600,600) ;
//  final_p->Draw() ;
//  c->Print("fig/fit.gif");
//  delete final_p;
//  delete c;

//  RooPlot* final_p2 = ws_->var("Mass")->frame(Bins(25));
//  global_data->plotOn(final_p2, Cut("etacat==etacat::etacat_3"));
//  ws_->pdf(pdf_name.c_str())->plotOn(final_p2, Slice(*ws_->cat("etacat"), "etacat_3"), ProjWData(projw_set, *global_data, false));
//  TCanvas* c2 = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",600,600) ;
//  final_p2->Draw() ;
//  c2->Print("fig/fit2.gif");
//  delete final_p2;
//  delete c2;
//  exit(0);

  return RFR;
}

void pdf_fitData::print() {
  cout << "printing" << endl;
  RooPlot *rp = ws_->var("Mass")->frame();
  global_data->plotOn(rp, Binning(40));

  if (!pee) {
    ws_->pdf("pdf_ext_total")->plotOn(rp, FillColor(kYellow), Range(range_.c_str()), LineWidth(3), VisualizeError(*RFR), MoveToBack());
    ws_->pdf("pdf_ext_total")->plotOn(rp, LineColor(kBlue), Range(range_.c_str()), LineWidth(3));
  }
  else {
//    ws_->pdf("pdf_ext_total")->plotOn(rp, FillColor(kYellow), Range(range_.c_str()), LineWidth(3), VisualizeError(*RFR), MoveToBack(), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE));
    ws_->pdf("pdf_ext_total")->plotOn(rp, LineColor(kBlue), Range(range_.c_str()), LineWidth(3), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE));
  }
  ws_->pdf("pdf_ext_total")->paramOn(rp, Layout(0.50, 0.9, 0.9));
  // components
  RooArgSet * set = ws_->pdf("pdf_ext_total")->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (!pee) {
        size_t found;
        found = name.find("pdf_bs");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kRed),          LineStyle(1), DrawOption("F"), FillColor(kRed), FillStyle(3001), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_bd");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kViolet - 4),   LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_comb");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kBlue - 5),   LineStyle(2), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_semi");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kGreen - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_peak");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kCyan - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
      }
      else {
        size_t found;
        found = name.find("pdf_bs");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kRed),          LineStyle(1), DrawOption("F"), FillColor(kRed), FillStyle(3001), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_bd");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kViolet - 4),   LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_comb");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kBlue - 5),   LineStyle(2), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_semi");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kGreen - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_peak");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kCyan - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
      }
    }
  }
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  c->Print((get_address("data_fitData_", "pdf_ext_total") + ".gif").c_str());
  c->Print((get_address("data_fitData_", "pdf_ext_total") + ".pdf").c_str());
  delete rp;
  delete c;
}

void pdf_fitData::print_each_channel(string var, string output, RooWorkspace* ws, RooDataSet* rds_) {
  if (ws == 0) ws = (RooWorkspace*)ws_->Clone();
  if (rds_ == 0) rds_ = (RooDataSet*)global_data->Clone();
  cout << red_color_bold << "printing" << default_console_color << endl;

  for (int i = 0; i < channels; i++) {
  	if (years_ == "0" && i > 1) break;
  	if (years_ == "1" && i < 2) continue;
    for (int j = 0; j < bdt_index_max(i); j++) {
      /// all texts
      string cut;
      string title;
      RooArgSet slice_set;
      RooArgSet projw_set;

      ws->cat("etacat")->setIndex(i);
      ws->cat("bdtcat")->setIndex(j);
      ws->cat("allcat")->setIndex(super_index(i, j));

      if (!simul_bdt_ && !simul_all_) {
        slice_set.add(*ws->cat("etacat"));
        projw_set.add(*ws->cat("etacat"));
        cut = Form("etacat==etacat::etacat_%d", i);
        title = Form("Candidate invariant mass for etacat %d", i);
      }
      else if (simul_bdt_ && !simul_all_) {
        slice_set.add(*ws->cat("etacat"));
        projw_set.add(*ws->cat("etacat"));
        cut = Form("etacat==etacat::etacat_%d&&bdtcat==bdtcat::bdtcat_%d", i, j);
        title = Form("Candidate invariant mass for etacat %d and bdtcat %d", i, j);
        slice_set.add(*ws->cat("bdtcat"));
        projw_set.add(*ws->cat("bdtcat"));
      }
      else if (!simul_bdt_ && simul_all_) {
        slice_set.add(*ws->cat("allcat"));
        projw_set.add(*ws->cat("allcat"));
        cut = Form("allcat==allcat::allcat_%d", super_index(i, j));
        title = Form("Candidate invariant mass for etacat %d and bdtcat %d", i, j);
      }
      if (pee) projw_set.add(*ws->var("MassRes"));

      RooPlot* final_p = ws->var(var.c_str())->frame((var == "bdt") ? Bins(20) : Bins(25), Title(title.c_str()), (var == "bdt") ? Range(0.1, 0.5) : RooCmdArg::none());

      rds_->plotOn(final_p, Cut(cut.c_str()));

      ws->pdf(pdfname.c_str())->plotOn(final_p, Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_), LineColor(kBlue), LineWidth(3));

      if (BF_ > 0) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_semi", i, j)), DrawOption("F"), FillColor(kGreen - 3), FillStyle(3001), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
      ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_bs", i, j)), DrawOption("F"), FillColor(kRed),        FillStyle(3001), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
      ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_bd", i, j)), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3001), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
      if (BF_ > 0) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_peak", i, j)), DrawOption("F"), FillColor(kBlack), FillStyle(3001), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));

      ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_comb", i, j)), LineColor(kBlue - 1),   LineStyle(2), LineWidth(3), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
      if (output == "") {
        ws->pdf(pdfname.c_str())->plotOn(final_p, VisualizeError(*RFR, 1, 1), FillColor(kYellow), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_), MoveToBack());
        if (BF_ > 0) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_semi", i, j)), LineColor(kBlack),  LineStyle(1), DrawOption("L"), LineWidth(2), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
        ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_bs", i, j)), LineColor(kBlack),        LineStyle(1), DrawOption("L"), LineWidth(2), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
        ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_bd", i, j)), LineColor(kBlack), LineStyle(1), DrawOption("L"), LineWidth(2), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
        if (BF_ > 0) ws->pdf(pdfname.c_str())->plotOn(final_p, Components(name("pdf_peak", i, j)), LineColor(kBlack),  LineStyle(1), DrawOption("L"), LineWidth(2), Slice(slice_set), ProjWData(projw_set, *rds_, bdt_fit_));
      }
      rds_->plotOn(final_p, Cut(cut.c_str()));
      final_p->SetMinimum(0);

      TCanvas* final_c = new TCanvas("final_c", "final_c", 600, 600);
      final_p->Draw();

      if (output == "") {
        // legend
        RooArgSet* vars =  ws->pdf(pdfname.c_str())->getVariables();
        RooRealVar* N_bs;
        if (BF_ == 0) N_bs = (RooRealVar*)vars->find(pdf_analysis::name("N_bs", i, j));
        else if (BF_ < 3) N_bs = (RooRealVar*)vars->find("BF_bs");
        else if (BF_ == 3) N_bs = (RooRealVar*)vars->find("BF_bsbd");
        RooRealVar* N_bd;
        if (!bd_constr_ && !SM_ && BF_ < 2) {
          N_bd = (RooRealVar*)vars->find(pdf_analysis::name("N_bd", i, j));
        }
        else if (bd_constr_) N_bd = (RooRealVar*)vars->find("Bd_over_Bs");
        else if (BF_ > 1) N_bd = (RooRealVar*)vars->find("BF_bd");
        RooRealVar* N_comb = (RooRealVar*)vars->find(pdf_analysis::name("N_comb", i, j));
        vector <string> fitresult_tex_vec;
        if (BF_ == 0) {
          ostringstream fitresult_tex;
          fitresult_tex << setprecision(2) << fixed << "N(B_{s}) = " << N_bs->getVal() << " ^{+" << getErrorHigh(N_bs) << "}_{" << getErrorLow(N_bs) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
        }
        else {
          ostringstream fitresult_tex;
          string title("BF(B^{0}_{s})");
          string name("BF_bs");
          if (BF_ == 3) {
            title = "BF(B^{0}_{s})/BF(B^{0})";
            name = "BF_bsbd";
          }
          fitresult_tex << setprecision(2) << scientific << title << " = " << N_bs->getVal() << " ^{+" << getErrorHigh(N_bs) << "}_{" << getErrorLow(N_bs) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
          ostringstream fitresult_tex2;
          fitresult_tex2 << "(N(B^{0}_{s}) = " << setprecision(2) << fixed << ws->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
          Double_t BF_bs_val = ws->var(name.c_str())->getVal();
          Double_t N_bs_ = ws->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
          ws->var(name.c_str())->setVal(BF_bs_val + getErrorHigh(N_bs));
          Double_t N_bs_up = ws->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
          Double_t N_bs_error_up = N_bs_up - N_bs_;
          fitresult_tex2 << " ^{+" << N_bs_error_up;
          ws->var(name.c_str())->setVal(BF_bs_val);
          ws->var(name.c_str())->setVal(BF_bs_val + getErrorLow(N_bs));
          Double_t N_bs_down = ws->function(pdf_analysis::name("N_bs_formula", i, j))->getVal();
          Double_t N_bs_error_down = N_bs_ - N_bs_down;
          fitresult_tex2 << "}_{-" << N_bs_error_down << "}" << ")";
          ws->var(name.c_str())->setVal(BF_bs_val);
          fitresult_tex_vec.push_back(fitresult_tex2.str());
        }
        if (!bd_constr_ && !SM_ && BF_ < 2) {
          ostringstream fitresult_tex;
          fitresult_tex << setprecision(2) << fixed << "N(B_{d}) = " << N_bd->getVal() << " ^{+" << getErrorHigh(N_bd) << "}_{" << getErrorLow(N_bd) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
        }
        else if (bd_constr_) {
          ostringstream fitresult_tex;
          fitresult_tex << setprecision(2) << fixed << "N(B_{d}) / N(B_{s}) = " << N_bd->getVal() << " ^{+" << getErrorHigh(N_bd) << "}_{" << getErrorLow(N_bd) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
        }
        else if (BF_ > 1) {
          ostringstream fitresult_tex;
          fitresult_tex << setprecision(2) << scientific << "BF(B^{0}) = " << N_bd->getVal() << " ^{+" << getErrorHigh(N_bd) << "}_{" << getErrorLow(N_bd) << "}";
          fitresult_tex_vec.push_back(fitresult_tex.str());
          ostringstream fitresult_tex2;
          fitresult_tex2 << "(N(B^{0}) = " << setprecision(2) << fixed << ws->function(name("N_bd_formula", i, j))->getVal();
          Double_t BF_bd_val = ws->var("BF_bd")->getVal();
          Double_t N_bd_ = ws->function(name("N_bd_formula", i, j))->getVal();
          ws->var("BF_bd")->setVal(BF_bd_val + getErrorHigh(N_bd));
          Double_t N_bd_up = ws->function(name("N_bd_formula", i, j))->getVal();
          Double_t N_bd_error_up = N_bd_up - N_bd_;
          fitresult_tex2 << " ^{+" << N_bd_error_up;
          ws->var("BF_bd")->setVal(BF_bd_val);
          ws->var("BF_bd")->setVal(BF_bd_val + getErrorLow(N_bd));
          Double_t N_bd_down = ws->function(name("N_bd_formula", i, j))->getVal();
          Double_t N_bd_error_down = N_bd_ - N_bd_down;
          fitresult_tex2 << "}_{-" << N_bd_error_down << "}" << ")";
          ws->var("BF_bd")->setVal(BF_bd_val);
          fitresult_tex_vec.push_back(fitresult_tex2.str());
        }
        ostringstream fitresult_tex;
        fitresult_tex << setprecision(2) << fixed << "N(comb. bkg) = " << N_comb->getVal() << " ^{+" << getErrorHigh(N_comb) << "}_{" << getErrorLow(N_comb) << "}";
        fitresult_tex_vec.push_back(fitresult_tex.str());
        if (BF_ > 0) {
        	RooRealVar* N_semi = (RooRealVar*)vars->find(pdf_analysis::name("N_semi", i, j));
        	RooRealVar* N_peak = (RooRealVar*)vars->find(pdf_analysis::name("N_peak", i, j));
        	ostringstream fitresult_tex2;
        	fitresult_tex2 << setprecision(2) << fixed << "N(semi bkg) = " << N_semi->getVal() << " ^{+" << getErrorHigh(N_semi) << "}_{" << getErrorLow(N_semi) << "}";
        	fitresult_tex_vec.push_back(fitresult_tex2.str());
        	ostringstream fitresult_tex3;
        	fitresult_tex3 << setprecision(2) << fixed << "N(peak bkg) = " << N_peak->getVal() << " ^{+" << getErrorHigh(N_peak) << "}_{" << getErrorLow(N_peak) << "}";
        	fitresult_tex_vec.push_back(fitresult_tex3.str());
        }

        TPaveText* fitresults = new TPaveText(0.57, 0.54, 0.89, 0.89, "NDCR");
        for (unsigned int jj = 0; jj < fitresult_tex_vec.size(); jj++) {
          fitresults->AddText(fitresult_tex_vec[jj].c_str());
        }
        fitresults->SetFillColor(0);
        fitresults->SetShadowColor(0);
        fitresults->SetTextSize(0.03);
        fitresults->SetTextAlign(11);
        fitresults->SetLineColor(0);
        fitresults->Draw();
      }

      channel = i;
      channel_bdt = j;
      string output_name;
      if (output == "") output_name = get_address("data_" + var, "", true);
      else output_name = get_address(pdfname, output, true);
      final_c->Print( (output_name + ".gif").c_str() );
      final_c->Print( (output_name + ".pdf").c_str() );
      delete final_p;
      delete final_c;
    }
  }
}

void pdf_fitData::FillRooDataSet(RooDataSet* dataset, bool cut_b, vector <double> cut_, string cuts, TTree* tree, int offset) {
  int events = 0;
  if (!strcmp(tree->GetName(), "SgData_bdt")) {
    TTree* reduced_tree = tree->CopyTree(cuts.c_str());
    Double_t m1eta_t, m2eta_t, m_t, eta_B_t, bdt_t, me_t;
    Int_t evt_t;
    reduced_tree->SetBranchAddress("m1eta", &m1eta_t);
    reduced_tree->SetBranchAddress("m2eta", &m2eta_t);
    reduced_tree->SetBranchAddress("m", &m_t);
    reduced_tree->SetBranchAddress("me", &me_t);
    reduced_tree->SetBranchAddress("eta", &eta_B_t);
    reduced_tree->SetBranchAddress("bdt", &bdt_t);
    reduced_tree->SetBranchAddress("evt", &evt_t);
    for (int i = 0; i < reduced_tree->GetEntries(); i++) {
      reduced_tree->GetEntry(i);
      if (m_t > 4.9 && m_t < 5.9) {
      	if (me_t < 0.0 || me_t > 0.1) continue; //skip wrong mass scale
        events++;
        Mass->setVal(m_t);
        bdt->setVal(bdt_t);
        MassRes->setVal(me_t);
        Int_t eta_channel = -1;
        if (fabs(m1eta_t)<1.4 && fabs(m2eta_t)<1.4) {
          eta_channel = 0 + offset*2;
          if (cut_b && bdt_t < cut_[eta_channel]) continue;
          channels_cat->setIndex(eta_channel);
        }
        else {
          eta_channel = 1 + offset*2;
          if (cut_b && bdt_t < cut_[eta_channel]) continue;
          channels_cat->setIndex(eta_channel);
        }
        int bdt_channel = bdt_index(eta_channel, bdt_t);
        if (simul_bdt_ || simul_all_) {
          if (bdt_channel == -1) continue; /// bdt < 0.1
          bdt_cat->setIndex(bdt_channel);
        }
        if (simul_all_) all_cat->setIndex(super_index(eta_channel, bdt_channel));

        RooArgSet varlist_tmp(*Mass, *MassRes, *bdt, *channels_cat);
        if (simul_bdt_ || simul_all_) varlist_tmp.add(*bdt_cat);
        if (simul_all_) varlist_tmp.add(*all_cat);
        dataset->add(varlist_tmp);
      }
    }
  }
  else {
    cout << "tree name is not SgData_bdt" << endl;
    exit(EXIT_FAILURE);
  }
  cout << "total events = " << events << endl;
}

void pdf_fitData::define_dataset() {
  cout << red_color_bold << "defining dataset" << default_console_color << endl;
  RooArgSet varlist(*Mass, *MassRes, *bdt, *channels_cat);
  if (simul_bdt_ || simul_all_) varlist.add(*bdt_cat);
  if (simul_all_) varlist.add(*all_cat);
  global_data = new RooDataSet("global_data", "global_data", varlist);
}

void pdf_fitData::make_dataset(bool cut_b, vector <double> cut_, string cuts, TTree* tree, int offset) {
  cout << red_color_bold << "making dataset" << default_console_color << endl;

  if (!random) FillRooDataSet(global_data, cut_b, cut_, cuts, tree, offset);

  else {
    if (syst) randomize_constraints(ws_);

    if (!simul_) {
      ws_->var("N_bs")->setVal(estimate_bs[ch_i_]);
      if (!SM_ && !bd_constr_) ws_->var("N_bd")->setVal(estimate_bd[ch_i_]);
      else if (bd_constr_) ws_->var("bd_over_bs")->setVal(estimate_bd[ch_i_]/estimate_bd[ch_i_]);
      ws_->var("N_semi")->setVal(estimate_semi[ch_i_]);
      ws_->var("N_comb")->setVal(estimate_comb[ch_i_]);
      global_data = ws_->pdf("pdf_ext_total")->generate(RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->var("bdt"), *ws_->cat("etacat")), Extended());
    }
    else if (simul_ && !simul_bdt_ && !simul_all_) {
      for (int i = 0; i < channels; i++) {
        ws_->var(name("N_bs", i))->setVal(estimate_bs[i]);
        if (!SM_ && !bd_constr_) ws_->var(name("N_bd", i))->setVal(estimate_bd[i]);
        else if (bd_constr_) ws_->var("bd_over_bs")->setVal(estimate_bd[i]/estimate_bd[i]);
        ws_->var(name("N_semi", i))->setVal(estimate_semi[i]);
        ws_->var(name("N_comb", i))->setVal(estimate_comb[i]);
      }
      if (!asimov_) global_data = ws_->pdf("pdf_ext_simul")->generate(RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->var("bdt"), *ws_->cat("etacat")), Extended());
      else global_datahist = ws_->pdf("pdf_ext_simul")->generateBinned(RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->var("bdt"), *ws_->cat("etacat")), ExpectedData());
    }
    else if (simul_ && (simul_bdt_ || simul_all_)) {
      RooArgSet set(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->cat("etacat"), *ws_->cat("bdtcat"));
      if (simul_all_) set.add(*ws_->cat("allcat"));
      global_data = new RooDataSet("global_data", "global_data", set);
      for (int i = 0; i < channels; i++) {
      	if (years_ == "0" && i > 1) break;
      	if (years_ == "1" && i < 2) continue;
        for (int j = 0; j < bdt_index_max(i); j++) {
          ws_->var(name("N_bs", i, j))->setVal(estimate2D_bs[i][j]);
          if (!SM_ && !bd_constr_) ws_->var(name("N_bd", i, j))->setVal(estimate2D_bd[i][j]);
          else if (bd_constr_) ws_->var("bd_over_bs")->setVal(estimate2D_bd[i][j]/estimate2D_bd[i][j]);
          if (!rare_constr_) ws_->var(name("N_semi", i, j))->setVal(estimate2D_semi[i][j]);
          ws_->var(name("N_comb", i, j))->setVal(estimate2D_comb[i][j]);
          RooDataSet* data_i = ws_->pdf(name("pdf_ext_total", i, j))->generate(RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes")), Extended());
          channels_cat->setIndex(i);
          bdt_cat->setIndex(j);
          data_i->addColumn(*channels_cat);
          data_i->addColumn(*bdt_cat);
          if (simul_all_) {
            all_cat->setIndex(super_index(i, j));
            data_i->addColumn(*all_cat);
          }
          global_data->append(*data_i);
        }
      }
    }
  }
  global_data->SetName("global_data");
  cout << " entries = " <<  global_data->sumEntries() << endl;
}

void pdf_fitData::make_pdf_input(string root_s) {
  ws_file_input = new TFile(root_s.c_str());
  if (!ws_file_input) {cout << root_s.c_str() << " does not exist" << endl; exit(EXIT_FAILURE);}
  ws_input = (RooWorkspace*)ws_file_input->Get("ws");
  if (!ws_input) {cout << "ws does not exist" << endl; exit(EXIT_FAILURE);}
  cout << "ws file: " << root_s << endl;
}

void pdf_fitData::make_pdf() {
  cout << red_color_bold << "making pdf" << default_console_color << endl;
  if (random) {
    if (simul_ && !simul_bdt_ && !simul_all_) {
      for (int i = 0; i < channels; i++) {
      	ws_input->var(name("N_bs", i))->setVal(estimate_bs[i]);
        ws_input->var(name("N_bd", i))->setVal(estimate_bd[i]);
        ws_input->var(name("N_semi", i))->setVal(estimate_semi[i]);
        ws_input->var(name("N_comb", i))->setVal(estimate_comb[i]);
      }
    }
    else if (!simul_) {
      ws_input->var("N_bs")->setVal(estimate_bs[ch_i_]);
      ws_input->var("N_bd")->setVal(estimate_bd[ch_i_]);
      ws_input->var("N_semi")->setVal(estimate_semi[ch_i_]);
      ws_input->var("N_comb")->setVal(estimate_comb[ch_i_]);
    }
    else {
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < bdt_index_max(i); j++) {
          ws_input->var(name("N_bs", i, j))->setVal(estimate2D_bs[i][j]);
          ws_input->var(name("N_bd", i, j))->setVal(estimate2D_bd[i][j]);
          ws_input->var(name("N_semi", i, j))->setVal(estimate2D_semi[i][j]);
          ws_input->var(name("N_comb", i, j))->setVal(estimate2D_comb[i][j]);
        }
      }
    }
    if (bd_constr_) ws_input->var("Bd_over_Bs")->setVal(estimate_bd[0] / estimate_bs[0]);
  }

  set_ws(ws_input);
}

void pdf_fitData::set_final_pdf() {
	cout << red_color_bold << "setting final pdf" << default_console_color << endl;
	define_perchannel_pdf();
	if (simul_) define_simul();
	ws_->Print();
	cout << "done" << endl;
}

void pdf_fitData::define_perchannel_pdf () {
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      if (BF_ > 0) {
        define_constraints(i, j);
      }
      define_total_extended(i, j);
    }
  }
}

void pdf_fitData::define_constraints(int i, int j) {

  if (i == 0 && j == 0) {
    RooGaussian fs_over_fu_gau("fs_over_fu_gau", "fs_over_fu_gau", *ws_->var("fs_over_fu"), RooConst(fs_over_fu_val), RooConst(fs_over_fu_err));
    ws_->import(fs_over_fu_gau);

    RooGaussian one_over_BRBR_gau("one_over_BRBR_gau", "one_over_BRBR_gau", *ws_->var("one_over_BRBR"), RooConst(one_over_BRBR_val), RooConst(one_over_BRBR_err));
    ws_->import(one_over_BRBR_gau);
  }

  RooGaussian N_bu_gau(name("N_bu_gau", i, j), "N_bu_gau", *ws_->var(name("N_bu", i, j)), RooConst(N_bu_val[i][j]), RooConst(N_bu_err[i][j]));
  ws_->import(N_bu_gau);

  RooGaussian effratio_gau_bs(name("effratio_gau_bs", i, j), "effratio_gau_bs", *ws_->var(name("effratio_bs", i, j)), RooConst(effratio_bs_val[i][j]), RooConst(effratio_bs_err[i][j]));
  RooGaussian effratio_gau_bd(name("effratio_gau_bd", i, j), "effratio_gau_bd", *ws_->var(name("effratio_bd", i, j)), RooConst(effratio_bd_val[i][j]), RooConst(effratio_bd_err[i][j]));
  ws_->import(effratio_gau_bs);
  ws_->import(effratio_gau_bd);

  if (rare_constr_) {
  	RooGaussian N_peak_gau(name("N_peak_gau", i, j), "N_peak_gau", *ws_->var(name("N_peak", i, j)), RooConst(ws_->var(name("N_peak", i, j))->getVal()), RooConst(ws_->var(name("N_peak", i, j))->getError()));
  	ws_->import(N_peak_gau);
  	RooGaussian N_semi_gau(name("N_semi_gau", i, j), "N_semi_gau", *ws_->var(name("N_semi", i, j)), RooConst(ws_->var(name("N_semi", i, j))->getVal()), RooConst(ws_->var(name("N_semi", i, j))->getError()));
  	ws_->import(N_semi_gau);
  }
  int eta = -1;
  if (!simul_all_) eta = i / 2;
  else {
  	vector <int> indexes(get_EtaBdt_bins(i));
  	eta = indexes[0] / 2;
  }
  RooGaussian Mean_gau_bs(name("Mean_gau_bs", i, j), "Mean_gau_bs", *ws_->var(name("Mean_bs", i, j)), RooConst(ws_->var(name("Mean_bs", i, j))->getVal()), RooConst(mass_scale_sys[eta]));
  RooGaussian Mean_gau_bd(name("Mean_gau_bd", i, j), "Mean_gau_bd", *ws_->var(name("Mean_bd", i, j)), RooConst(ws_->var(name("Mean_bd", i, j))->getVal()), RooConst(mass_scale_sys[eta]));
  ws_->import(Mean_gau_bs);
  ws_->import(Mean_gau_bd);

  RooGaussian Enne_gau_bs(name("Enne_gau_bs", i, j), "Enne_gau_bs", *ws_->var(name("Enne_bs", i, j)), RooConst(ws_->var(name("Enne_bs", i, j))->getVal()), RooConst(ws_->var(name("Enne_bs", i, j))->getError()));
  RooGaussian Alpha_gau_bs(name("Alpha_gau_bs", i, j), "Alpha_gau_bs", *ws_->var(name("Alpha_bs", i, j)), RooConst(ws_->var(name("Alpha_bs", i, j))->getVal()), RooConst(ws_->var(name("Alpha_bs", i, j))->getError()));
  RooGaussian Enne_gau_bd(name("Enne_gau_bd", i, j), "Enne_gau_bd", *ws_->var(name("Enne_bd", i, j)), RooConst(ws_->var(name("Enne_bd", i, j))->getVal()), RooConst(ws_->var(name("Enne_bd", i, j))->getError()));
  RooGaussian Alpha_gau_bd(name("Alpha_gau_bd", i, j), "Alpha_gau_bd", *ws_->var(name("Alpha_bd", i, j)), RooConst(ws_->var(name("Alpha_bd", i, j))->getVal()), RooConst(ws_->var(name("Alpha_bd", i, j))->getError()));
  ws_->import(Enne_gau_bs);
  ws_->import(Alpha_gau_bs);
  ws_->import(Enne_gau_bd);
  ws_->import(Alpha_gau_bd);

}

void pdf_fitData::define_total_extended(int i, int j) {
  RooArgList pdf_list(*ws_->pdf(name("pdf_bs", i, j)), *ws_->pdf(name("pdf_bd", i, j)), *ws_->pdf(name("pdf_comb", i, j)));
  if (BF_ > 0) {
  	pdf_list.add(*ws_->pdf(name("pdf_semi", i, j)));
  	pdf_list.add(*ws_->pdf(name("pdf_peak", i, j)));
  }
  RooArgList N_list("varlist");

  if (BF_ > 0) N_list.add(*ws_->function(name("N_bs_formula", i, j)));
  else N_list.add(*ws_->var(name("N_bs", i, j)));
  if ((SM_ || bd_constr_) && BF_ < 2) N_list.add(*ws_->function(name("N_bd_constr", i, j)));
  else if (BF_ > 1) N_list.add(*ws_->function(name("N_bd_formula", i, j)));
  else N_list.add(*ws_->var(name("N_bd", i, j)));
  N_list.add(*ws_->var(name("N_comb", i, j)));
  if (BF_ > 0) {
  	N_list.add(*ws_->var(name("N_semi", i, j)));
  	N_list.add(*ws_->var(name("N_peak", i, j)));
  }

  RooAddPdf pdf_ext_sum(name("pdf_ext_sum", i, j), "pdf_ext_sum", pdf_list, N_list);
  RooArgList constraints_list(*ws_->pdf(name("N_bu_gau", i, j)), *ws_->pdf("fs_over_fu_gau"), *ws_->pdf(name("effratio_gau_bs", i, j)), *ws_->pdf("one_over_BRBR_gau"));
  if (rare_constr_) {
  	constraints_list.add(*ws_->pdf(name("N_peak_gau", i, j)));
  	constraints_list.add(*ws_->pdf(name("N_semi_gau", i, j)));
  }
  if (BF_ > 1) {
  	constraints_list.add(*ws_->pdf(name("effratio_gau_bd", i, j)));
  }
  constraints_list.add(*ws_->pdf(name("Mean_gau_bs", i, j)));
  constraints_list.add(*ws_->pdf(name("Mean_gau_bd", i, j)));

  constraints_list.add(*ws_->pdf(name("Alpha_gau_bs", i, j)));
  constraints_list.add(*ws_->pdf(name("Enne_gau_bs", i, j)));
  constraints_list.add(*ws_->pdf(name("Alpha_gau_bd", i, j)));
  constraints_list.add(*ws_->pdf(name("Enne_gau_bd", i, j)));
  RooProdPdf constraints_pdfs(name("pdf_constraints", i, j), "pdf_constraints", constraints_list);
  RooProdPdf pdf_ext_total(name("pdf_ext_total", i, j), "pdf_ext_total", RooArgList(pdf_ext_sum, constraints_pdfs));
  ws_->import(pdf_ext_total);

}

void pdf_fitData::define_simul() {
  if (!simul_bdt_ && !simul_all_) {
    RooSimultaneous pdf_sim("pdf_ext_simul", "simultaneous pdf", *ws_->cat("etacat"));
    RooSimultaneous pdf_sim_noconstr("pdf_ext_simul_noconstr", "simultaneous pdf without constraints", *ws_->cat("etacat"));
    for (int i = 0; i < channels; i++) {
    	if (years_ == "0" && i > 1) break;
    	if (years_ == "1" && i < 2) continue;
      pdf_sim.addPdf(*ws_->pdf(name("pdf_ext_total", i)), name("etacat", i));
      pdf_sim_noconstr.addPdf(*ws_->pdf(name("pdf_ext_sum", i)), name("etacat", i));
    }
    ws_->import(pdf_sim);
    ws_->import(pdf_sim_noconstr);
    pdf_sim.graphVizTree("pdf_ext_simul.dot");
    pdf_sim_noconstr.graphVizTree("pdf_sim_noconstr.dot");
  }

  if (!simul_bdt_ && simul_all_) {
    RooSimultaneous pdf_sim("pdf_ext_simul", "simultaneous pdf", *ws_->cat("allcat"));
    RooSimultaneous pdf_sim_noconstr("pdf_ext_simul_noconstr", "simultaneous pdf without constraints", *ws_->cat("allcat"));
    for (int i = 0; i < channels_all; i++) {
      vector <int> indexes(get_EtaBdt_bins(i));
      if (years_=="0" && indexes[0] > 1) continue;
      if (years_=="1" && indexes[0] < 2) continue;
      pdf_sim.addPdf(*ws_->pdf(name("pdf_ext_total", indexes[0], indexes[1])), Form("allcat_%d", i));
      pdf_sim_noconstr.addPdf(*ws_->pdf(name("pdf_ext_sum", indexes[0], indexes[1])), Form("allcat_%d", i));
    }
    ws_->import(pdf_sim);
    pdf_sim.graphVizTree("pdf_ext_simul_all.dot");
    ws_->import(pdf_sim_noconstr);
  }

  if (simul_bdt_ && !simul_all_) {
    RooSuperCategory* rsc = dynamic_cast<RooSuperCategory*> (ws_->obj("super_cat"));
    RooSimultaneous pdf_sim("pdf_ext_simul", "simultaneous pdf", *rsc);
    RooSimultaneous pdf_sim_noconstr("pdf_ext_simul_noconstr", "simultaneous pdf without constraints", *rsc);
    for (int i = 0; i < channels; i++) {
    	if (years_=="0" && i > 1) break;
    	if (years_=="1" && i < 2) continue;
      for (int j = 0; j < bdt_index_max(i); j++) {
        RooArgSet icl = rsc->inputCatList();
        RooCategory* eta_c = (RooCategory*)icl.find("etacat");
        RooCategory* bdt_c = (RooCategory*)icl.find("bdtcat");
        eta_c->setIndex(i);
        bdt_c->setIndex(j);
        cout << rsc->getLabel() << " (" << rsc->getIndex() << ") " <<  i << " " << j << endl;
        pdf_sim.addPdf(*ws_->pdf(Form("pdf_ext_total_%d_%d", i, j)), rsc->getLabel());
        pdf_sim_noconstr.addPdf(*ws_->pdf(Form("pdf_ext_sum_%d_%d", i, j)), rsc->getLabel());
      }
    }
    ws_->import(pdf_sim);
    pdf_sim.graphVizTree("pdf_ext_simulBdt.dot");
    if (BF_ > 0) {
      ws_->import(pdf_sim_noconstr);
    }
  }

  if (simul_bdt_ && simul_all_) {
    cout << "simul_bdt_ can't be with simul_all_" << endl;
    exit(1);
  }
}

void pdf_fitData::set_syst() {

	ws_->var("fs_over_fu")->setConstant(!syst);
	ws_->var("one_over_BRBR")->setConstant(!syst);
	for (int i = 0; i < channels; i++) {
		for (int j = 0; j < bdt_index_max(i); j++) {
			ws_->var(name("N_bu", i, j))->setConstant(!syst);
			ws_->var(name("effratio_bs", i, j))->setConstant(!syst);
			if (rare_constr_) {
				ws_->var(name("N_peak", i, j))->setConstant(!syst);
				ws_->var(name("N_semi", i, j))->setConstant(!syst);
      }
			else {
				ws_->var(name("N_peak", i, j))->setConstant(true);
				ws_->var(name("N_semi", i, j))->setConstant(true);
			}
			if (BF_ > 1) ws_->var(name("effratio_bd", i, j))->setConstant(!syst);
			ws_->var(name("Mean_bs", i, j))->setConstant(!syst);
			ws_->var(name("Mean_bd", i, j))->setConstant(!syst);
			ws_->var(name("Enne_bs", i, j))->setConstant(!syst);
			ws_->var(name("Alpha_bs", i, j))->setConstant(!syst);
			ws_->var(name("Enne_bd", i, j))->setConstant(!syst);
			ws_->var(name("Alpha_bd", i, j))->setConstant(!syst);

//			if (BF_ == 0) {
//				ws_->var(name("Mean_bs", i, j))->setConstant(false);
//				ws_->var(name("Mean_bd", i, j))->setConstant(false);
//			}
		}
	}
	string constraints("");
	if (syst) {
		constraints += "fs_over_fu,one_over_BRBR";
		for (int i = 0; i < channels; i++) {
			if (years_=="0" && i > 1) break;
			if (years_=="1" && i < 2) continue;
			for (int j = 0; j < bdt_index_max(i); j++) {
				constraints += "," + (string)name("N_bu", i, j);
				constraints += "," + (string)name("effratio_bs", i, j);
				if (rare_constr_) {
					constraints += "," + (string)name("N_peak", i, j);
					constraints += "," + (string)name("N_semi", i, j);
				}
				if (BF_ > 1) constraints += "," + (string)name("effratio_bd", i, j);
				constraints += "," + (string)name("Mean_bs", i, j);
				constraints += "," + (string)name("Mean_bd", i, j);

				constraints += "," + (string)name("Enne_bs", i, j);
				constraints += "," + (string)name("Alpha_bs", i, j);
				constraints += "," + (string)name("Enne_bd", i, j);
				constraints += "," + (string)name("Alpha_bd", i, j);
			}
		}
	}
	ws_->defineSet("constr", constraints.c_str());
	cout << "constrain set: ";
	ws_->set("constr")->Print();
}

void pdf_fitData::save() {
  ostringstream output_ws;
  output_ws << "./output/ws_fitData_" << meth_;
  if (simul_) output_ws << "_simul" << channels;
  else output_ws << "_" << ch_s_;
  if (SM_)             output_ws << "_SM";
  else if (bd_constr_) output_ws << "_BdConst";
  if (bdt_fit_) output_ws << "_2D";
  if (pee) output_ws << "_PEE";
  output_ws << ".root";
  ws_->SaveAs(output_ws.str().c_str());
}

void pdf_fitData::significance() {
	cout << red_color_bold << "significance" << default_console_color << endl;
//  ProfileLikelihoodTestStat::SetAlwaysReuseNLL(true);
//  RatioOfProfiledLikelihoodsTestStat::SetAlwaysReuseNLL(true);
	if (sign < 0) return;
//	make_prior();
  if (sign == 0) sig_hand();
  else {
    make_models();
    if (sign == 1) sig_plhc();
    if (sign == 2) sig_plhts();
    if (sign == 3) sig_hybrid_plhts();
    if (sign == 4) sig_hybrid_roplhts();
  }
}

Double_t pdf_fitData::sig_hand() {
  /// by hand
//  RooProdPdf* modelpdf = new RooProdPdf("model_pdf", "model_pdf", *ws_->pdf(pdfname.c_str()), *ws_->pdf("prior") );
//	ws_->import(*modelpdf);
	RooAbsPdf* modelpdf = ws_->pdf(pdfname.c_str());
  Double_t minNLL = fit_pdf(true, modelpdf->GetName())->minNll();
  RooRealVar *arg_var = 0;
  RooArgSet *all_vars = ws_->pdf(modelpdf->GetName())->getVariables();
  TIterator* vars_it = all_vars->createIterator();
  size_t found;
  string alt_name("N_bs");
  if (BF_ > 0) alt_name = "BF_bs";
  if (Bd && BF_ == 0) alt_name = "N_bd";
  if (Bd && BF_ > 0) alt_name = "BF_bd";
  double null = 0.;
  if (SMIsNull && BF_ > 1) {
    if (Bd) null = Bd2MuMu_SM_BF_val;
    else null = Bs2MuMu_SM_BF_val;
  }
  else if (SMIsNull) {
    cout << "SMIsNull works only with BF = 2" << endl;
    exit(1);
  }
  while ( (arg_var = (RooRealVar*)vars_it->Next())) {
    string name(arg_var->GetName());
    found = name.find(alt_name);
    if (found != string::npos) {
      arg_var->setVal(null);
      arg_var->setConstant(1);
    }
    if (SM_ || bd_constr_){
      found = name.find("Bd_over_Bs");
      if (found!=string::npos) {
        arg_var->setVal(0);
        arg_var->setConstant(1);
      }
    }
  }
  Double_t newNLL = fit_pdf(true, modelpdf->GetName())->minNll();

  TIterator* vars_after = all_vars->createIterator();
  while ( (arg_var = (RooRealVar*)vars_after->Next())) {
    string name(arg_var->GetName());
    found = name.find(alt_name);
    if (found != string::npos) arg_var->setConstant(0);
    if (SM_ || bd_constr_) {
      found = name.find("Bd_over_Bs");
      if (found!=string::npos) arg_var->setConstant(0);
    }
  }

  Double_t deltaLL = newNLL - minNLL;
  Double_t signif = deltaLL>0 ? sqrt(2*deltaLL) : -sqrt(-2*deltaLL) ;
  if (verbosity > 0) {
    cout << "H1 minNLL = " << minNLL << endl;
    cout << "H0 minNLL = " << newNLL << endl;
    cout << "significance (by hand) = " << signif << endl << endl;
  }
  return signif;
}

void pdf_fitData::sig_plhc() {
  ModelConfig model;
  model.SetName("model");
  RooArgSet poi;
  RooArgSet CO;
  if (pee) {
    CO.add(*ws_->var("MassRes"));
  }
  model.SetWorkspace(*ws_);
//  RooProdPdf* modelpdf = new RooProdPdf("model_pdf", "model_pdf", *ws_->pdf(pdfname.c_str()), *ws_->pdf("prior") );
  model.SetPdf(*ws_->pdf(pdfname.c_str()));
//  model.SetPdf(*modelpdf);

  string alt_name("N_bs");
  if (BF_ > 0) alt_name = "BF_bs";
  if (Bd && BF_ == 0) alt_name = "N_bd";
  if (Bd && BF_ > 0) alt_name = "BF_bd";
  double null = 0.;
  if (SMIsNull && BF_ > 1) {
    if (Bd) null = Bd2MuMu_SM_BF_val;
    else null = Bs2MuMu_SM_BF_val;
  }
  else if (SMIsNull) {
    cout << "SMIsNull works only with BF = 2" << endl;
    exit(1);
  }

  if (BF_ == 0) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < bdt_index_max(i); j++) {
        poi.add(*ws_->var(name(alt_name.c_str(), i, j)));
        poi.setRealValue(name(alt_name.c_str(), i, j), 0);
      }
    }
  }
  else {
    poi.add(*ws_->var(alt_name.c_str()));
    poi.setRealValue(alt_name.c_str(), null);
  }
  if (bd_constr_) {
    poi.add(*ws_->var("Bd_over_Bs"));
    poi.setRealValue("Bd_over_Bs", 0);
  }

  ProfileLikelihoodCalculator plc;
  plc.SetData(*ws_->data("global_data"));
  plc.SetModel(model);
  if (pee) plc.SetConditionalObservables(CO);
  plc.SetNullParameters(poi);
  HypoTestResult* htr = plc.GetHypoTest();
  cout << "ProfileLikelihoodCalculator: The p-value for the null is " << htr->NullPValue() << "; The significance for the null is " << htr->Significance() << endl;
  model.SetSnapshot(poi);
  ws_->import(model);
}

void pdf_fitData::make_models() {
  vector <vector <double> > N_bs(channels, vector <double> (channels_bdt) );
  vector <vector <double> > N_bd(channels, vector <double> (channels_bdt) );
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      if (BF_ == 0) N_bs[i][j] = ws_->var(name("N_bs", i, j))->getVal();
      else N_bs[i][j] = ws_->function(name("N_bs_formula", i, j))->getVal();
      if (!SM_ && !bd_constr_ && BF_ < 2) N_bd[i][j] = ws_->var(name("N_bd", i, j))->getVal();
      else if (BF_ > 1) N_bd[i][j] = ws_->function(name("N_bd_formula", i, j))->getVal();
    }
  }

  /// obs
  string observables("Mass");
  if (bdt_fit_) observables += ",bdt";
  if (simul_) {
  	if (!simul_bdt_ && !simul_all_) observables += ",etacat";
  	else if (simul_bdt_) observables += ",bdtcat";
  	else if (simul_all_) observables += ",allcat";
  }
  ws_->defineSet("obs", observables.c_str());

  string alt_name("N_bs");
  if (BF_ > 0) alt_name = "BF_bs";
  if (Bd && BF_ == 0) alt_name = "N_bd";
  if (Bd && BF_ > 0) alt_name = "BF_bd";

  /// poi
  ostringstream name_poi;
  if (BF_ == 0) {
    if (simul_) {
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < bdt_index_max(i); j++) {
          if (i != 0 || j != 0) name_poi << ",";
          name_poi << name(alt_name.c_str(), i, j);
        }
      }
    }
    else {
      name_poi << alt_name;
    }
  }
  else {
    name_poi << alt_name;
  }
  if (bd_constr_) name_poi << ",Bd_over_Bs";
  ws_->defineSet("poi", name_poi.str().c_str());

  /// nui
  RooArgSet nuisanceParams;
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      nuisanceParams.add(*ws_->var(name("N_comb", i, j)));
      nuisanceParams.add(*ws_->var(name("N_semi", i, j)));
      if (rare_constr_) nuisanceParams.add(*ws_->var(name("N_peak", i, j)));
      nuisanceParams.add(*ws_->var(name("exp_comb", i, j)));
      if (!SM_ && !bd_constr_ && BF_ < 2 && !Bd) nuisanceParams.add(*ws_->var(name("N_bd", i, j)));
      else if (Bd && BF_ < 1) nuisanceParams.add(*ws_->var(name("N_bs", i, j)));
      if (BF_ > 0 && syst) {
        nuisanceParams.add(*ws_->var(name("effratio_bs", i, j)));
        nuisanceParams.add(*ws_->var(name("N_bu", i, j)));
        if (BF_ > 1) {
          nuisanceParams.add(*ws_->var(name("effratio_bd", i, j)));
        }
      }
    }
  }
  if (BF_ > 0) {
    nuisanceParams.add(*ws_->var("fs_over_fu"));
    nuisanceParams.add(*ws_->var("one_over_BRBR"));
  }
  if (BF_ > 1 && !Bd) nuisanceParams.add(*ws_->var("BF_bd"));
  else if (BF_ > 1 && Bd) nuisanceParams.add(*ws_->var("BF_bs"));
  ws_->defineSet("nui", nuisanceParams);

  double null = 0.;
  if (SMIsNull && BF_ > 1) {
    if (Bd) null = Bd2MuMu_SM_BF_val;
    else null = Bs2MuMu_SM_BF_val;
  }
  else if (SMIsNull) {
    cout << "SMIsNull works only with BF = 2" << endl;
    exit(1);
  }

  RooArgSet CO;
  if (pee) {
    CO.add(*ws_->var("MassRes"));
    ws_->defineSet("CO", "MassRes");
  }

// RooProdPdf* modelpdf = new RooProdPdf("model_pdf", "model_pdf", *ws_->pdf(pdfname.c_str()), *ws_->pdf("prior") );
//  ws_->import(*modelpdf);

  ModelConfig* H1 = new ModelConfig("H1", "background + signal hypothesis", ws_);
  if (pee) H1->SetConditionalObservables(*ws_->set("CO"));
  H1->SetPdf(*ws_->pdf(pdfname.c_str()));
  H1->SetParametersOfInterest(*ws_->set("poi"));
  H1->SetObservables(*ws_->set("obs"));
  H1->SetNuisanceParameters(*ws_->set("nui"));
  if (BF_ > 0 && syst) H1->SetConstraintParameters(*ws_->set("constr"));
  H1->SetSnapshot(*ws_->set("poi"));

  ModelConfig* H0 = new ModelConfig("H0", "null hypothesis", ws_);
  if (pee) H0->SetConditionalObservables(*ws_->set("CO"));
  H0->SetPdf(*ws_->pdf(pdfname.c_str()));
  H0->SetParametersOfInterest(*ws_->set("poi"));
  H0->SetObservables(*ws_->set("obs"));
  H0->SetNuisanceParameters(*ws_->set("nui"));
  if (BF_ > 0 && syst) H0->SetConstraintParameters(*ws_->set("constr"));
  if (BF_ == 0) {
    if (simul_) {
      for (int i = 0; i < channels; i++) {
        for (int j = 0; j < bdt_index_max(i); j++) {
          ws_->var(name(alt_name.c_str(), i, j))->setVal(0.0);
        }
      }
    }
    else {
      ws_->var(alt_name.c_str())->setVal(0.0);
    }
  }
  else {
    ws_->var(alt_name.c_str())->setVal(null);
  }
  if (bd_constr_) {
    ws_->var("Bd_over_Bs")->setVal(0.0);
  }
  H0->SetSnapshot(*ws_->set("poi"));

  ws_->import(*H0);
  ws_->import(*H1);

  RooAbsPdf* nuisPdf = RooStats::MakeNuisancePdf(*H0,"nui_pdf");
  cout <<"nuisance pdf = " << endl;
  nuisPdf->Print();
  ws_->import(*nuisPdf);
}

void pdf_fitData::sig_plhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  ws_->Print();

  ProfileLikelihoodTestStat pl_ts(*ws_->pdf(pdfname.c_str()));
//  pl_ts.SetPrintLevel(2);
  pl_ts.SetOneSidedDiscovery(true);
  if (pee) pl_ts.SetConditionalObservables(*ws_->set("CO"));

  ProofConfig* pc = NULL;
  pc = new ProofConfig(*ws_, proof, Form("workers=%d", proof), kTRUE); // machine with 4 cores

  ToyMCSampler *mcSampler_pl = new ToyMCSampler(pl_ts, NExp);
  if (pc && proof > 1) mcSampler_pl->SetProofConfig(pc);

  FrequentistCalculator frequCalc(*ws_->data("global_data"), *H1,*H0, mcSampler_pl); // null = bModel interpreted as signal, alt = s+b interpreted as bkg
  frequCalc.SetToys(NExp, NExp/4.);
  HypoTestResult *htr_pl = frequCalc.GetHypoTest();
  plot_hypotest(htr_pl);
  cout << "ProfileLikelihoodTestStat + frequentist: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << " \\pm " << htr_pl->SignificanceError() << endl;
}

void pdf_fitData::sig_hybrid_plhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  ws_->Print();

  ProfileLikelihoodTestStat pl_ts(*ws_->pdf(pdfname.c_str()));
//  pl_ts.SetPrintLevel(2);
  pl_ts.SetOneSidedDiscovery(true);
  if (pee) pl_ts.SetConditionalObservables(*ws_->set("CO"));

  ProofConfig* pc = NULL;
  pc = new ProofConfig(*ws_, proof, Form("workers=%d", proof), kTRUE); // machine with "proof" cores

  ToyMCSampler *mcSampler_pl = new ToyMCSampler(pl_ts, NExp);
  if (pc && proof > 1) mcSampler_pl->SetProofConfig(pc);

//  RooSimultaneous *sim  = dynamic_cast<RooSimultaneous *>(ws_->pdf("pdf_ext_simul"));
//  RooAbsCategoryLValue *cat = (RooAbsCategoryLValue *) sim->indexCat().clone(sim->indexCat().GetName());
//  for (int ic = 0, nc = cat->numBins((const char *)0); ic < nc; ++ic) {
//    cat->setBin(ic);
//    cout << cat->getLabel() << "";
//    RooAbsPdf* pdf_i = sim->getPdf(cat->getLabel());
//    cout << pdf_i->GetName() << endl;
//  }

  HybridCalculator hibrCalc(*ws_->data("global_data"), *H1, *H0, mcSampler_pl);
  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("nui_pdf"));
  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("nui_pdf"));
//  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("prior"));
//  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("prior"));
  hibrCalc.SetToys(NExp, NExp/4);

  HypoTestResult *htr_pl = hibrCalc.GetHypoTest();
  plot_hypotest(htr_pl);
  cout << "ProfileLikelihoodTestStat + hybrid: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << " \\pm " << htr_pl->SignificanceError() << endl;
}

void pdf_fitData::sig_hybrid_roplhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  ws_->Print();

  RatioOfProfiledLikelihoodsTestStat ropl_ts(*H1->GetPdf(),*H0->GetPdf(), ws_->set("poi"));
  ropl_ts.SetSubtractMLE(false);
  if (pee) ropl_ts.SetConditionalObservables(*ws_->set("CO"));

  ProofConfig* pc = NULL;
  pc = new ProofConfig(*ws_, proof, Form("workers=%d", proof), kTRUE); // machine with 4 cores

  ToyMCSampler *mcSampler_pl = new ToyMCSampler(ropl_ts, NExp);
  if(pc && proof > 1) mcSampler_pl->SetProofConfig(pc);

  HybridCalculator hibrCalc(*ws_->data("global_data"), *H1, *H0, mcSampler_pl);
  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("nui_pdf"));
  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("nui_pdf"));
//  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("prior"));
//  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("prior"));
  hibrCalc.SetToys(NExp, NExp/4.);

  HypoTestResult *htr_pl = hibrCalc.GetHypoTest();
  plot_hypotest(htr_pl);
  cout << "RatioOfProfiledLikelihoodsTestStat + hybrid: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << " \\pm " << htr_pl->SignificanceError() << endl;
}

void pdf_fitData::plot_hypotest(HypoTestResult *hts) {
  hts->Print();
  TCanvas c_plot("c_plot", "c_plot", 600, 600);
  HypoTestPlot * plot = new HypoTestPlot(*hts, 100);
  plot->SetLogYaxis(true);
  plot->Draw();
  c_plot.Print( (get_address(Form("hypotest%d", sign), Bd ? "Bd" : "") + ".gif").c_str());
  c_plot.Print( (get_address(Form("hypotest%d", sign), Bd ? "Bd" : "") + ".pdf").c_str());
  c_plot.Print( (get_address(Form("hypotest%d", sign), Bd ? "Bd" : "") + ".root").c_str());
  c_plot.Print( (get_address(Form("hypotest%d", sign), Bd ? "Bd" : "") + ".C").c_str());
}

void pdf_fitData::make_prior() {
  cout << "making prior" << endl;
  vector <vector <RooGaussian*> > prior_bd(channels, vector <RooGaussian*> (channels_bdt));
  vector <vector <RooGaussian*> > prior_semi(channels, vector <RooGaussian*> (channels_bdt));
  vector <vector <RooGaussian*> > prior_comb(channels, vector <RooGaussian*> (channels_bdt));
  vector <vector <RooGaussian*> > prior_exp(channels, vector <RooGaussian*> (channels_bdt));

  vector <vector <RooGamma*> > prior_comb2(channels, vector <RooGamma*> (channels_bdt));

  RooArgList prior_list("prior_list");

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      if (!SM_ && !bd_constr_ && BF_ < 2 && !Bd) prior_bd[i][j] = new RooGaussian(name("prior_bd", i, j), name("prior_bd", i, j), *ws_->var(name("N_bd", i, j)), RooConst(ws_->var(name("N_bd", i, j))->getVal()), RooConst(ws_->var(name("N_bd", i, j))->getError()));
      else if (BF_ < 2 && Bd) prior_bd[i][j] = new RooGaussian(name("prior_bs", i, j), name("prior_bs", i, j), *ws_->var(name("N_bs", i, j)), RooConst(ws_->var(name("N_bs", i, j))->getVal()), RooConst(ws_->var(name("N_bs", i, j))->getError()));
      prior_semi[i][j] = new RooGaussian(name("prior_semi", i, j), name("prior_semi", i, j), *ws_->var(name("N_semi", i, j)), RooConst(ws_->var(name("N_semi", i, j))->getVal()), RooConst(ws_->var(name("N_semi", i, j))->getError()));
      prior_comb[i][j] = new RooGaussian(name("prior_comb", i, j), name("prior_comb", i, j), *ws_->var(name("N_comb", i, j)), RooConst(ws_->var(name("N_comb", i, j))->getVal()), RooConst(ws_->var(name("N_comb", i, j))->getError()));
      prior_exp[i][j] = new RooGaussian(name("prior_exp", i, j), name("prior_exp", i, j), *ws_->var(name("exp_comb", i, j)), RooConst(ws_->var(name("exp_comb", i, j))->getVal()), RooConst(ws_->var(name("exp_comb", i, j))->getError()));
      prior_list.add(*prior_bd[i][j]);
      if (!rare_constr_) {
      	prior_list.add(*prior_semi[i][j]);
      }
      prior_list.add(*prior_comb[i][j]);
      prior_list.add(*prior_exp[i][j]);

//      prior_comb2[i][j] = new RooGamma(name("prior_comb", i, j), name("prior_comb", i, j), *ws_->var(name("N_comb", i, j)), RooConst(ws_->var(name("N_comb", i, j))->getVal() + 1), RooConst(1.), RooConst(0.));
//      prior_list.add(*prior_comb2[i][j]);
    }
  }
  if (BF_ > 1 && !Bd) {
    RooGaussian* prior_bf_bd = new RooGaussian("prior_bf_bd", "prior_bf_bd", *ws_->var("BF_bd"), RooConst(ws_->var("BF_bd")->getVal()), RooConst(ws_->var("BF_bd")->getError()));
    prior_list.add(*prior_bf_bd);
  }
  else if (BF_ > 1 && Bd) {
    RooGaussian* prior_bf_bs = new RooGaussian("prior_bf_bs", "prior_bf_bs", *ws_->var("BF_bs"), RooConst(ws_->var("BF_bs")->getVal()), RooConst(ws_->var("BF_bs")->getError()));
    prior_list.add(*prior_bf_bs);
  }

  RooProdPdf prior("prior", "prior", prior_list);
  ws_->import(prior);
}

void pdf_fitData::setnewlumi() {
  if (BF_ > 0) {
    for (int i = 0; i < channels; i++) {
      for (int j = 0; j < bdt_index_max(i); j++) {
        double old_val = ws_->var(name("N_bu", i, j))->getVal();
        ws_->var(name("N_bu", i, j))->setVal(old_val * lumi);
        cout << "channel " << i << " " << j << "; Bs expected: " << ws_->function(name("N_bs_formula", i, j))->getVal() << "; Bd expected: " << ws_->function(name("N_bd_formula", i, j))->getVal() << endl;
      }
    }
  }
}

void pdf_fitData::randomize_constraints(RooWorkspace* ws) {
   //generating random constrains
  if (BF_ == 0) {
    cout << "no BF!" << endl;
    return;
  }
  if (!randomsyst) return;

  RooDataSet* fs_over_fu_ds = ws->pdf("fs_over_fu_gau")->generate(RooArgSet(*ws->var("fs_over_fu")), 1);
  ws->var("fs_over_fu")->setVal(fs_over_fu_ds->get(0)->getRealValue("fs_over_fu"));

  RooDataSet* one_over_BRBR_ds = ws->pdf("one_over_BRBR_gau")->generate(RooArgSet(*ws->var("one_over_BRBR")), 1);
  ws->var("one_over_BRBR")->setVal(one_over_BRBR_ds->get(0)->getRealValue("one_over_BRBR"));

  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      RooDataSet* N_bu_ds = ws->pdf(name("N_bu_gau", i, j))->generate(RooArgSet(*ws->var(name("N_bu", i, j))), 1);
      RooDataSet* effratio_bs_ds = ws->pdf(name("effratio_gau_bs", i, j))->generate(RooArgSet(*ws->var(name("effratio_bs", i, j))), 1);
      if (rare_constr_) {
      	RooDataSet* N_peak_ds = ws->pdf(name("N_peak_gau", i, j))->generate(RooArgSet(*ws->var(name("N_peak", i, j))), 1);
      	RooDataSet* N_semi_ds = ws->pdf(name("N_semi_gau", i, j))->generate(RooArgSet(*ws->var(name("N_semi", i, j))), 1);
      	ws->var(name("N_peak", i, j))->setVal(N_peak_ds->get(0)->getRealValue(name("N_peak", i, j)));
      	ws->var(name("N_semi", i, j))->setVal(N_semi_ds->get(0)->getRealValue(name("N_semi", i, j)));
      }
      ws->var(name("N_bu", i, j))->setVal(N_bu_ds->get(0)->getRealValue(name("N_bu", i, j)));
      ws->var(name("effratio_bs", i, j))->setVal(effratio_bs_ds->get(0)->getRealValue(name("effratio_bs", i, j)));
      if (BF_ > 1) {
        RooDataSet* effratio_bd_ds = ws->pdf(name("effratio_gau_bd", i, j))->generate(RooArgSet(*ws->var(name("effratio_bd", i, j))), 1);
        ws->var(name("effratio_bd", i, j))->setVal(effratio_bd_ds->get(0)->getRealValue(name("effratio_bd", i, j)));
      }
    }
  }
}

void pdf_fitData::extract_N_inRanges() {
  if (BF_ < 2) return;
  cout << red_color_bold << "extracting events in ranges..." << default_console_color << endl;
  string full_output = "output/yields.tex";
  FILE* file_out = fopen(full_output.c_str(), "w");

  string title_i[5] = {"$N_{B_s^0}$", "$N_{B^0}$", "$N_{\\mathrm{peak}}$", "$N_{\\mathrm{semi}}$", "$N_{\\mathrm{comb}}$"};
  string N_i[5] = {"N_bs_formula", "N_bd_formula", "N_peak", "N_semi", "N_comb"};
  string pdf_i[5] = {"pdf_bs", "pdf_bd", "pdf_peak", "pdf_semi", "pdf_comb"};

  fprintf(file_out, "\\begin{table}\n");
  fprintf(file_out, "\\centering\n");
  fprintf(file_out, "\\caption{Final invariant mass yields evaluated with the UML}\n");
  fprintf(file_out, "\\label{tab:UMLfinalyields}\n");
  fprintf(file_out, "\\begin{tabular}{|l|c|c|c|c|c|}\n");
  fprintf(file_out, "\\hline \n");

  for (int i = 0; i < channels; i++) {
    fprintf(file_out, "\\hline \n");
    fprintf(file_out, " \\multicolumn{6}{|c|}{Channel: %s %s} \\\\ \n", i%2==0? "barrel" : "endcap", i < 2 ? "2011" : "2012");
    fprintf(file_out, "\\hline \n");
    fprintf(file_out, "Variable  & low SB & $B^0$ window & $B_s^0$ window & high SB & \\textbf{all} \\\\ \n");
    fprintf(file_out, "\\hline \n");
    for (int j = 0; j < bdt_index_max(i); j++) {
      vector < Double_t> total(5, 0.);
      for (int l = 0; l < 5; l++) {
        fprintf(file_out, "%s ", title_i[l].c_str());
        for (unsigned int k = 0; k < massrange_names.size(); k++) {
          RooAbsReal* rar = ws_->pdf(name(pdf_i[l].c_str(), i, j))->createIntegral(RooArgSet(*ws_->var("Mass")), RooArgSet(*ws_->var("Mass")), massrange_names[k].c_str());
          Double_t N;
          if (l < 2) {
            RooAbsReal* rrv = ws_->function(name(N_i[l].c_str(), i, j));
            N = rrv->getVal() * rar->getVal();
          }
          else {
            RooRealVar* rrv = ws_->var(name(N_i[l].c_str(), i, j));
            N = rrv->getVal() * rar->getVal();
          }
          fprintf(file_out, "& %.2f ", N);
          total[k] += N;
        }
        Double_t N_all;
        if (l < 2) {
          N_all = ws_->function(name(N_i[l].c_str(), i, j))->getVal();
        }
        else {
          N_all = ws_->var(name(N_i[l].c_str(), i, j))->getVal();
        }
        fprintf(file_out, "& %.2f ", N_all);
        fprintf(file_out, " \\\\ \n");
        total[4] += N_all;
      }
      fprintf(file_out, "\\hline \n");
      fprintf(file_out, "$N_{\\mathrm{all}}$ ");
      for (unsigned int k = 0; k < massrange_names.size() + 1; k++) {
        fprintf(file_out, "& %.2f ", total[k]);
      }
      fprintf(file_out, " \\\\ \n");
    }
    fprintf(file_out, "\\hline \n");
  }
  fprintf(file_out, "\\end{tabular} \n");
  fprintf(file_out, "\\end{table} \n");
  fclose(file_out);
//  system(Form("cat %s", full_output.c_str()));
  cout << "tex file saved in " << full_output << endl;
}

void pdf_fitData::profile_NLL() {
  if (BF_ < 2) return;

//  RooNLLVar nll0("nll0","nll0",ws_->pdf(pdfname.c_str()),*global_data) ;

  // Construct unbinned likelihood
  RooAbsReal* nll = ws_->pdf(pdfname.c_str())->createNLL(*global_data, NumCPU(2), Extended(), pee ? ConditionalObservables(*ws_->var("MassRes")) : RooCmdArg::none(), syst ? Constrain(*ws_->set("constr")) : RooCmdArg::none()) ;
  // Minimize likelihood w.r.t all parameters before making plots
  RooMinuit(*nll).migrad() ;

  string var_alt("BF_bs");
  if (Bd) var_alt = "BF_bd";
  RooAbsReal* pll = nll->createProfile(*ws_->var(var_alt.c_str())) ;
  RooFormulaVar* double_pll = new RooFormulaVar("double_pll", "double_pll", "2*@0", RooArgList(*pll));
  // Plot the profile likelihood
  RooPlot* frame = ws_->var(var_alt.c_str())->frame(Bins(20), Range(0, Bd ? 1.4e-9 : 4e-9), Title(Form("profileLL in %s", var_alt.c_str()))) ;

  double_pll->plotOn(frame, LineColor(kRed)) ;
//  nll0.plotOn(frame, ShiftToZero()) ;

  TCanvas *c = new TCanvas("c","c",600, 600);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->SetYTitle("- 2 ln L");
  frame->Draw();
  c->Print((get_address("profileLL", var_alt, false) + ".gif").c_str());
  c->Print((get_address("profileLL", var_alt, false) + ".pdf").c_str());
}

void pdf_fitData::hack_ws(string frozen_ws_address) {
  cout << "hacking 2011 shape with 2012 shape from " << frozen_ws_address << endl;
  TFile * frozen_f = new TFile(frozen_ws_address.c_str());
  RooWorkspace * frozen_ws = (RooWorkspace*)frozen_f->Get("ws");
  for (int i = 2; i < 4; i++) {
    ws_->var(name("C0_semi", i-2))->setVal(frozen_ws->var(name("C0_semi", i))->getVal());
    ws_->var(name("C1_semi", i-2))->setVal(frozen_ws->var(name("C1_semi", i))->getVal());
    ws_->var(name("C2_semi", i-2))->setVal(frozen_ws->var(name("C2_semi", i))->getVal());
    ws_->var(name("C3_semi", i-2))->setVal(frozen_ws->var(name("C3_semi", i))->getVal());
    ws_->var(name("tau_semi", i-2))->setVal(frozen_ws->var(name("tau_semi", i))->getVal());
  }
}

void pdf_fitData::reset_minmax() {
  for (int i = 0; i < channels; i++) {
    for (int j = 0; j < bdt_index_max(i); j++) {
      ws_->var(name("N_peak", i, j))->setMax(ws_->var(name("N_peak", i, j))->getVal() + 10 * ws_->var(name("N_peak", i, j))->getError());
      ws_->var(name("N_semi", i, j))->setMax(ws_->var(name("N_semi", i, j))->getVal() + 10 * ws_->var(name("N_semi", i, j))->getError());
      ws_->var(name("N_comb", i, j))->setMax(ws_->var(name("N_comb", i, j))->getVal() + 10 * ws_->var(name("N_comb", i, j))->getError());
    }
  }
}
