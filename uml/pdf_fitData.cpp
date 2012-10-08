#include "pdf_fitData.h"

pdf_fitData::pdf_fitData(bool print, int inputs, string input_estimates, string meth, string range, bool SM, bool bd_constr, TTree *input_tree, bool simul, bool pee_, bool bdt_fit, string ch_s, int sig): pdf_analysis(print, meth, ch_s, range, SM, bd_constr, simul, pee_, bdt_fit) {
  cout << "fitData constructor" << endl;
  channels = inputs;
  simul_ = simul;
  input_estimates_ = input_estimates;
  if (simul_) channels = inputs;
  estimate_channel.resize(channels);
  estimate_bs.resize(channels);
  estimate_bd.resize(channels);
  estimate_rare.resize(channels);
  estimate_comb.resize(channels);

  pee = pee_;
  simul_ = simul;
  bdt_fit_ = bdt_fit;

  parse_estimate();
  if (input_tree) {
    tree = input_tree;
    random = false;
  }
  else {
    cout << "no input tree, making random distribution" << endl;
    random = true;
  }
  RooRandom::randomGenerator()->SetSeed(0);

  ws_file_input.resize(channels);
  ws_input.resize(channels);

  sign = sig;

  eff_bd.resize(channels);
  eff_bs.resize(channels);
  eff_bu.resize(channels);
  N_bu.resize(channels);
  BF_bd.resize(channels);
  BF_bs.resize(channels);
}

pdf_fitData::~pdf_fitData() {
  cout << "pdf_fitData destructor" << endl;
}

void pdf_fitData::parse_estimate(){
  char buffer[1024];
  char cutName[128];
  float cut;
  FILE *estimate_file = fopen(input_estimates_.c_str(), "r");
  cout << "event estimates in " << input_estimates_ << " :" << endl;
  while (fgets(buffer, sizeof(buffer), estimate_file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '#') continue;
    sscanf(buffer, "%s %f", cutName, &cut);
    if (!parse(cutName, cut)) {
      cout << "==> Error parsing variable " << cutName << endl;
      exit(EXIT_FAILURE);
    }
  }
  if (estimate_file) fclose(estimate_file);
  return;
}

bool pdf_fitData::parse(char *cutName, float cut) {
  if (simul_) {
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
      sprintf(test_cut, "rare_%d", i);
      if (!strcmp(cutName, test_cut)) {
        estimate_rare[i] = (int)cut;
        cout << "rare[" << i <<"]: " << estimate_rare[i] << endl;
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
  else {
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
    sprintf(test_cut, "rare_%d", i);
    if (!strcmp(cutName, test_cut)) {
      estimate_rare[0] = (int)cut;
      cout << "rare[" << i <<"]: " << estimate_rare[0] << endl;
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
  return true;
}

void pdf_fitData::fit_pdf(bool do_not_import) {
  cout << "making fit" << endl;
  if (simul_) {
    if (verbosity > 0) cout << "fitting " << global_data->GetName() << " in range " << range_ << " with pdf_ext_simul:" << endl;
    ws_->pdf("pdf_ext_simul")->Print();

    if (!pee) RFR = ws_->pdf("pdf_ext_simul")->fitTo(*global_data, Extended(1), Save(1), Minos());
    else RFR = ws_->pdf("pdf_ext_simul")->fitTo(*global_data, Extended(1), Save(1), Minos(), ConditionalObservables(*ws_->var("MassRes")));
  }
  else {
    RooAbsData* subdata = global_data->reduce(Form("channels==channels::channel_%d", channel));
    global_data = (RooDataSet*)subdata;
    if (verbosity > 0) cout << "fitting " << global_data->GetName() << " in range " << range_ << " with pdf_ext_total:" << endl;
    ws_->pdf("pdf_ext_total")->Print();

    if (!pee) RFR = ws_->pdf("pdf_ext_total")->fitTo(*global_data, Extended(1), Save(1), Minos());
    else RFR = ws_->pdf("pdf_ext_total")->fitTo(*global_data, Extended(1), Save(1), Minos(), ConditionalObservables(*ws_->var("MassRes")));
  }
  if (!do_not_import) ws_->import(*global_data);
  if (verbosity > 0) RFR->Print();
}

void pdf_fitData::print() {
  cout << "printing" << endl;
  RooPlot *rp = ws_->var("Mass")->frame();
  //RooAbsData* subdata_res = global_data->reduce(Form("channels==channels::channel_%d", channel));
  global_data->plotOn(rp, Binning(40));

  if (!pee) {
    ws_->pdf("pdf_ext_total")->plotOn(rp, FillColor(kYellow), Range(range_.c_str()), LineWidth(3), VisualizeError(*RFR), MoveToBack());
    ws_->pdf("pdf_ext_total")->plotOn(rp, LineColor(kBlue), Range(range_.c_str()), LineWidth(3));
  }
  else {
//    RooAbsData* subdata_res = ws_->pdf()
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
        found = name.find("pdf_rare");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kGreen - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
      }
      else {
        size_t found;
        found = name.find("pdf_bs");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kRed),          LineStyle(1), DrawOption("F"), FillColor(kRed), FillStyle(3001), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_bd");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kViolet - 4),   LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_comb");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kBlue - 5),   LineStyle(2), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
        found = name.find("pdf_rare");
        if (found!=string::npos) ws_->pdf("pdf_ext_total")->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), ProjWData(RooArgSet(*ws_->var("MassRes")), *global_data, kFALSE), LineColor(kGreen - 7), LineStyle(1), LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()));
      }
    }
  }
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  string address = "fig/data_fitData_" + (string)ws_->pdf("pdf_ext_total")->GetName() + "_" + meth_ + "_" + ch_s_;
  if (SM_) address += "_SM";
  if (pee)  address += "_PEE";
  c->Print( (address + ".gif").c_str());
  c->Print( (address + ".pdf").c_str());
  delete rp;
  delete c;
}

void pdf_fitData::print_each_channel() {
  cout <<"printing"<< endl;
  for (int i = 0; i < channels; i++) {

    RooPlot* final_p = ws_->var("Mass")->frame(Bins(40), Title(Form("Candidate invariant mass for channel %d", i)));
    global_data->plotOn(final_p, Cut(Form("channels==channels::channel_%d", i)));
    if (!pee) {
      ws_->pdf("pdf_ext_simul")->plotOn(final_p, VisualizeError(*RFR, 1), FillColor(kYellow), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(*ws_->cat("channels"), *global_data), Range(range_.c_str()), MoveToBack()) ;
      ws_->pdf("pdf_ext_simul")->plotOn(final_p, Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(*ws_->cat("channels"), *global_data), LineColor(kBlue), Range(range_.c_str()), LineWidth(3));
    }
    else {
      ws_->pdf("pdf_ext_simul")->plotOn(final_p, Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(RooArgSet(*ws_->cat("channels"), *ws_->var("MassRes")), *global_data, kFALSE), LineColor(kBlue), Range(range_.c_str()), LineWidth(3));
    }
    //ws_->pdf("pdf_ext_simul")->paramOn(final_p, Layout(0.30, 0.95, 0.95), Format("NEAU")/*, Parameters(*param)*/);

    // components
    RooArgSet* set = ws_->pdf("pdf_ext_simul")->getComponents();
    TIterator* it = set->createIterator();
    TObject* var_Obj = 0;
    while((var_Obj = it->Next())){
      string name = var_Obj->GetName();
      if (name != ws_->pdf("pdf_ext_simul")->GetName()) {
        if (!pee) {
          size_t found;
          found = name.find(Form("pdf_bs_%d", i));
          if (found != string::npos) ws_->pdf("pdf_ext_simul")->plotOn(final_p, Components(name.c_str()), LineColor(kRed),        LineStyle(1), DrawOption("F"), FillColor(kRed),        FillStyle(3001)/*, LineWidth(3)*/, Range(range_.c_str()), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(*ws_->cat("channels"), *global_data));
          found = name.find(Form("pdf_bd_%d", i));
          if (found != string::npos) ws_->pdf("pdf_ext_simul")->plotOn(final_p, Components(name.c_str()), LineColor(kViolet - 4), LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144)/*, LineWidth(3)*/, Range(range_.c_str()), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(*ws_->cat("channels"), *global_data));
          found = name.find(Form("pdf_comb_%d", i));
          if (found != string::npos) ws_->pdf("pdf_ext_simul")->plotOn(final_p, Components(name.c_str()), LineColor(kBlue - 5),   LineStyle(2)/*, DrawOption("F"), FillColor(kBlue - 5), FillStyle(3001)*/, LineWidth(3), Range(range_.c_str()), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(*ws_->cat("channels"), *global_data));
          found = name.find(Form("pdf_rare_%d", i));
          if (found != string::npos) ws_->pdf("pdf_ext_simul")->plotOn(final_p, Components(name.c_str()), LineColor(kGreen - 7),  LineStyle(1)/*, DrawOption("F"), FillColor(kGreen - 7), FillStyle(3001)*/, LineWidth(2), Range(range_.c_str()), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(*ws_->cat("channels"), *global_data));
        }
        else {
          size_t found;
          found = name.find(Form("pdf_bs_%d", i));
          if (found != string::npos) ws_->pdf("pdf_ext_simul")->plotOn(final_p, Components(name.c_str()), LineColor(kRed),        LineStyle(1), DrawOption("F"), FillColor(kRed),        FillStyle(3001)/*, LineWidth(3)*/, Range(range_.c_str()), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(RooArgSet(*ws_->cat("channels"), *ws_->var("MassRes")), *global_data, kFALSE));
          found = name.find(Form("pdf_bd_%d", i));
          if (found != string::npos) ws_->pdf("pdf_ext_simul")->plotOn(final_p, Components(name.c_str()), LineColor(kViolet - 4), LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144)/*, LineWidth(3)*/, Range(range_.c_str()), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(RooArgSet(*ws_->cat("channels"), *ws_->var("MassRes")), *global_data, kFALSE));
          found = name.find(Form("pdf_comb_%d", i));
          if (found != string::npos) ws_->pdf("pdf_ext_simul")->plotOn(final_p, Components(name.c_str()), LineColor(kBlue - 5),   LineStyle(2)/*, DrawOption("F"), FillColor(kBlue - 5), FillStyle(3001)*/, LineWidth(3), Range(range_.c_str()), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(RooArgSet(*ws_->cat("channels"), *ws_->var("MassRes")), *global_data, kFALSE));
          found = name.find(Form("pdf_rare_%d", i));
          if (found != string::npos) ws_->pdf("pdf_ext_simul")->plotOn(final_p, Components(name.c_str()), LineColor(kGreen - 7),  LineStyle(1)/*, DrawOption("F"), FillColor(kGreen - 7), FillStyle(3001)*/, LineWidth(2), Range(range_.c_str()), Slice(*ws_->cat("channels"), Form("channel_%d", i)), ProjWData(RooArgSet(*ws_->cat("channels"), *ws_->var("MassRes")), *global_data, kFALSE));
        }
      }
    }
    //global_data->plotOn(final_p, Cut( Form("channels==channels::channel_%d", i)));

    TCanvas* final_c = new TCanvas("final_c", "final_c", 600, 600);
    final_p->Draw();

    RooArgSet* vars =  ws_->pdf("pdf_ext_simul")->getVariables();
    RooRealVar* N_bs = (RooRealVar*)vars->find(Form("N_bs_%d", i));
    RooRealVar* N_bd;
    if (!bd_constr_ && !SM_) {
      N_bd = (RooRealVar*)vars->find(Form("N_bd_%d", i));
    }
    else if (bd_constr_) N_bd = (RooRealVar*)vars->find("Bd_over_Bs");
    else N_bd = 0;
    RooRealVar* N_comb = (RooRealVar*)vars->find(Form("N_comb_%d", i));
    RooRealVar* N_rare = (RooRealVar*)vars->find(Form("N_rare_%d", i));
    ostringstream fitresult_tex[4];
    fitresult_tex[0] << setprecision(2) << fixed << "N(B_{s}) = " << N_bs->getVal() << " ^{+" << getErrorHigh(N_bs) << "}_{" << getErrorLow(N_bs) << "}";
    if (!bd_constr_ && !SM_) fitresult_tex[1] <<setprecision(2) << fixed << "N(B_{d}) = " << N_bd->getVal() << " ^{+" << getErrorHigh(N_bd) << "}_{" << getErrorLow(N_bd) << "}";
    else if (bd_constr_) fitresult_tex[1] << setprecision(2) << fixed << "N(B_{d}) / N(B_{s}) = " << N_bd->getVal() << " ^{+" << getErrorHigh(N_bd) << "}_{" << getErrorLow(N_bd) << "}";
    else fitresult_tex[1] << setprecision(2) << fixed << "";
    fitresult_tex[2] << setprecision(2) << fixed<< "N(comb. bkg) = " << N_comb->getVal() << " ^{+" << getErrorHigh(N_comb) << "}_{" << getErrorLow(N_comb) << "}";
    fitresult_tex[3] << setprecision(2) << fixed<< "N(rare bkg) = " << N_rare->getVal() << " ^{+" << getErrorHigh(N_rare) << "}_{" << getErrorLow(N_rare) << "}";

    TPaveText* fitresults = new TPaveText(0.57, 0.66, 0.9, 0.9, "NDCR");
    for (int j = 0; j < 4; j++) {
      if (SM_ && j == 1) continue;
      TText* text =  fitresults->AddText(fitresult_tex[j].str().c_str());
    }
    fitresults->SetFillColor(0);
    fitresults->SetShadowColor(0);
//    fitresults->SetLineColor(0);
    fitresults->SetTextSize(0.03);
    fitresults->SetTextAlign(11);
    fitresults->Draw();

    ostringstream output;
    output << "fig/data_simul_" << meth_ << "_" << i;
    if (SM_) output << "_SM";
    if (bd_constr_) output << "_bdConst";
    if (pee)  output << "_PEE";
    final_c->Print( (output.str() + ".gif").c_str() );
    final_c->Print( (output.str() + ".pdf").c_str() );
    delete final_p;
    delete final_c;

  }
}

void pdf_fitData::FillRooDataSet(RooDataSet* dataset, bool cut_b, vector <double> cut_, TF1 *MassRes_f, string cuts) {
  int events = 0;
  if (!strcmp(tree->GetName(), "bdt") || !strcmp(tree->GetName(), "cnc")) {
    TTree* reduced_tree = tree->CopyTree(cuts.c_str());
    Double_t m1eta_t, m2eta_t, m_t, eta_B_t, bdt_t;
    reduced_tree->SetBranchAddress("m1eta", &m1eta_t);
    reduced_tree->SetBranchAddress("m2eta", &m2eta_t);
    reduced_tree->SetBranchAddress("m", &m_t);
    reduced_tree->SetBranchAddress("eta", &eta_B_t);
    reduced_tree->SetBranchAddress("bdt", &bdt_t);
    for (int i = 0; i < tree->GetEntries(); i++){
      reduced_tree->GetEntry(i);
      if (m_t > 4.9 && m_t < 5.9) {
        events++;
        Mass->setVal(m_t);
        eta->setVal(eta_B_t);
        m1eta->setVal(m1eta_t);
        m2eta->setVal(m2eta_t);
        bdt->setVal(bdt_t);
        MassRes->setVal(MassRes_f->Eval(eta_B_t));
        if (fabs(m1eta_t)<1.4 && fabs(m2eta_t)<1.4) {
          if (cut_b && bdt_t < cut_[0]) continue;
          channels_cat->setIndex(0);
        }
        else {
          if (cut_b && bdt_t < cut_[1]) continue;
          channels_cat->setIndex(1);
        }
        if (bdt_t < 0.1) bdt_cat->setIndex(0);
        else if (bdt_t < 0.18) bdt_cat->setIndex(1);
        else bdt_cat->setIndex(2);
        RooArgSet varlist_tmp(*Mass, *MassRes, *eta, *m1eta, *m2eta, *bdt, *channels_cat, *bdt_cat);
        dataset->add(varlist_tmp);
      }
    }
  }
  else {
    cout << "tree name is not bdt nor cnc" << endl;
    exit(1);
  }
  cout << "total events = " << events << endl;
}

void pdf_fitData::define_channels() {
  cout << "making channels" << endl;
  for (int i = 0; i < channels; i++) {
    ostringstream channelName;
    channelName << "channel_" << i;
    channels_cat->defineType(channelName.str().c_str(), i);
  }
}

void pdf_fitData::make_dataset(bool cut_b, vector <double> cut_, TF1 *MassRes_f, string cuts) {
  cout << "making dataset" << endl;
  RooArgList varlist(*Mass, *MassRes, *eta, *m1eta, *m2eta, *bdt);
  varlist.add(*channels_cat);
  global_data = new RooDataSet("global_data", "global_data", varlist);
  if (!random) FillRooDataSet(global_data, cut_b, cut_, MassRes_f, cuts);
  else {
    //RooRandom::randomGenerator()->SetSeed(12345);
    if (!simul_) {
      ws_->var("N_bs")->setVal(estimate_bs[ch_i_]);
      if (!SM_ && !bd_constr_) ws_->var("N_bd")->setVal(estimate_bd[ch_i_]);
      else if (bd_constr_) ws_->var("bd_over_bs")->setVal(estimate_bd[ch_i_]/estimate_bd[ch_i_]);
      ws_->var("N_rare")->setVal(estimate_rare[ch_i_]);
      ws_->var("N_comb")->setVal(estimate_comb[ch_i_]);
      global_data = ws_->pdf("pdf_ext_total")->generate(RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->var("bdt"), *ws_->cat("channels")));
    }
    else {
      for (int i = 0; i < channels; i++) {
        ws_->var(name("N_bs", i))->setVal(estimate_bs[i]);
        if (!SM_ && !bd_constr_) ws_->var(name("N_bd", i))->setVal(estimate_bd[i]);
        else if (bd_constr_) ws_->var("bd_over_bs")->setVal(estimate_bd[i]/estimate_bd[i]);
        ws_->var(name("N_rare", i))->setVal(estimate_rare[i]);
        ws_->var(name("N_comb", i))->setVal(estimate_comb[i]);
      }
      global_data = ws_->pdf("pdf_ext_simul")->generate(RooArgSet(*ws_->var("Mass"), *ws_->var("MassRes"), *ws_->var("bdt"), *ws_->cat("channels")));
    }
  }
  global_data->SetName("global_data");
  //ws_->import(*global_data);
}

void pdf_fitData::changeName(RooWorkspace *ws, int str) {
  ostringstream chan;
  chan << str;
  RooArgSet all_vars(ws->allVars());
  TIterator* vars_it = all_vars.createIterator();
  RooRealVar *arg_var = 0;
  while ( (arg_var = (RooRealVar*)vars_it->Next())) {
    if (  (strcmp( arg_var->GetName(), "Mass")) && strcmp( arg_var->GetName(), "Bd_over_Bs")) {
      arg_var->SetName( Form("%s_%s", arg_var->GetName(), chan.str().c_str()));
    }
  }
  delete vars_it;
  RooArgSet all_pdf(ws->allPdfs());
  TIterator* pdf_it = all_pdf.createIterator();
  RooAbsPdf *pdf_arg = 0;
  while ( (pdf_arg = (RooAbsPdf*)pdf_it->Next())) {
    pdf_arg->SetName( Form("%s_%s", pdf_arg->GetName(), chan.str().c_str()));
  }
  delete pdf_it;
  RooArgSet all_fun(ws->allFunctions());
  TIterator* fun_it = all_fun.createIterator();
  RooFormulaVar *fun_arg = 0;
  while ( (fun_arg = (RooFormulaVar*)fun_it->Next())) {
    if ( !(strcmp( fun_arg->GetName(), "N_bd_constr"))) {
      fun_arg->SetName( Form("%s_%s", fun_arg->GetName(), chan.str().c_str()));
    }
  }
  delete fun_it;
  return;
}

void pdf_fitData::make_pdf_input() {
  cout << "making inputs of pdf" << endl;
  string root_s = "output/ws_pdf_" + meth_ + "_";
  if (simul_) root_s = "output/ws_simul_" + meth_;
  string tail_s;
  if (SM_) tail_s = "_SM";
  else if (bd_constr_) tail_s = "_BdConst";
  if (bdt_fit_) tail_s += "_2D";
  if (pee) tail_s += "_PEE";
  tail_s += ".root";
  if (simul_) {
    root_s += tail_s;
    ws_file_input[0] = new TFile(root_s.c_str());
    if (!ws_file_input[0]) {cout << root_s.c_str() << " does not exist" << endl; exit(EXIT_FAILURE);}
    ws_input[0] = (RooWorkspace*)ws_file_input[0]->Get("ws");
    if (!ws_input[0]) {cout << "ws does not exist" << endl; exit(EXIT_FAILURE);}
    cout << "ws file: " << root_s << endl;
  }
  else {
    ostringstream input_oss;
    input_oss << root_s << ch_s_ << tail_s;
    ws_file_input[0] = new TFile(input_oss.str().c_str());
    if (!ws_file_input[0]) {cout << input_oss.str().c_str() << " does not exist" << endl; exit(EXIT_FAILURE);}
    ws_input[0] = (RooWorkspace*)ws_file_input[0]->Get("ws");
    if (!ws_input[0]) {cout << "ws does not exist" << endl; exit(EXIT_FAILURE);}
    cout << "ws file: " << input_oss.str() << endl;
  }

}

void pdf_fitData::make_pdf() {
  cout << "making pdf" << endl;
  if (random) {
    if (simul_) {
      for (int i = 0; i < channels; i++) {
        ws_input[0]->var(name("N_bs", i))->setVal(estimate_bs[i]);
        ws_input[0]->var(name("N_bd", i))->setVal(estimate_bd[i]);
        ws_input[0]->var(name("N_rare", i))->setVal(estimate_rare[i]);
        ws_input[0]->var(name("N_comb", i))->setVal(estimate_comb[i]);
      }
    }
    else {
      int i = atoi(ch_s_.c_str());
      ws_input[0]->var("N_bs")->setVal(estimate_bs[i]);
      ws_input[0]->var("N_bd")->setVal(estimate_bd[i]);
      ws_input[0]->var("N_rare")->setVal(estimate_rare[i]);
      ws_input[0]->var("N_comb")->setVal(estimate_comb[i]);
    }
    if (bd_constr_) ws_input[0]->var("Bd_over_Bs")->setVal(estimate_bd[0] / estimate_bs[0]);
  }
  ws_ = ws_input[0];
  ws_->Print();
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
  if (sign == 0) sig_hand();
  else if (sign == 1) sig_plhc();
  else if (sign >= 1) {
    make_models();
    if (sign == 2) sig_plhts();
    if (sign == 3) sig_hybrid_plhts();
    if (sign == 4) sig_hybrid_roplhts();
  }
}

Double_t pdf_fitData::sig_hand() {
  /// by hand
  Double_t minNLL = RFR->minNll();
  RooRealVar *arg_var = 0;
  RooArgSet *all_vars;
  if (simul_) all_vars = ws_->pdf("pdf_ext_simul")->getVariables();
  else all_vars = ws_->pdf("pdf_ext_total")->getVariables();
  TIterator* vars_it = all_vars->createIterator();
  size_t found;
  while ( (arg_var = (RooRealVar*)vars_it->Next())) {
    string name(arg_var->GetName());
    found = name.find("N_bs");
    if (found!=string::npos) {
      arg_var->setVal(0.0);
      arg_var->setConstant(1);
    }
//    found = name.find("N_bd");
//    if (found!=string::npos) {
//      arg_var->setVal(0.0);
//      arg_var->setConstant(1);
//    }
    if (SM_ || bd_constr_){
      found = name.find("Bd_over_Bs");
      if (found!=string::npos) {
        arg_var->setVal(0.0);
        arg_var->setConstant(1);
      }
    }
  }
  fit_pdf(true);

  TIterator* vars_after = all_vars->createIterator();
  while ( (arg_var = (RooRealVar*)vars_after->Next())) {
    string name(arg_var->GetName());
    found = name.find("N_bs");
    if (found!=string::npos) arg_var->setConstant(0);
    found = name.find("N_bd");
    if (found!=string::npos) arg_var->setConstant(0);
    if (SM_ || bd_constr_) {
      found = name.find("Bd_over_Bs");
      if (found!=string::npos) arg_var->setConstant(0);
    }
  }

  Double_t newNLL = RFR->minNll();
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
  if (simul_) {
    model.SetPdf(*ws_->pdf("pdf_ext_simul"));
    for (int i = 0; i < channels; i++) {
      poi.add(*ws_->var(name("N_bs", i)));
      poi.setRealValue(name("N_bs", i), 0);
//      if (!bd_constr_ && !SM_) {
//        poi.add(*ws_->var(name("N_bd", i)));
//        poi.setRealValue(name("N_bd", i), 0);
//      }
    }
  }
  else {
    model.SetPdf(*ws_->pdf("pdf_ext_total"));
    poi.add(*ws_->var("N_bs"));
    poi.setRealValue("N_bs", 0);
//    if (!bd_constr_ && !SM_) {
//      poi.add(*ws_->var("N_bd"));
//      poi.setRealValue("N_bd", 0);
//    }
  }
  if (bd_constr_) {
    poi.add(*ws_->var("Bd_over_Bs"));
    poi.setRealValue("Bd_over_Bs", 0);
  }

  ProfileLikelihoodCalculator plc;
  plc.SetData(*ws_->data("global_data"));
  plc.SetModel(model);
  plc.SetConditionalObservables(CO);
  plc.SetNullParameters(poi);
  HypoTestResult* htr = plc.GetHypoTest();
  cout << "ProfileLikelihoodCalculator: The p-value for the null is " << htr->NullPValue() << "; The significance for the null is " << htr->Significance() << endl;
  model.SetSnapshot(poi);
  ws_->import(model);
}

void pdf_fitData::make_models() {
  vector <double> N_bs(channels);
  vector <double> N_bd(channels);
  for (int i = 0; i < channels; i++) {
    N_bs[i] = ws_->var(name("N_bs", i))->getVal();
    if (!SM_ && !bd_constr_) N_bd[i] = ws_->var(name("N_bd", i))->getVal();
  }

  if (simul_) ws_->defineSet("obs", "Mass,channels");
  else ws_->defineSet("obs", "Mass");
  ostringstream name_poi;
  if (simul_) {
    for (int i = 0; i < channels; i++) {
      if (i != 0) name_poi << ",";
      name_poi << "N_bs_" << i;
//      if (!bd_constr_ && !SM_) {
//        name_poi << ",N_bd_" << i;
//      }
    }
  }
  else {
    name_poi << "N_bs";
//    if (!bd_constr_ && !SM_) name_poi << ",N_bd";
  }
  if (bd_constr_) name_poi << ",Bd_over_Bs";
  ws_->defineSet("poi", name_poi.str().c_str());

  ModelConfig* H0 = new ModelConfig("H0", "null hypothesis", ws_);
  RooArgSet CO;
  if (pee) {
    CO.add(*ws_->var("MassRes"));
    H0->SetConditionalObservables(CO);
    ws_->defineSet("CO", "MassRes");
  }
  if (simul_) H0->SetPdf(*ws_->pdf("pdf_ext_simul"));
  else H0->SetPdf(*ws_->pdf("pdf_ext_total"));
  H0->SetParametersOfInterest(*ws_->set("poi"));
  H0->SetObservables(*ws_->set("obs"));
  if (simul_) {
    for (int i = 0; i < channels; i++) {
      ostringstream name_oss;
      name_oss << "N_bs_" << i;
      ws_->var(name_oss.str().c_str())->setVal(0.0);
//      if (!bd_constr_ && !SM_) {
//        ws_->var(name("N_bd", i))->setVal(0.0);
//      }
    }
  }
  else {
    ws_->var("N_bs")->setVal(0.0);
//    if (!bd_constr_ && !SM_) ws_->var("N_bd")->setVal(0.0);
  }
  if (bd_constr_) {
    ws_->var("Bd_over_Bs")->setVal(0.0);
  }
  H0->SetSnapshot(*ws_->set("poi"));

  ModelConfig* H1 = new ModelConfig("H1", "background + signal hypothesis", ws_);
  if (pee) {
    H1->SetConditionalObservables(CO);
  }
  if (simul_) H1->SetPdf(*ws_->pdf("pdf_ext_simul"));
  else H1->SetPdf(*ws_->pdf("pdf_ext_total"));
  H1->SetParametersOfInterest(*ws_->set("poi"));
  H1->SetObservables(*ws_->set("obs"));
  parse_estimate();
  if (simul_) {
    for (int i = 0; i < channels; i++) {
      ws_->var(name("N_bs", i))->setVal(N_bs[i]);
//      if (!bd_constr_ && !SM_) {
//        ws_->var(name("N_bd", i))->setVal(N_bd[i]);
//      }
    }
  }
  else {
    ws_->var("N_bs")->setVal(N_bs[0]);
//    ws_->var("N_bd")->setVal(N_bd[0]);
  }
  if (bd_constr_) {
    int index = simul_ ? 0 : atoi(ch_s_.c_str());
    double ratio = (double) estimate_bd[index] / estimate_bs[index];
    ws_->var("Bd_over_Bs")->setVal(ratio);
  }
  H1->SetSnapshot(*ws_->set("poi"));

  ws_->import(*H0);
  ws_->import(*H1);
}

void pdf_fitData::sig_plhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  ws_->Print();

  string name_of_pdf = "pdf_ext_simul";
  if (!simul_) name_of_pdf = "pdf_ext_total";
  ProfileLikelihoodTestStat pl_ts(*ws_->pdf(name_of_pdf.c_str()));
  pl_ts.SetOneSidedDiscovery(true);
  if (pee) pl_ts.SetConditionalObservables(*ws_->set("CO"));
  ToyMCSampler *mcSampler_pl = new ToyMCSampler(pl_ts, 500);
  FrequentistCalculator frequCalc(*ws_->data("global_data"), *H1,*H0, mcSampler_pl); // null = bModel interpreted as signal, alt = s+b interpreted as bkg
  HypoTestResult *htr_pl = frequCalc.GetHypoTest();
  htr_pl->Print();
//  HypoTestPlot *plot = new HypoTestPlot(*htr_pl);
//  TCanvas* c_hypotest = new TCanvas("c_hypotest", "c_hypotest", 600, 600);
//  plot->Draw();
//  c_hypotest->Print("fig/ProfileLikelihoodTestStat.gif");
//  delete plot;
//  delete c_hypotest;
  cout << "ProfileLikelihoodTestStat + frequentist: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << endl;
}

void pdf_fitData::sig_hybrid_plhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  make_prior();
  ws_->Print();

  string name_of_pdf = "pdf_ext_simul";
  if (!simul_) name_of_pdf = "pdf_ext_total";
  ProfileLikelihoodTestStat pl_ts(*ws_->pdf(name_of_pdf.c_str()));
  pl_ts.SetOneSidedDiscovery(true);
  if (pee) pl_ts.SetConditionalObservables(*ws_->set("CO"));
  ToyMCSampler *mcSampler_pl = new ToyMCSampler(pl_ts, 500);

  HybridCalculator hibrCalc(*ws_->data("global_data"), *H1, *H0, mcSampler_pl);
  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("prior"));
  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("prior"));

  HypoTestResult *htr_pl = hibrCalc.GetHypoTest();
  htr_pl->Print();
  cout << "ProfileLikelihoodTestStat + hybrid: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << endl;
}

void pdf_fitData::sig_hybrid_roplhts() {

  RooStats::ModelConfig *H0 = dynamic_cast<ModelConfig*> (ws_->obj("H0"));
  RooStats::ModelConfig *H1 = dynamic_cast<ModelConfig*> (ws_->obj("H1"));

  make_prior();
  ws_->Print();

  string name_of_pdf = "pdf_ext_simul";
  if (!simul_) name_of_pdf = "pdf_ext_total";
  RatioOfProfiledLikelihoodsTestStat ropl_ts(*H1->GetPdf(),*H0->GetPdf(), ws_->set("poi"));
  ropl_ts.SetSubtractMLE(false);
  if (pee) ropl_ts.SetConditionalObservables(*ws_->set("CO"));
  ToyMCSampler *mcSampler_pl = new ToyMCSampler(ropl_ts, 500);

  HybridCalculator hibrCalc(*ws_->data("global_data"), *H1, *H0, mcSampler_pl);
  hibrCalc.ForcePriorNuisanceAlt(*ws_->pdf("prior"));
  hibrCalc.ForcePriorNuisanceNull(*ws_->pdf("prior"));

  HypoTestResult *htr_pl = hibrCalc.GetHypoTest();
  htr_pl->Print();
  cout << "RatioOfProfiledLikelihoodsTestStat + hybrid: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << endl;
}

void pdf_fitData::make_prior() {
  vector <RooGaussian*> prior_bd(channels);
  vector <RooGaussian*> prior_rare(channels);
  vector <RooGaussian*> prior_comb(channels);

//  vector <RooGamma*> prior_bd1(channels)

  RooArgList prior_list("prior_list");

  for (int i = 0; i < channels; i++) {
    prior_bd[i] = new RooGaussian(name("prior_bd", i), name("prior_bd", i), *ws_->var(name("N_bd", i)), RooConst(ws_->var(name("N_bd", i))->getVal()), RooConst(ws_->var(name("N_bd", i))->getError()));
    prior_rare[i] = new RooGaussian(name("prior_rare", i), name("prior_rare", i), *ws_->var(name("N_rare", i)), RooConst(ws_->var(name("N_rare", i))->getVal()), RooConst(ws_->var(name("N_rare", i))->getError()));
    prior_comb[i] = new RooGaussian(name("prior_comb", i), name("prior_comb", i), *ws_->var(name("N_comb", i)), RooConst(ws_->var(name("N_comb", i))->getVal()), RooConst(ws_->var(name("N_comb", i))->getError()));
    //prior_bd1[i] = new RooGamma(name("prior_bd", i), name("prior_bd", i), *ws_->var(name("N_bd", i)), RooConst(1), )

    prior_list.add(*prior_bd[i]);
    prior_list.add(*prior_rare[i]);
    prior_list.add(*prior_comb[i]);
  }

  RooProdPdf prior("prior", "prior", prior_list);
  ws_->import(prior);
}

void pdf_fitData::BF(string eff_filename, string numbers_filename) {
  parse_external_numbers(numbers_filename);
  parse_efficiency_numbers(eff_filename);

  /// well, all this assuming uncorrelated varables...
  cout << "=========== Branching Fractions ==========" << endl;
  int ii = -1;
  for (int i = 0; i < channels; i++) {
    BF_bs.at(i).first = ws_->var(name("N_bs", i))->getVal() / N_bu.at(i).first / fs_over_fu.first * (eff_bu.at(i).first / eff_bs.at(i).first) * Bu2JpsiK_BF.first * Jpsi2MuMu_BF.first;
    BF_bs.at(i).second = BF_bs.at(i).first * sqrt(pow(ws_->var(name("N_bs", i))->getError() / ws_->var(name("N_bs", i))->getVal(), 2) +
                                                  pow(N_bu.at(i).second / N_bu.at(i).first, 2) +
                                                  pow(fs_over_fu.second / fs_over_fu.first, 2) +
                                                  pow(eff_bu.at(i).second / eff_bu.at(i).first, 2) +
                                                  pow(eff_bs.at(i).second / eff_bs.at(i).first, 2) +
                                                  pow(Bu2JpsiK_BF.second / Bu2JpsiK_BF.first, 2) +
                                                  pow(Jpsi2MuMu_BF.second / Jpsi2MuMu_BF.first, 2));

    BF_bd.at(i).first = ws_->var(name("N_bd", i))->getVal() / N_bu.at(i).first * (eff_bu.at(i).first / eff_bd.at(i).first) * Bu2JpsiK_BF.first * Jpsi2MuMu_BF.first;
    BF_bd.at(i).second = BF_bd.at(i).first * sqrt(pow(ws_->var(name("N_bd", i))->getError() / ws_->var(name("N_bd", i))->getVal(), 2) +
                                                  pow(N_bu.at(i).second / N_bu.at(i).first, 2) +
                                                  pow(eff_bu.at(i).second / eff_bu.at(i).first, 2) +
                                                  pow(eff_bd.at(i).second / eff_bd.at(i).first, 2) +
                                                  pow(Bu2JpsiK_BF.second / Bu2JpsiK_BF.first, 2) +
                                                  pow(Jpsi2MuMu_BF.second / Jpsi2MuMu_BF.first, 2));

    if (simul_) ii = i;
    else ii = ch_i_;
    cout << "channel " << ii << ":" << endl;
    cout << "Bs2MuMu BF = " << BF_bs.at(i).first << " \\pm " << BF_bs.at(i).second << "  channel " << ii <<  endl;
    cout << "Bd2MuMu BF = " << BF_bd.at(i).first << " \\pm " << BF_bd.at(i).second << "  channel " << ii << endl;
  }

  double bs_num = 0, bd_num = 0;
  double bs_den = 0, bd_den = 0;
  for (int i = 0; i < channels; i++) {
    bs_num += BF_bs.at(i).first / pow(BF_bs.at(i).second, 2);
    bs_den += 1. / pow(BF_bs.at(i).second, 2);
    bd_num += BF_bd.at(i).first / pow(BF_bd.at(i).second, 2);
    bd_den += 1. / pow(BF_bd.at(i).second, 2);
  }
  double bs_final_bf = bs_num / bs_den;
  double bs_error_bf = sqrt(1. / bs_den);
  double bd_final_bf = bd_num / bd_den;
  double bd_error_bf = sqrt(1. / bs_den);

  cout << "============= final numbers =============" << endl;
  cout << "Bs2MuMu BF = " << bs_final_bf << " \\pm " << bs_error_bf << endl;
  cout << "Bd2MuMu BF = " << bd_final_bf << " \\pm " << bd_error_bf << endl;
  cout << "============= ============= =============" << endl;
}

void pdf_fitData::parse_efficiency_numbers(string filename) {
  cout << "parsing " << filename << endl;

  FILE *file = fopen(filename.c_str(), "r");
  if (!file) {cout << "file " << filename << " does not exist"; exit(1);}

  char buffer[1024];
  char left[1024];
  double number;

  vector < pair<string, string> > end_bd(channels);
  vector < pair<string, string> > end_bs(channels);
  vector < pair<string, string> > end_bu(channels);
  vector < pair<string, string> > end_Nbu(channels);
  int ii = -1;
  for (int i = 0; i < channels; i++) {
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
    sscanf(buffer, "%s   {\\ensuremath{{%lf } } }", left, &number);
    string left_s(left);
    for (int i = 0; i < channels; i++) {
      size_t found;
      found = left_s.find(end_bd.at(i).first);
      if (found != string::npos) eff_bd.at(i).first = number;
      found = left_s.find(end_bd.at(i).second);
      if (found != string::npos) eff_bd.at(i).second = number;
      found = left_s.find(end_bs.at(i).first);
      if (found != string::npos) eff_bs.at(i).first = number;
      found = left_s.find(end_bs.at(i).second);
      if (found != string::npos) eff_bs.at(i).second = number;
      found = left_s.find(end_bu.at(i).first);
      if (found != string::npos) eff_bu.at(i).first = number;
      found = left_s.find(end_bu.at(i).second);
      if (found != string::npos) eff_bu.at(i).second = number;
      found = left_s.find(end_Nbu.at(i).first);
      if (found != string::npos) N_bu.at(i).first = number;
      found = left_s.find(end_Nbu.at(i).second);
      if (found != string::npos) N_bu.at(i).second = number;
    }
  }

  for (int i = 0; i < channels; i++) {
    if (simul_) ii = i;
    else ii = ch_i_;
    cout << "channel " << ii << endl;
    cout << "bd eff = " << eff_bd.at(i).first << " \\pm " << eff_bd.at(i).second << endl;
    cout << "bs eff = " << eff_bs.at(i).first << " \\pm " << eff_bs.at(i).second << endl;
    cout << "bu eff = " << eff_bu.at(i).first << " \\pm " << eff_bu.at(i).second << endl;
    cout << "N bu   = " << N_bu.at(i).first << " \\pm " << N_bu.at(i).second << endl;
    cout << endl;
  }
}

void pdf_fitData::parse_external_numbers(string filename) {
  cout << "parsing " << filename << endl;
  FILE *file = fopen(filename.c_str(), "r");
  if (!file) {cout << "file " << filename << " does not exist"; exit(1);}
  char buffer[1024];
  while (fgets(buffer, sizeof(buffer), file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    if (buffer[0] == '\045') continue;
    sscanf(buffer, "fsfu\t%lf\t%lf", &fs_over_fu.first, &fs_over_fu.second);
    sscanf(buffer, "BuToKpsiK\t%lf\t%lf", &Bu2JpsiK_BF.first, &Bu2JpsiK_BF.second);
    sscanf(buffer, "JpsiToMuMu\t%lf\t%lf", &Jpsi2MuMu_BF.first, &Jpsi2MuMu_BF.second);
  }
  cout << "fs/fu        " <<  fs_over_fu.first << " \\pm " << fs_over_fu.second << endl;
  cout << "Bu2JpsiK_BF  " <<  Bu2JpsiK_BF.first << " \\pm " << Bu2JpsiK_BF.second << endl;
  cout << "Jpsi2MuMu_BF " <<  Jpsi2MuMu_BF.first << " \\pm " << Jpsi2MuMu_BF.second << endl;
}
