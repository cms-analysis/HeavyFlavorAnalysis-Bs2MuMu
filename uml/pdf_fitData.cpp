#include "pdf_fitData.h"

pdf_fitData::pdf_fitData(bool print, int inputs, string input_estimates, string meth, string range, bool SM, bool bd_constr, TTree *input_tree, bool simul, string ch_s): pdf_analysis(print, meth, ch_s, range, SM, bd_constr) {
  cout << "fitData constructor" << endl;
  channels = inputs;
  simul_ = simul;
  input_estimates_ = input_estimates;
  if (simul_) channels = inputs;
  else channels = 1;
  estimate_channel.resize(channels);
  estimate_bs.resize(channels);
  estimate_bd.resize(channels);
  estimate_rare.resize(channels);
  estimate_comb.resize(channels);
  parse_estimate();
  if (input_tree) {
    tree = input_tree;
    random = false;
  }
  else {
    cout << "no input tree, making random distribution" << endl;
    random = true;
  }
  channel = new RooCategory("channel", "channel categories");
  simul_pdf = new RooSimultaneous("simul_pdf", "simultaneous fit for all channels", *channel);

  ws_file_input.resize(channels);
  ws_input.resize(channels);

      ws_ = new RooWorkspace("ws", "simultaneous fit workspace results");

}

void pdf_fitData::parse_estimate(){
  char buffer[1024];
  char cutName[128];
  float cut;
  int ok;
  FILE *estimate_file = fopen(input_estimates_.c_str(), "r");
  cout << "event estimates in " << input_estimates_ << " :" << endl;
  int counter = 0;
  while (fgets(buffer, sizeof(buffer), estimate_file)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
    counter++;
    // skip comment line
    if (buffer[0] == '#') continue;
    ok = sscanf(buffer, "%s %f", cutName, &cut);
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

void pdf_fitData::fit_pdf() {
  cout << "making simultaneous fit" << endl;
  if (simul_) RFR = simul_pdf->fitTo(*global_data, Extended(1), Save(1)/*, Verbose(true), PrintLevel(3), PrintEvalErrors(20)*/);
  else RFR = total_pdf_i[0]->fitTo(*global_data, Extended(1), Save(1)/*, Verbose(true), PrintLevel(3), PrintEvalErrors(20)*/);
  RFR->Print();
}

void pdf_fitData::fit_pdf (string pdf, RooAbsData* data, bool extended, bool sumw2error, bool hesse) {
  rds_ = data;
  pdf_name = "pdf_" + pdf;
  if (extended) pdf_name = "pdf_ext_" + pdf;
  string rdh_name = rds_->GetName();
  cout << "fitting " << rdh_name << " in range " << range_ << " with " << pdf_name << ":" << endl;
  ws_->pdf( pdf_name.c_str())->Print();

  RFR = ws_->pdf( pdf_name.c_str())->fitTo(*rds_, Extended(extended), Save(1) ,SumW2Error(sumw2error)/*, Range(range_.c_str()), SumCoefRange(range_.c_str())*/, Hesse(hesse));
  if (print_) print();
  RFR->Print();
}

void pdf_fitData::print() {
  cout << "printing" << endl;
  RooPlot *rp = ws_->var("Mass")->frame();

  global_data->plotOn(rp, Binning(40));

//  ws_->pdf(pdf_name.c_str())->plotOn(rp, FillColor(kYellow), Range(range_.c_str()), LineWidth(3), VisualizeError(*RFR));
  total_pdf_i[0]->plotOn(rp, LineColor(kBlue), Range(range_.c_str()), LineWidth(3));
  global_data->plotOn(rp, Binning(40));
  total_pdf_i[0]->paramOn(rp, Layout(0.50, 0.9, 0.9));
  // components
  RooArgSet * set = total_pdf_i[0]->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (name == "pdf_bs") total_pdf_i[0]->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kRed),          LineStyle(1), DrawOption("F"), FillColor(kRed), FillStyle(3001), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
      if (name == "pdf_bd") total_pdf_i[0]->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kViolet - 4),   LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144), LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()));
      double Ncomb = ws_->var("N_comb")->getVal();
      double Nrare = ws_->var("N_rare")->getVal();
      if (name == "pdf_comb") total_pdf_i[0]->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kBlue - 5),   LineStyle(2)/*, DrawOption("F"), FillColor(kBlue - 5), FillStyle(3001)*/, LineWidth(3), Range(range_.c_str()), NormRange(range_.c_str()), RooFit::Normalization(Ncomb, 2));
      if (name == "pdf_rare") total_pdf_i[0]->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kGreen - 7), LineStyle(1)/*, DrawOption("F"), FillColor(kGreen - 7), FillStyle(3001)*/, LineWidth(2), Range(range_.c_str()), NormRange(range_.c_str()), RooFit::Normalization(Nrare, 2));
    }
  }
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  string address = "fig/data_fitData_" + pdf_name + "_" + meth_ + "_" + ch_s_;
  if (SM_) address += "_SM";
  c->Print( (address + ".gif").c_str());
  c->Print( (address + ".pdf").c_str());
  delete rp;
  delete c;
}

void pdf_fitData::print_each_channel() {

  for (int i = 0; i < channels; i++) {
    RooPlot* final_p = ws_->var("Mass")->frame(Bins(20));

    global_data->plotOn(final_p, Cut( Form("channel==channel::channel_%d", i)), Binning(40));

    simul_pdf->plotOn(final_p, Slice(*channel, Form("channel_%d", i)), ProjWData(*channel, *global_data), LineColor(kBlue), Range(range_.c_str()), LineWidth(3));
    simul_pdf->paramOn(final_p, Layout(0.50, 0.9, 0.9));

    // components
    RooArgSet* set = simul_pdf->getComponents();
    TIterator* it = set->createIterator();
    TObject* var_Obj = 0;
    while((var_Obj = it->Next())){
      string name = var_Obj->GetName();
      if (name != simul_pdf->GetName()) {
        size_t found;
        found = name.find("pdf_bs");
        if (found != string::npos) simul_pdf->plotOn(final_p, Components(var_Obj->GetName()), LineColor(kRed),        LineStyle(1), DrawOption("F"), FillColor(kRed),        FillStyle(3001), LineWidth(3), Range(range_.c_str()), Slice(*channel, Form("channel_%d", i)), ProjWData(*channel, *global_data));
        found = name.find("pdf_bd");
        if (found != string::npos) simul_pdf->plotOn(final_p, Components(var_Obj->GetName()), LineColor(kViolet - 4), LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144), LineWidth(3), Range(range_.c_str()), Slice(*channel, Form("channel_%d", i)), ProjWData(*channel, *global_data));
        found = name.find("pdf_comb");
        if (found != string::npos) simul_pdf->plotOn(final_p, Components(var_Obj->GetName()), LineColor(kBlue - 5),   LineStyle(2)/*, DrawOption("F"), FillColor(kBlue - 5), FillStyle(3001)*/, LineWidth(3), Range(range_.c_str()), Slice(*channel, Form("channel_%d", i)), ProjWData(*channel, *global_data));
        found = name.find("pdf_rare");
        if (found != string::npos) simul_pdf->plotOn(final_p, Components(var_Obj->GetName()), LineColor(kGreen - 7),  LineStyle(1)/*, DrawOption("F"), FillColor(kGreen - 7), FillStyle(3001)*/, LineWidth(2), Range(range_.c_str()), Slice(*channel, Form("channel_%d", i)), ProjWData(*channel, *global_data));
      }
    }

    TCanvas* final_c = new TCanvas("final_c", "final_c", 600, 600);
    final_p->Draw();
    ostringstream output;
    output << "fig/data_simul_" << meth_ << "_" << i;
    if (SM_) output << "_SM";
    if (bd_constr_) output << "_bdConst";
    final_c->Print( (output.str() + ".gif").c_str() );
    final_c->Print( (output.str() + ".pdf").c_str() );
    delete final_p;
    delete final_c;
  }
}

void pdf_fitData::FillRooDataSet(TTree* tree, RooDataSet* dataset, RooRealVar *Mass, int ch_i) {
  double m1eta, m2eta, m;
  tree->SetBranchAddress("m1eta", &m1eta);
  tree->SetBranchAddress("m2eta", &m2eta);
  tree->SetBranchAddress("m", &m);
  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    Mass->setVal(m);
    if (ch_i == -1) {
      dataset->add(*Mass);
    }
    else if (ch_i == 0) {
      if (fabs(m1eta) < 1.4 && fabs(m2eta) < 1.4) {
        dataset->add(*Mass);
      }
    }
    else if (ch_i == 1) {
      if (fabs(m1eta) > 1.4 || fabs(m2eta) > 1.4) {
        dataset->add(*Mass);
      }
    }
    else {
      cout << "wrong channel: " << ch_i << endl;
      exit(EXIT_FAILURE);
    }
  }
}

void pdf_fitData::make_dataset() {
  cout << "making datasets" << endl;
  if (simul_) {
    map<string, RooDataSet* > data_map;
    for (int i = 0; i < channels; i++) {
      ostringstream channelName;
      channelName << "channel_" << i;
      RooDataSet* data_i;
      if (!random) {
        data_i = new RooDataSet( Form("data_%i", i), Form("data %i", i), *ws_->var("Mass"));
        FillRooDataSet(tree, data_i, ws_->var("Mass"), i);
      }
      else {
        RooRandom::randomGenerator()->SetSeed(0);
        int events = estimate_bs[i] + estimate_bd[i] + estimate_rare[i] + estimate_comb[i];
        data_i = total_pdf_i[i]->generate(*ws_->var("Mass"), events);
      }
      data_map.insert(make_pair(channelName.str(), data_i));
      channel->defineType(channelName.str().c_str(), i);
    }
//    cout << " >>>>>>>>>>>> "<< endl;
//    channel->Print();
    global_data = new RooDataSet("global_data", "data for all channels", *ws_->var("Mass"), Index(*channel), Import(data_map));
  }
  else {
    int i = atoi(ch_s_.c_str());
    RooDataSet* data_i;
    if (!random) {
      data_i = new RooDataSet("data", "data", *ws_->var("Mass"));
      FillRooDataSet(tree, data_i, ws_->var("Mass"), i);
    }
    else {
      RooRandom::randomGenerator()->SetSeed(0);
      int events = estimate_bs[i] + estimate_bd[i] + estimate_rare[i] + estimate_comb[i];
      data_i = total_pdf_i[i]->generate(*ws_->var("Mass"), events);
    }
    global_data = data_i;
  }
  ws_->import(*global_data);
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
  string root_s = "output/fit_ws_bdt_";
  string tail_s;
  if (SM_)             tail_s = "_SM.root";
  else if (bd_constr_) tail_s = "_BdConst.root";
  else                 tail_s = ".root";
  if (simul_) {
    for (int i = 0; i < channels; i++) {  // load each ws, containg pdf
      ostringstream input_oss;
      input_oss << root_s << i << tail_s;
      ws_file_input[i] = new TFile(input_oss.str().c_str());
      if (!ws_file_input[i]) {cout << input_oss.str().c_str() << " does not exist" << endl; exit(EXIT_FAILURE);}
      ws_input[i] = (RooWorkspace*)ws_file_input[i]->Get("ws");
      if (!ws_input[i]) {cout << "ws does not exist" << endl; exit(EXIT_FAILURE);}
    }
  }
  else {
    ostringstream input_oss;
    input_oss << root_s << ch_s_ << tail_s;
    ws_file_input[0] = new TFile(input_oss.str().c_str());
    if (!ws_file_input[0]) {cout << input_oss.str().c_str() << " does not exist" << endl; exit(EXIT_FAILURE);}
    ws_input[0] = (RooWorkspace*)ws_file_input[0]->Get("ws");
    if (!ws_input[0]) {cout << "ws does not exist" << endl; exit(EXIT_FAILURE);}
  }
}

void pdf_fitData::make_pdf() {
  cout << "making simultaneous pdf" << endl;
  total_pdf_i.resize(channels);
  if (simul_) {
    for (int i = 0; i < channels; i++) {
      if (random) {
        ws_input[i]->var("N_bs")->setVal(estimate_bs[i]);
        ws_input[i]->var("N_bd")->setVal(estimate_bd[i]);
        ws_input[i]->var("N_rare")->setVal(estimate_rare[i]);
        ws_input[i]->var("N_comb")->setVal(estimate_comb[i]);
        if (bd_constr_) ws_input[i]->var("Bd_over_Bs")->setVal(estimate_bd[i] / estimate_bs[i]);
      }
      changeName(ws_input[i], i);
      total_pdf_i[i] = (RooAbsPdf*)ws_input[i]->pdf(Form("pdf_ext_total_%d", i));
      simul_pdf->addPdf(*total_pdf_i[i], Form("channel_%d", i));
    }

    ws_->import(*simul_pdf);
    ws_->Print();
  }
  else {
    if (random) {
      int i = atoi(ch_s_.c_str());
      ws_input[0]->var("N_bs")->setVal(estimate_bs[i]);
      ws_input[0]->var("N_bd")->setVal(estimate_bd[i]);
      ws_input[0]->var("N_rare")->setVal(estimate_rare[i]);
      ws_input[0]->var("N_comb")->setVal(estimate_comb[i]);
      if (bd_constr_) ws_input[0]->var("Bd_over_Bs")->setVal(estimate_bd[i] / estimate_bs[i]);
    }
    total_pdf_i[0] = (RooAbsPdf*)ws_input[0]->pdf("pdf_ext_total");
    ws_ = new RooWorkspace("ws", "fit workspace results");
    ws_->import(*total_pdf_i[0]);
    ws_->Print();
  }
}

void pdf_fitData::save() {
  ostringstream output_ws;
  output_ws << "./output/ws_" << meth_;
  if (simul_) output_ws << "_simul" << channels;
  else output_ws << "_" << ch_s_;
  if (SM_)             output_ws << "_SM";
  else if (bd_constr_) output_ws << "_BdConst";
  output_ws << ".root";
  ws_->SaveAs(output_ws.str().c_str());

}

void pdf_fitData::significance(int meth) {

  if (meth == 0) sig_hand();
  else if (meth == 1) sig_plhc();
  else if (meth == 2) sig_plhts();
  else return;
}

void pdf_fitData::sig_hand() {
  /// by hand
  Double_t minNLL = RFR->minNll();
  TIterator* vars_it;
  RooRealVar *arg_var = 0;
  RooArgSet *all_vars;
  if (simul_) all_vars = simul_pdf->getVariables();
  else all_vars = total_pdf_i[0]->getVariables();
  vars_it = all_vars->createIterator();
  size_t found;
  while ( (arg_var = (RooRealVar*)vars_it->Next())) {
    string name(arg_var->GetName());
    found = name.find("N_bs");
    if (found!=string::npos) {
      arg_var->setVal(0.0);
      arg_var->setConstant(1);
    }
  }
  fit_pdf();
  Double_t newNLL = RFR->minNll();
  cout << "H1 minNLL = " << minNLL << endl;
  cout << "H0 minNLL = " << newNLL << endl;
  Double_t deltaLL = newNLL - minNLL;
  Double_t signif = deltaLL>0 ? sqrt(2*deltaLL) : -sqrt(-2*deltaLL) ;
  cout << "significance (by hand) = " << signif << endl << endl;
}

void pdf_fitData::sig_plhc() {
  using namespace RooStats;
  TIterator* vars_it;
  RooRealVar *arg_var = 0;
  RooArgSet *all_vars;
  if (simul_) all_vars = simul_pdf->getVariables();
  else all_vars = total_pdf_i[0]->getVariables();
  vars_it = all_vars->createIterator();
  while ( (arg_var = (RooRealVar*)vars_it->Next())) {
    string name(arg_var->GetName());
    size_t found;
    found = name.find("N_bs");
    if (found!=string::npos) {
      arg_var->setConstant(0);
    }
  }
  RooWorkspace* ws1 = new RooWorkspace("ws1", "ws1");
  ModelConfig model;
  RooArgSet poi;
  if (simul_) {
    ws1->import(*simul_pdf);
    model.SetWorkspace(*ws1);
    model.SetPdf(*ws1->pdf("simul_pdf"));
    //model.SetPdf(*simul_ws->pdf("pdf_ext_total_0"));
    for (int i = 0; i < channels; i++) {
      ostringstream name;
      name << "N_bs_" << i;
      poi.add(*ws1->var(name.str().c_str()));
      poi.setRealValue(name.str().c_str(), 0);
    }
  }
  else {
    ws1->import(*total_pdf_i[0]);
    model.SetWorkspace(*ws1);
    model.SetPdf(*ws1->pdf("pdf_ext_total"));
    //model.SetPdf(*simul_ws->pdf("pdf_ext_total_0"));
    poi.add(*ws1->var("N_bs"));
    poi.setRealValue("N_bs", 0);
  }
  ProfileLikelihoodCalculator plc;
  plc.SetData(*global_data);
  plc.SetModel(model);
  plc.SetNullParameters(poi);
  HypoTestResult* htr = plc.GetHypoTest();
  cout << "ProfileLikelihoodCalculator: The p-value for the null is " << htr->NullPValue() << "; The significance for the null is " << htr->Significance() << endl;
}

void pdf_fitData::sig_plhts() {
  using namespace RooStats;
  RooWorkspace* ws2 = new RooWorkspace("ws2", "ws2");
  if (simul_) ws2->import(*simul_pdf);
  else ws2->import(*total_pdf_i[0]);
  ws2->import(*global_data);
  if (simul_) ws2->defineSet("obs", "Mass,channel");
  else ws2->defineSet("obs", "Mass");
  ostringstream name_poi;
  if (simul_) {
    for (int i = 0; i < channels; i++) {
      if (i != 0) name_poi << ",";
      name_poi << "N_bs_" << i;
    }
  }
  else name_poi << "N_bs";
  ws2->defineSet("poi", name_poi.str().c_str());

  ModelConfig* H0 = new ModelConfig("H0", "background only hypothesis", ws2);
  if (simul_) H0->SetPdf(*ws2->pdf("simul_pdf"));
  else H0->SetPdf(*ws2->pdf("pdf_ext_total"));
  H0->SetParametersOfInterest(*ws2->set("poi"));
  H0->SetObservables(*ws2->set("obs"));
  if (simul_) {
    for (int i = 0; i < channels; i++) {
      ostringstream name;
      name << "N_bs_" << i;
      ws2->var(name.str().c_str())->setVal(0.0);
    }
  }
  else ws2->var("N_bs")->setVal(0.0);
  H0->SetSnapshot(*ws2->set("poi"));

  ModelConfig* H1 = new ModelConfig("H1", "background + signal hypothesis", ws2);
  if (simul_) H1->SetPdf(*ws2->pdf("simul_pdf"));
  else H1->SetPdf(*ws2->pdf("pdf_ext_total"));
  H1->SetParametersOfInterest(*ws2->set("poi"));
  H1->SetObservables(*ws2->set("obs"));
  parse_estimate();
  if (simul_) {
    for (int i = 0; i < channels; i++) {
      ostringstream name;
      name << "N_bs_" << i;
      ws2->var(name.str().c_str())->setVal(estimate_bs[i]);
    }
  }
  else {
    ws2->var("N_bs")->setVal(estimate_bs[atoi(ch_s_.c_str())]);
  }
  if (bd_constr_) ws2->var("Bd_over_Bs")->setVal(0.1);
  H1->SetSnapshot(*ws2->set("poi"));

  ws2->import(*H0);
  ws2->import(*H1);
  ws2->Print();

  string name_of_pdf = "simul_pdf";
  if (!simul_) name_of_pdf = "pdf_ext_total";
  ProfileLikelihoodTestStat pl_ts(*ws2->pdf(name_of_pdf.c_str()));
//  pl_ts.SetPrintLevel(3);
//  pl_ts.EnableDetailedOutput(true, true);
//  pl_ts.SetOneSidedDiscovery(true);

  ToyMCSampler *mcSampler_pl = new ToyMCSampler(pl_ts, 1000);
  FrequentistCalculator frequCalc(*global_data , *H1,*H0, mcSampler_pl); // null = bModel interpreted as signal, alt = s+b interpreted as bkg
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  HypoTestResult *htr_pl = frequCalc.GetHypoTest();
  RooMsgService::instance().cleanup();
  htr_pl->Print();

  cout << "ProfileLikelihoodTestStat + frequentist: The p-value for the null is " << htr_pl->NullPValue() << "; The significance for the null is " << htr_pl->Significance() << endl;

}
