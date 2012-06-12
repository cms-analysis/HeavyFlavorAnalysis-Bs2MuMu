#include "pdf_fitData.h"

pdf_fitData::pdf_fitData(bool print, string meth, string ch_s, string range, bool SM, bool bd_constr, TTree *input_tree, int inputs): pdf_analysis(print, meth, ch_s, range, SM, bd_constr) {
  tree = input_tree;

  channels = inputs;
  channel = new RooCategory("channel", "channel categories");
  simul_pdf = new RooSimultaneous("simul_pdf", "simultaneous fit for all channels", *channel);

  ws_file_input.resize(channels);
  ws_input.resize(channels);

}

void pdf_fitData::print(string output) {

  RooPlot *rp = ws_->var("Mass")->frame();

  rds_->plotOn(rp, Binning(40));

  ws_->pdf(pdf_name.c_str())->plotOn(rp, LineColor(kBlue), Range(range_.c_str()), LineWidth(3));
  ws_->pdf(pdf_name.c_str())->paramOn(rp, Layout(0.50, 0.9, 0.9));

  // components
  RooArgSet * set = ws_->pdf(pdf_name.c_str())->getComponents();
  TIterator* it = set->createIterator();
  TObject* var_Obj = 0;
  while((var_Obj = it->Next())){
    string name = var_Obj->GetName();
    if (name != pdf_name) {
      if (name == "pdf_bs") ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kRed),          LineStyle(1), DrawOption("F"), FillColor(kRed), FillStyle(3001), LineWidth(3), Range(range_.c_str()));
      if (name == "pdf_bd") ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kViolet - 4),   LineStyle(1), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3144), LineWidth(3), Range(range_.c_str()));
      if (name == "pdf_comb") ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kBlue - 5),   LineStyle(2)/*, DrawOption("F"), FillColor(kBlue - 5), FillStyle(3001)*/, LineWidth(3), Range(range_.c_str()));
      if (name == "pdf_rare") ws_->pdf(pdf_name.c_str())->plotOn(rp, Components(*ws_->pdf(var_Obj->GetName())), LineColor(kGreen - 7), LineStyle(1)/*, DrawOption("F"), FillColor(kGreen - 7), FillStyle(3001)*/, LineWidth(2), Range(range_.c_str()));
    }
  }

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  rp->Draw();
  string address = "fig/data_" + pdf_name + "_" + meth_ + "_" + ch_s_ + output;
  if (SM_) address += "_SM";
  c->Print( (address + ".gif").c_str());
  c->Print( (address + ".pdf").c_str());
  delete rp;
  delete c;
  return;
}

void pdf_fitData::print_each_channel() {

  for (int i = 0; i < channels; i++) {
    RooPlot* final_p = ws_->var("Mass")->frame(Bins(20));

    global_data->plotOn(final_p, Cut( Form("channel==channel::channel_%d", i)));

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
    output << "fig/data_simul_" << i;
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
    if (ch_i == 0) {
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
  return;
}

void pdf_fitData::make_dataset() {

  map<string, RooDataSet* > data_map;

  for (int i = 0; i < channels; i++) {
    ostringstream channelName;
    channelName << "channel_" << i;
    RooDataSet* data_i = new RooDataSet( Form("data_%i", i), "data i", *ws_->var("Mass"));
    FillRooDataSet(tree, data_i, ws_->var("Mass"), i);
    data_map.insert(make_pair(channelName.str(), data_i));
    channel->defineType(Form("channel_%d", i));
  }

  global_data = new RooDataSet("global_data", "data for all channels", *ws_->var("Mass"), Index(*channel), Import(data_map));

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
  return;
}

void pdf_fitData::make_pdf_input() {
  string root_s = "output/fit_ws_bdt_";
  string tail_s;
  if (SM_) tail_s = "_SM.root";
  else if (bd_constr_) tail_s = "_BdConst.root";
  else tail_s = ".root";
  for (int i = 0; i < channels; i++) {  // load each ws, containg pdf
    ostringstream input_oss;
    input_oss << root_s << i << tail_s;
    ws_file_input[i] = new TFile(input_oss.str().c_str());
    if (!ws_file_input[i]) {cout << input_oss.str().c_str() << " does not exist" << endl; exit(EXIT_FAILURE);}
    ws_input[i] = (RooWorkspace*)ws_file_input[i]->Get("ws");
    if (!ws_input[i]) {cout << "ws does not exist" << endl; exit(EXIT_FAILURE);}
    ws_input[i]->Print();
  }
}

void pdf_fitData::make_pdf() {

  vector <RooAbsPdf*> total_pdf_i(channels);
  for (int i = 0; i < channels; i++) {
    changeName(ws_input[i], i);
    total_pdf_i[i] = (RooAbsPdf*)ws_input[i]->pdf(Form("pdf_ext_total_%d", i))->Clone();
    simul_pdf->addPdf(*total_pdf_i[i], Form("channel_%d", i));
  }
}
