/*
 * File:   main_toyMC.cc
 * Author: lucamartini
 *
 * Created on 24 aprile 2012, 18.25
 */
#include "pdf_analysis.h"
#include "pdf_fitData.h"

#include <stdio.h>
#include <stdlib.h>

#include "CommonFun.h"

#include "TFile.h"
#include "TTree.h"

#include "RooDataSet.h"
#include "RooWorkspace.h"


/*
 *
 */

void changeName (RooWorkspace *ws, int str) {
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

void FillRooDataSet(TTree* tree, RooDataSet* dataset, RooRealVar *Mass, int ch_i) {
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

int main(int argc, char** argv) {

  parse_options(argc, argv);
  if (!input || (SM && bd_const)) help();

  TFile* input_data_f = new TFile(input_name.c_str());
  TTree* data_t = (TTree*)input_data_f->Get(meth.c_str());
  if (!data_t) {cout << "no tree called " << meth.c_str() << endl; return (EXIT_FAILURE);}

  if (simul) {
    pdf_fitData* fitdata = new pdf_fitData(false, "bdt", "0", "all", SM, bd_const, data_t, inputs);
    fitdata->make_dataset();
    fitdata->make_pdf_input();
    fitdata->make_pdf();
    fitdata->simul_pdf->fitTo(*fitdata->global_data, Extended(1));
    fitdata->print_each_channel();
    return (EXIT_SUCCESS);
  }

  //*******************
  if (!simul) { /// ancient for cross checking

    vector < TFile*> ws_file_input(inputs);
    vector < RooWorkspace*> ws_input(inputs);
    string root_s = "output/fit_ws_bdt_";
    string tail_s;
    if (SM) tail_s = "_SM.root";
    else if (bd_const) {cout << "can't constrain bd with only one channel; use also -simul # option" << endl; exit(EXIT_FAILURE);}
    else tail_s = ".root";
    for (int i = 0; i < inputs; i++) {  // load each ws, containg pdf
      ostringstream input_oss;
      input_oss << root_s << i << tail_s;
      ws_file_input[i] = new TFile(input_oss.str().c_str());
      if (!ws_file_input[i]) {cout << input_oss.str().c_str() << " does not exist" << endl; return (EXIT_FAILURE);}
      ws_input[i] = (RooWorkspace*)ws_file_input[i]->Get("ws");
      if (!ws_input[i]) {cout << "ws does not exist" << endl; return (EXIT_FAILURE);}
      ws_input[i]->Print();
    }
    RooDataSet* RealData = new RooDataSet("RealData", "RealData", *ws_input[ch_i]->var("Mass"));
    FillRooDataSet(data_t, RealData, ws_input[0]->var("Mass"), ch_i);
    pdf_fitData* fitdata = new pdf_fitData(false, meth.c_str(), ch_s.c_str());
    fitdata->set_ws(ws_input[ch_i]);
    fitdata->SM_ = SM;
    fitdata->fit_pdf("total", RealData, true);
    if (print) fitdata->print();
    return (EXIT_SUCCESS);
  }

  return (EXIT_SUCCESS);

}


