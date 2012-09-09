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

void FillRooDataSet(TTree* tree, RooDataSet* dataset, RooRealVar *Mass, int ch_i, string cuts_f) {
  int events = 0;
  if (!strcmp(tree->GetName(), "bdt") || !strcmp(tree->GetName(), "cnc")) {
    double m1eta, m2eta, m;
    tree->SetBranchAddress("m1eta", &m1eta);
    tree->SetBranchAddress("m2eta", &m2eta);
    tree->SetBranchAddress("m", &m);
    for (int i = 0; i < tree->GetEntries(); i++){
      tree->GetEntry(i);
      if (m > 4.9 && m < 5.9) {
        Mass->setVal(m);
        if (ch_i == -1) {
          dataset->add(*Mass);
          events++;
        }
        else if (ch_i == 0) {
          if (fabs(m1eta) < 1.4 && fabs(m2eta) < 1.4) {
            dataset->add(*Mass);
            events++;
          }
        }
        else if (ch_i == 1) {
          if (fabs(m1eta) > 1.4 || fabs(m2eta) > 1.4) {
            dataset->add(*Mass);
            events++;
          }
        }
        else {
          cout << "wrong channel: " << ch_i << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if (!strcmp(tree->GetName(), "T")) {
    float mlp_0_cut = -1, mlp_1_cut = -1;
    if (cuts_f.compare("no")) {
      FILE *cut_file = fopen(cuts_f.c_str(), "r");
      char buffer[1024];
      char cutName[1024];
      float cut;
      while (fgets(buffer, sizeof(buffer), cut_file)) {
        if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
        if (buffer[0] == '#') continue;
        sscanf(buffer, "%s %f", cutName, &cut);
        if (!strcmp(cutName, "NN_0")) mlp_0_cut = cut;
        if (!strcmp(cutName, "NN_1")) mlp_1_cut = cut;
      }
      cout << "NN_0 = " <<  mlp_0_cut << "  ";
      cout << "NN_1 = " <<  mlp_1_cut << endl;
    }
    else {
      mlp_0_cut = 1.0036;
      mlp_1_cut = 1.0041;
    }
    if (mlp_0_cut == -1 || mlp_1_cut == -1) { cout << "wrong parsing cut file" << endl; exit(EXIT_FAILURE);}

    Float_t mass, eta, mlp_0, mlp_1;
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("mlp_0", &mlp_0);
    tree->SetBranchAddress("mlp_1", &mlp_1);
    for (int i = 0; i < tree->GetEntries(); i++){
      tree->GetEntry(i);
      if (mass > 4.9 && mass < 5.9) {
        Mass->setVal(mass);
        if (ch_i == -1) {
          if (mlp_0 > mlp_0_cut && mlp_0 < 2.0) {
            dataset->add(*Mass);
            events++;
          }
        }
        else if (ch_i == 0) {
          if (fabs(eta) < 1.4) {
            if (mlp_0 > mlp_0_cut && mlp_0 < 2.0) {
              dataset->add(*Mass);
              events++;
            }
          }
        }
        else if (ch_i == 1) {
          if (fabs(eta) > 1.4) {
            if (mlp_1 > mlp_1_cut && mlp_1 < 2.0) {
              dataset->add(*Mass);
              events++;
            }
          }
        }
        else {
          cout << "wrong channel: " << ch_i << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  cout << "total events = " << events << endl;
}

int main(int argc, char** argv) {

  parse_options(argc, argv);
  if (SM && bd_const) help();
  if (simul && channel) help();

  TFile* input_data_f = new TFile(input_name.c_str());
  TTree* data_t = 0;
  if (input) {
    data_t = (TTree*)input_data_f->Get(tree_name.c_str());
    if (!data_t) {cout << "no tree called " << tree_name.c_str() << endl; return EXIT_FAILURE;}
  }

  pdf_fitData* fitdata = new pdf_fitData(false,  inputs, input_estimates, cuts_f, meth, "all", SM, bd_const, data_t, simul, ch_s);
  if (pee) fitdata->pee = true;
  fitdata->initialize();
  fitdata->make_pdf_input();
  fitdata->make_pdf();
  if (strcmp(rare_f.c_str(),"no")) fitdata->set_rare_normalization(rare_f, true);

  fitdata->make_dataset();
  fitdata->fit_pdf();
  if (print) {
    if (simul) fitdata->print_each_channel();
    else fitdata->print();
  }
  fitdata->significance(sig_meth);
  fitdata->save();

  return EXIT_SUCCESS;
}


