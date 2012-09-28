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

  pdf_fitData* fitdata = new pdf_fitData(false,  inputs, input_estimates, cuts_f, meth, "all", SM, bd_const, data_t, simul, pee, bdt_fit, ch_s, sig_meth);
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
  fitdata->significance();
  fitdata->save();
  return EXIT_SUCCESS;
}


