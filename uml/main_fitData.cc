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

  pdf_fitData* fitdata = new pdf_fitData(false, inputs, inputs_bdt, input_estimates, meth, "all", BF, SM, bd_const, simul, simul_bdt, pee, bdt_fit, ch_s, sig_meth, asimov);
  fitdata->initialize();
  fitdata->make_pdf_input();
  fitdata->make_pdf();
  fitdata->setnewlumi();
  //if (strcmp(rare_f.c_str(),"no")) fitdata->set_rare_normalization(rare_f, true);

  vector <double> cuts_v(inputs, -10);
  if (cuts_f_b) cuts_v = cut_bdt_file();

  vector <string> inputnamess(0);
  FILE *file = fopen(input_name.c_str(), "r");
  if (!file) {
    cout << "no file name, making random" << endl;
    fitdata->random = true;
    fitdata->make_dataset(cuts_f_b, cuts_v, cuts, 0, 0);
  }
  else {
    fitdata->random = false;
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), file)) {
      if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
      if (buffer[0] == '#') continue;
      char inputname[1024];
      sscanf(buffer, "%s", inputname);
      string inputname_s(inputname);
      inputnamess.push_back(inputname_s);
    }
    for (unsigned int i = 0; i < inputnamess.size(); i++) {
      TFile* input_data_f = new TFile(inputnamess[i].c_str());
      TTree* data_t = 0;
      if (input_data_f) {
        data_t = (TTree*)input_data_f->Get(tree_name.c_str());
        if (!data_t) {cout << "no tree called " << tree_name.c_str() << endl; return EXIT_FAILURE;}
      }
      fitdata->make_dataset(cuts_f_b, cuts_v, cuts, data_t, i);
    }
  }

  if (!asimov) {
    fitdata->fit_pdf();
    if (print) {
      if (simul || simul_bdt) fitdata->print_each_channel();
      else fitdata->print();
    }
    //  fitdata->BF("./input/anaBmm.plotResults.default-11.tex", "input/external_numbers.txt");
  }
  //if (systname!="") fitdata->parse_systematics(systname);
  fitdata->significance();
  fitdata->save();
  return EXIT_SUCCESS;
}


