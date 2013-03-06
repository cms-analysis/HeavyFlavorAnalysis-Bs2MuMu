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

  pdf_fitData* fitdata = new pdf_fitData(false, input_estimates, "all", BF, SM, bd_const, inputs, (!simul_all) ? inputs_bdt : 1, inputs_all, pee, bdt_fit, ch_s, sig_meth, asimov, syst, randomsyst, rare_constr, NExp, Bd, years_opt);
  fitdata->make_pdf_input(input_ws);
  fitdata->make_pdf();
  fitdata->set_final_pdf();
  if (hack_semi2011) fitdata->hack_ws("output/frozen/ws_simul4_bdt_BF2_PEE.root");
  fitdata->setnewlumi();
  fitdata->set_syst();
  fitdata->define_dataset();

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
      cout << "opening " << inputname << endl;
      inputnamess.push_back(inputname_s);
    }
    int years, years_i;
    if (years_opt == "all") years = 2;
    else {
    	years = 1;
    	years_i = atoi(years_opt.c_str());
    }
    if (inputnamess.size() == 1 && years == 2) {
    	cout << "years = " << years << " but there is only one file in " << input_name << endl;
    	return 1;
    }
    for (unsigned int i = 0; i < inputnamess.size(); i++) {
    	int y = (years == 2) ? i : years_i;
      TFile* input_data_f = new TFile(inputnamess[i].c_str(), "UPDATE");
      TTree* data_t = 0;
      if (input_data_f) {
        data_t = (TTree*)input_data_f->Get(tree_name.c_str());
        if (!data_t) {cout << "no tree called " << tree_name.c_str() << endl; return EXIT_FAILURE;}
      }
      fitdata->make_dataset(cuts_f_b, cuts_v, cuts, data_t, y);
    }
  }

  if (!asimov) {
    fitdata->fit_pdf();
    if (print) {
      if (simul) {
        fitdata->print_each_channel();
        if (bdt_fit) fitdata->print_each_channel("bdt");
      }
      else fitdata->print();
    }
    fitdata->extract_N_inRanges();
  }
  fitdata->proof = proof;
  if (SMIsNull) fitdata->SMIsNull = true;
  fitdata->reset_minmax();
  fitdata->significance();
  if (LLprofile) fitdata->profile_NLL();
  //fitdata->save();
  return EXIT_SUCCESS;
}


