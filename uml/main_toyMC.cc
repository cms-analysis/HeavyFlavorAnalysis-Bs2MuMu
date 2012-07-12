/* 
 * File:   main_toyMC.cc
 * Author: lucamartini
 *
 * Created on 24 aprile 2012, 18.25
 */
#include "pdf_analysis.h"
#include "pdf_toyMC.h"

#include <stdio.h>
#include <stdlib.h>

#include "CommonFun.h"

#include "TFile.h"

#include "RooWorkspace.h"

/*
 * 
 */

int main(int argc, char** argv) {

  parse_options(argc, argv);
  if (!input || !estimate || !pdf) help();
  if (!(!bias_s.compare("no") || !bias_s.compare("c+") || !bias_s.compare("c-") || !bias_s.compare("p+") || !bias_s.compare("p-"))) { cout << "I don't understand what to bias: please enter -bias c+, c-, p+, p-" << endl; exit(EXIT_SUCCESS);}

  TFile* input_f = new TFile(input_name.c_str());
  RooWorkspace* ws = (RooWorkspace*)input_f->Get("ws");
  ws->Print();

  parse_input(input_name);
  if (!pdf_test_b) pdf_test = pdf_toy;
  pdf_toyMC toy1(print, inputs, input_estimates, meth, "all", SM, bd_const, 0, bias_s, simul, ch_s);
  toy1.set_ws(ws);
  toy1.make_pdf_input();
  toy1.make_pdf();
  toy1.make_dataset();
  if (roomcs) toy1.mcstudy(NExp, pdf_toy);
  if (!roomcs) toy1.generate(NExp, pdf_toy, pdf_test);
  
  delete input_f;

  return EXIT_SUCCESS;
}

