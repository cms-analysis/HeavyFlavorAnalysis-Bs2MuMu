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
  if (!input || !estimate || !pdf || !method || !channel) help();
  
  TFile* input_f = new TFile(input_name.c_str());
  RooWorkspace* ws = (RooWorkspace*)input_f->Get("ws");
  
  pdf_toyMC toy1(input_estimates, false, meth, ch_s);
  toy1.set_ws(ws);
 // toy1.unset_constant();
  toy1.generate(NExp, pdf_toy);
  
  return (EXIT_SUCCESS);
}

